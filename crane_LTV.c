#include "osqp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_timing.h"


//parts to modify: A,B,Q,R,P,q; max/min ; NX NU DT N ; initial and desired point



/* ---------------------------------------------------------
 * Problem Constants & Parameters
 * --------------------------------------------------------- */
#define NX 6    // Number of States (CHANGE if system dimension changes)
#define NU 2    // Number of Controls (CHANGE if input dimension changes)
#define N 100   // Prediction Horizon steps
#define ND  6    // Disturbance dimension (matches NX here)
#define NH  3    // Number of nonlinear inequality constraints h(x)
#define DT 0.020 // Sampling time (s)
#define Tsim 20.0// Full simulation time (s)
#define R_BAND 30 // Truncation bandwidth for condensed band (controls coupling span)



/* ---------------------------------------------------------
 * Sparse Condensed (SC) OSQP dimensions
 * --------------------------------------------------------- */
// Decision variables: virtual inputs z only
OSQPInt n_vars = N * NU;
// Constraints: states + inputs + nonlinear per step
OSQPInt n_cons = (N * NX) + (N * NU) + (N * NH);

typedef struct {
    OSQPInt *i;
    OSQPInt *p;
    OSQPFloat *x;
    OSQPInt nz;
    OSQPInt n_cols;
} CscBuilder; // lightweight CSC assembler

static void mat_copy(OSQPFloat *dst, const OSQPFloat *src, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) dst[i] = src[i];
}

static void mat_eye(OSQPFloat *m, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) m[i] = 0.0;
    for (OSQPInt i = 0; i < (rows < cols ? rows : cols); i++) m[i * cols + i] = 1.0;
}

static void mat_mul(const OSQPFloat *A, const OSQPFloat *B, OSQPFloat *C, OSQPInt m, OSQPInt k, OSQPInt n) {
    for (OSQPInt i = 0; i < m; i++) {
        for (OSQPInt j = 0; j < n; j++) {
            OSQPFloat acc = 0.0f;
            for (OSQPInt t = 0; t < k; t++) acc += A[i * k + t] * B[t * n + j];
            C[i * n + j] = acc;
        }
    }
}

static void mat_add(OSQPFloat *A, const OSQPFloat *B, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] += B[i];
}

static void mat_scale(OSQPFloat *A, OSQPFloat s, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] *= s;
}

static void mat_set(OSQPFloat *A, OSQPFloat v, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] = v;
}

static void mat_add_scaled(OSQPFloat *A, const OSQPFloat *B, OSQPFloat s, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] += s * B[i];
}

static void mat_power(const OSQPFloat *A, OSQPInt n, OSQPFloat *out, OSQPInt dim) {
    // out = A^n (n >= 1) using repeated squaring (small dim so simple loop)
    OSQPFloat temp[dim * dim];
    mat_copy(out, A, dim, dim);
    for (OSQPInt p = 1; p < n; p++) {
        mat_mul(out, A, temp, dim, dim, dim);
        mat_copy(out, temp, dim, dim);
    }
}

static void mat_transpose(const OSQPFloat *A, OSQPFloat *AT, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows; i++) {
        for (OSQPInt j = 0; j < cols; j++) AT[j * rows + i] = A[i * cols + j];
    }
}

// Euclidean 2-norm helper
static OSQPFloat vec_norm(const OSQPFloat *v, OSQPInt len) {
    OSQPFloat acc = 0.0f;
    for (OSQPInt i = 0; i < len; i++) acc += v[i] * v[i];
    return sqrtf(acc);
}

// Nonlinear crane dynamics f(x,u); REPLACE when switching to another plant
static void crane_dynamics(const OSQPFloat x[NX], const OSQPFloat u[NU], OSQPFloat f_out[NX]) {
    const OSQPFloat g = 9.81f;
    f_out[0] = x[1];
    f_out[1] = u[0];
    f_out[2] = x[3];
    f_out[3] = u[1];
    f_out[4] = x[5];
    f_out[5] = -((g * sinf(x[4]) + cosf(x[4]) * u[0] + 2.0f * x[3] * x[5]) / x[2]);
}

/* Pseudo-inverse via normal equations: pinv(B) = (B^T B + eps I)^{-1} B^T
static void pseudo_inverse(const OSQPFloat *B, OSQPFloat *pinv, OSQPInt rows, OSQPInt cols) {
    // Allocate temporary buffers dynamically to handle arbitrary cols
    OSQPFloat *BtB = malloc(cols * cols * sizeof(OSQPFloat));
    mat_set(BtB, 0.0, cols, cols);
    for (OSQPInt i = 0; i < cols; i++) {
        for (OSQPInt j = 0; j < cols; j++) {
            OSQPFloat acc = 0.0;
            for (OSQPInt k = 0; k < rows; k++) acc += B[k * cols + i] * B[k * cols + j];
            BtB[i * cols + j] = acc;
        }
    }
    OSQPFloat eps = 1e-6f;
    for (OSQPInt i = 0; i < cols; i++) BtB[i * cols + i] += eps;

    // Gauss-Jordan inversion
    OSQPFloat *aug = malloc(cols * 2 * cols * sizeof(OSQPFloat));
    for (OSQPInt i = 0; i < cols; i++) {
        for (OSQPInt j = 0; j < cols; j++) aug[i * 2 * cols + j] = BtB[i * cols + j];
        for (OSQPInt j = 0; j < cols; j++) aug[i * 2 * cols + cols + j] = (i == j) ? 1.0 : 0.0;
    }
    for (OSQPInt i = 0; i < cols; i++) {
        OSQPFloat pivot = aug[i * 2 * cols + i];
        if (fabs(pivot) < 1e-12) pivot = (pivot >= 0 ? 1e-12 : -1e-12);
        for (OSQPInt j = 0; j < 2 * cols; j++) aug[i * 2 * cols + j] /= pivot;
        for (OSQPInt k = 0; k < cols; k++) {
            if (k == i) continue;
            OSQPFloat factor = aug[k * 2 * cols + i];
            for (OSQPInt j = 0; j < 2 * cols; j++) aug[k * 2 * cols + j] -= factor * aug[i * 2 * cols + j];
        }
    }

    OSQPFloat *invBtB = malloc(cols * cols * sizeof(OSQPFloat));
    for (OSQPInt i = 0; i < cols; i++) {
        for (OSQPInt j = 0; j < cols; j++) invBtB[i * cols + j] = aug[i * 2 * cols + cols + j];
    }

    OSQPFloat *Bt = malloc(cols * rows * sizeof(OSQPFloat));
    mat_transpose(B, Bt, rows, cols);
    mat_mul(invBtB, Bt, pinv, cols, cols, rows);

    free(BtB); free(aug); free(invBtB); free(Bt);
}
*/
//update sub x, i(do once for every nonzero entry )
void csc_set(CscBuilder *m, OSQPInt row, OSQPInt col, OSQPFloat val) {
    m->x[m->nz] = val;//store value(x)
    m->i[m->nz] = row;//store row(i)
    m->nz++;//nz plus 1
}

//update sub p(only do when all rows looped through once)
void csc_col_done(CscBuilder *m) {
    m->p[m->n_cols + 1] = m->nz;//store column
    m->n_cols++;//move to next column
}

// Compute stabilizing gain via hardcoded LQR and precompute SC blocks with truncated band
void compute_sc_blocks(OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX],
                       OSQPFloat K[NU][NX], OSQPFloat AK[NX][NX],
                       OSQPFloat P_blocks[], OSQPFloat M_blocks[], OSQPFloat H_blocks[]) {

    // Hardcoded stabilizing discrete-time LQR around target x=[2,0,2,0,0,0] (u = K*x + z)
    // If system matrices or target change, update these gains or recompute offline.
    OSQPFloat K_lqr[NU][NX] = {
        {-2.7510233, -5.6501608, 0.0, 0.0, 17.3512501, 3.4993693},
        {0.0, 0.0, -4.7604404, -4.8105959, 0.0, 0.0}
    };

    for (OSQPInt i = 0; i < NU; i++) {
        for (OSQPInt j = 0; j < NX; j++) K[i][j] = K_lqr[i][j];
    }

    // 1) Closed-loop matrix AK = Ad + Bd*K
    for (OSQPInt i = 0; i < NX; i++) {
        for (OSQPInt j = 0; j < NX; j++) {
            OSQPFloat acc = 0.0;
            for (OSQPInt c = 0; c < NU; c++) acc += Bd[i][c] * K[c][j];
            AK[i][j] = Ad[i][j] + acc;
        }
    }

    // 2) Precompute A_K powers up to R_BAND
    OSQPFloat AK_powers[R_BAND + 1][NX * NX];
    mat_eye(AK_powers[0], NX, NX);
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        mat_mul(AK_powers[i - 1], (OSQPFloat *)AK, AK_powers[i], NX, NX, NX);
    }

    // 3) P_i blocks: P_i = AK^i * Bd
    OSQPFloat Bd_flat[NX * NU];
    for (OSQPInt i = 0; i < NX; i++) {
        for (OSQPInt j = 0; j < NU; j++) Bd_flat[i * NU + j] = Bd[i][j];
    }
    for (OSQPInt i = 0; i < R_BAND; i++) {
        mat_mul(AK_powers[i + 1], Bd_flat, &P_blocks[i * NX * NU], NX, NX, NU);
    }

    // 4) M_i blocks: M_0 = I, M_i = K * AK^{i-1} * Bd for i>=1
    mat_set(M_blocks, 0.0, (R_BAND + 1) * NU, NU);
    for (OSQPInt i = 0; i < NU; i++) M_blocks[i * NU + i] = 1.0; // M_0
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        OSQPFloat KAK[NU * NX];
        mat_mul((OSQPFloat *)K, AK_powers[i - 1], KAK, NU, NX, NX);
        mat_mul(KAK, Bd_flat, &M_blocks[i * NU * NU], NU, NX, NU);
    }

    // 5) H_i blocks: jacobian_H * AK^{i-1} * Bd
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        OSQPFloat CAK[NH * NX];
        mat_mul((OSQPFloat *)jacobian_H, AK_powers[i - 1], CAK, NH, NX, NX);
        mat_mul(CAK, Bd_flat, &H_blocks[(i - 1) * NH * NU], NH, NX, NU);
    }
}
// Linearize crane dynamics around (x_current, u_applied) and discretize via Euler
// IMPORTANT: replace with your system's Jacobians and drift when changing plants.
void linearization( OSQPFloat x_current[NX], OSQPFloat u_applied[NU], 
                                OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat d_lin[NX]) {//Ad Bd can be place holder
    
    // 1. initializations
    OSQPFloat g=9.81;
    

    //precalculate variables
    OSQPFloat t1 = (2*x_current[3]*x_current[5] + u_applied[0]*cos(x_current[4]) + g*sin(x_current[4]))/(x_current[2]*x_current[2]);
    OSQPFloat t2 = -(2*x_current[5])/x_current[2];
    OSQPFloat t3 = -(g*cos(x_current[4]) - u_applied[0]*sin(x_current[4]))/x_current[2];
    OSQPFloat t4 = -(2*x_current[3])/x_current[2];
    OSQPFloat t5 = -cos(x_current[4])/x_current[2];

    // Assign to Ac (df/dx) using the optimized variables
    OSQPFloat Ac[NX][NX] = {
                                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, 
                                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 
                                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}, 
                                {0.0, 0.0, t1, t2, t3, t4}};

    // Assign to Bc (df/du)

    OSQPFloat Bc[NX][NU] = {
                                {0.0, 0.0}, 
                                {1.0, 0.0},
                                {0.0, 0.0}, 
                                {0.0, 1.0},
                                {0.0, 0.0},
                                {t5, 0.0}};
    

    // fc = continuous nonlinear dynamics f(x,u) at linearization point
    OSQPFloat fc[NX];
    fc[0] = x_current[1];
    fc[1] = u_applied[0];
    fc[2] = x_current[3];
    fc[3] = u_applied[1];
    fc[4] = x_current[5];
    fc[5] = -((9.81 * sin(x_current[4]) + cos(x_current[4]) * u_applied[0] + 2 * x_current[3] * x_current[5]) / x_current[2]);

    // Discretize Ad, Bd, and calculate the d_lin offset
    for (int i = 0; i < NX; i++) {
        OSQPFloat Ac_x = 0.0;
        OSQPFloat Bc_u = 0.0;
        for (int j = 0; j < NX; j++) {
            Ad[i][j] = (i == j ? 1.0 : 0.0) + Ac[i][j] * DT;
            Ac_x += Ac[i][j] * x_current[j];
        }
        for (int j = 0; j < NU; j++) {
            Bd[i][j] = Bc[i][j] * DT;
            Bc_u += Bc[i][j] * u_applied[j];
        }
        // This is the missing physics vector
        d_lin[i] = DT * (fc[i] - Ac_x - Bc_u);
    }
    /*                           
    // 3. Discretization (Euler: Ad = I + Ac*dt, Bd = Bc*dt)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NX; j++) {
            Ad[i][j] = (i == j ? 1.0 : 0.0) + Ac[i][j] * DT;
        }
        for (int j = 0; j < NU; j++) {
            Bd[i][j] = Bc[i][j] * DT;
        }
    }*/
}
// Evaluate nonlinear sway constraints h(x) and their Jacobian
// Replace with your problem-specific constraints and gradients
void calc_h_gradient(const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_val[NH]) {//look at grampc dhdx_vect function
   
    // Parameters from pSys in GRAMPC
    OSQPFloat p0 = 0.2;
    OSQPFloat p1 = 1.25;
    OSQPFloat p2 = 0.3;

    // Pre-calculate temp terms
    OSQPFloat x_pos = x_current[0] + x_current[2] * sin(x_current[4]); // Horizontal load position
    OSQPFloat temp = p0 * x_pos;

    // 1. Calculate h(x) value for the RHS bound
    h_val[0] = x_current[2] * cos(x_current[4]) - p0 * (x_pos * x_pos) - p1;
    h_val[1] =x_current[5]-p2;
    h_val[2] =-x_current[5]-p2;
    // 2. Calculate Gradient (Partial derivatives)
    for(int j=0; j<NH; j++) {
        for(int i=0; i<NX; i++) jacobian_H[j][i] = 0.0;
    }
    jacobian_H[0][0]=-2.0*temp;// dh/dx0
    jacobian_H[0][2] = cos(x_current[4]) + temp * sin(x_current[4]);// dh/dx2
    jacobian_H[0][4] = -x_current[2] * sin(x_current[4]) -2.0*temp * x_current[2] * cos(x_current[4]); // dh/dx4 
    jacobian_H[1][5] = 1.0;
    jacobian_H[2][5] = -1.0;
}

// Assemble banded Hessian P and gradient q in CSC form
// recompute_hessian==0: skip Hessian rebuild, only refresh gradient q
OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX],
                       const OSQPFloat K[NU][NX],
                       const OSQPFloat P_blocks[], const OSQPFloat M_blocks[],
                       const OSQPFloat e_free[][NX], OSQPFloat q[], OSQPCscMatrix *P_existing,
                       OSQPInt r_local, OSQPInt recompute_hessian) {

    OSQPInt setup = (P_existing == NULL);
    // Each column couples with at most r_local*NU previous rows + itself
    OSQPInt P_nnz_est = n_vars * (r_local * NU + NU);

    OSQPCscMatrix *P;
    CscBuilder P_build = {0};
    OSQPInt val_idx = 0;
    if (setup) {
        P = malloc(sizeof(OSQPCscMatrix));
        P->x = malloc(P_nnz_est * sizeof(OSQPFloat));
        P->i = malloc(P_nnz_est * sizeof(OSQPInt));
        P->p = malloc((n_vars + 1) * sizeof(OSQPInt));
        P_build.x = P->x; P_build.i = P->i; P_build.p = P->p;
        P_build.nz = 0; P_build.n_cols = 0; P_build.p[0] = 0;
    } else {
        P = P_existing;
    }

    // zero q before filling
    for (OSQPInt i = 0; i < n_vars; i++) q[i] = 0.0f;

    // dense accumulator for Hessian (symmetric)
    static OSQPFloat H_dense[N * NU][N * NU];
    if (recompute_hessian) {
        for (OSQPInt i = 0; i < n_vars; i++) {
            for (OSQPInt j = 0; j < n_vars; j++) H_dense[i][j] = 0.0f;
        }
    }

    // Precompute v_free = K * e_free for each step (baseline around target)
    OSQPFloat v_free[N][NU];
    for (OSQPInt k = 0; k < N; k++) {
        for (OSQPInt i = 0; i < NU; i++) {
            OSQPFloat acc = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) acc += K[i][j] * e_free[k][j];
            v_free[k][i] = acc;
        }
    }

    // Build dense matrix in band and gradient
    for (OSQPInt k1 = 0; k1 < N; k1++) {
        for (OSQPInt j1 = 0; j1 < NU; j1++) {
            OSQPInt col1 = k1 * NU + j1;

            // Gradient
            for (OSQPInt t = k1; t <= N && (t - k1) <= r_local; t++) {
                OSQPInt idx1 = t - k1;
                if (idx1 >= r_local) break;

                const OSQPFloat *weight = (t == N) ? P_diag : Q_diag;
                for (OSQPInt s = 0; s < NX; s++) {
                    q[col1] += P_blocks[idx1 * NX * NU + s * NU + j1] * weight[s] * e_free[t][s];
                }
                if (t < N) {
                    for (OSQPInt u_i = 0; u_i < NU; u_i++) {
                        q[col1] += M_blocks[idx1 * NU * NU + u_i * NU + j1] * R_diag[u_i] * v_free[t][u_i];
                    }
                }
            }

            if (recompute_hessian) {
                // Hessian upper triangle within band
                for (OSQPInt k2 = k1; k2 < N && (k2 - k1) < r_local; k2++) {
                    for (OSQPInt j2 = 0; j2 < NU; j2++) {
                        OSQPInt col2 = k2 * NU + j2;
                        if (col2 < col1) continue;
                        OSQPFloat h_val = 0.0f;
                        for (OSQPInt t = k2; t <= N; t++) {
                            OSQPInt idx1 = t - k1;
                            OSQPInt idx2 = t - k2;
                            if (idx1 >= r_local || idx2 >= r_local) break;

                            const OSQPFloat *weight = (t == N) ? P_diag : Q_diag;
                            for (OSQPInt s = 0; s < NX; s++) {
                                h_val += P_blocks[idx1 * NX * NU + s * NU + j1] * weight[s] * P_blocks[idx2 * NX * NU + s * NU + j2];
                            }
                            if (t < N) {
                                for (OSQPInt u_i = 0; u_i < NU; u_i++) {
                                    h_val += M_blocks[idx1 * NU * NU + u_i * NU + j1] * R_diag[u_i] * M_blocks[idx2 * NU * NU + u_i * NU + j2];
                                }
                            }
                        }
                        H_dense[col1][col2] += h_val;
                        if (col1 != col2) H_dense[col2][col1] += h_val;
                    }
                }
            }
        }
    }

    // CSC Extraction: upper triangular, banded
    if (recompute_hessian) {
        if (setup) {
            for (OSQPInt col = 0; col < n_vars; col++) {
                OSQPInt row_start = (col >= r_local * NU) ? (col - r_local * NU) : 0;
                for (OSQPInt row = row_start; row <= col; row++) {
                    OSQPFloat v = H_dense[row][col];
                    csc_set(&P_build, row, col, v);
                }
                csc_col_done(&P_build);
            }
            P->m = n_vars; P->n = n_vars; P->nz = -1; P->nzmax = P_build.nz;
        } else {
            for (OSQPInt col = 0; col < n_vars; col++) {
                for (OSQPInt idx = P->p[col]; idx < P->p[col + 1]; idx++) {
                    OSQPInt row = P->i[idx];
                    P->x[val_idx++] = H_dense[row][col];
                }
            }
        }
    }
    return P;
}

// Assemble constraint matrix A (states, inputs, nonlinear) in CSC
OSQPCscMatrix* build_A(OSQPCscMatrix* A,
                       const OSQPFloat P_blocks[], const OSQPFloat M_blocks[], const OSQPFloat H_blocks[],
                       OSQPInt r_local) {

    OSQPInt setup = (A == NULL);
    CscBuilder A_build = {0};
    OSQPInt A_nnz_est = N * NU * (r_local * NX + (r_local + 1) * NU + r_local * NH);

    if (setup) {
        A = malloc(sizeof(OSQPCscMatrix));
        A->x = malloc(A_nnz_est * sizeof(OSQPFloat));
        A->i = malloc(A_nnz_est * sizeof(OSQPInt));
        A->p = malloc((n_vars + 1) * sizeof(OSQPInt));
        A_build.x = A->x; A_build.i = A->i; A_build.p = A->p;
        A_build.nz = 0; A_build.n_cols = 0; A_build.p[0] = 0;
    } else {
        A_build.x = A->x; A_build.i = A->i; A_build.p = A->p;
        A_build.nz = 0; A_build.n_cols = 0;
    }

    OSQPInt val_idx = 0;

    for (OSQPInt k = 0; k < N; k++) {
        for (OSQPInt j = 0; j < NU; j++) {
            OSQPInt current_col = k * NU + j;

            // 1. State constraints (top stack)
            for (OSQPInt i = 0; i < r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (k + i) * NX;
                const OSQPFloat *P_i = &P_blocks[i * NX * NU + j];
                for (OSQPInt r = 0; r < NX; r++) {
                    OSQPFloat v = P_blocks[i * NX * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v);
                    else { A->x[val_idx] = v; val_idx++; }
                }
            }

            // 2. Input constraints (middle stack)
            for (OSQPInt i = 0; i <= r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (N * NX) + (k + i) * NU;
                const OSQPFloat *M_i = &M_blocks[i * NU * NU + j];
                for (OSQPInt r = 0; r < NU; r++) {
                    OSQPFloat v = M_blocks[i * NU * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v);
                    else { A->x[val_idx] = v; val_idx++; }
                }
            }

            // 3. Nonlinear constraints (bottom stack)
            for (OSQPInt i = 1; i <= r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (N * NX) + (N * NU) + (k + i) * NH;
                const OSQPFloat *H_i = &H_blocks[(i - 1) * NH * NU + j];
                for (OSQPInt r = 0; r < NH; r++) {
                    OSQPFloat v = H_blocks[(i - 1) * NH * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v);
                    else { A->x[val_idx] = v; val_idx++; }
                }
            }

            if (setup) csc_col_done(&A_build);
        }
    }

    if (setup) {
        A->m = n_cons; A->n = n_vars; A->nz = -1; A->nzmax = A_build.nz;
    }
    return A;
};

// Build lower/upper bounds for state, input, and nonlinear constraints
void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX],
                  const OSQPFloat e_free[][NX], const OSQPFloat x_target[NX], const OSQPFloat u_target[NU], const OSQPFloat K[NU][NX],
                  const OSQPFloat jacobian_H[NH][NX], const OSQPFloat h_vals[NH],
                  const OSQPFloat x_min[NX], const OSQPFloat x_max[NX],
                  const OSQPFloat u_min[NU], const OSQPFloat u_max[NU],
                  OSQPInt r_local) {

    OSQPInt setup = (*l_ptr == NULL);
    if (setup) {
        *l_ptr = malloc(n_cons * sizeof(OSQPFloat));
        *u_ptr = malloc(n_cons * sizeof(OSQPFloat));
    }

    OSQPFloat *l = *l_ptr;
    OSQPFloat *u = *u_ptr;

    // 1) State bounds rows
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt i = 0; i < NX; i++) {
            OSQPFloat x_pred = e_free[n + 1][i] + x_target[i];
            l[n * NX + i] = x_min[i] - x_pred;
            u[n * NX + i] = x_max[i] - x_pred;
        }
    }

    // 2) Input bounds rows
    OSQPInt u_offset_rows = N * NX;
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt i = 0; i < NU; i++) {
            OSQPFloat v_free_val = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) v_free_val += K[i][j] * e_free[n][j];
            OSQPFloat u_pred = u_target[i] + v_free_val;
            l[u_offset_rows + n * NU + i] = u_min[i] - u_pred;
            u[u_offset_rows + n * NU + i] = u_max[i] - u_pred;
        }
    }

    // 3) Nonlinear bounds rows
    OSQPInt h_offset_rows = N * NX + N * NU;
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt ih = 0; ih < NH; ih++) {
            OSQPFloat acc = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) acc += jacobian_H[ih][j] * (e_free[n + 1][j] + x_target[j]);
            l[h_offset_rows + n * NH + ih] = -OSQP_INFTY;
            u[h_offset_rows + n * NH + ih] = OSQP_INFTY; // temporarily relax nonlinear constraints for feasibility check
        }
    }
}


typeTime solver_start, solver_end, iter_start, iter_end;
OSQPFloat solver_total_ms = 0.0f;
OSQPFloat iter_total_ms = 0.0f;
OSQPFloat solver_avg_ms = 0.0f;
OSQPFloat iter_avg_ms = 0.0f;



int main(void) {


    // System Dynamics (Discretized Euler): x_{k+1} = Ad*x_k + Bd*u_k + d_lin
    // NOTE: Ad, Bd, d_lin are updated by linearization(); change linearization() for new plants.
    OSQPFloat Ad[NX][NX] = {0};
    OSQPFloat Bd[NX][NU] = {0};
    OSQPFloat jacobian_H[NH][NX] = {0};
    OSQPFloat h_val[NH] = {0};
    OSQPFloat d_lin[NX]={0};
    OSQPFloat x_lin_pred[NX] = {0};
    const OSQPFloat eta_pri = 0.09; // Primal CMoN threshold (tunable)
    OSQPInt linearized_times = 1;
    // Cost Weights (Derived from your pCost)
    // Q (State cost), R (Control cost), P (Terminal cost)
    //1/2( xQx + uRu + x'Px' ) where x' is desired position
    const OSQPFloat Q_diag[NX] = {1.0*DT, 2.0*DT, 2.0*DT, 1.0*DT, 1.0*DT, 4.0*DT}; //for fast response
    const OSQPFloat R_diag[NU] = {0.05*DT,0.05*DT}; //penalize large input       
    const OSQPFloat P_diag[NX] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};//for damping fater to desired point(optimize?)

    // Desired Setpoints (Target)
    // CHANGE x_current / x_target / u_target for new problems
    OSQPFloat x_current[NX] = { -3.0, 0.0, 2.0, 0.0, 0.0, 0.0 };//initial point
    OSQPFloat x_target[NX] = { 3.0, 0.0, 3.0, 0.0, 0.0, 0.0 };//target point
    OSQPFloat u_applied[NU] = {0.0, 0.0};
    OSQPFloat u_target[NU] = {0.0, 0.0};
    linearization(x_current, u_applied, Ad, Bd, d_lin);
    calc_h_gradient(x_current, jacobian_H, h_val);

    // upper lower bounds
    const OSQPFloat x_min[NX] = {-OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY};
    const OSQPFloat x_max[NX] = {OSQP_INFTY, OSQP_INFTY, OSQP_INFTY, OSQP_INFTY, OSQP_INFTY, OSQP_INFTY};
    const OSQPFloat u_min[NU] = {-2.0, -2.0};
    const OSQPFloat u_max[NU] = { 2.0,  2.0};

    //estimate possible disturbance
    OSQPFloat d_est[ND] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    // SC precomputation containers
    OSQPFloat K[NU][NX] = {0};
    OSQPFloat AK[NX][NX] = {0};
    OSQPFloat P_blocks[NU * NX * R_BAND];
    OSQPFloat M_blocks[NU * NU * (R_BAND + 1)];
    OSQPFloat H_blocks[NH * NU * R_BAND];
    OSQPInt r_local = R_BAND;

    // Free response buffer (absolute and error coordinates)
    OSQPFloat x_free[N + 1][NX];
    OSQPFloat e_free[N + 1][NX];
    OSQPFloat v_free[N][NU];

    // Hessian/gradient containers
    OSQPCscMatrix *P = NULL;
    OSQPFloat *q = malloc(n_vars * sizeof(OSQPFloat));

    // Constraint matrix & bounds
    OSQPCscMatrix *A = NULL;
    OSQPFloat *l = NULL, *u = NULL;

    // Initial SC blocks and matrices
    compute_sc_blocks(Ad, Bd, jacobian_H, K, AK, P_blocks, M_blocks, H_blocks);
    // Steady-state tracking residual d_err = Ad*x_t + Bd*u_t + d_lin - x_t
    OSQPFloat d_err[NX];
    for (OSQPInt i = 0; i < NX; i++) {
        OSQPFloat ax = 0.0f, bu = 0.0f;
        for (OSQPInt j = 0; j < NX; j++) ax += Ad[i][j] * x_target[j];
        for (OSQPInt j = 0; j < NU; j++) bu += Bd[i][j] * u_target[j];
        d_err[i] = ax + bu + d_lin[i] - x_target[i];
    }

    // Free response (with corrected affine shift for non-zero target)
    for (OSQPInt i = 0; i < NX; i++) x_free[0][i] = x_current[i];
    for (OSQPInt k = 0; k < N; k++) {
        for (OSQPInt i = 0; i < NX; i++) {
            // 1. A_K * x_free
            OSQPFloat acc = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) acc += AK[i][j] * x_free[k][j];

            // 2. The missing shift term: Bd * (u_target - K * x_target)
            OSQPFloat shift_term = 0.0f;
            for (OSQPInt c = 0; c < NU; c++) {
                OSQPFloat K_xtarg = 0.0f;
                for (OSQPInt j = 0; j < NX; j++) K_xtarg += K[c][j] * x_target[j];
                shift_term += Bd[i][c] * (u_target[c] - K_xtarg);
            }

            // 3. Total update
            x_free[k + 1][i] = acc + d_lin[i] + shift_term;
        }
    }
    // Convert to error coordinates
    for (OSQPInt k = 0; k <= N; k++) {
        for (OSQPInt i = 0; i < NX; i++) e_free[k][i] = x_free[k][i] - x_target[i];
    }
    P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, NULL, r_local, 1);
    A = build_A(NULL, P_blocks, M_blocks, H_blocks, r_local);
    build_bounds(&l, &u, x_current, e_free, x_target, u_target, K, jacobian_H, h_val, x_min, x_max, u_min, u_max, r_local);

    //setup for osqp solver
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);

    OSQPInt opt_time=0;
    if(opt_time==0){//balanced default
        settings->alpha = 1.6;      // higher relaxation for faster convergence
        settings->scaling= 0;      // stronger scaling to improve conditioning
        settings->verbose = 0;
        settings->adaptive_rho = 1;
        settings->check_termination = 1;
        settings->max_iter = 10000;//make sure is not reached for 100% solved
        settings->warm_starting = 1;//default for grampc
        settings->eps_abs = 1e-3;   // slightly looser tolerances to aid feasibility
        settings->eps_rel = 1e-3;
    }
    

    OSQPSolver *solver;
    OSQPInt setup_status = osqp_setup(&solver, P, q, A, l, u, n_cons, n_vars, settings);
    if (setup_status != 0 || solver == NULL) {
        printf("OSQP setup failed (code %lld)\n", (long long)setup_status);
        return 1;
    }
    
    

    //MPC loop
    int MAX_STEPS = Tsim/DT; 

    
    //write track file
    FILE *f_out = fopen("../mpc_data.csv", "w");
    if (f_out == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    // Write Header (CSV format)
    fprintf(f_out, "step");
    for(int i=0; i<NX; i++) fprintf(f_out, ",x%d", i);
    for(int i=0; i<NU; i++) fprintf(f_out, ",u%d", i);
    fprintf(f_out, "\n");

    printf("Starting MPC Loop (NX=%d, NU=%d)...\n", NX, NU);
    printf("Step |");
    for(int i=0; i<NX; i++) printf("   x%d   ", i);
    printf("|");
    for(int i=0; i<NU; i++) printf("   u%d   ", i);
    printf("\n");
    
    printf("-----|");
    for(int i=0; i<NX; i++) printf("-------");
    OSQPFloat log_x[NX] = {0};
    OSQPFloat log_u[NU] = {0};
    printf("|");
    for(int i=0; i<NU; i++) printf("-------");
    printf("\n");

    OSQPFloat x_next[NX];
    OSQPFloat z_warm[n_vars];//for warm start shifts
    OSQPInt steps_completed = 0;


    for (int step = 0; step < MAX_STEPS; step++) {

        timer_now(&iter_start);

    // -----------------------------------------------------
    // Phase 1: Solve & Apply (OSQP solve using current matrices)
    // -----------------------------------------------------
        timer_now(&solver_start);
        osqp_solve(solver);
        timer_now(&solver_end);
        OSQPFloat solver_dur_ms = timer_diff_ms(&solver_start, &solver_end);
        solver_total_ms += solver_dur_ms;

        if (solver->info->status_val != OSQP_SOLVED &&
            solver->info->status_val != OSQP_SOLVED_INACCURATE &&
            solver->info->status_val != OSQP_MAX_ITER_REACHED) {
            printf("Solver failed! status=%d (%s)\n", (int)solver->info->status_val, solver->info->status);
            break;
        }

        // Extract z0 and recover physical input u0 = u_target + K*(x_current - x_target) + z0
        OSQPFloat z0[NU];
        for (int i = 0; i < NU; i++) z0[i] = solver->solution->x[i];
        for (int i = 0; i < NU; i++) {
            OSQPFloat acc = u_target[i] + z0[i];
            for (int j = 0; j < NX; j++) acc += K[i][j] * (x_current[j] - x_target[j]);
            u_applied[i] = acc;
        }

        for (int i = 0; i < NX; i++) log_x[i] = x_current[i];
        for (int i = 0; i < NU; i++) log_u[i] = u_applied[i];

    // -----------------------------------------------------
    // Phase 2: Prediction vs Reality (linear prediction vs RK2 ground truth)
    // -----------------------------------------------------
        for (int i = 0; i < NX; i++) {
            OSQPFloat ax = 0.0f, bu = 0.0f;
            for (int j = 0; j < NX; j++) ax += Ad[i][j] * x_current[j];
            for (int j = 0; j < NU; j++) bu += Bd[i][j] * u_applied[j];
            x_lin_pred[i] = ax + bu + d_lin[i];
        }

    // Heun (RK2) plant simulation using crane_dynamics(); replace if plant changes
        OSQPFloat k1[NX], k2[NX], x_temp[NX];
        crane_dynamics(x_current, u_applied, k1);
        for (int i = 0; i < NX; i++) x_temp[i] = x_current[i] + DT * k1[i];
        crane_dynamics(x_temp, u_applied, k2);
        for (int i = 0; i < NX; i++) x_next[i] = x_current[i] + 0.5f * DT * (k1[i] + k2[i]);

        for (int i = 0; i < NX; i++) x_current[i] = x_next[i];

    // -----------------------------------------------------
    // Phase 3: CMoN Gatekeeper
    // kappa = ||x_actual - x_lin_pred|| / (||x_lin_pred - x_prev|| + eps)
    // If kappa > eta_pri -> rebuild Jacobians & Hessian, else reuse cached matrices
    // -----------------------------------------------------
        OSQPFloat diff_num[NX];
        OSQPFloat diff_den[NX];
        for (int i = 0; i < NX; i++) {
            diff_num[i] = x_current[i] - x_lin_pred[i];
            diff_den[i] = x_lin_pred[i] - log_x[i];
        }
        OSQPFloat num = vec_norm(diff_num, NX);
        OSQPFloat den = vec_norm(diff_den, NX);
        OSQPFloat kappa = num / (den + 1e-8f);
        OSQPInt recompute_mats = (kappa > eta_pri);

        if (recompute_mats) {
            linearization(x_current, u_applied, Ad, Bd, d_lin);
            calc_h_gradient(x_current, jacobian_H, h_val);
            compute_sc_blocks(Ad, Bd, jacobian_H, K, AK, P_blocks, M_blocks, H_blocks);
            P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, P, r_local, 1);
            A = build_A(A, P_blocks, M_blocks, H_blocks, r_local);
            osqp_update_data_mat(solver, P->x, NULL, P->p[n_vars], A->x, NULL, A->p[n_vars]);
            linearized_times++;
        }

    // -----------------------------------------------------
    // Phase 4: Universal vector update (always update q, bounds, warm start)
    // -----------------------------------------------------
        OSQPFloat d_err[NX];
        for (OSQPInt i = 0; i < NX; i++) {
            OSQPFloat ax = 0.0f, bu = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) ax += Ad[i][j] * x_target[j];
            for (OSQPInt j = 0; j < NU; j++) bu += Bd[i][j] * u_target[j];
            d_err[i] = ax + bu + d_lin[i] - x_target[i];
        }

        for (OSQPInt i = 0; i < NX; i++) x_free[0][i] = x_current[i];
        for (OSQPInt k = 0; k < N; k++) {
            for (OSQPInt i = 0; i < NX; i++) {
                OSQPFloat acc = 0.0f;
                for (OSQPInt j = 0; j < NX; j++) acc += AK[i][j] * x_free[k][j];

                OSQPFloat shift_term = 0.0f;
                for (OSQPInt c = 0; c < NU; c++) {
                    OSQPFloat K_xtarg = 0.0f;
                    for (OSQPInt j = 0; j < NX; j++) K_xtarg += K[c][j] * x_target[j];
                    shift_term += Bd[i][c] * (u_target[c] - K_xtarg);
                }

                x_free[k + 1][i] = acc + d_lin[i] + shift_term;
            }
        }
        for (OSQPInt k = 0; k <= N; k++) {
            for (OSQPInt i = 0; i < NX; i++) e_free[k][i] = x_free[k][i] - x_target[i];
        }

        P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, P, r_local, recompute_mats);
        build_bounds(&l, &u, x_current, e_free, x_target, u_target, K, jacobian_H, h_val, x_min, x_max, u_min, u_max, r_local);
        osqp_update_data_vec(solver, q, l, u);

        // Warm start: shift z* by NU, pad zeros
        for (int idx = 0; idx < n_vars - NU; idx++) z_warm[idx] = solver->solution->x[idx + NU];
        for (int idx = n_vars - NU; idx < n_vars; idx++) z_warm[idx] = 0.0f;
        osqp_warm_start(solver, z_warm, NULL);

        timer_now(&iter_end);
        OSQPFloat iter_dur_ms = timer_diff_ms(&iter_start, &iter_end);
        iter_total_ms += iter_dur_ms;
        steps_completed++;

        fprintf(f_out, "%d", step);
        for (int i = 0; i < NX; i++) fprintf(f_out, ",%f", log_x[i]);
        for (int i = 0; i < NU; i++) fprintf(f_out, ",%f", log_u[i]);
        fprintf(f_out, "\n");

        printf("%d|", step);
        for(int i=0; i<NX; i++) {
            printf(" %6.5f", log_x[i]);
        }
        printf("|");
        for(int i=0; i<NU; i++) {
            printf(" %6.5f", log_u[i]);
        }
        printf("\n");

    }

    if (steps_completed > 0) {
        solver_avg_ms = solver_total_ms / steps_completed;
        iter_avg_ms   = iter_total_ms   / steps_completed;
    }

    printf("Solver time   -> total: %.4f ms, avg: %.4f ms\n", solver_total_ms, solver_avg_ms);
    printf("Iteration time-> total: %.4f ms, avg: %.4f ms\n", iter_total_ms, iter_avg_ms);
    printf("Linearizations performed: %lld out of %lld steps\n", linearized_times, steps_completed);


    fclose(f_out);


    osqp_cleanup(solver);
    free(P->x); free(P->i); free(P->p);
    free(A->x); free(A->i); free(A->p);
    free(q); free(l); free(u);
    free(settings);

    return 0;
}