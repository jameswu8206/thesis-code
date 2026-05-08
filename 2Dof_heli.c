#include "osqp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_timing.h"

/* =========================================================
 * FORMULATION SELECTOR
 * ---------------------------------------------------------
 * Change this to switch the mathematical formulation of the QP.
 * 0 : STANDARD CONDENSED (Dense, states eliminated)
 * 1 : NON-CONDENSED (Sparse, states and inputs as variables)
 * 2 : SPARSE CONDENSED (SC, virtual inputs, LQR banded)
 * ========================================================= */
#define OPT_CONDENSED         0
#define OPT_NON_CONDENSED     1
#define OPT_SPARSE_CONDENSED  2

#define FORMULATION 2  // <--- SET YOUR METHOD HERE

/* ---------------------------------------------------------
 * Problem Constants & Parameters
 * --------------------------------------------------------- */
#define NX 4    // 4 States: Pitch, Yaw, Pitch_dot, Yaw_dot
#define NU 2    // 2 Inputs: Pitch voltage, Yaw voltage
#define N 40    // Adjust based on your benchmark scope (e.g., 10, 20, 40)
#define ND 1    
#define NH 0    // No nonlinear constraints needed; handled natively via box bounds
#define DT 0.0250 // 50ms sampling time
#define Tsim 5.0
#define R_BAND 4

/* ---------------------------------------------------------
 * Formulation-Specific Dimensions
 * --------------------------------------------------------- */
#if FORMULATION == OPT_CONDENSED
    #define N_VARS (N * NU)
    #define N_CONS (N * NU + N * NX + N * NH)
#elif FORMULATION == OPT_NON_CONDENSED
    #define N_VARS ((N + 1) * NX + N * NU)
    #define N_CONS ((N + 1) * NX + (N * NH) + ((N + 1) * NX + N * NU))
#elif FORMULATION == OPT_SPARSE_CONDENSED
    #define N_VARS (N * NU)
    #define N_CONS ((N * NX) + (N * NU) + (N * NH))
#endif

// Expose dimensions as variables for OSQP setup
OSQPInt n_vars = N_VARS;
OSQPInt n_cons = N_CONS;


/* =========================================================
 * UNIVERSAL HELPER FUNCTIONS (Shared across all methods)
 * ========================================================= */
typedef struct {
    OSQPInt *i;
    OSQPInt *p;
    OSQPFloat *x;
    OSQPInt nz;
    OSQPInt n_cols;
} CscBuilder;

void csc_set(CscBuilder *m, OSQPInt row, OSQPInt col, OSQPFloat val) {
    m->x[m->nz] = val;
    m->i[m->nz] = row;
    m->nz++;
}

void csc_col_done(CscBuilder *m) {
    m->p[m->n_cols + 1] = m->nz;
    m->n_cols++;
}

static OSQPFloat vec_norm(const OSQPFloat *v, OSQPInt len) {
    OSQPFloat acc = 0.0f;
    for (OSQPInt i = 0; i < len; i++) acc += v[i] * v[i];
    return sqrtf(acc);
}

// Basic Matrix Math Helpers
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
static void mat_set(OSQPFloat *A, OSQPFloat v, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] = v;
}



// --- Helicopter Physical Constants (Replace with your exact lab parameters) ---
const OSQPFloat m   = 1.3872;  // Total mass
const OSQPFloat l   = 0.186;   // Center of mass length
const OSQPFloat g   = 9.81;    // Gravity
const OSQPFloat Jp  = 0.0384;  // Jeq_p: Pitch equivalent inertia
const OSQPFloat Jy  = 0.0432;  // Jeq_y: Yaw equivalent inertia
const OSQPFloat Kpp = 0.204;   // Torque constant pitch
const OSQPFloat Kpy = 0.0068;  // Cross-torque pitch
const OSQPFloat Kyp = 0.0219;  // Cross-torque yaw
const OSQPFloat Kyy = 0.072;   // Torque constant yaw
const OSQPFloat Bp  = 0.800;   // Damping pitch
const OSQPFloat By  = 0.318;   // Damping yaw
static void dynamics(const OSQPFloat x[NX], const OSQPFloat u[NU], OSQPFloat f_out[NX]) {
    // x[0] = Pitch (theta), x[1] = Yaw (psi)
    // x[2] = Pitch_dot,     x[3] = Yaw_dot
    
    
    OSQPFloat sx = sin(x[0]);
    OSQPFloat cx = cos(x[0]);

    OSQPFloat Dp = Jp + m * l * l;
    OSQPFloat Dy = Jy + m * l * l * (cx * cx);

    f_out[0] = x[2];
    f_out[1] = x[3];
    f_out[2] = (Kpp * u[0] + Kpy * u[1] - Bp * x[2] - m * g * l * cx - m * l * l * sx * cx * x[3] * x[3]) / Dp;
    f_out[3] = (Kyp * u[0] + Kyy * u[1] - By * x[3] + 2.0 * m * l * l * sx * cx * x[2] * x[3]) / Dy;
}

void calc_h_gradient(const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_val[NH]) {
    // NH is 0. Left empty as constraints are handled via box bounds.
}

void linearization(OSQPFloat x_current[NX], OSQPFloat u_applied[NU], OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat d_lin[NX]) {
    
    OSQPFloat Ac[NX][NX] = {0};
    OSQPFloat Bc[NX][NU] = {0};

    OSQPFloat sx = sin(x_current[0]);
    OSQPFloat cx = cos(x_current[0]);
    OSQPFloat s2x = 2.0 * sx * cx; 
    OSQPFloat c2x = cx * cx - sx * sx;

    OSQPFloat Dp = Jp + m * l * l;
    OSQPFloat Dy = Jy + m * l * l * (cx * cx);
    
    // Numerator for f4 (Yaw dot dynamics)
    OSQPFloat Ny = Kyp * u_applied[0] + Kyy * u_applied[1] - By * x_current[3] + 2.0 * m * l * l * sx * cx * x_current[2] * x_current[3];

    // Derivatives of Dy and Ny w.r.t x[0]
    OSQPFloat dDy_dx0 = -2.0 * m * l * l * sx * cx;
    OSQPFloat dNy_dx0 = 2.0 * m * l * l * c2x * x_current[2] * x_current[3];

    // --- Compute Continuous Jacobian Ac = df/dx ---
    Ac[0][2] = 1.0;
    Ac[1][3] = 1.0;
    
    // df3/dx
    Ac[2][0] = (m * g * l * sx - m * l * l * c2x * x_current[3] * x_current[3]) / Dp;
    Ac[2][2] = -Bp / Dp;
    Ac[2][3] = (-2.0 * m * l * l * sx * cx * x_current[3]) / Dp;

    // df4/dx (Quotient Rule)
    Ac[3][0] = (dNy_dx0 * Dy - Ny * dDy_dx0) / (Dy * Dy);
    Ac[3][2] = (2.0 * m * l * l * sx * cx * x_current[3]) / Dy;
    Ac[3][3] = (-By + 2.0 * m * l * l * sx * cx * x_current[2]) / Dy;

    // --- Compute Continuous Input Jacobian Bc = df/du ---
    Bc[2][0] = Kpp / Dp;
    Bc[2][1] = Kpy / Dp;
    Bc[3][0] = Kyp / Dy;
    Bc[3][1] = Kyy / Dy;

    // --- Discretization (Forward Euler) & Affine Term ---
    OSQPFloat fc[NX];
    dynamics(x_current, u_applied, fc);

    for (int i = 0; i < NX; i++) {
        OSQPFloat Ac_x = 0.0, Bc_u = 0.0;
        for (int j = 0; j < NX; j++) {
            Ad[i][j] = (i == j ? 1.0 : 0.0) + Ac[i][j] * DT;
            Ac_x += Ac[i][j] * x_current[j];
        }
        for (int j = 0; j < NU; j++) {
            Bd[i][j] = Bc[i][j] * DT;
            Bc_u += Bc[i][j] * u_applied[j];
        }
        // Compute linearization offset error (d_lin)
        d_lin[i] = DT * (fc[i] - Ac_x - Bc_u);
    }
}

/* =========================================================
 * FORMULATION 0: STANDARD CONDENSED
 * ========================================================= */
#if FORMULATION == OPT_CONDENSED

void build_AB(const OSQPFloat Ad[NX][NX], const OSQPFloat Bd[NX][NU], OSQPFloat A_pow[N + 1][NX][NX], OSQPFloat AB_mat[N][NX][NU]) {
    for (int r = 0; r < NX; r++) for (int c = 0; c < NX; c++) A_pow[0][r][c] = (r == c) ? 1.0 : 0.0;
    for (int k = 1; k <= N; k++) {
        for (int r = 0; r < NX; r++) {
            for (int c = 0; c < NX; c++) {
                A_pow[k][r][c] = 0.0;
                for (int m = 0; m < NX; m++) A_pow[k][r][c] += Ad[r][m] * A_pow[k - 1][m][c];
            }
        }
        for (int r = 0; r < NX; r++) {
            for (int c = 0; c < NU; c++) {
                AB_mat[k - 1][r][c] = 0.0f;
                for (int m = 0; m < NX; m++) AB_mat[k - 1][r][c] += A_pow[k - 1][r][m] * Bd[m][c];
            }
        }
    }
}

OSQPCscMatrix* build_P(OSQPCscMatrix* P, const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){
    OSQPInt setup = (P == NULL);
    CscBuilder P_build = {0};
    if (setup) {
        OSQPInt P_nnz = (N_VARS + 1) * N_VARS / 2;
        P = malloc(sizeof(OSQPCscMatrix));
        P->x = malloc(P_nnz * sizeof(OSQPFloat));
        P->i = malloc(P_nnz * sizeof(OSQPInt));
        P->p = malloc((N_VARS + 1) * sizeof(OSQPInt));
        P->nzmax = P_nnz;
    }
    P_build.x = P->x; P_build.i = P->i; P_build.p = P->p; P_build.nz = 0; P_build.n_cols = 0; P_build.p[0] = 0;

    for (int j = 0; j < N_VARS; j++) { 
        OSQPInt step_j = j / NU; OSQPInt ctrl_j = j % NU;           
        for (int i = 0; i <= j; i++) { 
            OSQPInt step_i = i / NU; OSQPInt ctrl_i = i % NU;
            OSQPFloat val = (i == j) ? R_diag[ctrl_i] : 0.0;
            OSQPInt start_k = (step_i > step_j ? step_i : step_j) + 1;
            for (int k = start_k; k <= N; k++) {
                const OSQPFloat* Q_k = (k == N) ? P_diag : Q_diag;
                for (int s = 0; s < NX; s++) {
                    val += AB_mat[k - step_i - 1][s][ctrl_i] * Q_k[s] * AB_mat[k - step_j - 1][s][ctrl_j];
                }
            }
            val *= 2.0;
            if (setup) csc_set(&P_build, i, j, val);
            else { P_build.x[P_build.nz] = val; P_build.nz++; }
        }
        if (setup) csc_col_done(&P_build); else P_build.n_cols++;
    }
    if (setup) { P->m = N_VARS; P->n = N_VARS; P->nz = -1; }
    return P;
}

OSQPFloat* build_q(OSQPFloat *q, const OSQPFloat Ad[NX][NX], const OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat x_curr[NX], const OSQPFloat x_target[NX], const OSQPFloat A_pow[N+1][NX][NX], const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat Q_diag[NX], const OSQPFloat P_diag[NX]){
    OSQPInt setup = (q == NULL);
    if (setup) q = calloc(N_VARS, sizeof(OSQPFloat));
    else for (int i = 0; i < N_VARS; i++) q[i] = 0.0;

    OSQPFloat x_free[N + 1][NX];
    for (int r = 0; r < NX; r++) x_free[0][r] = x_curr[r];
    for (int k = 1; k <= N; k++) {
        for (int r = 0; r < NX; r++) {
            OSQPFloat ax = 0.0f;
            for (int c = 0; c < NX; c++) ax += Ad[r][c] * x_free[k - 1][c];
            x_free[k][r] = ax + d_lin[r] + Edist[r];
        }
    }

    for (int j = 0; j < N_VARS; j++) {
        OSQPInt step_j = j / NU; OSQPInt ctrl_j = j % NU;           
        for (int k = step_j + 1; k <= N; k++) {
            const OSQPFloat* Q_k = (k == N) ? P_diag : Q_diag;
            OSQPFloat error_k[NX];
            for (int r = 0; r < NX; r++) error_k[r] = x_free[k][r] - x_target[r];
            OSQPInt p_idx = k - step_j - 1; 
            for (int r = 0; r < NX; r++) q[j] += 2.0 * error_k[r] * Q_k[r] * AB_mat[p_idx][r][ctrl_j];
        }
    }
    return q;
}

OSQPCscMatrix* build_A(OSQPCscMatrix* A, const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat jacobian_H[NH][NX]){
    OSQPInt setup=(A==NULL);
    CscBuilder A_build = {0};
    if(setup){
        OSQPInt A_nnz = (N * (N + 1) / 2) * (NX * NU) + N_VARS + (NH * NU * N * (N + 1) / 2);
        A = malloc(sizeof(OSQPCscMatrix));
        A->x = malloc(A_nnz * sizeof(OSQPFloat));
        A->i = malloc(A_nnz * sizeof(OSQPInt));
        A->p = malloc((N_VARS + 1) * sizeof(OSQPInt));
        A->nzmax = A_nnz; 
    }
    A_build.x = A->x; A_build.i = A->i; A_build.p = A->p; A_build.nz = 0; A_build.n_cols = 0; A_build.p[0] = 0;

    for (int j = 0; j < N_VARS; j++) { 
        int step_j = j / NU; int ctrl_j = j % NU;           
        csc_set(&A_build, j, j, 1.0);
        for (int k = step_j + 1; k <= N; k++) {
            int p_idx = k - step_j - 1;
            for (int r = 0; r < NX; r++) {
                OSQPInt row_idx = N_VARS + (k - 1) * NX + r;
                if(setup) csc_set(&A_build, row_idx, j, AB_mat[p_idx][r][ctrl_j]);
                else { A->x[A_build.nz] = AB_mat[p_idx][r][ctrl_j]; A_build.nz++;}
            }

            // Nonlinear inequality rows: grad(x_k) * delta_u
            OSQPInt h_row_offset = N_VARS + N * NX + (k - 1) * NH;
            for (int ih = 0; ih < NH; ih++) {
                OSQPFloat coeff = 0.0f;
                for (int s = 0; s < NX; s++) coeff += jacobian_H[ih][s] * AB_mat[p_idx][s][ctrl_j];
                if (setup) csc_set(&A_build, h_row_offset + ih, j, coeff);
                else { A->x[A_build.nz] = coeff; A_build.nz++; }
            }
        }
        if(setup) csc_col_done(&A_build);
    }
    if(setup){ A->m = N_CONS; A->n = N_VARS; A->nz = -1;}
    return A;
}

void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat Ad[NX][NX], const OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat jacobian_H[NH][NX], const OSQPFloat h_val[NH], const OSQPFloat x_current[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPFloat A_pow[N+1][NX][NX]){
    OSQPInt setup=(*l_ptr==NULL);
    if(setup){ *l_ptr = malloc(N_CONS * sizeof(OSQPFloat)); *u_ptr = malloc(N_CONS * sizeof(OSQPFloat)); }
    OSQPFloat *l = *l_ptr, *u = *u_ptr;

    for (int i = 0; i < N * NU; i++) { l[i] = u_min[i % NU]; u[i] = u_max[i % NU]; }

    OSQPFloat x_free[N + 1][NX];
    for (int r = 0; r < NX; r++) x_free[0][r] = x_current[r];
    for (int k = 1; k <= N; k++) {
        for (int r = 0; r < NX; r++) {
            OSQPFloat ax = 0.0f;
            for (int c = 0; c < NX; c++) ax += Ad[r][c] * x_free[k - 1][c];
            x_free[k][r] = ax + d_lin[r] + Edist[r];
        }
    }

    int b_idx = N * NU; 
    for (int k = 1; k <= N; k++) {
        for (int i = 0; i < NX; i++) {
            l[b_idx] = x_min[i] - x_free[k][i]; u[b_idx] = x_max[i] - x_free[k][i]; b_idx++;
        }
    }

    for (int k = 1; k <= N; k++) {
        OSQPFloat grad_xfree[NH];
        for (int ih = 0; ih < NH; ih++) {
            grad_xfree[ih] = 0.0f;
            for (int s = 0; s < NX; s++) grad_xfree[ih] += jacobian_H[ih][s] * x_free[k][s];
        }
        OSQPFloat grad_xcurr[NH];
        for (int ih = 0; ih < NH; ih++) {
            grad_xcurr[ih] = 0.0f;
            for (int s = 0; s < NX; s++) grad_xcurr[ih] += jacobian_H[ih][s] * x_current[s];
        }
        for (int ih = 0; ih < NH; ih++) {
            l[b_idx] = -OSQP_INFTY;
            u[b_idx] = -h_val[ih] + grad_xcurr[ih] - grad_xfree[ih];
            b_idx++;
        }
    }
}


/* =========================================================
 * FORMULATION 1: NON-CONDENSED
 * ========================================================= */
#elif FORMULATION == OPT_NON_CONDENSED

OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){

    //build big P the costs
    OSQPInt P_nnz = n_vars;//each x and u has cost
    OSQPFloat *P_x = malloc(P_nnz * sizeof(OSQPFloat));
    OSQPInt *P_i = malloc(P_nnz * sizeof(OSQPInt));
    OSQPInt *P_p = malloc((n_vars + 1) * sizeof(OSQPInt));

    CscBuilder P_build = {P_i, P_p, P_x, 0, 0};
    P_p[0] = 0;
    //filling costs(big P)
    for (int k = 0; k < N; k++) {
        // x_k cost
        for (int i = 0; i < NX; i++) {
            csc_set(&P_build, P_build.n_cols, P_build.n_cols, Q_diag[i]);
            csc_col_done(&P_build);
        }
        // u_k cost
        for (int i = 0; i < NU; i++) {
            csc_set(&P_build, P_build.n_cols, P_build.n_cols, R_diag[i]);
            csc_col_done(&P_build);
        }
    }
    // Terminal x_N cost
    for (int i = 0; i < NX; i++) {
        csc_set(&P_build, P_build.n_cols, P_build.n_cols, P_diag[i]);
        csc_col_done(&P_build);
    }

    //copy filled big P
    OSQPCscMatrix* P = malloc(sizeof(OSQPCscMatrix));
    P->m = n_vars; P->n = n_vars; P->nz = -1; P->nzmax = P_nnz;
    P->x = P_x; P->i = P_i; P->p = P_p;
    return P;
}

OSQPFloat* build_q(const OSQPFloat* x_target, const OSQPFloat* u_target, const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){

    OSQPFloat *q = calloc(n_vars, sizeof(OSQPFloat));
    int idx = 0;
    for (int k = 0; k < N; k++) {
        
        // State Linear Cost: q = -Q * x_target
        for (int i = 0; i < NX; i++) {
            q[idx] = -Q_diag[i] * x_target[i]; 
            idx++;
        }

        // Control Linear Cost: q = -R * u_target
        for (int i = 0; i < NU; i++) {
            q[idx] = -R_diag[i] * u_target[i];
            idx++;
        }
    }
    for (int i = 0; i < NX; i++) {
        q[idx] = -P_diag[i] * x_target[i];
        idx++;
    }
    return q;

}

OSQPCscMatrix* build_A(OSQPCscMatrix* A,OSQPFloat Ad[NX][NX],OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX] ){
    //build big A
    // Init state: NX
    // Dynamics: N * (NX*NX for A + NX*NU for B + NX for -I) (Ad + Bd + -I)
    // Bounds: n_vars (Identity diagonal)
    OSQPInt setup=(A==NULL);
    CscBuilder A_build = {0};

    if(setup){
        //(N * NH * NX) to account for nonlinear constraints
        OSQPInt A_nnz_est = NX + N * (NX * NX + NX * NU + NX) + n_vars+ (N * NH * NX);
    
        A = malloc(sizeof(OSQPCscMatrix));
        
        // Allocate the Arrays
        A->x = malloc(A_nnz_est * sizeof(OSQPFloat));
        A->i = malloc(A_nnz_est * sizeof(OSQPInt));
        A->p = malloc((n_vars + 1) * sizeof(OSQPInt));
        
        // Initialize the Builder
        A_build.x = A->x; 
        A_build.i = A->i; 
        A_build.p = A->p;
        A_build.nz = 0;
        A_build.n_cols = 0;
        A_build.p[0] = 0;


    }
    


    OSQPInt var_idx = 0;//current column
    OSQPInt val_idx = 0;


    //filling big A
    for (int k = 0; k <= N; k++) { // Note: include xN
        
        // only fill A if k<N
        for (int i = 0; i < NX; i++) {
            
            // 1. The "Self" Term (-I)
            // Row Index: k * NX + i
            if (setup) csc_set(&A_build, k*NX + i, var_idx, -1.0);
            else  A->x[val_idx] = -1.0;
            val_idx++;

            // 2. The "Next" Term (A)
            // Row Index: (k + 1) * NX + r
            if (k < N) { 
                for (int r = 0; r < NX; r++) {
                    
                    if (setup) csc_set(&A_build, (k+1)*NX + r, var_idx, Ad[r][i]);
                        else  A->x[val_idx] = Ad[r][i];
                        val_idx++;
                        //+1 since there is I above it
                    
                }
            }
            
            // 3. UNIVERSAL NONLINEAR SECTION (States)
            // Reserves NH rows for EVERY state column at time k
            if (k > 0) {
                for (int ih = 0; ih < NH; ih++) {
                    OSQPInt obs_row = (N + 1) * NX + (k - 1) * NH + ih;
    
                    if (setup) {
                        csc_set(&A_build, obs_row, var_idx, 0.0);
                    } else {
                        // You fill the specific gradient for constraint ih
                        A->x[val_idx] = jacobian_H[ih][i]; 
                    } // Gradient dh/dx_i goes here
                    val_idx++;
                }
            }

            // 4. Box Bounds (I)
            OSQPInt bound_row =(N + 1) * NX + (N * NH) + var_idx; 
            if (setup) csc_set(&A_build, bound_row, var_idx, 1.0);
            else A->x[val_idx] = 1.0;
            val_idx++;

            if (setup) csc_col_done(&A_build);
            var_idx++;
        }

        // Control u[k] k<N
        if (k < N) {
            for (int i = 0; i < NU; i++) {
                // 1. Dynamics Input (B)
                for (int r = 0; r < NX; r++) {
                    
                    if (setup) csc_set(&A_build, (k+1)*NX + r, var_idx, Bd[r][i]);
                    else  A->x[val_idx] = Bd[r][i];
                    val_idx++;
                }
                // 2. Bounds(I)
                if (setup) csc_set(&A_build, (N+1)*NX + var_idx, var_idx, 1.0);
                else  A->x[val_idx] = 1.0;
                val_idx++;

                if (setup) csc_col_done(&A_build);
                var_idx++;
            }
        }
    }

    //copy filled big A
    if(setup){
        A->m = n_cons; A->n = n_vars; A->nz = -1; A->nzmax = A_build.nz;
    }
    
    return A;
};

void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_vals[NH],OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU]){

    OSQPInt setup=(*l_ptr==NULL);

    if(setup){
        *l_ptr = malloc(n_cons * sizeof(OSQPFloat));
        *u_ptr = malloc(n_cons * sizeof(OSQPFloat));
    }

    OSQPFloat *l = *l_ptr;
    OSQPFloat *u = *u_ptr;

    //update x0 starting point for each MPC loop
    for (int i = 0; i < NX; i++) {
        l[i] = -x_current[i];
        u[i] = -x_current[i];
    }
    //Dynamics Equality (Rows NX to (N+1)*NX - 1)
    // Ax + Bu - x_next = -Ed
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < NX; i++) {
            l[NX + k * NX + i] = -d_lin[i] - Edist[i];
            u[NX + k * NX + i] = -d_lin[i] - Edist[i];
        }
    } 
    // --- PART 2: NEW NONLINEAR SECTION (Inequality) ---
    OSQPInt safety_start = (N + 1) * NX;
    for (int k = 1; k <= N; k++) {
        for (int ih = 0; ih < NH; ih++) {
            OSQPInt row = safety_start + (k - 1) * NH + ih;
            
            // h(x) <= 0 becomes l_obs <= grad*x <= u_obs
            l[row] = -OSQP_INFTY; 
            
            // Calculate u = grad' * x_current - h(x_current)
            OSQPFloat grad_dot_x = 0.0;
            for (int i = 0; i < NX; i++) {
                grad_dot_x += jacobian_H[ih][i] * x_current[i];
            }
            u[row] = grad_dot_x - h_vals[ih];
            
        }
    }
    
    if(setup){

        //Variable Bounds lower<x,u<upper (Rows (N+1)*NX to end)
        int b_idx = (N+1)*NX+ (N * NH);
        for (int k = 0; k < N; k++) {
            // x bounds
            for (int i = 0; i < NX; i++) {
                l[b_idx] = x_min[i]; u[b_idx] = x_max[i]; b_idx++;
            }
            // u bounds
            for (int i = 0; i < NU; i++) {
                l[b_idx] = u_min[i]; u[b_idx] = u_max[i]; b_idx++;
            }
        }
        // x_N bounds
        for (int i = 0; i < NX; i++) {
            l[b_idx] = x_min[i]; u[b_idx] = x_max[i]; b_idx++;
        }

    }
    
}


/* =========================================================
 * FORMULATION 2: SPARSE CONDENSED
 * ========================================================= */
#elif FORMULATION == OPT_SPARSE_CONDENSED


void compute_sc_blocks(OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX], OSQPFloat K[NU][NX], OSQPFloat AK[NX][NX], OSQPFloat P_blocks[], OSQPFloat M_blocks[], OSQPFloat H_blocks[]) {
    
    /* 
     * AGGRESSIVE LQR GAIN FOR 2-DOF HELICOPTER (Ts = 0.05s)
     * Tuned to damp the prediction matrices without causing numeric overflow in OSQP.
     * K[0] -> Pitch Voltage
     * K[1] -> Yaw Voltage
     */
    //OSQPFloat K_lqr[NU][NX] = {{-156.38, -1.70, -12.27, 0.08},{19.49, -147.39, 1.76, -19.39}};
    OSQPFloat K_lqr[NU][NX] = {{-160, -2.0, -13.0, 0.1},{20.0, -150.5, 2.0, -20.0}};

    // 1. Assign K and compute closed-loop matrix AK = Ad + Bd*K
    for (OSQPInt i = 0; i < NU; i++) {
        for (OSQPInt j = 0; j < NX; j++) {
            K[i][j] = K_lqr[i][j];
        }
    }

    for (OSQPInt i = 0; i < NX; i++) {
        for (OSQPInt j = 0; j < NX; j++) {
            OSQPFloat acc = 0.0;
            for (OSQPInt c = 0; c < NU; c++) {
                acc += Bd[i][c] * K[c][j];
            }
            AK[i][j] = Ad[i][j] + acc;
        }
    }
    
    // 2. Compute powers of AK
    OSQPFloat AK_powers[R_BAND + 1][NX * NX];
    mat_eye(AK_powers[0], NX, NX);
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        mat_mul(AK_powers[i - 1], (OSQPFloat *)AK, AK_powers[i], NX, NX, NX);
    }
    
    // 3. Flatten Bd
    OSQPFloat Bd_flat[NX * NU];
    for (OSQPInt i = 0; i < NX; i++) {
        for (OSQPInt j = 0; j < NU; j++) {
            Bd_flat[i * NU + j] = Bd[i][j];
        }
    }
    
    // 4. Compute P_blocks
    for (OSQPInt i = 0; i < R_BAND; i++) {
        mat_mul(AK_powers[i], Bd_flat, &P_blocks[i * NX * NU], NX, NX, NU);    
    }
    
    // 5. Compute M_blocks
    mat_set(M_blocks, 0.0, (R_BAND + 1) * NU, NU);
    for (OSQPInt i = 0; i < NU; i++) {
        M_blocks[i * NU + i] = 1.0;
    }
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        OSQPFloat KAK[NU * NX];
        mat_mul((OSQPFloat *)K, AK_powers[i - 1], KAK, NU, NX, NX);
        mat_mul(KAK, Bd_flat, &M_blocks[i * NU * NU], NU, NX, NU);
    }
    
    // 6. Compute H_blocks
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        OSQPFloat CAK[NH * NX];
        mat_mul((OSQPFloat *)jacobian_H, AK_powers[i - 1], CAK, NH, NX, NX);
        mat_mul(CAK, Bd_flat, &H_blocks[(i - 1) * NH * NU], NH, NX, NU);
    }
}
OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX], const OSQPFloat K[NU][NX], const OSQPFloat P_blocks[], const OSQPFloat M_blocks[], const OSQPFloat e_free[][NX], OSQPFloat q[], OSQPCscMatrix *P_existing, OSQPInt r_local, OSQPInt recompute_hessian) {
    OSQPInt setup = (P_existing == NULL);
    OSQPInt P_nnz_est = N_VARS * (r_local * NU + NU);
    OSQPCscMatrix *P; CscBuilder P_build = {0}; OSQPInt val_idx = 0;
    if (setup) {
        P = malloc(sizeof(OSQPCscMatrix)); P->x = malloc(P_nnz_est * sizeof(OSQPFloat)); P->i = malloc(P_nnz_est * sizeof(OSQPInt)); P->p = malloc((N_VARS + 1) * sizeof(OSQPInt));
        P_build.x = P->x; P_build.i = P->i; P_build.p = P->p; P_build.nz = 0; P_build.n_cols = 0; P_build.p[0] = 0;
    } else { P = P_existing; }

    for (OSQPInt i = 0; i < N_VARS; i++) q[i] = 0.0f;
    static OSQPFloat H_dense[N * NU][N * NU];
    if (recompute_hessian) for (OSQPInt i = 0; i < N_VARS; i++) for (OSQPInt j = 0; j < N_VARS; j++) H_dense[i][j] = 0.0f;

    OSQPFloat v_free[N][NU];
    for (OSQPInt k = 0; k < N; k++) {
        for (OSQPInt i = 0; i < NU; i++) {
            OSQPFloat acc = 0.0f; for (OSQPInt j = 0; j < NX; j++) acc += K[i][j] * e_free[k][j]; v_free[k][i] = acc;
        }
    }

    for (OSQPInt k1 = 0; k1 < N; k1++) {
        for (OSQPInt j1 = 0; j1 < NU; j1++) {
            OSQPInt col1 = k1 * NU + j1;
            for (OSQPInt t = k1; t <= N && (t - k1) <= r_local; t++) {
                OSQPInt idx1 = t - k1; if (idx1 >= r_local) break;
                const OSQPFloat *weight = (t == N) ? P_diag : Q_diag;
                for (OSQPInt s = 0; s < NX; s++) q[col1] += P_blocks[idx1 * NX * NU + s * NU + j1] * weight[s] * e_free[t][s];
                if (t < N) for (OSQPInt u_i = 0; u_i < NU; u_i++) q[col1] += M_blocks[idx1 * NU * NU + u_i * NU + j1] * R_diag[u_i] * v_free[t][u_i];
            }
            if (recompute_hessian) {
                for (OSQPInt k2 = k1; k2 < N && (k2 - k1) < r_local; k2++) {
                    for (OSQPInt j2 = 0; j2 < NU; j2++) {
                        OSQPInt col2 = k2 * NU + j2; if (col2 < col1) continue;
                        OSQPFloat h_val = 0.0f;
                        for (OSQPInt t = k2; t <= N; t++) {
                            OSQPInt idx1 = t - k1, idx2 = t - k2; if (idx1 >= r_local || idx2 >= r_local) break;
                            const OSQPFloat *w = (t == N) ? P_diag : Q_diag;
                            for (OSQPInt s = 0; s < NX; s++) h_val += P_blocks[idx1 * NX * NU + s * NU + j1] * w[s] * P_blocks[idx2 * NX * NU + s * NU + j2];
                            if (t < N) for (OSQPInt u_i = 0; u_i < NU; u_i++) h_val += M_blocks[idx1 * NU * NU + u_i * NU + j1] * R_diag[u_i] * M_blocks[idx2 * NU * NU + u_i * NU + j2];
                        }
                        H_dense[col1][col2] += h_val; if (col1 != col2) H_dense[col2][col1] += h_val;
                    }
                }
            }
        }
    }
    if (recompute_hessian) {
        if (setup) {
            for (OSQPInt col = 0; col < N_VARS; col++) {
                OSQPInt row_start = (col >= r_local * NU) ? (col - r_local * NU) : 0;
                for (OSQPInt row = row_start; row <= col; row++) csc_set(&P_build, row, col, H_dense[row][col]);
                csc_col_done(&P_build);
            }
            P->m = N_VARS; P->n = N_VARS; P->nz = -1; P->nzmax = P_build.nz;
        } else {
            for (OSQPInt col = 0; col < N_VARS; col++) for (OSQPInt idx = P->p[col]; idx < P->p[col + 1]; idx++) P->x[val_idx++] = H_dense[P->i[idx]][col];
        }
    }
    return P;
}

OSQPCscMatrix* build_A(OSQPCscMatrix* A, const OSQPFloat P_blocks[], const OSQPFloat M_blocks[], const OSQPFloat H_blocks[], OSQPInt r_local) {
    OSQPInt setup = (A == NULL); CscBuilder A_build = {0};
    OSQPInt A_nnz_est = N * NU * (r_local * NX + (r_local + 1) * NU + r_local * NH);
    if (setup) {
        A = malloc(sizeof(OSQPCscMatrix)); A->x = malloc(A_nnz_est * sizeof(OSQPFloat)); A->i = malloc(A_nnz_est * sizeof(OSQPInt)); A->p = malloc((N_VARS + 1) * sizeof(OSQPInt));
        A_build.x = A->x; A_build.i = A->i; A_build.p = A->p; A_build.nz = 0; A_build.n_cols = 0; A_build.p[0] = 0;
    } else { A_build.x = A->x; A_build.i = A->i; A_build.p = A->p; A_build.nz = 0; A_build.n_cols = 0; }
    OSQPInt val_idx = 0;
    for (OSQPInt k = 0; k < N; k++) {
        for (OSQPInt j = 0; j < NU; j++) {
            OSQPInt current_col = k * NU + j;
            for (OSQPInt i = 0; i < r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (k + i) * NX;
                for (OSQPInt r = 0; r < NX; r++) {
                    OSQPFloat v = P_blocks[i * NX * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v); else { A->x[val_idx++] = v; }
                }
            }
            for (OSQPInt i = 0; i <= r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (N * NX) + (k + i) * NU;
                for (OSQPInt r = 0; r < NU; r++) {
                    OSQPFloat v = M_blocks[i * NU * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v); else { A->x[val_idx++] = v; }
                }
            }
            for (OSQPInt i = 1; i <= r_local && (k + i) < N; i++) {
                OSQPInt row_offset = (N * NX) + (N * NU) + (k + i) * NH;
                for (OSQPInt r = 0; r < NH; r++) {
                    OSQPFloat v = H_blocks[(i - 1) * NH * NU + r * NU + j];
                    if (setup) csc_set(&A_build, row_offset + r, current_col, v); else { A->x[val_idx++] = v; }
                }
            }
            if (setup) csc_col_done(&A_build);
        }
    }
    if (setup) { A->m = N_CONS; A->n = N_VARS; A->nz = -1; A->nzmax = A_build.nz; }
    return A;
}

void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], const OSQPFloat e_free[][NX], const OSQPFloat x_target[NX], const OSQPFloat u_target[NU], const OSQPFloat K[NU][NX], const OSQPFloat jacobian_H[NH][NX], const OSQPFloat h_vals[NH], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPInt r_local) {
    if (*l_ptr == NULL) { *l_ptr = malloc(N_CONS * sizeof(OSQPFloat)); *u_ptr = malloc(N_CONS * sizeof(OSQPFloat)); }
    OSQPFloat *l = *l_ptr, *u = *u_ptr;
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt i = 0; i < NX; i++) {
            OSQPFloat x_pred = e_free[n + 1][i] + x_target[i];
            l[n * NX + i] = x_min[i] - x_pred; u[n * NX + i] = x_max[i] - x_pred;
        }
    }
    OSQPInt u_off = N * NX;
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt i = 0; i < NU; i++) {
            OSQPFloat v_free_val = 0.0f; for (OSQPInt j = 0; j < NX; j++) v_free_val += K[i][j] * e_free[n][j];
            OSQPFloat u_pred = u_target[i] + v_free_val;
            l[u_off + n * NU + i] = u_min[i] - u_pred; u[u_off + n * NU + i] = u_max[i] - u_pred;
        }
    }
    OSQPInt h_off = N * NX + N * NU;
    for (OSQPInt n = 0; n < N; n++) {
        for (OSQPInt ih = 0; ih < NH; ih++) {
            l[h_off + n * NH + ih] = -OSQP_INFTY; u[h_off + n * NH + ih] = OSQP_INFTY;
        }
    }
}
#endif


/* =========================================================
 * MAIN NMPC PIPELINE
 * ========================================================= */
int main(void) {

    // --- Dynamic Physics Variables ---
    OSQPFloat Ad[NX][NX] = {0};
    OSQPFloat Bd[NX][NU] = {0};
    OSQPFloat jacobian_H[NH][NX] = {};
    OSQPFloat h_val[NH] = {};
    OSQPFloat d_lin[NX] = {0};
    OSQPFloat x_lin_pred[NX] = {0};
    OSQPInt linearized_times = 1;
    OSQPFloat d_est[ND] = {0.0};
    const OSQPFloat E[NX][ND] = {{0.00277189},
                                {0.55371331}};
    
    
    // --- Cost Weights ---
    const OSQPFloat Q_diag[NX] = {1000.0 * DT, 500.0 * DT, 0.0, 0.0}; 
    const OSQPFloat R_diag[NU] = {0.001 * DT, 0.01 * DT};      
    const OSQPFloat P_diag[NX] = {1000.0, 500.0, 0.0, 0.0}; 

    // --- Setpoints & Bounds ---
    OSQPFloat x_current[NX] = {-0.7, 0.0, 0.0, 0.0}; // Initial point
    OSQPFloat x_target[NX]  = { 0.0, 0.0, 0.0, 0.0}; // Target point (Hover)
    OSQPFloat u_applied[NU] = {0.0, 0.0};
    OSQPFloat u_target[NU]  = {0.0, 0.0};

    // State Bounds: Pitch constrained to +- 40 degrees (0.698 rad)
    const OSQPFloat x_min[NX] = {-1.0, -100.0, -100.0, -100.0};
    const OSQPFloat x_max[NX] = { 1.0,  100.0,  100.0,  100.0};
    
    // Input Bounds: Pitch +-24V, Yaw +-15V
    const OSQPFloat u_min[NU] = {-24.0, -15.0};
    const OSQPFloat u_max[NU] = { 24.0,  15.0};


    linearization(x_current, u_applied, Ad, Bd, d_lin);
    calc_h_gradient(x_current, jacobian_H, h_val);
    OSQPFloat Edist[NX] = {0};
    for(int i=0; i<NX; i++) {
        for(int j=0; j<ND; j++) {
            Edist[i] += E[i][j] * d_est[j];
        }
    }
    // --- Formulation Specific Pointers & Precomputation arrays ---
    OSQPCscMatrix *P = NULL;
    OSQPFloat *q = malloc(n_vars * sizeof(OSQPFloat));

    // Constraint matrix & bounds
    OSQPCscMatrix *A = NULL;
    OSQPFloat *l = NULL, *u = NULL;

    #if FORMULATION == OPT_CONDENSED
        OSQPFloat A_pow[N + 1][NX][NX];
        OSQPFloat AB_mat[N][NX][NU];
    #elif FORMULATION == OPT_SPARSE_CONDENSED
        OSQPFloat K[NU][NX] = {0};
        OSQPFloat AK[NX][NX] = {0};
        OSQPFloat P_blocks[NU * NX * R_BAND];
        OSQPFloat M_blocks[NU * NU * (R_BAND + 1)];
        OSQPFloat H_blocks[NH * NU * R_BAND];
        OSQPInt r_local = R_BAND;
        OSQPFloat x_free[N + 1][NX];
        OSQPFloat e_free[N + 1][NX];
    #endif



    const OSQPFloat eta_pri = 1000.0; // CMoN Gatekeeper Tolerance


    #if FORMULATION == OPT_CONDENSED
        build_AB(Ad, Bd, A_pow, AB_mat);
        P = build_P(NULL, AB_mat, Q_diag, R_diag, P_diag);
        q = build_q(NULL, Ad, d_lin, Edist, x_current, x_target, A_pow, AB_mat, Q_diag, P_diag);
        A = build_A(NULL, AB_mat, jacobian_H);
        build_bounds(&l, &u, Ad, d_lin, Edist, jacobian_H, h_val, x_current, x_min, x_max, u_min, u_max, A_pow);

    #elif FORMULATION == OPT_NON_CONDENSED
        P = build_P(Q_diag, R_diag, P_diag); // Constant
        q = build_q(x_target, u_target, Q_diag, R_diag, P_diag); // Constant for fixed target
        A = build_A(NULL, Ad, Bd, jacobian_H);
        build_bounds(&l, &u, x_current, jacobian_H, h_val, d_lin, Edist, x_min, x_max, u_min, u_max);
    #elif FORMULATION == OPT_SPARSE_CONDENSED
        
        compute_sc_blocks(Ad, Bd, jacobian_H, K, AK, P_blocks, M_blocks, H_blocks);
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

                x_free[k + 1][i] = acc + d_lin[i] + Edist[i] + shift_term;
            }
        }
        // Convert to error coordinates
        for (OSQPInt k = 0; k <= N; k++) {
            for (OSQPInt i = 0; i < NX; i++) e_free[k][i] = x_free[k][i] - x_target[i];
        }
        P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, NULL, r_local, 1);
        A = build_A(NULL, P_blocks, M_blocks, H_blocks, r_local);
        build_bounds(&l, &u, x_current, e_free, x_target, u_target, K, jacobian_H, h_val, x_min, x_max, u_min, u_max, r_local);
    #endif

    // Steady-state tracking residual d_err = Ad*x_t + Bd*u_t + d_lin - x_t
    OSQPFloat d_err[NX];
    for (OSQPInt i = 0; i < NX; i++) {
        OSQPFloat ax = 0.0f, bu = 0.0f;
        for (OSQPInt j = 0; j < NX; j++) ax += Ad[i][j] * x_target[j];
        for (OSQPInt j = 0; j < NU; j++) bu += Bd[i][j] * u_target[j];
        d_err[i] = ax + bu + d_lin[i] - x_target[i];
    }

    

    // --- Timing accumulators ---
    typeTime solver_start, solver_end, iter_start, iter_end;
    OSQPFloat solver_total_ms = 0.0f;
    OSQPFloat iter_total_ms   = 0.0f;
    OSQPFloat solver_avg_ms   = 0.0f;
    OSQPFloat iter_avg_ms     = 0.0f;

    // --- OSQP Setup ---
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);
    #if FORMULATION == OPT_CONDENSED
        settings->adaptive_rho = 1.0;
        settings->check_termination = 1;//exit immediately when hit accurarcy target
        settings->eps_abs = 1e-5; 
        settings->eps_rel = 1e-5;
        settings->alpha = 1.2;
        //settings->adaptive_rho_interval = 1;
        settings->warm_starting = 1;
        settings->max_iter = 50;
        settings->verbose = 0;
        settings->scaling = 1;//no scaling since stacking already condensed
        settings->polishing = 1;
    #elif FORMULATION == OPT_NON_CONDENSED
        settings->alpha = 1.2;//too loose will fail for stabilize to unstable equilibrium
        //settings->scaling= 1;
        settings->verbose = 0;
        settings->adaptive_rho = 1;
        settings->check_termination = 1;
        settings->max_iter = 50;//make sure is not reached for 100% solved
        settings->warm_starting = 1;//default for grampc
        settings->eps_abs = 1e-5;
        settings->eps_rel = 1e-5;
        settings->polishing = 1;
    #elif FORMULATION == OPT_SPARSE_CONDENSED
        settings->alpha = 1.2;      // higher relaxation for faster convergence
        settings->scaling= 1;      // stronger scaling to improve conditioning
        settings->verbose = 0;
        settings->adaptive_rho = 1;
        settings->check_termination = 1;
        settings->max_iter = 50;//make sure is not reached for 100% solved
        settings->warm_starting = 1;//default for grampc
        settings->eps_abs = 1e-5;   // slightly looser tolerances to aid feasibility
        settings->eps_rel = 1e-5;
        settings->polishing = 1;
    
    #endif

    OSQPSolver *solver;
    OSQPInt setup_status = osqp_setup(&solver, P, q, A, l, u, n_cons, n_vars, settings);
    if (setup_status != 0 || solver == NULL) {
        printf("OSQP setup failed (code %lld)\n", (long long)setup_status);
        return 1;
    }

    // --- Simulation Loop ---
    int MAX_STEPS = Tsim / DT; 
    
    OSQPFloat x_next[NX];
    OSQPFloat z_warm[n_vars];//for warm start shifts
    OSQPInt steps_completed = 0;

    // Logging setup (CSV + console header)
    FILE *f_out = fopen("../mpc_data.csv", "w");
    if (f_out == NULL) {
        printf("Error opening ../mpc_data.csv for write\n");
        return 1;
    }
    fprintf(f_out, "step");
    for (int i = 0; i < NX; i++) fprintf(f_out, ",x%d", i);
    for (int i = 0; i < NU; i++) fprintf(f_out, ",u%d", i);
    fprintf(f_out, "\n");

    printf("Starting MPC Loop (NX=%d, NU=%d, FORMULATION=%d)\n", NX, NU, FORMULATION);
    printf("Step |");
    for (int i = 0; i < NX; i++) printf("   x%d   ", i);
    printf("|");
    for (int i = 0; i < NU; i++) printf("   u%d   ", i);
    printf("\n");
    printf("-----|");
    for (int i = 0; i < NX; i++) printf("-------");
    OSQPFloat log_x[NX] = {0};
    OSQPFloat log_u[NU] = {0};
    printf("|");
    for (int i = 0; i < NU; i++) printf("-------");
    printf("\n");

    for (int step = 0; step < MAX_STEPS; step++) {

        /*if(step%400==0&&step>0){
            x_current[1] += 0.1;//disturbance at step 50, should trigger CMoN update
        }*/

        

    // -----------------------------------------------------
    // Phase 1: Solve & Apply (OSQP solve using current matrices)
    // -----------------------------------------------------
        timer_now(&iter_start);
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

        // 2. Extract control applied based on formulation
        #if FORMULATION == OPT_CONDENSED
            for (int i = 0; i < NU; i++) u_applied[i] = solver->solution->x[i];
        #elif FORMULATION == OPT_NON_CONDENSED
            for (int i = 0; i < NU; i++) u_applied[i] = solver->solution->x[NX + i];
        #elif FORMULATION == OPT_SPARSE_CONDENSED
            OSQPFloat z0[NU];
            for (int i = 0; i < NU; i++) z0[i] = solver->solution->x[i];
            for (int i = 0; i < NU; i++) {
                OSQPFloat acc = u_target[i] + solver->solution->x[i];
                for (int j = 0; j < NX; j++) acc += K[i][j] * (x_current[j] - x_target[j]);
                u_applied[i] = acc;
            }
        #endif

        // 3. RK2 Simulation vs Linear Prediction (CMoN Base)
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

        // Heun (RK2) plant simulation using dynamics(); replace if plant changes
            OSQPFloat k1[NX], k2[NX], x_temp[NX];
            dynamics(x_current, u_applied, k1);
            for (int i = 0; i < NX; i++) x_temp[i] = x_current[i] + DT * k1[i];
            dynamics(x_temp, u_applied, k2);
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

        // 5. Matrix and Vector Updates
        if (recompute_mats) {
            linearization(x_current, u_applied, Ad, Bd, d_lin);
            calc_h_gradient(x_current, jacobian_H, h_val);
            for(int i=0; i<NX; i++) { Edist[i] = 0.0f; for(int j=0; j<ND; j++) Edist[i] += E[i][j] * d_est[j]; }

            #if FORMULATION == OPT_CONDENSED
                build_AB(Ad, Bd, A_pow, AB_mat);
                build_P(P, AB_mat, Q_diag, R_diag, P_diag);
                build_A(A, AB_mat, jacobian_H);
                osqp_update_data_mat(solver, P->x, NULL, P->nzmax, A->x, NULL, A->nzmax);
            #elif FORMULATION == OPT_NON_CONDENSED
                // P is strictly diagonal Q and R, it never changes. Only A needs updating.
                build_A(A, Ad, Bd, jacobian_H);
                osqp_update_data_mat(solver, NULL, NULL, 0, A->x, NULL, A->nzmax);
            #elif FORMULATION == OPT_SPARSE_CONDENSED
                compute_sc_blocks(Ad, Bd, jacobian_H, K, AK, P_blocks, M_blocks, H_blocks);
                build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, P, r_local, 1);
                build_A(A, P_blocks, M_blocks, H_blocks, r_local);
                osqp_update_data_mat(solver, P->x, NULL, P->nzmax, A->x, NULL, A->nzmax);
            #endif
            linearized_times++;
        }

        OSQPFloat d_err[NX];
        for (OSQPInt i = 0; i < NX; i++) {
            OSQPFloat ax = 0.0f, bu = 0.0f;
            for (OSQPInt j = 0; j < NX; j++) ax += Ad[i][j] * x_target[j];
            for (OSQPInt j = 0; j < NU; j++) bu += Bd[i][j] * u_target[j];
            d_err[i] = ax + bu + d_lin[i] + Edist[i] - x_target[i];
        }

        // Always update vectors
        #if FORMULATION == OPT_CONDENSED
            build_q(q, Ad, d_lin, Edist, x_current, x_target, A_pow, AB_mat, Q_diag, P_diag);
            build_bounds(&l, &u, Ad, d_lin, Edist, jacobian_H, h_val, x_current, x_min, x_max, u_min, u_max, A_pow);
        #elif FORMULATION == OPT_NON_CONDENSED
            // q is constant (-Q * x_target). Only bounds change based on current state.
            build_bounds(&l, &u, x_current, jacobian_H, h_val, d_lin, Edist, x_min, x_max, u_min, u_max);
        #elif FORMULATION == OPT_SPARSE_CONDENSED
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

                    x_free[k + 1][i] = acc + d_lin[i] + Edist[i] + shift_term;
                }
            }
            for (OSQPInt k = 0; k <= N; k++) {
                for (OSQPInt i = 0; i < NX; i++) e_free[k][i] = x_free[k][i] - x_target[i];
            }
            P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, P, r_local, recompute_mats);
            build_bounds(&l, &u, x_current, e_free, x_target, u_target, K, jacobian_H, h_val, x_min, x_max, u_min, u_max, r_local);
        #endif
        
        osqp_update_data_vec(solver, q, l, u);

        // 6a. Log results (state before propagation, applied control)
        fprintf(f_out, "%d", step);
        for (int i = 0; i < NX; i++) fprintf(f_out, ",%f", log_x[i]);
        for (int i = 0; i < NU; i++) fprintf(f_out, ",%f", u_applied[i]);
        fprintf(f_out, "\n");
        printf("%d|", step);
        for (int i = 0; i < NX; i++) printf(" %6.5f", log_x[i]);
        printf("|");
        for (int i = 0; i < NU; i++) printf(" %6.5f", u_applied[i]);
        printf("\n");

        // 6. Warm Start
        if(recompute_mats==0) {
            #if FORMULATION == OPT_CONDENSED || FORMULATION == OPT_SPARSE_CONDENSED
                for (int idx = 0; idx < n_vars - NU; idx++) z_warm[idx] = solver->solution->x[idx + NU];
                for (int idx = n_vars - NU; idx < n_vars; idx++) z_warm[idx] = 0.0f;
                osqp_warm_start(solver, z_warm, NULL);
            #elif FORMULATION == OPT_NON_CONDENSED
                OSQPInt stride = NX + NU;
                for (int k = 0; k < N - 1; k++) {
                    for (int i = 0; i < NX; i++) z_warm[k * stride + i] = solver->solution->x[(k + 1) * stride + i];
                    for (int i = 0; i < NU; i++) z_warm[k * stride + NX + i] = solver->solution->x[(k + 1) * stride + NX + i];
                }
                for (int i = 0; i < NX; i++) z_warm[(N - 1) * stride + i] = solver->solution->x[N * stride + i];
                for (int i = 0; i < NU; i++) z_warm[(N - 1) * stride + NX + i] = solver->solution->x[(N - 1) * stride + NX + i];
                for (int i = 0; i < NX; i++) z_warm[N * stride + i] = solver->solution->x[N * stride + i];
                osqp_warm_start(solver, z_warm, NULL);
            #endif
        }

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

    // Calculate final true cost of objective function
    OSQPFloat final_cost = 0.0;
    // The final state is x_current after the loop
    // In OSQP, Q and R diagonals are pre-multiplied by DT. We must divide by DT
    // to compare fairly against GRAMPC's raw penalty scale.
    for (int i = 0; i < NX; i++) {
        OSQPFloat err_x = x_current[i] - x_target[i];
        final_cost += err_x * (Q_diag[i] / DT) * err_x + err_x * P_diag[i] * err_x;
    }
    for (int i = 0; i < NU; i++) {
        OSQPFloat err_u = u_applied[i] - u_target[i];
        final_cost += err_u * (R_diag[i] / DT) * err_u;
    }

    printf("Final End Cost (Raw)   -> %.3f\n", (double)final_cost);
    printf("Final End Cost (log10) -> %.3f\n", (final_cost > 0) ? log10((double)final_cost) : -99.0);
    printf("%.3f/%.3f/%.3f\n\n", solver_avg_ms,iter_avg_ms,(final_cost > 0) ? log10((double)final_cost) : -99.0);

    fclose(f_out);


    osqp_cleanup(solver);
    free(P->x); free(P->i); free(P->p);
    free(A->x); free(A->i); free(A->p);
    free(q); free(l); free(u);
    free(settings);

    return 0;
}
