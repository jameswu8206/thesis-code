#include "inverted_pendulum_formulation_helpers.h"
#include "mpc_pipeline_helpers.h"

#include <stdlib.h>

extern OSQPInt n_vars;
extern OSQPInt n_cons;

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
    (void)A_pow;
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
    (void)A_pow;
}

#elif FORMULATION == OPT_NON_CONDENSED

OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){

    OSQPInt P_nnz = n_vars;
    OSQPFloat *P_x = malloc(P_nnz * sizeof(OSQPFloat));
    OSQPInt *P_i = malloc(P_nnz * sizeof(OSQPInt));
    OSQPInt *P_p = malloc((n_vars + 1) * sizeof(OSQPInt));

    CscBuilder P_build = {P_i, P_p, P_x, 0, 0};
    P_p[0] = 0;
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < NX; i++) {
            csc_set(&P_build, P_build.n_cols, P_build.n_cols, Q_diag[i]);
            csc_col_done(&P_build);
        }
        for (int i = 0; i < NU; i++) {
            csc_set(&P_build, P_build.n_cols, P_build.n_cols, R_diag[i]);
            csc_col_done(&P_build);
        }
    }
    for (int i = 0; i < NX; i++) {
        csc_set(&P_build, P_build.n_cols, P_build.n_cols, P_diag[i]);
        csc_col_done(&P_build);
    }

    OSQPCscMatrix* P = malloc(sizeof(OSQPCscMatrix));
    P->m = n_vars; P->n = n_vars; P->nz = -1; P->nzmax = P_nnz;
    P->x = P_x; P->i = P_i; P->p = P_p;
    return P;
}

OSQPFloat* build_q(const OSQPFloat* x_target, const OSQPFloat* u_target, const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){

    OSQPFloat *q = calloc(n_vars, sizeof(OSQPFloat));
    int idx = 0;
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < NX; i++) {
            q[idx] = -Q_diag[i] * x_target[i];
            idx++;
        }

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
    OSQPInt setup=(A==NULL);
    CscBuilder A_build = {0};

    if(setup){
        OSQPInt A_nnz_est = NX + N * (NX * NX + NX * NU + NX) + n_vars+ (N * NH * NX);

        A = malloc(sizeof(OSQPCscMatrix));

        A->x = malloc(A_nnz_est * sizeof(OSQPFloat));
        A->i = malloc(A_nnz_est * sizeof(OSQPInt));
        A->p = malloc((n_vars + 1) * sizeof(OSQPInt));

        A_build.x = A->x;
        A_build.i = A->i;
        A_build.p = A->p;
        A_build.nz = 0;
        A_build.n_cols = 0;
        A_build.p[0] = 0;


    }



    OSQPInt var_idx = 0;
    OSQPInt val_idx = 0;


    for (int k = 0; k <= N; k++) {

        for (int i = 0; i < NX; i++) {

            if (setup) csc_set(&A_build, k*NX + i, var_idx, -1.0);
            else  A->x[val_idx] = -1.0;
            val_idx++;

            if (k < N) {
                for (int r = 0; r < NX; r++) {

                    if (setup) csc_set(&A_build, (k+1)*NX + r, var_idx, Ad[r][i]);
                        else  A->x[val_idx] = Ad[r][i];
                        val_idx++;

                }
            }

            if (k > 0) {
                for (int ih = 0; ih < NH; ih++) {
                    OSQPInt obs_row = (N + 1) * NX + (k - 1) * NH + ih;

                    if (setup) {
                        csc_set(&A_build, obs_row, var_idx, 0.0);
                    } else {
                        A->x[val_idx] = jacobian_H[ih][i];
                    }
                    val_idx++;
                }
            }

            OSQPInt bound_row =(N + 1) * NX + (N * NH) + var_idx;
            if (setup) csc_set(&A_build, bound_row, var_idx, 1.0);
            else A->x[val_idx] = 1.0;
            val_idx++;

            if (setup) csc_col_done(&A_build);
            var_idx++;
        }

        if (k < N) {
            for (int i = 0; i < NU; i++) {
                for (int r = 0; r < NX; r++) {

                    if (setup) csc_set(&A_build, (k+1)*NX + r, var_idx, Bd[r][i]);
                    else  A->x[val_idx] = Bd[r][i];
                    val_idx++;
                }
                if (setup) csc_set(&A_build, (N+1)*NX + var_idx, var_idx, 1.0);
                else  A->x[val_idx] = 1.0;
                val_idx++;

                if (setup) csc_col_done(&A_build);
                var_idx++;
            }
        }
    }

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

    for (int i = 0; i < NX; i++) {
        l[i] = -x_current[i];
        u[i] = -x_current[i];
    }
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < NX; i++) {
            l[NX + k * NX + i] = -d_lin[i] - Edist[i];
            u[NX + k * NX + i] = -d_lin[i] - Edist[i];
        }
    }
    OSQPInt safety_start = (N + 1) * NX;
    for (int k = 1; k <= N; k++) {
        for (int ih = 0; ih < NH; ih++) {
            OSQPInt row = safety_start + (k - 1) * NH + ih;

            l[row] = -OSQP_INFTY;

            OSQPFloat grad_dot_x = 0.0;
            for (int i = 0; i < NX; i++) {
                grad_dot_x += jacobian_H[ih][i] * x_current[i];
            }
            u[row] = grad_dot_x - h_vals[ih];

        }
    }

    if(setup){

        int b_idx = (N+1)*NX+ (N * NH);
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < NX; i++) {
                l[b_idx] = x_min[i]; u[b_idx] = x_max[i]; b_idx++;
            }
            for (int i = 0; i < NU; i++) {
                l[b_idx] = u_min[i]; u[b_idx] = u_max[i]; b_idx++;
            }
        }
        for (int i = 0; i < NX; i++) {
            l[b_idx] = x_min[i]; u[b_idx] = x_max[i]; b_idx++;
        }

    }

}

#elif FORMULATION == OPT_SPARSE_CONDENSED

void compute_sc_blocks(OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX], OSQPFloat K[NU][NX], OSQPFloat AK[NX][NX], OSQPFloat P_blocks[], OSQPFloat M_blocks[], OSQPFloat H_blocks[]) {
    
    // 1. Automatically extract dt from the Ad matrix
    OSQPFloat dt = Ad[0][1];
    if (dt < 1e-6) dt = 1e-6; // Safety bounds to avoid division by zero
    
    // 2. Extract the dynamically mapped B-matrix coefficient
    OSQPFloat b2 = Bd[1][0];
    if (b2 > -1e-6 && b2 <= 0.0) b2 = -1e-6;
    if (b2 < 1e-6 && b2 > 0.0) b2 = 1e-6;

    // 3. Compute exact deadbeat gains directly from the discrete-time matrices
    // This mathematically forces Trace = 0 and Determinant = 0 for (Ad + Bd*K)
    OSQPFloat k0 = (-1.0 / dt - Ad[1][0]) / b2;
    OSQPFloat k1 = (-1.0 - Ad[1][1]) / b2;

    K[0][0] = k0;
    K[0][1] = k1;

    // 4. Proceed with standard Sparse Condensed matrix building
    for (OSQPInt i = 0; i < NX; i++) {
        for (OSQPInt j = 0; j < NX; j++) {
            OSQPFloat acc = 0.0;
            for (OSQPInt c = 0; c < NU; c++) acc += Bd[i][c] * K[c][j];
            AK[i][j] = Ad[i][j] + acc;
        }
    }
    
    OSQPFloat AK_powers[R_BAND + 1][NX * NX];
    mat_eye(AK_powers[0], NX, NX);
    for (OSQPInt i = 1; i <= R_BAND; i++) mat_mul(AK_powers[i - 1], (OSQPFloat *)AK, AK_powers[i], NX, NX, NX);
    
    OSQPFloat Bd_flat[NX * NU];
    for (OSQPInt i = 0; i < NX; i++) for (OSQPInt j = 0; j < NU; j++) Bd_flat[i * NU + j] = Bd[i][j];
    
    for (OSQPInt i = 0; i < R_BAND; i++) mat_mul(AK_powers[i + 1], Bd_flat, &P_blocks[i * NX * NU], NX, NX, NU);
    mat_set(M_blocks, 0.0, (R_BAND + 1) * NU, NU);
    for (OSQPInt i = 0; i < NU; i++) M_blocks[i * NU + i] = 1.0;
    
    for (OSQPInt i = 1; i <= R_BAND; i++) {
        OSQPFloat KAK[NU * NX];
        mat_mul((OSQPFloat *)K, AK_powers[i - 1], KAK, NU, NX, NX);
        mat_mul(KAK, Bd_flat, &M_blocks[i * NU * NU], NU, NX, NU);
    }
    
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
    (void)x_current;
    (void)jacobian_H;
    (void)h_vals;
    (void)r_local;
}

#endif
