#include "osqp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_timing.h"
#include "mpc_pipeline_config.h"
#include "mpc_pipeline_helpers.h"
#include "mpc_pipeline_formulation_helpers.h"

// Expose dimensions as variables for OSQP setup
OSQPInt n_vars = N_VARS;
OSQPInt n_cons = N_CONS;

/* =========================================================
 * MAIN NMPC PIPELINE
 * ========================================================= */
int main(void) {

    // --- Dynamic Physics Variables ---
    OSQPFloat Ad[NX][NX] = {0};
    OSQPFloat Bd[NX][NU] = {0};
    OSQPFloat jacobian_H[NH][NX] = {0};
    OSQPFloat h_val[NH] = {0};
    OSQPFloat d_lin[NX] = {0};
    OSQPFloat x_lin_pred[NX] = {0};
    OSQPInt linearized_times = 1;
    OSQPFloat d_est[ND] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    const OSQPFloat E[NX][NX] = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
    };

    // --- Cost Weights ---
    const OSQPFloat Q_diag[NX] = {1.0*DT, 2.0*DT, 2.0*DT, 1.0*DT, 1.0*DT, 4.0*DT};
    const OSQPFloat R_diag[NU] = {0.05*DT, 0.05*DT};
    const OSQPFloat P_diag[NX] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // --- Setpoints & Bounds ---
    OSQPFloat x_current[NX] = {-2.0, 0.0, 2.0, 0.0, 0.0, 0.0}; // initial point
    OSQPFloat x_target[NX]  = { 2.0, 0.0, 2.0, 0.0, 0.0, 0.0};
    if(ANGLE){
         x_current[4] = 0.25;
         x_current[5] = 0.25;
    }
    
    OSQPFloat u_applied[NU] = {0.0, 0.0};
    OSQPFloat u_target[NU]  = {0.0, 0.0};


    const OSQPFloat x_min[NX] = {-OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY, -OSQP_INFTY};
    const OSQPFloat x_max[NX] = { OSQP_INFTY,  OSQP_INFTY,  OSQP_INFTY,  OSQP_INFTY,  OSQP_INFTY,  OSQP_INFTY};
    const OSQPFloat u_min[NU] = {-2.0, -2.0};
    const OSQPFloat u_max[NU] = { 2.0,  2.0};

    // Bootstrap linearization at initial state
    linearization(x_current, u_applied, Ad, Bd, d_lin);
    calc_h_gradient(x_current, jacobian_H, h_val);
    
    OSQPFloat Edist[NX] = {0};
    for(int i = 0; i < NX; i++) {
        for(int j = 0; j < ND; j++) {
            Edist[i] += E[i][j] * d_est[j];
        }
    }

    // --- Formulation Specific Pointers & Precomputation arrays ---
    OSQPCscMatrix *P = NULL;
    OSQPFloat *q = malloc(n_vars * sizeof(OSQPFloat));
    OSQPCscMatrix *A = NULL;
    OSQPFloat *l = NULL, *u = NULL;

    #if FORMULATION == OPT_CONDENSED
        OSQPFloat A_pow[N + 1][NX][NX];
        OSQPFloat AB_mat[N][NX][NU];
        build_AB(Ad, Bd, A_pow, AB_mat);
        P = build_P(NULL, AB_mat, Q_diag, R_diag, P_diag);
        q = build_q(NULL, Ad, d_lin, Edist, x_current, x_target, A_pow, AB_mat, Q_diag, P_diag);
        A = build_A(NULL, AB_mat, jacobian_H);
        build_bounds(&l, &u, Ad, d_lin, Edist, jacobian_H, h_val, x_current, x_min, x_max, u_min, u_max, A_pow);

    #elif FORMULATION == OPT_NON_CONDENSED
        P = build_P(Q_diag, R_diag, P_diag); 
        q = build_q(x_target, u_target, Q_diag, R_diag, P_diag);
        A = build_A(NULL, Ad, Bd, jacobian_H);
        build_bounds(&l, &u, x_current, jacobian_H, h_val, d_lin, Edist, x_min, x_max, u_min, u_max);

    #elif FORMULATION == OPT_SPARSE_CONDENSED
        OSQPFloat K[NU][NX] = {0};
        OSQPFloat AK[NX][NX] = {0};
        OSQPFloat P_blocks[NU * NX * R_BAND];
        OSQPFloat M_blocks[NU * NU * (R_BAND + 1)];
        OSQPFloat H_blocks[NH * NU * R_BAND];
        OSQPInt r_local = R_BAND;
        OSQPFloat x_free[N + 1][NX];
        OSQPFloat e_free[N + 1][NX];
        
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
        for (OSQPInt k = 0; k <= N; k++) {
            for (OSQPInt i = 0; i < NX; i++) e_free[k][i] = x_free[k][i] - x_target[i];
        }
        P = build_P(Q_diag, R_diag, P_diag, K, P_blocks, M_blocks, e_free, q, NULL, r_local, 1);
        A = build_A(NULL, P_blocks, M_blocks, H_blocks, r_local);
        build_bounds(&l, &u, x_current, e_free, x_target, u_target, K, jacobian_H, h_val, x_min, x_max, u_min, u_max, r_local);
    #endif

    // --- Tuning Parameters ---
    // --- Tuning Parameters ---
    OSQPFloat eps_abs;
    OSQPFloat eps_rel;

    if((FORMULATION == OPT_NON_CONDENSED && MAXITER == 3) ||(FORMULATION == OPT_CONDENSED && N>=40)|| ANGLE){
            eps_abs = 0.05; // Update these if you intended them to be different
            eps_rel = 0.05;
    }
    else{
            eps_abs = 100; // Update these if you intended them to be different
            eps_rel = 100;
    }
    
    

    // --- Timing accumulators ---
    typeTime solver_start, solver_end, iter_start, iter_end;
    OSQPFloat solver_total_ms = 0.0f;
    OSQPFloat iter_total_ms = 0.0f;
    OSQPFloat solver_avg_ms = 0.0f;
    OSQPFloat iter_avg_ms = 0.0f;

    // --- OSQP Setup ---
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);
    
    settings->alpha = 1.6; 
    settings->adaptive_rho = 1;
    settings->warm_starting = 1;
    settings->max_iter = MAXITER;
    settings->verbose = 0;
    settings->eps_abs = 1e-3;
    settings->eps_rel = 1e-3;
    #if FORMULATION == OPT_CONDENSED
        settings->check_termination = 1;
        settings->scaling = 0;
    #elif FORMULATION == OPT_NON_CONDENSED
        settings->check_termination = 10;
        settings->scaling = 2;
    #elif FORMULATION == OPT_SPARSE_CONDENSED
        settings->check_termination = 1;
        //settings->scaling = N;
        //settings->polishing = 1;
    #endif

    OSQPSolver *solver;
    OSQPInt setup_status = osqp_setup(&solver, P, q, A, l, u, n_cons, n_vars, settings);
    if (setup_status != 0 || solver == NULL) {
        printf("OSQP setup failed (code %lld)\n", (long long)setup_status);
        return 1;
    }

    // --- Simulation Loop Setup ---
    int MAX_STEPS = Tsim / DT;
    OSQPFloat x_next[NX];
    OSQPFloat z_warm[n_vars];
    OSQPFloat *zero_x = (OSQPFloat *)calloc(n_vars, sizeof(OSQPFloat));
    OSQPFloat *zero_y = (OSQPFloat *)calloc(n_cons, sizeof(OSQPFloat));
    OSQPInt steps_completed = 0;

    // Logging setup
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
    
    OSQPFloat log_x[NX] = {0};
    OSQPFloat log_u[NU] = {0};

    // --- Main Simulation Loop ---
    for (int step = 0; step < MAX_STEPS; step++) {

        // 1) Solve using currently cached QP matrices
        timer_now(&iter_start);
        timer_now(&solver_start);
        osqp_solve(solver);
        timer_now(&solver_end);
        
        solver_total_ms += timer_diff_ms(&solver_start, &solver_end);
        OSQPInt status = solver->info->status_val;

        if (status != OSQP_SOLVED && status != OSQP_SOLVED_INACCURATE && status != OSQP_MAX_ITER_REACHED) {
            printf("Solver failed! status=%d (%s)\n", (int)status, solver->info->status);
            osqp_warm_start(solver, zero_x, zero_y);
            break;
        }

        // Extract control action
        #if FORMULATION == OPT_CONDENSED
            for (int i = 0; i < NU; i++) u_applied[i] = solver->solution->x[i];
        #elif FORMULATION == OPT_NON_CONDENSED
            for (int i = 0; i < NU; i++) u_applied[i] = solver->solution->x[NX + i];
        #elif FORMULATION == OPT_SPARSE_CONDENSED
            for (int i = 0; i < NU; i++) {
                OSQPFloat acc = u_target[i] + solver->solution->x[i];
                for (int j = 0; j < NX; j++) acc += K[i][j] * (x_current[j] - x_target[j]);
                u_applied[i] = acc;
            }
        #endif

        // Keep pre-propagation values for logging
        for (int i = 0; i < NX; i++) log_x[i] = x_current[i];
        for (int i = 0; i < NU; i++) log_u[i] = u_applied[i];

        // 2) Prediction vs plant reality
        for (int i = 0; i < NX; i++) {
            OSQPFloat ax = 0.0f, bu = 0.0f;
            for (int j = 0; j < NX; j++) ax += Ad[i][j] * x_current[j];
            for (int j = 0; j < NU; j++) bu += Bd[i][j] * u_applied[j];
            x_lin_pred[i] = ax + bu + d_lin[i];
        }

        // Plant propagation
        OSQPFloat k1[NX], k2[NX], x_temp[NX];
        crane_dynamics(x_current, u_applied, k1);
        for (int i = 0; i < NX; i++) x_temp[i] = x_current[i] + DT * k1[i];
        crane_dynamics(x_temp, u_applied, k2);
        for (int i = 0; i < NX; i++) x_next[i] = x_current[i] + 0.5f * DT * (k1[i] + k2[i]);

        for (int i = 0; i < NX; i++) x_current[i] = x_next[i];

        // 3) Trigger Logic
        OSQPFloat diff_num[NX], diff_den[NX];
        for (int i = 0; i < NX; i++) {
            diff_num[i] = x_current[i] - x_lin_pred[i];
            diff_den[i] = x_lin_pred[i] - log_x[i];
        }
        
        OSQPFloat eta_pri = eps_abs * sqrt(NX) + eps_rel * vec_norm(x_current, NX);
        OSQPFloat num = vec_norm(diff_num, NX);
        OSQPFloat den = vec_norm(diff_den, NX);
        OSQPFloat kappa = num / (den + 1e-8f);
        
        OSQPInt recompute_mats = (kappa > eta_pri);

        // 4) Matrix refresh path
        if (recompute_mats) {
            linearization(x_current, u_applied, Ad, Bd, d_lin);
            calc_h_gradient(x_current, jacobian_H, h_val);
            for(int i = 0; i < NX; i++) { 
                Edist[i] = 0.0f; 
                for(int j = 0; j < ND; j++) Edist[i] += E[i][j] * d_est[j]; 
            }

            #if FORMULATION == OPT_CONDENSED
                build_AB(Ad, Bd, A_pow, AB_mat);
                build_P(P, AB_mat, Q_diag, R_diag, P_diag);
                build_A(A, AB_mat, jacobian_H);
                osqp_update_data_mat(solver, P->x, NULL, P->nzmax, A->x, NULL, A->nzmax);
            #elif FORMULATION == OPT_NON_CONDENSED
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

        // 5) Vector update path
        #if FORMULATION == OPT_CONDENSED
            build_q(q, Ad, d_lin, Edist, x_current, x_target, A_pow, AB_mat, Q_diag, P_diag);
            build_bounds(&l, &u, Ad, d_lin, Edist, jacobian_H, h_val, x_current, x_min, x_max, u_min, u_max, A_pow);
        #elif FORMULATION == OPT_NON_CONDENSED
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

        // 6) Warm-start policy
        if (recompute_mats) {
            osqp_warm_start(solver, NULL, zero_y); // Zero duals on re-linearization to prevent instability
        } else {
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

        // Final Logging for Step
        fprintf(f_out, "%d", step);
        for (int i = 0; i < NX; i++) fprintf(f_out, ",%f", log_x[i]);
        for (int i = 0; i < NU; i++) fprintf(f_out, ",%f", log_u[i]);
        fprintf(f_out, "\n");

        timer_now(&iter_end);
        iter_total_ms += timer_diff_ms(&iter_start, &iter_end);
        steps_completed++;
    }
    
    if (steps_completed > 0) {
        solver_avg_ms = solver_total_ms / steps_completed;
        iter_avg_ms = iter_total_ms / steps_completed;
    }

    // Calculate final true cost of objective function
    OSQPFloat final_cost = 0.0f;
    for (int i = 0; i < NX; i++) {
        OSQPFloat err_x = x_current[i] - x_target[i];
        final_cost += err_x * (Q_diag[i] / DT) * err_x + err_x * P_diag[i] * err_x;
    }
    for (int i = 0; i < NU; i++) {
        OSQPFloat err_u = u_applied[i] - u_target[i];
        final_cost += err_u * (R_diag[i] / DT) * err_u;
    }

    printf("Actual Tsim -> %.3f s\n", (double)(steps_completed * DT));
    printf("Solver time -> total: %.3f ms, avg: %.3f ms\n", solver_total_ms, solver_avg_ms);
    printf("Iteration time-> total: %.3f ms, avg: %.3f ms\n", iter_total_ms, iter_avg_ms);
    printf("Final End Cost (Raw) -> %.3e\n", (double)final_cost);
    printf("Final End Cost (log10) -> %.3f\n", (final_cost > 0) ? log10((double)final_cost) : -99.0);
    printf("Linearizations performed: %lld out of %lld steps\n\n", linearized_times, steps_completed);
    printf("%.3f/%.3f/%.3f\n\n", solver_avg_ms, iter_avg_ms, (final_cost > 0) ? log10((double)final_cost) : -99.0);

    fclose(f_out);

    osqp_cleanup(solver);
    free(P->x); free(P->i); free(P->p);
    free(A->x); free(A->i); free(A->p);
    free(q); free(l); free(u);
    free(zero_x); free(zero_y);
    free(settings);

    return 0;
}