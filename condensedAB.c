#include "osqp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_timing.h"


//parts to modify: A,B,Q,R,P,q; max/min ; NX NU DT N ; initial and desired point


/* ---------------------------------------------------------
 * Problem Constants & Parameters
 * --------------------------------------------------------- */
#define NX 6    // Number of States
#define NU 2    // Number of Controls
#define N  50   // Prediction Horizon steps
#define ND  6
#define DT 0.01 // Sampling time (s)
#define Tsim 30.0// default 11



// Total variables: (N+1)*NX states + N*NU controls
// Ordered as: [u0~uN-1]
OSQPInt n_vars = NU*N;
    
// Total constraints: 
// 1. inequality of u: N * NU 
// 2. states equality(x): NX*N
OSQPInt n_cons = NU*N+NX*N;

typedef struct {
    OSQPInt *i;
    OSQPInt *p;
    OSQPFloat *x;
    OSQPInt nz;
    OSQPInt n_cols;
} CscBuilder;

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

void build_AB(const OSQPFloat Ad[NX][NX], const OSQPFloat Bd[NX][NU], 
                               OSQPFloat A_pow[N + 1][NX][NX], OSQPFloat AB_mat[N][NX][NU]) {

    // 1. Initialize A^0 = Identity Matrix
    for (int r = 0; r < NX; r++) {
        for (int c = 0; c < NX; c++) {
            A_pow[0][r][c] = (r == c) ? 1.0 : 0.0;
        }
    }

    // 2. Recursive Generation Loop
    for (int k = 1; k <= N; k++) {
        
        // --- Calculate A^k = Ad * A^(k-1) ---
        for (int r = 0; r < NX; r++) {
            for (int c = 0; c < NX; c++) {
                A_pow[k][r][c] = 0.0;
                for (int m = 0; m < NX; m++) {
                    A_pow[k][r][c] += Ad[r][m] * A_pow[k - 1][m][c];
                }
            }
        }

        // --- Calculate A^(k-1)B = A^(k-1) * Bd ---
        // Note: AB_mat[0] is B, AB_mat[1] is AB, etc.
        for (int r = 0; r < NX; r++) {
            for (int c = 0; c < NU; c++) {
                AB_mat[k - 1][r][c] = 0.0f;
                for (int m = 0; m < NX; m++) {
                    AB_mat[k - 1][r][c] += A_pow[k - 1][r][m] * Bd[m][c];
                }
            }
        }
    }
}


OSQPCscMatrix* build_P(OSQPCscMatrix* P,const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]){


    OSQPInt setup = (P == NULL);
    CscBuilder P_build = {0};

    if (setup) {
        OSQPInt P_nnz = (n_vars + 1) * n_vars / 2;
        P = malloc(sizeof(OSQPCscMatrix));
        P->x = malloc(P_nnz * sizeof(OSQPFloat));
        P->i = malloc(P_nnz * sizeof(OSQPInt));
        P->p = malloc((n_vars + 1) * sizeof(OSQPInt));
        P->nzmax = P_nnz;
    }

    // Initialize/Reset Builder pointers
    P_build.x = P->x;
    P_build.i = P->i;
    P_build.p = P->p;
    P_build.nz = 0;
    P_build.n_cols = 0;
    P_build.p[0] = 0;
    //filling costs(big P)
    for (int j = 0; j < n_vars; j++) { // n_vars columns
        OSQPInt step_j = j / NU;           // Which time step (0 to 29)
        OSQPInt ctrl_j = j % NU;           // Which control input (0 or 1)

        for (int i = 0; i <= j; i++) { // Row i (Upper triangle: i <= j)
            OSQPInt step_i = i / NU;
            OSQPInt ctrl_i = i % NU;
            
            OSQPFloat val = 0.0;

            // 1. Add Control Cost (R) - only diagonal
            if (i == j) val += R_diag[ctrl_i];

            // 2. Add State Cost (T' * Q * T)
            // A control input u_i only affects states at k > step_i
            // So we sum from the max step of i and j up to N
            OSQPInt start_k = (step_i > step_j ? step_i : step_j) + 1;
            
            for (int k = start_k; k <= N; k++) {
                // Use Terminal Cost P if k == N, otherwise Stage Cost Q
                const OSQPFloat* Q_k = (k == N) ? P_diag : Q_diag;
                
                // Matrix math: P_ij += sum_{s=1 to NX} ( (A^{k-step_i-1}B)_s * Q_s * (A^{k-step_j-1}B)_s )
                for (int s = 0; s < NX; s++) {
                    OSQPFloat term_i = AB_mat[k - step_i - 1][s][ctrl_i];
                    OSQPFloat term_j = AB_mat[k - step_j - 1][s][ctrl_j];
                    val += term_i * Q_k[s] * term_j;
                }
            }
            val *= 2.0;
            if (setup) {
                csc_set(&P_build, i, j, val);
            } else {
                // 2. If already setup, just overwrite the value at the current nz index
                P_build.x[P_build.nz] = val;
                P_build.nz++;
            }
        }
        if (setup) csc_col_done(&P_build);
        else P_build.n_cols++;
    }

    //copy filled big P
    if (setup) {
        P->m = n_vars; P->n = n_vars; P->nz = -1;
    }
    return P;
}

OSQPFloat* build_q(OSQPFloat *q, 
                        const OSQPFloat x_curr[NX], 
                        const OSQPFloat x_target[NX], 
                        const OSQPFloat A_pow[N+1][NX][NX], 
                        const OSQPFloat AB_mat[N][NX][NU], 
                        const OSQPFloat Q_diag[NX], 
                        const OSQPFloat P_diag[NX]){

    OSQPInt setup = (q == NULL);

    if (setup) {
        // Allocate only on the first call
        q = calloc(n_vars, sizeof(OSQPFloat));
    } else {
        // Reset the vector for the update
        for (int i = 0; i < n_vars; i++) q[i] = 0.0;
    }

    // Outer Loop: For each control input u_j (Column of T)
    for (int j = 0; j < n_vars; j++) {
        OSQPInt step_j = j / NU;           // Time step (0 to N-1)
        OSQPInt ctrl_j = j % NU;           // Control index (0 to NU-1)

        // Inner Loop: Sum the impact of u_j on every future state x_k
        for (int k = step_j + 1; k <= N; k++) {
            // Select Q (stage cost) or P (terminal cost)
            const OSQPFloat* Q_k = (k == N) ? P_diag : Q_diag;

            // 1. Calculate Predicted Error: (A^k * x_0 - x_target)
            OSQPFloat error_k[NX];
            for (int r = 0; r < NX; r++) {
                OSQPFloat x_free_rk = 0.0f;
                for (int c = 0; c < NX; c++) {
                    x_free_rk += A_pow[k][r][c] * x_curr[c];
                }
                error_k[r] = x_free_rk - x_target[r];
            }

            // 2. Project Error onto Control j: 2 * error' * Q * Waterfall_Column
            OSQPInt p_idx = k - step_j - 1; // Correct index for AB_mat
            for (int r = 0; r < NX; r++) {
                // q_j += 2 * (Error * Weight * Impulse_Response)
                q[j] += 2.0 * error_k[r] * Q_k[r] * AB_mat[p_idx][r][ctrl_j];
            }
        }
    }
    return q;

}

OSQPCscMatrix* build_A(OSQPCscMatrix* A, const OSQPFloat AB_mat[N][NX][NU]){
    //build big A
    // Init state: NX
    // Dynamics: N * (NX*NX for A + NX*NU for B + NX for -I) (Ad + Bd + -I)
    // Bounds: n_vars (Identity diagonal)
    OSQPInt setup=(A==NULL);

    CscBuilder A_build = {0};

    if(setup){
        OSQPInt A_nnz = (N * (N + 1) / 2) * (NX * NU) + n_vars;
        A = malloc(sizeof(OSQPCscMatrix));
        A->x = malloc(A_nnz * sizeof(OSQPFloat));
        A->i = malloc(A_nnz * sizeof(OSQPInt));
        A->p = malloc((n_vars + 1) * sizeof(OSQPInt));
        A->nzmax = A_nnz; // Added: nzmax initialization
    }
    
    A_build.x = A->x; 
    A_build.i = A->i; 
    A_build.p = A->p;
    A_build.nz = 0;
    A_build.n_cols = 0;
    A_build.p[0] = 0;

    //CscBuilder A_build = {A->i, A->p, A->x, 0, 0}; 
    //A_build.p[0] = 0;

    //filling big A
    for (int j = 0; j < n_vars; j++) { // Column j
        int step_j = j / NU;           // Which time step (0 to 29)
        int ctrl_j = j % NU;           // Which control index (0 or 1)

        // --- PART 1: IDENTITY ZONE (u constraints) ---
        // Puts a 1.0 on the diagonal of the top 60x60 block
        csc_set(&A_build, j, j, 1.0);

        // --- PART 2: WATERFALL ZONE (x constraints) ---
        // Control u_j affects states starting from x_{step_j + 1}
        for (int k = step_j + 1; k <= N; k++) {
            // Calculate the power index: (k - step_j - 1)
            // This pulls the correct A^n*B block from your pre-computed array
            int p_idx = k - step_j - 1;

            for (int r = 0; r < NX; r++) {
                float val = AB_mat[p_idx][r][ctrl_j];
                
                // Row mapping: 
                // Offset by n_vars (60) + (k-1)*NX + current state row r
                OSQPInt row_idx = n_vars + (k - 1) * NX + r;
                csc_set(&A_build, row_idx, j, val);
            }
        }
        csc_col_done(&A_build);
    }
    //copy filled big A
    if(setup){
        A->m = n_cons; A->n = n_vars; A->nz = -1; A->nzmax = A_build.nz;
    }
    
    return A;
};

void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPFloat A_pow[N+1][NX][NX], OSQPFloat d_est[ND], const OSQPFloat E[NX][NX]){

    OSQPInt setup=(*l_ptr==NULL);

    if(setup){
        *l_ptr = malloc(n_cons * sizeof(OSQPFloat));
        *u_ptr = malloc(n_cons * sizeof(OSQPFloat));
    }

    OSQPFloat *l = *l_ptr;
    OSQPFloat *u = *u_ptr;

    //for upper identities to restrict u bounds
    for (int i = 0; i < N * NU; i++) {
        l[i] = u_min[i % NU];
        u[i] = u_max[i % NU];
    }

    // --- PART 2: SHIFTED STATE BOUNDS (Rows N*NU to End) ---
    // Formula: x_min - A^k * x_0 <= Tau * U <= x_max - A^k * x_0
    int b_idx = N * NU; 
    for (int k = 1; k <= N; k++) {
        // 1. Calculate the Free Response for step k (A^k * x_curr)
        float x_free_k[NX];
        for (int r = 0; r < NX; r++) {
            x_free_k[r] = 0.0f;
            for (int c = 0; c < NX; c++) {
                x_free_k[r] += A_pow[k][r][c] * x_current[c];
            }
        }

        // 2. Shift the original x_min/x_max by this free response
        for (int i = 0; i < NX; i++) {
            l[b_idx] = x_min[i] - x_free_k[i];
            u[b_idx] = x_max[i] - x_free_k[i];
            b_idx++;
        }
    }
    
}


OSQPFloat Ad[NX][NX];
OSQPFloat Bd[NX][NU];
void linearization( OSQPFloat x_current[NX], OSQPFloat u_applied[NU], 
                                OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU]) {
    
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
                                {0.0, 0.0, t1, t2, t3, t4},};

    // Assign to Bc (df/du)

    OSQPFloat Bc[NX][NU] = {
                                {0.0, 0.0}, 
                                {1.0, 0.0},
                                {0.0, 0.0}, 
                                {0.0, 1.0},
                                {0.0, 0.0},
                                {t5, 0.0}};
    
    // ====================================================================
    // 3. Discretization (Euler: Ad = I + Ac*dt, Bd = Bc*dt)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NX; j++) {
            Ad[i][j] = (i == j ? 1.0 : 0.0) + Ac[i][j] * DT;
        }
        for (int j = 0; j < NU; j++) {
            Bd[i][j] = Bc[i][j] * DT;
        }
    }
}


typeTime start, end;
OSQPFloat *CPUtimeVec;
OSQPFloat CPUtime = 0;
OSQPFloat total_time;
OSQPFloat avg_time;



int main(void) {


    // System Dynamics (Discretized Euler): x_k+1 = Ad*x_k + Bd*u_k
    // Ad = I +Ac*DT
    // Bd = Bc*DT
    // E for disturbance is known
    /*
    const OSQPFloat Ad[NX][NX] = {
                                {1.0, 0.0, DT, 0.0, 0.0, 0.0}, 
                                {0.0, 1.0, 0.0, DT, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                {DT, 0.0, 0.0, 0.0, 1.0, 0.0}, 
                                {0.0, DT, 0.0, 0.0, 0.0, 1.0},};
    const OSQPFloat Bd[NX][NU] = {
                                {0.0, 0.0}, 
                                {0.01*DT, -0.01*DT},
                                {0.19*DT, 0.19*DT}, 
                                {1.32*DT, -1.32*DT},
                                {0.0, 0.0},
                                {0.0, 0.0}};*/
    const OSQPFloat E[NX][NX] = {
                                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
                                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 
                                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},};
    
    OSQPFloat A_pow[N + 1][NX][NX];
    OSQPFloat AB_mat[N][NX][NU];
    // Cost Weights (Derived from your pCost)
    // Q (State cost), R (Control cost), P (Terminal cost)
    //1/2( xQx + uRu + x'Px' ) where x' is desired position
    const OSQPFloat Q_diag[NX] = {1.0, 2.0, 2.0, 1.0, 1.0, 4.0}; //for fast response
    const OSQPFloat R_diag[NU] = {0.05,0.05}; //penalize large input       
    const OSQPFloat P_diag[NX] = {1.0, 2.0, 2.0, 1.0, 1.0, 4.0};//for damping fater to desired point(optimize?)

    // upper lower bounds
    const OSQPFloat x_min[NX] = {-10.0, -10.0, -10.0, -10.0, -10.0, -10.0};
    const OSQPFloat x_max[NX] = { 10.0,  10.0,  10.0,  10.0,  10.0,  10.0};
    const OSQPFloat u_min[NU] = {-2.0, -2.0};
    const OSQPFloat u_max[NU] = { 2.0,  2.0};

    //estimate possible disturbance
    OSQPFloat d_est[ND] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    
    // Desired Setpoints (Target)
    //Initial State Equality -x0=-x0
    OSQPFloat x_current[NX] = { -2.000,0.000,2.000,0.000,0.000,0.000 };//initial point
    OSQPFloat x_target[NX] = { 2.000,0.000,2.000,0.000,0.000,0.000 };//target point
    OSQPFloat u_target[NU] = {0.0, 0.0};
    OSQPFloat x_next[NX];
    OSQPFloat u_applied[NU]={0.0,0.0};
    //for open loop
    OSQPFloat x_ol_next[NX];
    linearization(x_current,u_applied,Ad,Bd);

    build_AB(Ad,Bd,A_pow,AB_mat);
    OSQPCscMatrix* P =NULL;
    P =build_P(P,AB_mat,Q_diag, R_diag, P_diag);

    OSQPFloat *q = NULL;
    q=build_q(q,x_current,x_target,A_pow,AB_mat,Q_diag,P_diag);

    OSQPCscMatrix* A = build_A(NULL,AB_mat);

    

    //fill the bounds l, u
    OSQPFloat *l=NULL, *u=NULL;
    build_bounds(&l, &u, x_current, x_min, x_max, u_min, u_max, A_pow, d_est, E);

    //for open loop(no controll u)
    OSQPFloat x_open_loop[NX];
    for(int i=0; i<NX; i++) x_open_loop[i] = x_current[i];

    //setup for osqp solver
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);

    OSQPInt opt_time=1;
    if(opt_time==1){//best timing but jagged behavior
        settings->adaptive_rho = 1;//must enable
        settings->check_termination = 1;
        settings->eps_abs = 1e-1; //loosing the tolerance
        settings->eps_rel = 1e-1;
        settings->time_limit = 0.008;//stop asap
        settings->alpha = 1.5;//relaxation
        settings->verbose = 0;
    }if(opt_time==0){//do nothing
        settings->alpha = 1.0;
        settings->verbose = 0;
    }if(opt_time==2){// smooth while still better in time dur to optimized P cost
        settings->adaptive_rho = 1;//must enable
        settings->check_termination = 5;
        settings->eps_abs = 1e-2; 
        settings->eps_rel = 1e-2;
        settings->time_limit = 0.01;
        settings->alpha = 1.5;
        settings->verbose = 0;
    }
    if(opt_time==3){//special for this case
        settings->adaptive_rho = 1;
        settings->check_termination = 1;//exit immediately when hit accurarcy target
        settings->eps_abs = 1e-1; 
        settings->eps_rel = 1e-1;
        settings->alpha = 1.6;
        settings->adaptive_rho_interval = 1;
        //settings->time_limit = 0.0025;//has to change together with eps
        settings->verbose = 0;
        settings->scaling = 0;//no scaling since stacking already condensed
    }
    








    OSQPSolver *solver;
    osqp_setup(&solver, P, q, A, l, u, n_cons, n_vars, settings);
    
    

    //MPC loop
    int MAX_STEPS = Tsim/DT; 
    CPUtimeVec = (OSQPFloat*)calloc(MAX_STEPS + 1, sizeof(*CPUtimeVec));
	if (CPUtimeVec == NULL) {
        fprintf(stderr, "Error: Allocating memory for computation time measurement failed\n");
        exit(EXIT_FAILURE);	
    }

    
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
    // for open loop
    //for(int i=0; i<NX; i++) fprintf(f_out, ",x_ol%d", i);
    fprintf(f_out, "\n");

    printf("Starting MPC Loop (NX=%d, NU=%d)...\n", NX, NU);
    printf("Step |");
    for(int i=0; i<NX; i++) printf("   x%d   ", i);
    printf("|");
    for(int i=0; i<NU; i++) printf("   u%d   ", i);
    printf("\n");
    
    printf("-----|");
    for(int i=0; i<NX; i++) printf("-------");
    printf("|");
    for(int i=0; i<NU; i++) printf("-------");
    printf("\n");

    


    for (int step = 0; step < MAX_STEPS; step++) {

        //timer_now(&start);
        /*/dynamic disturbances
        if(step<=300){
            for(int i=0; i<NX; i++){
                d_est[i]=0.0;
            }
        }if(step<=600 && step>300){
            for(int i=0; i<NX; i++){
                d_est[i]=0.0;
            }
        }if(step<=900 && step>600){
            for(int i=0; i<NX; i++){
                d_est[i]=0.0;
            }
        }*/

       if((step%1==0&&step<250)||(step%3==0&&step>250&&step<700)||(step%30==0&&step>700&step<1500)||(step%50==0&&step>1500)){
            linearization(x_current,u_applied,Ad,Bd);
            build_AB(Ad,Bd,A_pow,AB_mat);
            P =build_P(P,AB_mat,Q_diag, R_diag, P_diag);
            A = build_A(A,AB_mat);
            osqp_update_data_mat(solver,P->x,NULL,P->nzmax,A->x,NULL,A->nzmax);
       }
        
        build_bounds(&l, &u, x_current, x_min, x_max, u_min, u_max,A_pow, d_est, E);
        q=build_q(q,x_current,x_target,A_pow,AB_mat,Q_diag,P_diag);
        osqp_update_data_vec(solver, q, l, u);
        
        //osqp_warm_start_x(solver, solver->solution->x);
        //osqp_warm_start_y(solver, solver->solution->y);

        timer_now(&start);
        osqp_solve(solver);
        timer_now(&end);
        CPUtimeVec[step] = timer_diff_ms(&start, &end);

        
        if (solver->info->status_val != OSQP_SOLVED) {
            printf("Solver failed!\n");
            break;
        }

        //extract next input u1
        for (int i = 0; i < NU; i++) {
            u_applied[i] = solver->solution->x[i];
        }

        fprintf(f_out, "%d", step);
        for (int i = 0; i < NX; i++) fprintf(f_out, ",%f", x_current[i]);
        for (int i = 0; i < NU; i++) fprintf(f_out, ",%f", u_applied[i]);
        //for open loop
        //for (int i = 0; i < NX; i++) fprintf(f_out, ",%f", x_open_loop[i]);
        fprintf(f_out, "\n");

        printf("%d|", step);
        // Print all states
        for(int i=0; i<NX; i++) {
            printf(" %6.5f", x_current[i]);
        }
        printf("|");
        // Print all controls
        for(int i=0; i<NU; i++) {
            printf(" %6.5f", u_applied[i]);
        }
        printf("\n");

        //manual calculate the next point if u1 appied, since x in sol is ideal condition, by manually adding we can add disturbances
        for (int i = 0; i < NX; i++) {
            x_next[i] = 0.0;
            
            // + Ad * x term
            for (int j = 0; j < NX; j++) {
                if (Ad[i][j] != 0.0) { // Optimization: Skip zeros
                    x_next[i] += Ad[i][j] * x_current[j];
                }
            }
            
            // + Bd * u term
            for (int j = 0; j < NU; j++) {
                if (Bd[i][j] != 0.0) {
                    x_next[i] += Bd[i][j] * u_applied[j];
                }
            }
        }
        //for open loop
        for (int i = 0; i < NX; i++) {
            x_ol_next[i] = 0.0;
            // Only add Natural Physics (Ad), ignore Control (Bd)
            for (int j = 0; j < NX; j++) {
                if (Ad[i][j] != 0.0) x_ol_next[i] += Ad[i][j] * x_open_loop[j];
            }
        }
        // Update state for next loop
        for (int i = 0; i < NX; i++) {
            x_current[i] = x_next[i];
            x_open_loop[i] = x_ol_next[i];
        }
        
        /*

        //alternative, assuming model mathematically correct, directly extract x from sol
        for (int i = 0; i < NX; i++) {
        x_current[i] = solver->solution->x[NX + NU + i];
        }
        
        
        */

        /*add disturbance
        if (step == 60|| step == 300|| step == 600|| step == 900|| step == 1100) {
             printf("   >> SPLASH! Disturbance added << \n");
             x_current[2]-=0.5;
             x_current[3]+=0.5;

             //for open loop(making it diverge)
             x_open_loop[2] -= 0.5;
             x_open_loop[3] += 0.5;
             
        }*/

        //check terminate condition
        int converged = 1;
        for(int i=0; i<NX; i++) {
            if (fabs(x_current[i] - x_target[i]) > 0.01) {
                converged = 0;
                break;
            }
        }
        if (converged && step > 10) { // Ensure we ran at least a bit
            printf("Target Reached at step %d!\n", step);
            break;
        }
    }


    
    for (int i = 0; i <= MAX_STEPS; i++) {
		CPUtime = CPUtime + CPUtimeVec[i];
	}
    printf("total time: %.4f ms\n", CPUtime);
	CPUtime = CPUtime / (MAX_STEPS + 1);
    printf("avg time: %.4f ms\n",CPUtime);


    fclose(f_out);


    osqp_cleanup(solver);
    free(CPUtimeVec);
    free(P->x); free(P->i); free(P->p);
    free(A->x); free(A->i); free(A->p);
    free(q); free(l); free(u);
    free(settings);

    return 0;
}