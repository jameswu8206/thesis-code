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
#define N 100   // Prediction Horizon steps
#define ND  6
#define NH  3 //numbers of nonlinear constraints
#define DT 0.02 // Sampling time (s)
#define Tsim 20.0// full simulation time



// Total variables: (N+1)*NX states + N*NU controls
// Ordered as: [x0, u0, x1, u1, ... , xN]
OSQPInt n_vars = (N + 1) * NX + N * NU;
    
// 1. Dynamics & x0 Equality: (N+1)*NX rows
// 2. General Inequalities:   N*NH rows (The "General Case" placeholder)
// 3. Simple Box Bounds:      n_vars rows (x_min < x < x_max)
OSQPInt n_cons = (N + 1) * NX + (N * NH) + ((N + 1) * NX + N * NU);    

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
//need to modify according to problem
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
//need to modify according to problem
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

void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_vals[NH],OSQPFloat d_lin[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPFloat d_est[ND], const OSQPFloat E[NX][NX]){

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
    //calculate Ed first
    //update equality bounds
    OSQPFloat Edist[NX];
    for(int i=0; i<NX; i++) {
        Edist[i] = 0.0;
        for(int j=0; j<ND; j++) {
            Edist[i] += E[i][j] * d_est[j];
        }
    }
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < NX; i++) {
            l[NX + k * NX + i] =d_lin[i] -Edist[i];
            u[NX + k * NX + i] =d_lin[i] -Edist[i];
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


typeTime start, end;
OSQPFloat *CPUtimeVec;
OSQPFloat CPUtime = 0;
OSQPInt iter_used=0;
OSQPFloat total_time;
OSQPFloat avg_time;



int main(void) {


    // System Dynamics (Discretized Euler): x_k+1 = Ad*x_k + Bd*u_k
    // Ad = I +Ac*DT
    // Bd = Bc*DT
    // E for disturbance is known
    OSQPFloat Ad[NX][NX] = {0};
    OSQPFloat Bd[NX][NU] = {0};
    const OSQPFloat E[NX][NX] = {
                                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
                                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 
                                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},};
    OSQPFloat jacobian_H[NH][NX] = {0};
    OSQPFloat h_val[NH] = {0};
    OSQPFloat d_lin[NX]={0};

    // Cost Weights (Derived from your pCost)
    // Q (State cost), R (Control cost), P (Terminal cost)
    //1/2( xQx + uRu + x'Px' ) where x' is desired position
    const OSQPFloat Q_diag[NX] = {1.0*DT, 2.0*DT, 2.0*DT, 1.0*DT, 1.0*DT, 4.0*DT}; //for fast response
    const OSQPFloat R_diag[NU] = {0.05*DT,0.05*DT}; //penalize large input       
    const OSQPFloat P_diag[NX] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};//for damping fater to desired point(optimize?)

    // Desired Setpoints (Target)
    //Initial State Equality -x0=-x0
    OSQPFloat x_current[NX] = { -2.0, 0.0, 2.0, 0.0, 0.0, 0.0 };//initial point
    OSQPFloat x_target[NX] = { 2.0, 0.0, 2.0, 0.0, 0.0, 0.0 };//target point
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


    OSQPCscMatrix* P =build_P(Q_diag, R_diag, P_diag);

    
    OSQPFloat *q = build_q(x_target, u_target, Q_diag, R_diag, P_diag);

    OSQPCscMatrix* A = build_A(NULL, Ad, Bd, jacobian_H);

    

    //fill the bounds l, u
    OSQPFloat *l=NULL, *u=NULL;
    build_bounds(&l, &u, x_current, jacobian_H, h_val, d_lin, x_min, x_max, u_min, u_max, d_est, E);

    //for open loop(no controll u)
    OSQPFloat x_open_loop[NX];
    for(int i=0; i<NX; i++) x_open_loop[i] = x_current[i];

    //setup for osqp solver
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    osqp_set_default_settings(settings);

    OSQPInt opt_time=0;
    if(opt_time==1){//best timing but jagged behavior
        settings->adaptive_rho = 1;//must enable
        settings->check_termination = 1;
        settings->eps_abs = 1e-1;//loosing the tolerance
        settings->eps_rel = 1e-1;
        settings->time_limit = 0.008;//stop asap
        settings->alpha = 1.5;//relaxation
        settings->verbose = 0;
    }if(opt_time==0){//do nothing
        settings->alpha = 1.0;
        settings->scaling= 10;
        settings->verbose = 0;
        settings->adaptive_rho = 1;
        settings->check_termination = 5;
        settings->max_iter = 10000;//make sure is not reached for 100% solved
        settings->warm_starting = 1;//default for grampc
        settings->eps_abs = 5e-3;
        settings->eps_rel = 5e-3;
    }if(opt_time==2){// smooth while still better in time dur to optimized P cost
        settings->adaptive_rho = 1;//must enable
        settings->check_termination = 5;
        settings->eps_abs = 1e-2; 
        settings->eps_rel = 1e-2;
        settings->time_limit = 0.01;
        settings->alpha = 1.5;
        settings->verbose = 0;
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

    OSQPFloat x_next[NX];
    OSQPFloat x_warm[n_vars];//for warm start shifts
    //for open loop
    OSQPFloat x_ol_next[NX];


    for (int step = 0; step < MAX_STEPS; step++) {

        linearization(x_current, u_applied, Ad, Bd, d_lin);
        calc_h_gradient(x_current, jacobian_H, h_val);

       if(step%10==0){

        
        build_A(A, Ad, Bd, jacobian_H);
        osqp_update_data_mat(solver, NULL, NULL, 0,A->x, NULL,A->nzmax);
        
       }
        
        build_bounds(&l, &u, x_current, jacobian_H, h_val, d_lin, x_min, x_max, u_min, u_max, d_est, E);
        osqp_update_data_vec(solver, NULL, l, u);
        
        timer_now(&start);
        osqp_solve(solver);
        timer_now(&end);
        CPUtimeVec[step] = timer_diff_ms(&start, &end);

        
        /*
        if (solver->info->status_val != OSQP_SOLVED && 
            solver->info->status_val != OSQP_SOLVED_INACCURATE && 
            solver->info->status_val != OSQP_MAX_ITER_REACHED) {
            
            printf("Solver fundamentally failed! Status: %d\n", solver->info->status_val);
            break;
        }
    
        */
       if (solver->info->status_val != OSQP_SOLVED) {
            printf("Solver failed!\n");
            break;
        }
        // Accept standard solve, inaccurate solve, or time-limit reached.
        // Only fail on primal/dual infeasibility or actual setup errors.
        

        //extract next input u1
        for (int i = 0; i < NU; i++) {
            u_applied[i] = solver->solution->x[NX + i];
        }
        //warm start shift(default for grampc)
        // --- CORRECT INTERLEAVED WARM START SHIFT ---
        OSQPInt stride = NX + NU;
        // 1. Shift x_k and u_k backward by one time step
        for (int k = 0; k < N - 1; k++) {
            // Shift States
            for (int i = 0; i < NX; i++) {
                x_warm[k * stride + i] = solver->solution->x[(k + 1) * stride + i];
            }
            // Shift Controls
            for (int i = 0; i < NU; i++) {
                x_warm[k * stride + NX + i] = solver->solution->x[(k + 1) * stride + NX + i];
            }
        }
        // 2. Handle the horizon tail (Steady State Guesses)
        // Shift the terminal state x_N backward into x_{N-1}
        for (int i = 0; i < NX; i++) {
            x_warm[(N - 1) * stride + i] = solver->solution->x[N * stride + i];
        }
        // Duplicate the final control u_{N-1} to act as the new tail control
        for (int i = 0; i < NU; i++) {
            x_warm[(N - 1) * stride + NX + i] = solver->solution->x[(N - 1) * stride + NX + i];
        }
        // Duplicate the terminal state x_N to act as the new terminal state
        for (int i = 0; i < NX; i++) {
            x_warm[N * stride + i] = solver->solution->x[N * stride + i];
        }
        // 3. Inject into OSQP
        osqp_warm_start(solver, x_warm, NULL);
        // --------------------------------------------


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