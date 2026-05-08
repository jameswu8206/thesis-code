#ifndef MPC_PIPELINE_HELPERS_H
#define MPC_PIPELINE_HELPERS_H

#include "osqp.h"

/*
 * Pipeline helper dimensions are kept fixed to match current mpc_pipeline.c
 * behavior exactly (NX=6, NU=2, NH=3, DT=0.020).
 */

typedef struct {
    OSQPInt *i;
    OSQPInt *p;
    OSQPFloat *x;
    OSQPInt nz;
    OSQPInt n_cols;
} CscBuilder;

void csc_set(CscBuilder *m, OSQPInt row, OSQPInt col, OSQPFloat val);
void csc_col_done(CscBuilder *m);

OSQPFloat vec_norm(const OSQPFloat *v, OSQPInt len);

void mat_copy(OSQPFloat *dst, const OSQPFloat *src, OSQPInt rows, OSQPInt cols);
void mat_eye(OSQPFloat *m, OSQPInt rows, OSQPInt cols);
void mat_mul(const OSQPFloat *A, const OSQPFloat *B, OSQPFloat *C, OSQPInt m, OSQPInt k, OSQPInt n);
void mat_set(OSQPFloat *A, OSQPFloat v, OSQPInt rows, OSQPInt cols);

void crane_dynamics(const OSQPFloat x[6], const OSQPFloat u[2], OSQPFloat f_out[6]);
void calc_h_gradient(const OSQPFloat x_current[6], OSQPFloat jacobian_H[3][6], OSQPFloat h_val[3]);
void linearization(OSQPFloat x_current[6], OSQPFloat u_applied[2], OSQPFloat Ad[6][6], OSQPFloat Bd[6][2], OSQPFloat d_lin[6]);

void pend_dynamics(const OSQPFloat x[2], const OSQPFloat u[1], OSQPFloat f_out[2]);
void pend_calc_h_gradient(const OSQPFloat x_current[2], OSQPFloat *jacobian_H, OSQPFloat *h_val, int nh);
void pend_linearization(OSQPFloat x_current[2], OSQPFloat u_applied[1], OSQPFloat Ad[2][2], OSQPFloat Bd[2][1], OSQPFloat d_lin[2], OSQPFloat dt);

#endif
