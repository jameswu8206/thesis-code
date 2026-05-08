#ifndef INVERTED_PENDULUM_FORMULATION_HELPERS_H
#define INVERTED_PENDULUM_FORMULATION_HELPERS_H

#include "osqp.h"
#include "inverted_pendulum_config.h"

#if FORMULATION == OPT_CONDENSED
void build_AB(const OSQPFloat Ad[NX][NX], const OSQPFloat Bd[NX][NU], OSQPFloat A_pow[N + 1][NX][NX], OSQPFloat AB_mat[N][NX][NU]);
OSQPCscMatrix* build_P(OSQPCscMatrix* P, const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]);
OSQPFloat* build_q(OSQPFloat *q, const OSQPFloat Ad[NX][NX], const OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat x_curr[NX], const OSQPFloat x_target[NX], const OSQPFloat A_pow[N+1][NX][NX], const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat Q_diag[NX], const OSQPFloat P_diag[NX]);
OSQPCscMatrix* build_A(OSQPCscMatrix* A, const OSQPFloat AB_mat[N][NX][NU], const OSQPFloat jacobian_H[NH][NX]);
void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat Ad[NX][NX], const OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat jacobian_H[NH][NX], const OSQPFloat h_val[NH], const OSQPFloat x_current[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPFloat A_pow[N+1][NX][NX]);

#elif FORMULATION == OPT_NON_CONDENSED
OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]);
OSQPFloat* build_q(const OSQPFloat* x_target, const OSQPFloat* u_target, const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX]);
OSQPCscMatrix* build_A(OSQPCscMatrix* A,OSQPFloat Ad[NX][NX],OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX]);
void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_vals[NH],OSQPFloat d_lin[NX], const OSQPFloat Edist[NX], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU]);

#elif FORMULATION == OPT_SPARSE_CONDENSED
void compute_sc_blocks(OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat jacobian_H[NH][NX], OSQPFloat K[NU][NX], OSQPFloat AK[NX][NX], OSQPFloat P_blocks[], OSQPFloat M_blocks[], OSQPFloat H_blocks[]);
OSQPCscMatrix* build_P(const OSQPFloat Q_diag[NX], const OSQPFloat R_diag[NU], const OSQPFloat P_diag[NX], const OSQPFloat K[NU][NX], const OSQPFloat P_blocks[], const OSQPFloat M_blocks[], const OSQPFloat e_free[][NX], OSQPFloat q[], OSQPCscMatrix *P_existing, OSQPInt r_local, OSQPInt recompute_hessian);
OSQPCscMatrix* build_A(OSQPCscMatrix* A, const OSQPFloat P_blocks[], const OSQPFloat M_blocks[], const OSQPFloat H_blocks[], OSQPInt r_local);
void build_bounds(OSQPFloat **l_ptr, OSQPFloat **u_ptr, const OSQPFloat x_current[NX], const OSQPFloat e_free[][NX], const OSQPFloat x_target[NX], const OSQPFloat u_target[NU], const OSQPFloat K[NU][NX], const OSQPFloat jacobian_H[NH][NX], const OSQPFloat h_vals[NH], const OSQPFloat x_min[NX], const OSQPFloat x_max[NX], const OSQPFloat u_min[NU], const OSQPFloat u_max[NU], OSQPInt r_local);
#endif

#endif
