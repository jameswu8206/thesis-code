#include "mpc_pipeline_helpers.h"

#include <math.h>

#define NX 6
#define NU 2
#define NH 3
#define DT 0.020

void csc_set(CscBuilder *m, OSQPInt row, OSQPInt col, OSQPFloat val) {
    (void)col;
    m->x[m->nz] = val;
    m->i[m->nz] = row;
    m->nz++;
}

void csc_col_done(CscBuilder *m) {
    m->p[m->n_cols + 1] = m->nz;
    m->n_cols++;
}

OSQPFloat vec_norm(const OSQPFloat *v, OSQPInt len) {
    OSQPFloat acc = 0.0f;
    for (OSQPInt i = 0; i < len; i++) acc += v[i] * v[i];
    return sqrtf(acc);
}

void mat_copy(OSQPFloat *dst, const OSQPFloat *src, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) dst[i] = src[i];
}

void mat_eye(OSQPFloat *m, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) m[i] = 0.0;
    for (OSQPInt i = 0; i < (rows < cols ? rows : cols); i++) m[i * cols + i] = 1.0;
}

void mat_mul(const OSQPFloat *A, const OSQPFloat *B, OSQPFloat *C, OSQPInt m, OSQPInt k, OSQPInt n) {
    for (OSQPInt i = 0; i < m; i++) {
        for (OSQPInt j = 0; j < n; j++) {
            OSQPFloat acc = 0.0f;
            for (OSQPInt t = 0; t < k; t++) acc += A[i * k + t] * B[t * n + j];
            C[i * n + j] = acc;
        }
    }
}

void mat_set(OSQPFloat *A, OSQPFloat v, OSQPInt rows, OSQPInt cols) {
    for (OSQPInt i = 0; i < rows * cols; i++) A[i] = v;
}

void crane_dynamics(const OSQPFloat x[NX], const OSQPFloat u[NU], OSQPFloat f_out[NX]) {
    const OSQPFloat g = 9.81f;
    f_out[0] = x[1];
    f_out[1] = u[0];
    f_out[2] = x[3];
    f_out[3] = u[1];
    f_out[4] = x[5];
    f_out[5] = -((g * sinf(x[4]) + cosf(x[4]) * u[0] + 2.0f * x[3] * x[5]) / x[2]);
}

void calc_h_gradient(const OSQPFloat x_current[NX], OSQPFloat jacobian_H[NH][NX], OSQPFloat h_val[NH]) {
    OSQPFloat p0 = 0.2, p1 = 1.25, p2 = 0.3;
    OSQPFloat x_pos = x_current[0] + x_current[2] * sin(x_current[4]);
    OSQPFloat temp = p0 * x_pos;

    h_val[0] = x_current[2] * cos(x_current[4]) - p0 * (x_pos * x_pos) - p1;
    h_val[1] = x_current[5] - p2;
    h_val[2] = -x_current[5] - p2;

    for (int j = 0; j < NH; j++) for (int i = 0; i < NX; i++) jacobian_H[j][i] = 0.0;

    jacobian_H[0][0] = -2.0 * temp;
    jacobian_H[0][2] = cos(x_current[4]) + temp * sin(x_current[4]);
    jacobian_H[0][4] = -x_current[2] * sin(x_current[4]) - 2.0 * temp * x_current[2] * cos(x_current[4]);
    jacobian_H[1][5] = 1.0;
    jacobian_H[2][5] = -1.0;
}

void linearization(OSQPFloat x_current[NX], OSQPFloat u_applied[NU], OSQPFloat Ad[NX][NX], OSQPFloat Bd[NX][NU], OSQPFloat d_lin[NX]) {
    OSQPFloat g = 9.81;
    OSQPFloat t1 = (2 * x_current[3] * x_current[5] + u_applied[0] * cos(x_current[4]) + g * sin(x_current[4])) / (x_current[2] * x_current[2]);
    OSQPFloat t2 = -(2 * x_current[5]) / x_current[2];
    OSQPFloat t3 = -(g * cos(x_current[4]) - u_applied[0] * sin(x_current[4])) / x_current[2];
    OSQPFloat t4 = -(2 * x_current[3]) / x_current[2];
    OSQPFloat t5 = -cos(x_current[4]) / x_current[2];

    OSQPFloat Ac[NX][NX] = {
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
        {0.0, 0.0, t1, t2, t3, t4}};

    OSQPFloat Bc[NX][NU] = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}, {t5, 0.0}};

    OSQPFloat fc[NX];
    crane_dynamics(x_current, u_applied, fc);

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
        d_lin[i] = DT * (fc[i] - Ac_x - Bc_u);
    }
}

void pend_dynamics(const OSQPFloat x[2], const OSQPFloat u[1], OSQPFloat f_out[2]) {
    const OSQPFloat g = 9.81f;
    const OSQPFloat l = 0.3f;
    const OSQPFloat m = 0.2f;
    const OSQPFloat b = 0.01f;
    f_out[0] = x[1];
    f_out[1] = -(g/l)*x[0]-(b/(m*l*l))*x[0]+(1.0/(m*l*l))*u[0];
}

void pend_calc_h_gradient(const OSQPFloat x_current[2], OSQPFloat *jacobian_H, OSQPFloat *h_val, int nh) {
    if (nh > 0) {
        h_val[0] = 0;
    }
    for(int j=0; j<nh; j++) {
        for(int i=0; i<2; i++) {
            jacobian_H[j * 2 + i] = 0.0;
        }
    }
}

void pend_linearization(OSQPFloat x_current[2], OSQPFloat u_applied[1], OSQPFloat Ad[2][2], OSQPFloat Bd[2][1], OSQPFloat d_lin[2], OSQPFloat dt) {
    const OSQPFloat g = 9.81f;
    const OSQPFloat l = 0.3f;
    const OSQPFloat m = 0.2f;
    const OSQPFloat b = 0.01f;
    OSQPFloat Ac[2][2] = {
        {0.0, 1.0}, 
        {-g/l, -(b/(m*l*l))}
    };

    OSQPFloat Bc[2][1] = {
        {0.0}, 
        {1.0/(m*l*l)}
    };
    
    OSQPFloat fc[2];
    pend_dynamics(x_current, u_applied, fc);

    for (int i = 0; i < 2; i++) {
        OSQPFloat Ac_x = 0.0, Bc_u = 0.0;
        for (int j = 0; j < 2; j++) {
            Ad[i][j] = (i == j ? 1.0 : 0.0) + Ac[i][j] * dt;
            Ac_x += Ac[i][j] * x_current[j];
        }
        for (int j = 0; j < 1; j++) {
            Bd[i][j] = Bc[i][j] * dt;
            Bc_u += Bc[i][j] * u_applied[j];
        }
        d_lin[i] = dt * (fc[i] - Ac_x - Bc_u);
    }
}
