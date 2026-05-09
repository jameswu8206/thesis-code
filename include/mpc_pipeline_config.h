#ifndef MPC_PIPELINE_CONFIG_H
#define MPC_PIPELINE_CONFIG_H

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

#define FORMULATION 2 /* <--- SET YOUR METHOD HERE */

/* ---------------------------------------------------------
 * Problem Constants & Parameters
 * --------------------------------------------------------- */
#define NX 6
#define NU 2
#define N 50
#define ND 6
#define NH 3
#define DT 0.0400
#define Tsim 25.0
#define R_BAND N
#define MAXITER 50
#define TERMINATION_TOLERANCE 1e-4

/* ---------------------------------------------------------
 * CMoN Linearization Cadence (User tunable)
 * ---------------------------------------------------------
 * FORCE_LINEARIZATION_BURST:
 *   Minimum consecutive linearizations once triggered.
 * FORCE_LINEARIZATION_COOLDOWN:
 *   Number of iterations to wait before allowing retrigger.
 * FORCE_LINEARIZATION_EMERGENCY_KAPPA:
 *   If kappa exceeds this value, cooldown is overridden and
 *   linearization is forced immediately.
 *
 * Example: set either value to 3, 5, or any positive integer.
 * --------------------------------------------------------- */
#define FORCE_LINEARIZATION_BURST 0
#define FORCE_LINEARIZATION_COOLDOWN 0
#define FORCE_LINEARIZATION_EMERGENCY_KAPPA 10.0

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

#endif
