#ifndef INVERTED_PENDULUM_CONFIG_H
#define INVERTED_PENDULUM_CONFIG_H

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

#define FORMULATION 0 /* <--- SET YOUR METHOD HERE */

#define CT 0 /* Set to 1 to enable explicit check_termination calls after each solve. Only applies to non-condensed formulation. */

/* ---------------------------------------------------------
 * Problem Constants & Parameters
 * --------------------------------------------------------- */
#define NX 2
#define NU 1
#define N 120
#define ND 1
#define NH 0
#define DT 0.01667
#define Tsim 11.0
#define R_BAND 2
#define MAXITER 10000
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
#define FORCE_LINEARIZATION_EMERGENCY_KAPPA 10000.0

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
