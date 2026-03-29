# Crane LTV OSQP (Sparse-Condensed + Primal CMoN Gatekeeper)

This README explains how `crane_LTV.c` is structured, what each block does, how it maps to the underlying MPC/CMoN theory, and what to change when porting the template to another system (e.g., a pendulum).

## High-level flow
1. **Phase 0 (pre-loop):**
   - Define dimensions/weights/limits, perform the first `linearization()`, `compute_sc_blocks()`, `build_P()`, `build_A()`, `build_bounds()`, and call `osqp_setup()`.
   - Cache arrays (`log_x`, `x_lin_pred`, SC blocks, etc.) are allocated once.
2. **Phase 1 (solve & apply):** Solve OSQP, extract virtual input `z0`, recover physical control `u = u_target + K(x - x_target) + z0`.
3. **Phase 2 (prediction vs. reality):** Predict the next state with cached Jacobians (`x_lin_pred = Ad x + Bd u + d_lin`), simulate true plant one step with Heun/RK2 using `crane_dynamics` to get `x_current`.
4. **Phase 3 (Primal CMoN gatekeeper):**
   - Compute κ = ‖x_actual − x_lin_pred‖ / (‖x_lin_pred − x_prev‖ + ε).
   - If κ > `eta_pri`, rebuild Jacobians, SC blocks, Hessian, and constraint matrix; else reuse cached matrices.
5. **Phase 4 (always update vectors):** Recompute steady-state residual `d_err`, propagate free response `e_free`, refresh gradient `q`, bounds `l/u`, update OSQP vectors, shift warm start, log, repeat.

## Module / block guide
- **Dimension & timing macros (`NX`, `NU`, `N`, `DT`, `R_BAND`, etc.)**: Define state/input sizes, horizon, sampling time, and condensation band. **Change these first** when moving to a new system.
- **`crane_dynamics()`**: Nonlinear continuous-time dynamics f(x,u). Used by Heun (RK2) for “ground truth” simulation. **Replace with your plant**.
- **`linearization()`**: Builds discrete `Ad`, `Bd`, and affine offset `d_lin` via Euler around `(x_current, u_applied)`. **Replace Jacobians and drift for new plants**.
- **`calc_h_gradient()`**: Nonlinear inequality constraints h(x) and Jacobian. **Edit or remove to match your constraints**; adjust `NH` accordingly.
- **`compute_sc_blocks()`**: Pre-stabilized LQR gain `K`, closed-loop `AK`, and banded blocks (`P_blocks`, `M_blocks`, `H_blocks`) for partial condensation. Uses `R_BAND` to truncate coupling. **If you change K or system matrices, adapt here.**
- **`build_P()`**: Forms Hessian/gradient in condensed space. `recompute_hessian=0` skips Hessian rebuild, only refreshes `q`—critical for the gatekeeper speed-up.
- **`build_A()`**: Sparse constraint matrix stacking state/input/nonlinear constraints per stage.
- **`build_bounds()`**: Lower/upper bounds for states, inputs, nonlinear constraints in condensed coordinates.
- **Main loop phases**: Implement the CMoN logic and OSQP data updates per the above flow.

## Mathematical linkage
- **Condensed MPC:** Decision variables are virtual inputs `z` (dimension `N*NU`). State/input coupling is captured via precomputed SC blocks derived from powers of `AK = Ad + Bd*K` with bandwidth `R_BAND`.
- **Cost:** ½∑(xᵀQx + uᵀRu) + ½ x_Nᵀ P x_N. Q, R, P are diagonal here (see `Q_diag`, `R_diag`, `P_diag`).
- **Constraints:**
  - State box: x_min ≤ x_k ≤ x_max (represented as offsets relative to the target using `e_free`).
  - Input box: u_min ≤ u_k ≤ u_max (around `u_target` + K e_k).
  - Nonlinear h(x) ≤ 0 style via `calc_h_gradient` (currently sway limits); currently relaxed in bounds (±∞) but Jacobian is prepared.
- **Primal CMoN κ:**
  - Numerator: linearization error ‖x_actual − x_lin_pred‖.
  - Denominator: curvature reference ‖x_lin_pred − x_prev‖ + ε.
  - Gate: if κ > `eta_pri` ⇒ rebuild Jacobians, Hessian, A; else reuse cached matrices and only update vectors.

## How to adapt to a new system
1. **Update dimensions & timing:** `NX`, `NU`, `N`, `DT`, `R_BAND`, `NH` (and any limits). Allocate arrays accordingly if sizes change.
2. **Replace dynamics:**
   - `crane_dynamics()` → your f(x,u).
   - `linearization()` → your Ac, Bc, and affine drift; set `Ad = I + Ac*DT`, `Bd = Bc*DT`, `d_lin = DT*(f - Ac x - Bc u)`.
3. **Constraints:** Edit `calc_h_gradient()` (and `NH`) to reflect your path constraints or remove if none.
4. **Cost weights:** Tune `Q_diag`, `R_diag`, `P_diag` per performance needs.
5. **Stabilizing gain:** Update `K_lqr` in `compute_sc_blocks()` (or compute online) to match the new linearization point.
6. **Targets/initials:** Set `x_current`, `x_target`, `u_target`, and bounds (`x_min/x_max`, `u_min/u_max`).
7. **CMoN threshold:** Tune `eta_pri` for how aggressively to rebuild sensitivities.

## Usage
- Build from `mpc_example/build` with CMake-generated files: `cmake --build .`
- Run: `./ltv` (logs to stdout and `../mpc_data.csv`).
- Adjust parameters, rebuild, rerun.

## Files
- `crane_LTV.c` — full MPC + CMoN gatekeeper implementation (partial condensation).
- `mpc_data.csv` — log output (step, states, inputs).

## Notes
- Hessian reuse is key for speed; ensure `recompute_hessian` flag is respected when modifying `build_P`.
- If you change `R_BAND`, ensure arrays sized with `R_BAND` stay consistent.
- If you remove nonlinear constraints, set `NH=0` and strip related code (or keep zeroed Jacobian/∞ bounds).
