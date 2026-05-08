import numpy as np
import scipy.linalg

# Exact Analytical Ad and Bd for 2-DOF Helicopter at Hover (Ts = 0.05s)
Ad = np.array([
    [1.0, 0.0, 0.050, 0.0],
    [0.0, 1.0, 0.0,   0.050],
    [0.0, 0.0, 0.537, 0.0],
    [0.0, 0.0, 0.0,   0.826]
])

Bd = np.array([
    [0.0,     0.0],
    [0.0,     0.0],
    [0.118,   0.004],
    [0.012,   0.039]
])

# OSQP Penalty Weights (Adjust R to tune the size of K)
Q = np.diag([1000, 500, 0.1, 0.1])
R = np.diag([0.001, 0.01])   # <-- INCREASE THESE if OSQP gives "Non Convex" errors

# Solve Discrete Algebraic Riccati Equation
P = scipy.linalg.solve_discrete_are(Ad, Bd, Q, R)

# Compute Discrete LQR Gain
K = -np.linalg.inv(R + Bd.T @ P @ Bd) @ (Bd.T @ P @ Ad)

print("K_lqr_heli[2][4] = {")
print(f"    {{{K[0,0]:.2f}, {K[0,1]:.2f}, {K[0,2]:.2f}, {K[0,3]:.2f}}},")
print(f"    {{{K[1,0]:.2f}, {K[1,1]:.2f}, {K[1,2]:.2f}, {K[1,3]:.2f}}}")
print("};")

g, l, m, b = 9.81, 0.3, 0.2, 0.01
DT = 0.02  # Make sure this matches #define DT in your C code!

# --- Continuous Dynamics Jacobian (Ac, Bc) ---
# Extracted exactly from your C code linear physics: -(g/l)*x[0] - (b/(m*l*l))*x[1]
Ac = np.array([
    [0.0,  1.0],
    [-g/l, -b / (m * l**2)]
])

Bc = np.array([
    [0.0],
    [1.0 / (m * l**2)]
])

# --- Forward Euler Discretization ---
# Matches: Ad = I + Ac*dt, Bd = Bc*dt
Ad = np.eye(2) + Ac * DT
Bd = Bc * DT

# --- Gentle LQR Tuning ---
# We heavily penalize R to force the K values to be small, keeping OSQP's Hessian stable.
Q_tune = np.diag([1.0, 1.0])
R_tune = np.diag([100.0])

# --- Solve Discrete LQR ---
P = scipy.linalg.solve_discrete_are(Ad, Bd, Q_tune, R_tune)
K = -np.linalg.inv(R_tune + Bd.T @ P @ Bd) @ (Bd.T @ P @ Ad)

print(f"OSQPFloat K_lqr_pend[NU][NX] = {{{{ {K[0,0]:.4f}, {K[0,1]:.4f} }}}};")

g = 9.81
DT = 0.01          # Make sure this matches #define DT in your C code!
L_target = 2.0     # Target rope length (x[2])

# --- Continuous Dynamics Jacobian (Ac, Bc) ---
# Linearized around hover: x[2] = 2.0, x[4] = 0.0, u[0] = 0.0
Ac = np.zeros((6, 6))
Ac[0, 1] = 1.0
Ac[2, 3] = 1.0
Ac[4, 5] = 1.0
Ac[5, 4] = -g / L_target

Bc = np.zeros((6, 2))
Bc[1, 0] = 1.0
Bc[3, 1] = 1.0
Bc[5, 0] = -1.0 / L_target

# --- Forward Euler Discretization ---
# Matches: Ad = I + Ac*dt, Bd = Bc*dt
Ad = np.eye(6) + Ac * DT
Bd = Bc * DT

# --- Gentle LQR Tuning ---
# High R penalties keep the feedback gains low for ADMM stability
Q_tune = np.diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
R_tune = np.diag([100.0, 100.0])

# --- Solve Discrete LQR ---
P = scipy.linalg.solve_discrete_are(Ad, Bd, Q_tune, R_tune)
K = -np.linalg.inv(R_tune + Bd.T @ P @ Bd) @ (Bd.T @ P @ Ad)

print("OSQPFloat K_lqr_crane[NU][NX] = {")
print(f"    {{{K[0,0]:.4f}, {K[0,1]:.4f}, {K[0,2]:.4f}, {K[0,3]:.4f}, {K[0,4]:.4f}, {K[0,5]:.4f}}},")
print(f"    {{{K[1,0]:.4f}, {K[1,1]:.4f}, {K[1,2]:.4f}, {K[1,3]:.4f}, {K[1,4]:.4f}, {K[1,5]:.4f}}}")
print("};")