import numpy as np
import scipy.linalg

# 1. Define physical constants
DT = 0.01

# 2. Define Discretized System Matrices (Ad, Bd)
Ad = np.array([
    [1.0, 0.0, DT,  0.0, 0.0, 0.0], 
    [0.0, 1.0, 0.0, DT,  0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [DT,  0.0, 0.0, 0.0, 1.0, 0.0], 
    [0.0, DT,  0.0, 0.0, 0.0, 1.0]
])

Bd = np.array([
    [0.0,     0.0], 
    [0.01*DT, -0.01*DT],
    [0.19*DT, 0.19*DT], 
    [1.32*DT, -1.32*DT],
    [0.0,     0.0],
    [0.0,     0.0]
])

# 3. Define Running Costs
Q = np.diag([100.0, 100.0, 10.0, 10.0, 400.0, 200.0])
R = np.diag([0.44, 0.6])

# 4. Solve Discrete Algebraic Riccati Equation (DARE)
P_dense = scipy.linalg.solve_discrete_are(Ad, Bd, Q, R)

print("--- P_diag (If you refuse to update your C code) ---")
print(np.diag(P_dense))

print("\n--- P_dense (The mathematically correct matrix) ---")
np.set_printoptions(suppress=True, precision=4, linewidth=100)
print(P_dense)