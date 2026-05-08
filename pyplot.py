import pandas as pd
import matplotlib.pyplot as plt
import math

# ==========================================
# USER CONFIGURATION: Set your parameters here
# ==========================================
DT = 0.01  # <--- Change your time step (dt) here
X_TARGET = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0]  # <--- Change your target states here
U_TARGET = [0.0, 0.0]                       # <--- Change your target inputs here
# ==========================================

# 1. Read the Generic CSV
try:
    # This reads the header automatically (step, x0, x1, ... u0, u1...)
    data = pd.read_csv('mpc_data.csv')
    data['time'] = data['step'] * DT
except FileNotFoundError:
    print("Error: 'mpc_data.csv' not found. Run the C simulation first!")
    exit()

# 2. Identify Columns dynamically
# Find all columns that start with 'x' (States) and 'u' (Controls)
x_cols = [col for col in data.columns if col.startswith('x')]
u_cols = [col for col in data.columns if col.startswith('u')]

print(f"Found {len(x_cols)} States: {x_cols}")
print(f"Found {len(u_cols)} Controls: {u_cols}")
print(f"Total time: {max(data['time']):.4f}s")

# --- Steady-State Error Calculation ---
try:
    # Calculate Euclidean norm of the state error at the final step
    last_row = data.iloc[-1]
    
    # Check if dimensions match
    if len(x_cols) != len(X_TARGET):
        print(f"Warning: CSV has {len(x_cols)} states, but {len(X_TARGET)} targets provided.")
        
    e_ss_sq = sum((last_row[col] - xt)**2 for col, xt in zip(x_cols, X_TARGET))
    e_ss_norm = math.sqrt(e_ss_sq)
    
    print("\n--- Steady-State Analysis ---")
    print(f"Final States: {[round(last_row[c], 6) for c in x_cols]}")
    print(f"Target States: {X_TARGET}")
    print(f"State Steady-State Error (Euclidean Norm): {e_ss_norm:.2e}")
    
except Exception as e:
    print(f"Error during steady-state calculation: {e}")
    e_ss_norm = None
# -------------------------------------------------------------

# 3. Setup the Plot
# We creates 2 main subplots: Top for States, Bottom for Controls
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

# --- Plot 1: All States (x0, x1, ... xN) ---
for col in x_cols:
    # Use a solid line for states
    ax1.plot(data['time'], data[col], label=f'State {col}', linewidth=2)

ax1.set_ylabel('State Value')
ax1.set_title('MPC System States')
ax1.legend(loc='upper right', bbox_to_anchor=(1.15, 1)) # Move legend outside if crowded
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.axhline(y=0.0, color='black', linewidth=0.5) # Zero line reference

# Add Steady-State Error text box to the state plot
if e_ss_norm is not None:
    textstr = f"Steady-State Error (Norm): {e_ss_norm:.4e}"
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
    ax1.text(0.02, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
             verticalalignment='top', bbox=props)

# --- Plot 2: All Controls (u0, u1, ... uN) ---
for col in u_cols:
    # Use 'step' plot for controls (Digital truth!)
    ax2.step(data['time'], data[col], label=f'Control {col}', where='post', linewidth=2)

ax2.set_ylabel('Control Input')
ax2.set_xlabel('Time (s)')
ax2.set_title('Control Inputs')
ax2.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax2.grid(True, linestyle='--', alpha=0.6)

# 4. Show the plot
plt.tight_layout()
plt.show()