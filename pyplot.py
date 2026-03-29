import pandas as pd
import matplotlib.pyplot as plt

# 1. Read the Generic CSV
try:
    # This reads the header automatically (step, x0, x1, ... u0, u1...)
    data = pd.read_csv('mpc_data.csv')
except FileNotFoundError:
    print("Error: 'mpc_data.csv' not found. Run the C simulation first!")
    exit()

# 2. Identify Columns dynamically
# Find all columns that start with 'x' (States) and 'u' (Controls)
x_cols = [col for col in data.columns if col.startswith('x')]
u_cols = [col for col in data.columns if col.startswith('u')]

print(f"Found {len(x_cols)} States: {x_cols}")
print(f"Found {len(u_cols)} Controls: {u_cols}")

# 3. Setup the Plot
# We creates 2 main subplots: Top for States, Bottom for Controls
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

# --- Plot 1: All States (x0, x1, ... xN) ---
for col in x_cols:
    # Use a solid line for states
    ax1.plot(data['step'], data[col], label=f'State {col}', linewidth=2)

ax1.set_ylabel('State Value')
ax1.set_title('MPC System States')
ax1.legend(loc='upper right', bbox_to_anchor=(1.15, 1)) # Move legend outside if crowded
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.axhline(y=0.0, color='black', linewidth=0.5) # Zero line reference

# --- Plot 2: All Controls (u0, u1, ... uN) ---
for col in u_cols:
    # Use 'step' plot for controls (Digital truth!)
    ax2.step(data['step'], data[col], label=f'Control {col}', where='post', linewidth=2)

ax2.set_ylabel('Control Input')
ax2.set_xlabel('Simulation Step (k)')
ax2.set_title('Control Inputs')
ax2.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax2.grid(True, linestyle='--', alpha=0.6)

# 4. Show the plot
plt.tight_layout()
plt.show()