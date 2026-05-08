import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Load and clean the data
df = pd.read_csv('mpc_data.csv')
df = df.drop_duplicates(subset=['step']).reset_index(drop=True)

# ==========================================
# PLOT: Dual-Density Physical Strobe Plot
# ==========================================
fig, ax = plt.subplots(figsize=(12, 6))
ax.axhline(0, color='black', linewidth=2, label='Track')

total_steps = len(df)

# --- DUAL DENSITY SAMPLING LOGIC ---
# Define the split point (e.g., 50% of the way through the trajectory)
split_point = int(total_steps * 0.15) # Adjust 0.4 to capture the exact end of the heavy swing

# Define how dense you want each section
# For example, ~50 frames in the transient phase, ~10 frames in the settling phase
dense_step = max(1, split_point // 30)  
sparse_step = max(1, (total_steps - split_point) // 40)

# Create a combined list of indices to plot
indices_first_half = list(range(0, split_point, dense_step))
indices_second_half = list(range(split_point, total_steps, sparse_step))
plot_indices = indices_first_half + indices_second_half

# -----------------------------------

for i in plot_indices:
    p = df['x0'].iloc[i]
    angle = df['x4'].iloc[i]
    L = df['x2'].iloc[i] 
    
    # Calculate payload X, Y coordinates
    x_payload = p + L * np.sin(angle)
    y_payload = -L * np.cos(angle)
    
    # Fading gradient based on total progress
    intensity = i / total_steps
    
    # Keep the lines slightly transparent so the dense area doesn't turn solid black
    color_val = str(0.85 - 0.7 * intensity)
    line_width = 1.0 if i < split_point else 1.5 # Thinner lines in dense area
    
    ax.plot([p, x_payload], [0, y_payload], color=color_val, linewidth=line_width)
    ax.plot(p, 0, marker='s', color=color_val, markersize=5 if i < split_point else 7)

# Plot the Final State boldly
final_p = df['x0'].iloc[-1]
final_theta = df['x4'].iloc[-1]
final_L = df['x2'].iloc[-1]
final_x_payload = final_p + final_L * np.sin(final_theta)
final_y_payload = -final_L * np.cos(final_theta)

ax.plot([final_p, final_x_payload], [0, final_y_payload], color='red', linewidth=2.5, label='Final State')
ax.plot(final_p, 0, marker='s', color='red', markersize=9)

ax.set_aspect('equal') # CRITICAL for accurate physical angles
ax.set_xlabel('Horizontal Position [m]')
ax.set_ylabel('Vertical Position [m]')
ax.set_title('Overhead Crane Trajectory: Transient vs. Steady-State')
ax.grid(True, linestyle=':')
ax.legend(loc='lower right')

# plt.savefig('crane_strobe_dual_density.pdf', dpi=300)
plt.show()