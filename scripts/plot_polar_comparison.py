#!/usr/bin/env python3
"""Plot RustFoil vs XFOIL polar comparison."""

import matplotlib.pyplot as plt
import numpy as np

# Data from test_full_polar_comparison (latest run)
# (alpha, xfoil_cl, rustfoil_cl, xfoil_cd, rustfoil_cd)
data = [
    (-4.0, -0.4278, -0.4710, 0.007279, 0.005531),
    (-3.0, -0.3199, -0.3434, 0.006390, 0.005007),
    (-2.0, -0.2142, -0.2353, 0.005802, 0.004416),
    (-1.0, -0.1074, -0.1063, 0.005487, 0.005613),
    ( 0.0,  0.0000,  0.0176, 0.005404, 0.004440),
    ( 1.0,  0.1074,  0.1414, 0.005487, 0.004766),
    ( 2.0,  0.2142,  0.2365, 0.005802, 0.008474),
    ( 3.0,  0.3200,  0.3462, 0.006390, 0.005427),
    ( 4.0,  0.4278,  0.4826, 0.007279, 0.007943),
    ( 5.0,  0.5580,  0.5915, 0.008477, 0.009613),
    ( 6.0,  0.6948,  0.7043, 0.009726, 0.010194),
    ( 7.0,  0.8264,  0.8120, 0.010939, 0.008281),
    ( 8.0,  0.9099,  0.9324, 0.012110, 0.009156),
    ( 9.0,  0.9948,  1.0605, 0.013406, 0.009260),
    (10.0,  1.0809,  1.1287, 0.014977, 0.012993),
]

alpha = np.array([d[0] for d in data])
xfoil_cl = np.array([d[1] for d in data])
rustfoil_cl = np.array([d[2] for d in data])
xfoil_cd = np.array([d[3] for d in data])
rustfoil_cd = np.array([d[4] for d in data])

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# CL vs Alpha
ax1.plot(alpha, xfoil_cl, 'b-o', label='XFOIL', linewidth=2, markersize=6)
ax1.plot(alpha, rustfoil_cl, 'r--s', label='RustFoil', linewidth=2, markersize=6)
ax1.set_xlabel('Angle of Attack (deg)', fontsize=12)
ax1.set_ylabel('Lift Coefficient CL', fontsize=12)
ax1.set_title('CL vs Alpha - NACA 0012, Re=1M', fontsize=14)
ax1.legend(loc='lower right', fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-5, 11)

# CD vs Alpha
ax2.plot(alpha, xfoil_cd * 1000, 'b-o', label='XFOIL', linewidth=2, markersize=6)
ax2.plot(alpha, rustfoil_cd * 1000, 'r--s', label='RustFoil', linewidth=2, markersize=6)
ax2.set_xlabel('Angle of Attack (deg)', fontsize=12)
ax2.set_ylabel('Drag Coefficient CD (counts)', fontsize=12)
ax2.set_title('CD vs Alpha - NACA 0012, Re=1M', fontsize=14)
ax2.legend(loc='upper left', fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(-5, 11)

plt.tight_layout()
plt.savefig('scripts/polar_comparison.png', dpi=150, bbox_inches='tight')
# plt.show()  # Don't show in non-interactive mode

print("Plot saved to scripts/polar_comparison.png")
