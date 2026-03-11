#!/usr/bin/env python3
"""
Diagnose the factor-of-2 discrepancy between RustFoil gamma and XFOIL Ue.

This script:
1. Loads XFOIL's coordinate file
2. Computes influence coefficients using XFOIL's exact formulas
3. Solves for GAMU base solutions
4. Compares with RustFoil's values at matching x-positions
"""

import numpy as np
import json
import sys

def load_xfoil_coords(filename):
    """Load coordinates from XFOIL .dat file."""
    coords = []
    with open(filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 0 and not line[0].isdigit():
                continue  # Skip header
            parts = line.split()
            if len(parts) >= 2:
                coords.append((float(parts[0]), float(parts[1])))
    return np.array(coords)

def compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2):
    """XFOIL's PSIS/PSID formulas."""
    psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2)
    dxinv = 1.0 / (x1 - x2) if abs(x1 - x2) > 1e-20 else 0.0
    psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1**2 - x2**2)) * dxinv
    return psis, psid

def compute_dzdg(x, y):
    """Build influence coefficient matrix (XFOIL's GGCALC)."""
    n = len(x)
    QOPI = 0.25 / np.pi
    dzdg = np.zeros((n, n))
    
    for i in range(n):
        xi, yi = x[i], y[i]
        for jo in range(n):
            jp = (jo + 1) % n
            if jo == n - 1:  # Skip TE panel
                continue
            dx, dy = x[jp] - x[jo], y[jp] - y[jo]
            ds_sq = dx**2 + dy**2
            if ds_sq < 1e-24:
                continue
            dso = np.sqrt(ds_sq)
            sx, sy = dx / dso, dy / dso
            
            rx1, ry1 = xi - x[jo], yi - y[jo]
            rx2, ry2 = xi - x[jp], yi - y[jp]
            
            x1 = sx * rx1 + sy * ry1
            x2 = sx * rx2 + sy * ry2
            yy = sx * ry1 - sy * rx1
            
            rs1, rs2 = rx1**2 + ry1**2, rx2**2 + ry2**2
            
            g1 = np.log(rs1) if i != jo and rs1 > 1e-20 else 0.0
            t1 = np.arctan2(x1, yy) if i != jo and rs1 > 1e-20 else 0.0
            g2 = np.log(rs2) if i != jp and rs2 > 1e-20 else 0.0
            t2 = np.arctan2(x2, yy) if i != jp and rs2 > 1e-20 else 0.0
            
            psis, psid = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2)
            dzdg[i, jo] += QOPI * (psis - psid)
            dzdg[i, jp] += QOPI * (psis + psid)
    
    return dzdg

def solve_gamu(x, y, dzdg):
    """Solve for GAMU base solutions (alpha=0, alpha=90)."""
    n = len(x)
    A = np.zeros((n+1, n+1))
    rhs_0 = np.zeros(n+1)
    rhs_90 = np.zeros(n+1)
    
    for i in range(n):
        for j in range(n):
            A[i, j] = dzdg[i, j]
        A[i, n] = -1.0
        rhs_0[i] = -y[i]
        rhs_90[i] = x[i]
    
    A[n, 0] = 1.0
    A[n, n-1] = 1.0
    
    gamu_0 = np.linalg.solve(A, rhs_0)[:n]
    gamu_90 = np.linalg.solve(A, rhs_90)[:n]
    
    return gamu_0, gamu_90

def compute_gamma(gamu_0, gamu_90, alpha_deg):
    """Compute gamma at given alpha."""
    alpha = np.radians(alpha_deg)
    return np.cos(alpha) * gamu_0 + np.sin(alpha) * gamu_90

def find_stagnation(gamma):
    """Find stagnation point index."""
    for i in range(len(gamma) - 1):
        if gamma[i] >= 0 and gamma[i+1] < 0:
            return i
    return len(gamma) // 2

def main():
    # Load coordinates
    coords = load_xfoil_coords('testdata/naca0012.dat')
    x, y = coords[:, 0], coords[:, 1]
    n = len(x)
    print(f"Loaded {n} coordinates")
    print(f"First: ({x[0]:.6f}, {y[0]:.6f})")
    print(f"Last:  ({x[-1]:.6f}, {y[-1]:.6f})")
    
    # Compute influence coefficients
    print("\nComputing influence coefficients...")
    dzdg = compute_dzdg(x, y)
    
    # Solve for GAMU
    print("Solving for base solutions...")
    gamu_0, gamu_90 = solve_gamu(x, y, dzdg)
    
    # Compute gamma at alpha=4
    alpha_deg = 4.0
    gamma = compute_gamma(gamu_0, gamu_90, alpha_deg)
    
    # Find stagnation
    ist = find_stagnation(gamma)
    print(f"\nStagnation at panel {ist}: x={x[ist]:.6f}, y={y[ist]:.6f}, gamma={gamma[ist]:.6f}")
    
    # Print values near stagnation (IBL = 1, 2, 3, ...)
    print(f"\n=== Upper surface from stagnation (XFOIL IBL convention) ===")
    print(f"{'ibl':>4} {'panel':>6} {'x':>12} {'y':>12} {'|gamma|':>12}")
    for ibl, panel in enumerate(range(ist, max(ist-15, -1), -1)):
        if panel >= 0:
            print(f"{ibl+1:4} {panel:6} {x[panel]:12.6f} {y[panel]:12.6f} {abs(gamma[panel]):12.6f}")
    
    # Compare with XFOIL reference
    print("\n=== XFOIL Reference (from viscous_ref.json) ===")
    print("IBL  x          Ue(XFOIL)")
    xfoil_data = [
        (2, 0.0011044, 0.060676),
        (3, 0.0033612, 0.207621),
        (4, 0.0054999, 0.367744),
        (5, 0.0075362, 0.535957),
    ]
    for ibl, x_ref, ue_ref in xfoil_data:
        print(f"{ibl:3}  {x_ref:.6f}  {ue_ref:.6f}")
    
    # Find ratio at IBL=2 position
    print("\n=== Comparison at x=0.0011 ===")
    x_target = 0.0011
    # Interpolate
    for i in range(ist-1, 0, -1):
        if x[i] <= x_target <= x[i-1] or x[i] >= x_target >= x[i-1]:
            t = (x_target - x[i]) / (x[i-1] - x[i]) if abs(x[i-1] - x[i]) > 1e-12 else 0.5
            gamma_interp = abs(gamma[i]) + t * (abs(gamma[i-1]) - abs(gamma[i]))
            print(f"Python |gamma| at x={x_target}: {gamma_interp:.6f}")
            print(f"XFOIL Ue at x={x_target}:       0.060676")
            print(f"Ratio: {gamma_interp / 0.060676:.2f}")
            break

if __name__ == "__main__":
    main()
