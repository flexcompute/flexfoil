#!/usr/bin/env python3
"""
Verify paneling parity between XFOIL and RustFoil.

This script ensures that both solvers use identical panel geometries
before proceeding to viscous comparisons.

Usage:
    python scripts/verify_paneling_parity.py [--tolerance 1e-6]
"""

import argparse
import sys
from pathlib import Path
import numpy as np


def load_dat_file(filepath: Path) -> tuple[str, np.ndarray, np.ndarray]:
    """
    Load airfoil coordinates from .dat file.
    
    Returns (name, x_array, y_array).
    Handles both XFOIL format (no header or name header) and RustFoil format.
    """
    x_coords = []
    y_coords = []
    name = filepath.stem
    
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    x_coords.append(x)
                    y_coords.append(y)
                except ValueError:
                    # Header line (name)
                    if i == 0:
                        name = line
    
    return name, np.array(x_coords), np.array(y_coords)


def compare_coordinates(
    x1: np.ndarray, y1: np.ndarray,
    x2: np.ndarray, y2: np.ndarray,
    tolerance: float
) -> dict:
    """
    Compare two sets of coordinates.
    
    Returns dict with comparison metrics and pass/fail status.
    """
    if len(x1) != len(x2):
        return {
            'pass': False,
            'error': f'Point count mismatch: {len(x1)} vs {len(x2)}',
            'n_points_1': len(x1),
            'n_points_2': len(x2),
        }
    
    n = len(x1)
    dx = x1 - x2
    dy = y1 - y2
    dist = np.sqrt(dx**2 + dy**2)
    
    max_dx = np.max(np.abs(dx))
    max_dy = np.max(np.abs(dy))
    max_dist = np.max(dist)
    rms_dist = np.sqrt(np.mean(dist**2))
    
    # Find worst point
    worst_idx = np.argmax(dist)
    
    passed = max_dist < tolerance
    
    return {
        'pass': passed,
        'n_points': n,
        'max_dx': max_dx,
        'max_dy': max_dy,
        'max_dist': max_dist,
        'rms_dist': rms_dist,
        'worst_idx': worst_idx,
        'worst_x1': x1[worst_idx],
        'worst_y1': y1[worst_idx],
        'worst_x2': x2[worst_idx],
        'worst_y2': y2[worst_idx],
        'tolerance': tolerance,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Verify paneling parity between XFOIL and RustFoil"
    )
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=1e-5,
        help="Maximum allowed coordinate difference (default: 1e-5)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-point comparison"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    
    # Foils to verify
    foils = ["naca0012", "naca2412", "naca4412"]
    
    print("=" * 70)
    print(" PANELING PARITY VERIFICATION")
    print("=" * 70)
    print(f"  Tolerance: {args.tolerance:.2e}")
    print()
    
    all_passed = True
    results = {}
    
    for foil in foils:
        print(f"\n{'='*50}")
        print(f" {foil.upper()}")
        print("=" * 50)
        
        # Look for XFOIL paneled file
        xfoil_path = workspace / f"{foil}_xfoil_paneled.dat"
        if not xfoil_path.exists():
            xfoil_path = workspace / f"{foil}_paneled_xfoil.dat"
        
        # Look for RustFoil paneled file
        rustfoil_path = workspace / f"{foil}_rustfoil_paneled.dat"
        if not rustfoil_path.exists():
            rustfoil_path = workspace / f"{foil}_paneled_rustfoil.dat"
        
        if not xfoil_path.exists():
            print(f"  WARNING: XFOIL paneled file not found: {xfoil_path}")
            all_passed = False
            continue
        
        if not rustfoil_path.exists():
            print(f"  WARNING: RustFoil paneled file not found: {rustfoil_path}")
            all_passed = False
            continue
        
        print(f"  XFOIL file:    {xfoil_path.name}")
        print(f"  RustFoil file: {rustfoil_path.name}")
        
        # Load coordinates
        xf_name, xf_x, xf_y = load_dat_file(xfoil_path)
        rf_name, rf_x, rf_y = load_dat_file(rustfoil_path)
        
        print(f"  XFOIL points:    {len(xf_x)}")
        print(f"  RustFoil points: {len(rf_x)}")
        
        # Compare
        result = compare_coordinates(xf_x, xf_y, rf_x, rf_y, args.tolerance)
        results[foil] = result
        
        if 'error' in result:
            print(f"\n  ERROR: {result['error']}")
            all_passed = False
            continue
        
        status = "PASS" if result['pass'] else "FAIL"
        print(f"\n  Status: {status}")
        print(f"  Max X diff:  {result['max_dx']:.2e}")
        print(f"  Max Y diff:  {result['max_dy']:.2e}")
        print(f"  Max dist:    {result['max_dist']:.2e}")
        print(f"  RMS dist:    {result['rms_dist']:.2e}")
        print(f"  Worst point: idx={result['worst_idx']}")
        print(f"    XFOIL:    ({result['worst_x1']:.8f}, {result['worst_y1']:.8f})")
        print(f"    RustFoil: ({result['worst_x2']:.8f}, {result['worst_y2']:.8f})")
        
        if not result['pass']:
            all_passed = False
        
        if args.verbose and result['n_points'] <= 200:
            print(f"\n  Per-point comparison (first 20):")
            print(f"  {'Idx':>4} | {'X_XF':>12} | {'X_RF':>12} | {'dX':>10} | {'Y_XF':>12} | {'Y_RF':>12} | {'dY':>10}")
            print("  " + "-" * 85)
            for i in range(min(20, result['n_points'])):
                dx = abs(xf_x[i] - rf_x[i])
                dy = abs(xf_y[i] - rf_y[i])
                flag = " *" if dx > args.tolerance or dy > args.tolerance else ""
                print(f"  {i:>4} | {xf_x[i]:>12.8f} | {rf_x[i]:>12.8f} | {dx:>10.2e} | "
                      f"{xf_y[i]:>12.8f} | {rf_y[i]:>12.8f} | {dy:>10.2e}{flag}")
    
    # Summary
    print(f"\n{'='*70}")
    print(" SUMMARY")
    print("=" * 70)
    
    print(f"{'Foil':<12} | {'Points':>6} | {'Max Dist':>12} | {'RMS Dist':>12} | {'Status':<6}")
    print("-" * 60)
    
    for foil, result in results.items():
        if 'error' in result:
            print(f"{foil:<12} | ERROR: {result['error']}")
        else:
            status = "PASS" if result['pass'] else "FAIL"
            print(f"{foil:<12} | {result['n_points']:>6} | {result['max_dist']:>12.2e} | "
                  f"{result['rms_dist']:>12.2e} | {status:<6}")
    
    print()
    if all_passed:
        print("OVERALL: PASS - All foils have matching paneling")
        print("         Inviscid results can be compared directly.")
        return 0
    else:
        print("OVERALL: FAIL - Paneling mismatch detected")
        print("         Fix paneling before comparing viscous results.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
