#!/usr/bin/env python3
"""
Compare XFOIL and RustFoil Jacobian matrices (BLDIF output).

This script analyzes VS1, VS2 matrices from XFOIL's instrumented output
and checks for structural issues that could cause Newton solver problems.

Usage:
    python compare_jacobians.py xfoil_debug.json
    python compare_jacobians.py xfoil_debug.json --station 10
"""

import json
import argparse
import numpy as np


def load_bldif_events(json_path, iteration=None, side=None, station=None):
    """Load BLDIF events from XFOIL debug JSON."""
    with open(json_path) as f:
        data = json.load(f)
    
    events = [e for e in data['events'] if e['subroutine'] == 'BLDIF']
    
    if iteration is not None:
        events = [e for e in events if e.get('iteration') == iteration]
    if side is not None:
        events = [e for e in events if e.get('side') == side]
    if station is not None:
        events = [e for e in events if e.get('ibl') == station]
    
    return events


def analyze_jacobian(vs1, vs2, vsrez, label=""):
    """Analyze a single Jacobian for structural issues."""
    vs1 = np.array(vs1)
    vs2 = np.array(vs2)
    vsrez = np.array(vsrez)
    
    issues = []
    
    # Check for all-zero columns (should never happen in well-posed system)
    for col in range(vs1.shape[1]):
        if np.allclose(vs1[:3, col], 0) and np.allclose(vs2[:3, col], 0):
            issues.append(f"Column {col} is all zeros in both VS1 and VS2")
    
    # Check VS1 column 2 (δ* derivatives) - was previously missing in RustFoil
    if vs1.shape[1] > 2:
        col2_vs1 = vs1[:3, 2]
        if np.allclose(col2_vs1, 0, atol=1e-20):
            issues.append("VS1 column 2 (δ* derivatives) is all zeros!")
    
    # Check 3x3 VA block determinant (singularity check)
    if vs1.shape[0] >= 3 and vs1.shape[1] >= 3:
        va = vs1[:3, :3]
        try:
            det = np.linalg.det(va)
            cond = np.linalg.cond(va)
            if abs(det) < 1e-10:
                issues.append(f"Near-singular VA block! det={det:.2e}")
            if cond > 1e10:
                issues.append(f"Ill-conditioned VA block! cond={cond:.2e}")
        except:
            issues.append("Could not compute VA determinant/condition")
    
    # Check residual magnitudes
    max_res = np.max(np.abs(vsrez[:3]))
    if max_res > 1e3:
        issues.append(f"Large residual magnitude: {max_res:.2e}")
    
    return issues


def compare_jacobians_xfoil(events, max_display=10):
    """Analyze XFOIL Jacobians for common issues."""
    print("\n" + "=" * 70)
    print("XFOIL Jacobian Analysis")
    print("=" * 70)
    
    total = len(events)
    issues_count = 0
    
    for i, ev in enumerate(events[:max_display]):
        label = f"Side {ev.get('side', '?')}, Station {ev.get('ibl', '?')}, Type {ev.get('flow_type', '?')}"
        issues = analyze_jacobian(ev['VS1'], ev['VS2'], ev['VSREZ'], label)
        
        if issues:
            issues_count += 1
            print(f"\n{label}:")
            for issue in issues:
                print(f"  WARNING: {issue}")
        
        vs1 = np.array(ev['VS1'])
        vs2 = np.array(ev['VS2'])
        
        if i < 3:  # Show first 3 in detail
            print(f"\n{label}:")
            print(f"  VS1[0,:] (row 0): {vs1[0, :]}")
            print(f"  VS1[1,:] (row 1): {vs1[1, :]}")
            print(f"  VS1[2,:] (row 2): {vs1[2, :]}")
            det = np.linalg.det(vs1[:3, :3])
            print(f"  VA block det: {det:.4e}")
    
    # Summary statistics
    print(f"\n--- Summary ---")
    print(f"Total BLDIF events: {total}")
    print(f"Events with issues: {issues_count}")
    
    # Check convergence of VA determinants
    dets = []
    for ev in events:
        vs1 = np.array(ev['VS1'])
        if vs1.shape[0] >= 3 and vs1.shape[1] >= 3:
            dets.append(np.linalg.det(vs1[:3, :3]))
    
    if dets:
        dets = np.array(dets)
        print(f"VA determinant range: [{dets.min():.2e}, {dets.max():.2e}]")
        print(f"VA determinant mean: {dets.mean():.2e}")
        
        # Check for sign changes (could indicate issues)
        sign_changes = np.sum(np.diff(np.sign(dets)) != 0)
        if sign_changes > 0:
            print(f"WARNING: {sign_changes} sign changes in VA determinant!")


def check_jacobian_structure(events):
    """Check that Jacobian has expected non-zero structure."""
    print("\n" + "=" * 70)
    print("Jacobian Structure Analysis")
    print("=" * 70)
    
    # Expected non-zero pattern for laminar flow:
    # VS1[0,:] = [ampl1, theta1, delta_star1, u1, x1]
    # VS1[1,:] = momentum derivatives
    # VS1[2,:] = shape parameter derivatives
    
    for ev in events[:5]:
        if ev.get('flow_type') == 1:  # Laminar
            vs1 = np.array(ev['VS1'])
            vs2 = np.array(ev['VS2'])
            
            print(f"\nStation {ev.get('ibl')} (Laminar):")
            print("  VS1 non-zero pattern:")
            for row in range(3):
                nz_cols = np.where(np.abs(vs1[row, :]) > 1e-20)[0]
                print(f"    Row {row}: cols {list(nz_cols)}")
            
            print("  VS2 non-zero pattern:")
            for row in range(3):
                nz_cols = np.where(np.abs(vs2[row, :]) > 1e-20)[0]
                print(f"    Row {row}: cols {list(nz_cols)}")


def main():
    parser = argparse.ArgumentParser(description='Analyze XFOIL Jacobian matrices')
    parser.add_argument('xfoil_json', help='Path to xfoil_debug.json')
    parser.add_argument('--station', type=int, default=None,
                        help='Filter to specific station')
    parser.add_argument('--side', type=int, default=None,
                        help='Filter to specific side (1 or 2)')
    parser.add_argument('--iteration', type=int, default=1,
                        help='Iteration number to analyze (default: 1)')
    args = parser.parse_args()
    
    print(f"Loading {args.xfoil_json}...")
    events = load_bldif_events(args.xfoil_json, 
                               iteration=args.iteration,
                               side=args.side,
                               station=args.station)
    
    print(f"Found {len(events)} BLDIF events")
    
    compare_jacobians_xfoil(events)
    check_jacobian_structure(events)
    
    print("\n" + "=" * 70)
    print("Analysis Complete")
    print("=" * 70)


if __name__ == '__main__':
    main()
