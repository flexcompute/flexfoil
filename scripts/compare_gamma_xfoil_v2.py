#!/usr/bin/env python3
"""
Compare RustFoil inviscid gamma with XFOIL Ue values by station index.

Key insight: XFOIL's BLVAR 'x' is arc length from stagnation, not x-coordinate.
We should compare by station index, not x-position.
"""

import json

def load_rustfoil_gamma(filepath):
    """Load gamma from RustFoil debug output."""
    with open(filepath) as f:
        data = json.load(f)
    
    for event in data.get('events', []):
        if event.get('subroutine') == 'FULL_INVISCID':
            return event.get('gamma', [])
    return []

def load_xfoil_blvar(filepath):
    """Load BLVAR Ue values from XFOIL debug output."""
    with open(filepath) as f:
        data = json.load(f)
    
    # Get unique BLVAR entries by (side, ibl)
    blvar_map = {}
    
    for event in data.get('events', []):
        if event.get('subroutine') == 'BLVAR':
            side = event.get('side')
            ibl = event.get('ibl')
            inp = event.get('input', {})
            x = inp.get('x')  # This is arc length from stagnation
            u = inp.get('u')  # Edge velocity
            
            key = (side, ibl)
            if key not in blvar_map:
                blvar_map[key] = {'x': x, 'ue': u, 'side': side, 'ibl': ibl}
    
    return blvar_map

def load_xfoil_coords(filepath):
    """Load XFOIL paneled coordinates."""
    coords = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    coords.append((x, y))
                except ValueError:
                    continue
    return coords

def main():
    rustfoil_debug = 'debug_identical.json'
    xfoil_debug = 'Xfoil-instrumented/bin/xfoil_debug_a4.json'
    xfoil_coords = 'naca0012_xfoil_paneled.dat'
    
    print("=" * 70)
    print("RustFoil vs XFOIL Gamma/Ue Comparison BY STATION INDEX")
    print("=" * 70)
    print()
    
    coords = load_xfoil_coords(xfoil_coords)
    gamma = load_rustfoil_gamma(rustfoil_debug)
    blvar_map = load_xfoil_blvar(xfoil_debug)
    
    print(f"RustFoil: {len(gamma)} gamma values")
    print(f"XFOIL: {len(blvar_map)} unique BLVAR entries")
    print()
    
    # Find stagnation in gamma (sign change)
    stag_idx = None
    for i in range(len(gamma) - 1):
        if gamma[i] > 0 and gamma[i+1] < 0:
            stag_idx = i
            break
    
    print(f"RustFoil stagnation panel index: {stag_idx}")
    print(f"  gamma[{stag_idx}] = {gamma[stag_idx]:.6f} (just before sign change)")
    print(f"  gamma[{stag_idx+1}] = {gamma[stag_idx+1]:.6f} (just after sign change)")
    print()
    
    # XFOIL IBL indexing:
    # - IBL=1 is stagnation point (virtual, Ue≈0)
    # - IBL=2 is first actual node
    # - Upper surface: side=1, IBL increases away from stagnation
    # - Lower surface: side=2, IBL increases away from stagnation
    
    print("=" * 70)
    print("UPPER SURFACE COMPARISON (side=1)")
    print("=" * 70)
    print()
    print("XFOIL IBL=2 corresponds to first panel node after stagnation on upper")
    print("In RustFoil 160-node geometry with stag at ~85:")
    print("  Upper: indices 84, 83, 82, ... (decreasing toward TE)")
    print()
    
    print(f"{'XFOIL':>6} {'RF Idx':>6} {'x-coord':>10} {'XFOIL Ue':>12} {'RF Gamma':>12} {'Match?':>8}")
    print("-" * 60)
    
    # Upper surface: XFOIL IBL 2,3,4,... -> RustFoil indices stag_idx-1, stag_idx-2, ...
    for ibl in range(2, 20):
        key = (1, ibl)  # side=1 is upper
        if key in blvar_map:
            xfoil_ue = blvar_map[key]['ue']
            xfoil_arc = blvar_map[key]['x']
            
            # RustFoil index: stag_idx - (ibl - 1)
            # IBL=2 -> stag_idx - 1 = 84
            # IBL=3 -> stag_idx - 2 = 83
            rf_idx = stag_idx - (ibl - 1)
            
            if 0 <= rf_idx < len(gamma):
                rf_gamma = gamma[rf_idx]
                x_coord = coords[rf_idx][0]
                
                match = "YES" if abs(xfoil_ue - rf_gamma) < 0.0001 else "NO"
                diff = (rf_gamma - xfoil_ue) / max(abs(xfoil_ue), 0.001) * 100
                
                print(f"IBL={ibl:2d} {rf_idx:6d} {x_coord:10.6f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {match:>8} ({diff:+.2f}%)")
    
    print()
    print("=" * 70)
    print("LOWER SURFACE COMPARISON (side=2)")
    print("=" * 70)
    print()
    
    print(f"{'XFOIL':>6} {'RF Idx':>6} {'x-coord':>10} {'XFOIL Ue':>12} {'RF Gamma':>12} {'Match?':>8}")
    print("-" * 60)
    
    # Lower surface: XFOIL IBL 2,3,4,... -> RustFoil indices stag_idx+1, stag_idx+2, ...
    # But gamma sign is negative on lower, XFOIL Ue is positive (absolute)
    for ibl in range(2, 20):
        key = (2, ibl)  # side=2 is lower
        if key in blvar_map:
            xfoil_ue = blvar_map[key]['ue']
            xfoil_arc = blvar_map[key]['x']
            
            # RustFoil index: stag_idx + (ibl - 1)
            rf_idx = stag_idx + (ibl - 1)
            
            if 0 <= rf_idx < len(gamma):
                rf_gamma = abs(gamma[rf_idx])  # Take absolute value for lower surface
                x_coord = coords[rf_idx][0]
                
                match = "YES" if abs(xfoil_ue - rf_gamma) < 0.0001 else "NO"
                diff = (rf_gamma - xfoil_ue) / max(abs(xfoil_ue), 0.001) * 100
                
                print(f"IBL={ibl:2d} {rf_idx:6d} {x_coord:10.6f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {match:>8} ({diff:+.2f}%)")
    
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    
    # Check if all match
    upper_matches = 0
    lower_matches = 0
    upper_total = 0
    lower_total = 0
    
    for ibl in range(2, min(30, stag_idx)):
        key = (1, ibl)
        if key in blvar_map:
            upper_total += 1
            rf_idx = stag_idx - (ibl - 1)
            if 0 <= rf_idx < len(gamma):
                if abs(blvar_map[key]['ue'] - gamma[rf_idx]) < 0.001:
                    upper_matches += 1
    
    for ibl in range(2, 30):
        key = (2, ibl)
        if key in blvar_map:
            lower_total += 1
            rf_idx = stag_idx + (ibl - 1)
            if 0 <= rf_idx < len(gamma):
                if abs(blvar_map[key]['ue'] - abs(gamma[rf_idx])) < 0.001:
                    lower_matches += 1
    
    print(f"Upper surface: {upper_matches}/{upper_total} stations match exactly")
    print(f"Lower surface: {lower_matches}/{lower_total} stations match exactly")
    print()
    
    if upper_matches == upper_total and lower_matches == lower_total:
        print("✓ RustFoil gamma EXACTLY MATCHES XFOIL Ue when using identical geometry!")
        print()
        print("CONCLUSION: The 2× discrepancy seen earlier was due to REPANELING,")
        print("            not the influence coefficient calculation.")
    else:
        print("✗ Some discrepancies found - further investigation needed")
        print()
        print("CONCLUSION: Issue may be in influence coefficient calculation")

if __name__ == '__main__':
    main()
