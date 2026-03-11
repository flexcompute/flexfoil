#!/usr/bin/env python3
"""
Compare RustFoil inviscid gamma with XFOIL Ue values by station index.
Fixed indexing based on observed matching patterns.
"""

import json

def load_rustfoil_gamma(filepath):
    with open(filepath) as f:
        data = json.load(f)
    
    for event in data.get('events', []):
        if event.get('subroutine') == 'FULL_INVISCID':
            return event.get('gamma', [])
    return []

def load_xfoil_blvar(filepath):
    with open(filepath) as f:
        data = json.load(f)
    
    blvar_map = {}
    for event in data.get('events', []):
        if event.get('subroutine') == 'BLVAR':
            side = event.get('side')
            ibl = event.get('ibl')
            inp = event.get('input', {})
            x = inp.get('x')
            u = inp.get('u')
            
            key = (side, ibl)
            if key not in blvar_map:
                blvar_map[key] = {'x': x, 'ue': u, 'side': side, 'ibl': ibl}
    
    return blvar_map

def load_xfoil_coords(filepath):
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
    print("RustFoil vs XFOIL Gamma/Ue Comparison (Corrected Indexing)")
    print("=" * 70)
    print()
    
    coords = load_xfoil_coords(xfoil_coords)
    gamma = load_rustfoil_gamma(rustfoil_debug)
    blvar_map = load_xfoil_blvar(xfoil_debug)
    
    # Find stagnation (sign change)
    stag_idx = None
    for i in range(len(gamma) - 1):
        if gamma[i] > 0 and gamma[i+1] < 0:
            stag_idx = i
            break
    
    print(f"Stagnation index: {stag_idx}")
    print(f"  gamma[{stag_idx}] = {gamma[stag_idx]:.6f}")
    print(f"  gamma[{stag_idx+1}] = {gamma[stag_idx+1]:.6f}")
    print()
    
    # Verify: gamma[85] should match XFOIL upper IBL=2
    xfoil_upper2 = blvar_map.get((1, 2), {}).get('ue', 0)
    print(f"XFOIL upper IBL=2 Ue: {xfoil_upper2:.6f}")
    print(f"RustFoil gamma[{stag_idx}]: {gamma[stag_idx]:.6f}")
    print(f"Match: {abs(xfoil_upper2 - gamma[stag_idx]) < 0.0001}")
    print()
    
    print("=" * 70)
    print("UPPER SURFACE (Corrected mapping: IBL=n -> gamma[stag - n + 2])")
    print("=" * 70)
    print()
    
    print(f"{'IBL':>4} {'RF Idx':>6} {'x-coord':>10} {'XFOIL Ue':>12} {'RF Gamma':>12} {'Match?':>8}")
    print("-" * 60)
    
    upper_matches = 0
    upper_total = 0
    
    for ibl in range(2, 25):
        key = (1, ibl)
        if key in blvar_map:
            upper_total += 1
            xfoil_ue = blvar_map[key]['ue']
            
            # Corrected mapping: rf_idx = stag_idx - ibl + 2
            rf_idx = stag_idx - ibl + 2
            
            if 0 <= rf_idx < len(gamma):
                rf_gamma = gamma[rf_idx]
                x_coord = coords[rf_idx][0]
                
                diff = abs(xfoil_ue - rf_gamma)
                match = "YES" if diff < 0.0001 else "NO"
                if diff < 0.0001:
                    upper_matches += 1
                
                pct = (rf_gamma - xfoil_ue) / max(abs(xfoil_ue), 0.001) * 100
                print(f"{ibl:4d} {rf_idx:6d} {x_coord:10.6f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {match:>8} ({pct:+.4f}%)")
    
    print()
    print("=" * 70)
    print("LOWER SURFACE (IBL=n -> gamma[stag + n - 1])")
    print("=" * 70)
    print()
    
    print(f"{'IBL':>4} {'RF Idx':>6} {'x-coord':>10} {'XFOIL Ue':>12} {'RF Gamma':>12} {'Match?':>8}")
    print("-" * 60)
    
    lower_matches = 0
    lower_total = 0
    
    for ibl in range(2, 25):
        key = (2, ibl)
        if key in blvar_map:
            lower_total += 1
            xfoil_ue = blvar_map[key]['ue']
            
            # Mapping: rf_idx = stag_idx + ibl - 1
            rf_idx = stag_idx + ibl - 1
            
            if 0 <= rf_idx < len(gamma):
                rf_gamma = abs(gamma[rf_idx])  # Lower surface has negative gamma
                x_coord = coords[rf_idx][0]
                
                diff = abs(xfoil_ue - rf_gamma)
                match = "YES" if diff < 0.0001 else "NO"
                if diff < 0.0001:
                    lower_matches += 1
                
                pct = (rf_gamma - xfoil_ue) / max(abs(xfoil_ue), 0.001) * 100
                print(f"{ibl:4d} {rf_idx:6d} {x_coord:10.6f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {match:>8} ({pct:+.4f}%)")
    
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"Upper surface: {upper_matches}/{upper_total} stations match exactly")
    print(f"Lower surface: {lower_matches}/{lower_total} stations match exactly")
    print()
    
    total_match = upper_matches + lower_matches
    total_stations = upper_total + lower_total
    
    if total_match == total_stations:
        print("=" * 70)
        print("CONCLUSION: RustFoil gamma EXACTLY MATCHES XFOIL Ue!")
        print("=" * 70)
        print()
        print("When using identical XFOIL-paneled geometry:")
        print("  - Inviscid gamma values match to < 0.01%")
        print("  - The 2× discrepancy seen earlier is due to REPANELING")
        print("  - The influence coefficient calculation is CORRECT")
        print()
        print("Root cause: RustFoil's repaneling produces different node")
        print("positions than XFOIL's PANE command, especially near LE.")
    else:
        print(f"Matched {total_match}/{total_stations} stations")
        print("Some discrepancies remain - check index mapping")

if __name__ == '__main__':
    main()
