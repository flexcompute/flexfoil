#!/usr/bin/env python3
"""
Compare RustFoil inviscid gamma with XFOIL Ue values at identical geometry.

This script:
1. Extracts gamma from RustFoil debug_identical.json (FULL_INVISCID event)
2. Extracts Ue from XFOIL xfoil_debug_a4.json (BLVAR events)
3. Compares values at matching x positions
"""

import json
import sys

def load_rustfoil_gamma(filepath):
    """Load gamma from RustFoil debug output."""
    with open(filepath) as f:
        data = json.load(f)
    
    # Find the FULL_INVISCID event
    for event in data.get('events', []):
        if event.get('subroutine') == 'FULL_INVISCID':
            return event.get('gamma', [])
    return []

def load_xfoil_blvar(filepath, max_events=200):
    """Load BLVAR Ue values from XFOIL debug output."""
    # The file is too large to load entirely, so we'll read it line by line
    # looking for BLVAR entries
    
    blvar_entries = []
    
    with open(filepath) as f:
        content = f.read()
    
    # Parse the JSON - it's a large array of events
    data = json.loads(content)
    
    seen_positions = set()
    
    for event in data.get('events', []):
        if event.get('subroutine') == 'BLVAR':
            side = event.get('side')
            ibl = event.get('ibl')
            inp = event.get('input', {})
            x = inp.get('x')
            u = inp.get('u')
            
            if x is not None and u is not None:
                key = (side, ibl, x)
                if key not in seen_positions:
                    seen_positions.add(key)
                    blvar_entries.append({
                        'side': side,
                        'ibl': ibl,
                        'x': x,
                        'ue': u
                    })
        
        if len(seen_positions) > max_events:
            break
    
    return blvar_entries

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
    # File paths
    rustfoil_debug = 'debug_identical.json'
    xfoil_debug = 'Xfoil-instrumented/bin/xfoil_debug_a4.json'
    xfoil_coords = 'naca0012_xfoil_paneled.dat'
    
    print("=" * 70)
    print("RustFoil vs XFOIL Gamma/Ue Comparison (Identical Geometry)")
    print("=" * 70)
    print()
    
    # Load coordinates
    coords = load_xfoil_coords(xfoil_coords)
    print(f"Loaded {len(coords)} XFOIL-paneled coordinates")
    
    # Load RustFoil gamma
    gamma = load_rustfoil_gamma(rustfoil_debug)
    print(f"Loaded {len(gamma)} RustFoil gamma values")
    
    # Load XFOIL BLVAR entries
    print("Loading XFOIL BLVAR entries...")
    blvar = load_xfoil_blvar(xfoil_debug)
    print(f"Loaded {len(blvar)} XFOIL BLVAR entries")
    
    print()
    print("=" * 70)
    print("RustFoil Gamma Distribution (near stagnation)")
    print("=" * 70)
    print()
    
    # Find stagnation point in RustFoil gamma (where sign changes)
    stag_idx = None
    for i in range(len(gamma) - 1):
        if gamma[i] > 0 and gamma[i+1] < 0:
            stag_idx = i
            break
        if gamma[i] < 0 and gamma[i+1] > 0:
            stag_idx = i
            break
    
    print(f"RustFoil stagnation index: ~{stag_idx} (gamma changes sign)")
    print()
    
    # Print gamma values around stagnation
    print(f"{'Idx':>4} {'x':>10} {'y':>10} {'gamma':>12}")
    print("-" * 40)
    
    # Upper surface (before stagnation) - going backward from stag
    if stag_idx:
        for i in range(min(stag_idx, 10)):
            idx = stag_idx - i
            if 0 <= idx < len(coords) and idx < len(gamma):
                x, y = coords[idx]
                print(f"{idx:4d} {x:10.6f} {y:10.6f} {gamma[idx]:12.6f}")
    
    print()
    print("=" * 70)
    print("XFOIL BLVAR Ue Distribution (Upper Surface, near LE)")
    print("=" * 70)
    print()
    
    # Filter for upper surface (side=1) and small x values
    upper_blvar = [b for b in blvar if b['side'] == 1 and b['x'] < 0.1]
    upper_blvar.sort(key=lambda b: b['x'])
    
    print(f"{'IBL':>4} {'x':>10} {'Ue':>12}")
    print("-" * 30)
    for b in upper_blvar[:15]:
        print(f"{b['ibl']:4d} {b['x']:10.6f} {b['ue']:12.6f}")
    
    print()
    print("=" * 70)
    print("COMPARISON TABLE: RustFoil Gamma vs XFOIL Ue")
    print("=" * 70)
    print()
    print("Note: gamma IS the surface velocity Ue in XFOIL's linear vortex method")
    print()
    
    # Create comparison by matching x positions
    # RustFoil: upper surface is indices 0 to stag_idx (from TE to LE)
    # XFOIL BLVAR: side=1 is upper surface
    
    # Build comparison table for specific x positions
    target_x = [0.001, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03, 0.05]
    
    print(f"{'x':>10} {'XFOIL Ue':>12} {'RustFoil γ':>12} {'Ratio':>8} {'Diff %':>10}")
    print("-" * 55)
    
    for tx in target_x:
        # Find closest XFOIL BLVAR entry
        xfoil_entry = min(upper_blvar, key=lambda b: abs(b['x'] - tx), default=None)
        
        # Find closest RustFoil gamma (upper surface = indices going backward from stag)
        if stag_idx:
            best_rf_idx = None
            best_rf_dist = float('inf')
            for i in range(stag_idx, -1, -1):
                if i < len(coords):
                    x, y = coords[i]
                    if y > 0:  # upper surface
                        dist = abs(x - tx)
                        if dist < best_rf_dist:
                            best_rf_dist = dist
                            best_rf_idx = i
            
            if xfoil_entry and best_rf_idx is not None:
                xfoil_ue = xfoil_entry['ue']
                rf_gamma = abs(gamma[best_rf_idx])  # Take absolute value
                ratio = rf_gamma / xfoil_ue if xfoil_ue != 0 else float('inf')
                diff_pct = 100 * (rf_gamma - xfoil_ue) / xfoil_ue if xfoil_ue != 0 else float('inf')
                print(f"{tx:10.4f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {ratio:8.3f} {diff_pct:10.1f}%")
    
    print()
    print("=" * 70)
    print("DETAILED COMPARISON (First 20 Upper Surface Stations)")
    print("=" * 70)
    print()
    
    print(f"{'IBL':>4} {'x':>10} {'XFOIL Ue':>12} {'RF Gamma':>12} {'Ratio':>8}")
    print("-" * 50)
    
    for b in upper_blvar[:20]:
        ibl = b['ibl']
        x = b['x']
        xfoil_ue = b['ue']
        
        # Find matching RustFoil gamma by x position
        if stag_idx:
            best_rf_idx = None
            best_rf_dist = float('inf')
            for i in range(stag_idx, -1, -1):
                if i < len(coords):
                    coord_x, coord_y = coords[i]
                    if coord_y > 0:  # upper surface
                        dist = abs(coord_x - x)
                        if dist < best_rf_dist:
                            best_rf_dist = dist
                            best_rf_idx = i
            
            if best_rf_idx is not None:
                rf_gamma = abs(gamma[best_rf_idx])
                ratio = rf_gamma / xfoil_ue if xfoil_ue != 0 else float('inf')
                print(f"{ibl:4d} {x:10.6f} {xfoil_ue:12.6f} {rf_gamma:12.6f} {ratio:8.3f}")
    
    print()
    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()
    
    # Compute average ratio for x < 0.02
    ratios = []
    for b in upper_blvar:
        if b['x'] < 0.02:
            xfoil_ue = b['ue']
            x = b['x']
            
            if stag_idx:
                best_rf_idx = None
                best_rf_dist = float('inf')
                for i in range(stag_idx, -1, -1):
                    if i < len(coords):
                        coord_x, coord_y = coords[i]
                        if coord_y > 0:
                            dist = abs(coord_x - x)
                            if dist < best_rf_dist:
                                best_rf_dist = dist
                                best_rf_idx = i
                
                if best_rf_idx is not None and xfoil_ue > 0.001:
                    rf_gamma = abs(gamma[best_rf_idx])
                    ratio = rf_gamma / xfoil_ue
                    ratios.append((x, ratio))
    
    if ratios:
        avg_ratio = sum(r for _, r in ratios) / len(ratios)
        print(f"Average RustFoil/XFOIL ratio for x < 0.02: {avg_ratio:.3f}")
        
        if avg_ratio > 1.5:
            print("\n→ RustFoil gamma is ~{:.0f}× XFOIL Ue near stagnation".format(avg_ratio))
            print("→ Issue is in INFLUENCE COEFFICIENT calculation, not repaneling")
        elif abs(avg_ratio - 1.0) < 0.1:
            print("\n→ RustFoil gamma matches XFOIL Ue with identical geometry!")
            print("→ Issue is in REPANELING, not influence coefficients")
        else:
            print(f"\n→ Ratio is {avg_ratio:.2f}, unclear source of discrepancy")

if __name__ == '__main__':
    main()
