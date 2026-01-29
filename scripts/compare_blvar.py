#!/usr/bin/env python3
"""
Compare BLVAR events between XFOIL and RustFoil to find first divergence.
BLVAR contains: x, u (Ue), theta, delta_star, Hk, Cf, etc. at each station.
"""

import json
import sys
import numpy as np
from collections import defaultdict

def load_trace(path):
    with open(path) as f:
        return json.load(f).get('events', [])

def get_blvar_by_station(events):
    """Group BLVAR events by (iteration, side, ibl)."""
    blvar = defaultdict(list)
    for e in events:
        if e.get('subroutine') == 'BLVAR':
            key = (e.get('iteration', 0), e.get('side'), e.get('ibl'))
            blvar[key].append(e)
    return blvar

def extract_blvar_value(event, field):
    """Extract a value from BLVAR event (handles different structures)."""
    # RustFoil structure: input/output dicts
    if 'input' in event:
        if field in event['input']:
            return event['input'][field]
        if 'output' in event and field in event['output']:
            return event['output'][field]
    # XFOIL structure: flat
    if field in event:
        return event[field]
    # Try uppercase
    if field.upper() in event:
        return event[field.upper()]
    return None

def compare_blvar_field(xf_events, rf_events, field, tolerance=0.05):
    """Compare a specific field across all stations."""
    xf_blvar = get_blvar_by_station(xf_events)
    rf_blvar = get_blvar_by_station(rf_events)
    
    # Focus on iteration 1, upper surface (side=1)
    results = []
    
    for ibl in range(1, 100):
        key = (1, 1, ibl)  # iteration 1, upper surface
        
        if key not in xf_blvar or key not in rf_blvar:
            continue
        
        xf_val = extract_blvar_value(xf_blvar[key][0], field)
        rf_val = extract_blvar_value(rf_blvar[key][0], field)
        
        if xf_val is None or rf_val is None:
            continue
        
        if abs(xf_val) > 1e-10:
            rel_err = abs(xf_val - rf_val) / abs(xf_val)
        else:
            rel_err = abs(xf_val - rf_val)
        
        results.append({
            'ibl': ibl,
            'xfoil': xf_val,
            'rustfoil': rf_val,
            'rel_err': rel_err
        })
    
    return results

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_blvar.py <xfoil_trace> <rustfoil_trace>")
        sys.exit(1)
    
    print(f"Loading traces...")
    xf_events = load_trace(sys.argv[1])
    rf_events = load_trace(sys.argv[2])
    
    print(f"XFOIL: {len(xf_events)} events")
    print(f"RustFoil: {len(rf_events)} events")
    
    # Fields to compare (in order of calculation)
    fields = [
        ('x', 'x-coordinate'),
        ('u', 'Edge velocity Ue'),
        ('theta', 'Momentum thickness'),
        ('delta_star', 'Displacement thickness'),
        ('Hk', 'Kinematic shape factor'),
        ('Cf', 'Skin friction'),
        ('Rtheta', 'Reynolds number based on theta'),
    ]
    
    print("\n" + "="*100)
    print("BLVAR COMPARISON - Iteration 1, Upper Surface")
    print("="*100)
    
    for field, desc in fields:
        results = compare_blvar_field(xf_events, rf_events, field)
        
        if not results:
            print(f"\n{field} ({desc}): No data")
            continue
        
        errors = [r['rel_err'] for r in results]
        max_err_idx = np.argmax(errors)
        max_err = results[max_err_idx]
        mean_err = np.mean(errors)
        
        # Find first station with > 5% error
        first_divergence = None
        for r in results:
            if r['rel_err'] > 0.05:
                first_divergence = r
                break
        
        print(f"\n{field} ({desc}):")
        print(f"  Stations compared: {len(results)}")
        print(f"  Mean error: {mean_err:.2%}")
        print(f"  Max error: {max_err['rel_err']:.2%} at station {max_err['ibl']}")
        print(f"    XFOIL:    {max_err['xfoil']:.6e}")
        print(f"    RustFoil: {max_err['rustfoil']:.6e}")
        
        if first_divergence:
            print(f"  FIRST DIVERGENCE (>5%): Station {first_divergence['ibl']}")
            print(f"    XFOIL:    {first_divergence['xfoil']:.6e}")
            print(f"    RustFoil: {first_divergence['rustfoil']:.6e}")
            print(f"    Error:    {first_divergence['rel_err']:.2%}")
        
        # Show first 10 stations
        print(f"\n  First 10 stations:")
        print(f"  {'IBL':>4} {'XFOIL':>12} {'RustFoil':>12} {'Error':>10}")
        print(f"  {'-'*42}")
        for r in results[:10]:
            print(f"  {r['ibl']:4d} {r['xfoil']:12.6e} {r['rustfoil']:12.6e} {r['rel_err']:10.2%}")

if __name__ == '__main__':
    main()
