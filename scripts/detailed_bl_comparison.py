#!/usr/bin/env python3
"""
Detailed comparison of boundary layer quantities at specific stations.
Extracts Ue, theta, delta_star, H, Cf from XFOIL BLVAR events.
"""

import json
import subprocess
from pathlib import Path
import numpy as np

XFOIL_BIN = Path(__file__).parent.parent / "Xfoil-instrumented/bin"
RUSTFOIL = Path(__file__).parent.parent / "target/release/rustfoil"
AIRFOIL = Path(__file__).parent.parent / "testdata/naca0012.dat"

def get_xfoil_bl_data(alpha):
    """Extract BL quantities from XFOIL debug at specific x locations."""
    json_path = XFOIL_BIN / f"xfoil_debug_a{alpha}.json"
    if not json_path.exists():
        return None
    
    with open(json_path) as f:
        data = json.load(f)
    
    events = data.get('events', [])
    
    # Get final iteration for forces
    viscres = [e for e in events if e.get('subroutine') == 'VISCAL_RESULT']
    if not viscres:
        return None
    
    # Use iteration 1 for BL data (has both surfaces)
    # Later iterations only have wake data
    blvar_iter = 1
    
    # Get BLVAR events
    blvar_events = [e for e in events if e.get('subroutine') == 'BLVAR' 
                    and e.get('iteration') == blvar_iter]
    
    # Organize by surface
    upper = []
    lower = []
    
    seen_upper = set()
    seen_lower = set()
    
    for e in blvar_events:
        side = e.get('side', 0)
        ibl = e.get('ibl', 0)
        inp = e.get('input', {})
        out = e.get('output', {})
        
        key = (side, ibl)
        
        station = {
            'ibl': ibl,
            'x': inp.get('x', 0),
            'u': inp.get('u', 0),
            'theta': inp.get('theta', 0),
            'delta_star': inp.get('delta_star', 0),
            'H': out.get('H', 0),
            'Hk': out.get('Hk', 0),
            'Cf': out.get('Cf', 0),
            'Rtheta': out.get('Rtheta', 0),
        }
        
        if side == 1 and key not in seen_upper:
            upper.append(station)
            seen_upper.add(key)
        elif side == 2 and key not in seen_lower:
            lower.append(station)
            seen_lower.add(key)
    
    # Sort by x
    upper.sort(key=lambda s: s['x'])
    lower.sort(key=lambda s: s['x'])
    
    # Get forces
    final = viscres[-1]
    
    return {
        'alpha': alpha,
        'CL': final.get('CL', 0),
        'CD': final.get('CD', 0),
        'CD_f': final.get('CD_friction', 0),
        'CD_p': final.get('CD_viscous', 0),
        'upper': upper,
        'lower': lower,
        'iterations': final.get('iteration', 0),
    }

def interpolate_at_x(stations, target_x):
    """Interpolate BL quantities at a specific x location."""
    if not stations:
        return None
    
    # Find bracketing stations
    x_vals = [s['x'] for s in stations]
    
    if target_x < min(x_vals) or target_x > max(x_vals):
        return None
    
    # Linear interpolation
    for i in range(len(stations) - 1):
        if stations[i]['x'] <= target_x <= stations[i+1]['x']:
            x1, x2 = stations[i]['x'], stations[i+1]['x']
            t = (target_x - x1) / (x2 - x1) if x2 != x1 else 0.5
            
            result = {}
            for key in ['u', 'theta', 'delta_star', 'H', 'Hk', 'Cf']:
                v1 = stations[i].get(key, 0)
                v2 = stations[i+1].get(key, 0)
                result[key] = v1 + t * (v2 - v1)
            result['x'] = target_x
            return result
    
    return None

def main():
    alphas = [0, 2, 4, 6, 8]
    x_locations = [0.1, 0.3, 0.5, 0.7, 0.9]
    
    print("=" * 80)
    print("DETAILED BL COMPARISON: XFOIL values at key x-locations")
    print("=" * 80)
    
    # Collect data for analysis
    all_xfoil_data = {}
    
    for alpha in alphas:
        xf = get_xfoil_bl_data(alpha)
        if xf is None:
            print(f"\nNo XFOIL data for alpha={alpha}°")
            continue
        
        all_xfoil_data[alpha] = xf
        
        print(f"\n{'='*80}")
        print(f"ALPHA = {alpha}°   CL={xf['CL']:.4f}  CD={xf['CD']:.5f}")
        print(f"{'='*80}")
        
        # Upper surface
        print(f"\nUPPER SURFACE ({len(xf['upper'])} stations):")
        print(f"{'x':>8} | {'Ue':>8} | {'theta':>10} | {'d*':>10} | {'H':>6} | {'Hk':>6} | {'Cf':>10}")
        print("-" * 80)
        
        for x_loc in x_locations:
            s = interpolate_at_x(xf['upper'], x_loc)
            if s:
                print(f"{s['x']:8.4f} | {s['u']:8.4f} | {s['theta']:10.2e} | {s['delta_star']:10.2e} | {s['H']:6.3f} | {s['Hk']:6.3f} | {s['Cf']:10.2e}")
        
        # Also print TE values
        if xf['upper']:
            te = xf['upper'][-1]
            print(f"   TE    | {te['u']:8.4f} | {te['theta']:10.2e} | {te['delta_star']:10.2e} | {te['H']:6.3f} | {te['Hk']:6.3f} | {te['Cf']:10.2e}")
        
        # Lower surface
        print(f"\nLOWER SURFACE ({len(xf['lower'])} stations):")
        print(f"{'x':>8} | {'Ue':>8} | {'theta':>10} | {'d*':>10} | {'H':>6} | {'Hk':>6} | {'Cf':>10}")
        print("-" * 80)
        
        for x_loc in x_locations:
            s = interpolate_at_x(xf['lower'], x_loc)
            if s:
                print(f"{s['x']:8.4f} | {s['u']:8.4f} | {s['theta']:10.2e} | {s['delta_star']:10.2e} | {s['H']:6.3f} | {s['Hk']:6.3f} | {s['Cf']:10.2e}")
        
        if xf['lower']:
            te = xf['lower'][-1]
            print(f"   TE    | {te['u']:8.4f} | {te['theta']:10.2e} | {te['delta_star']:10.2e} | {te['H']:6.3f} | {te['Hk']:6.3f} | {te['Cf']:10.2e}")
    
    # Summary of key trends
    print("\n" + "=" * 80)
    print("KEY TRENDS ACROSS ALPHA")
    print("=" * 80)
    
    print(f"\n{'Alpha':>5} | {'CL':>8} | {'CD':>8} | {'CD_f':>8} | {'CD_p':>8} | {'CD_p/CD':>8}")
    print("-" * 60)
    for alpha in alphas:
        if alpha in all_xfoil_data:
            xf = all_xfoil_data[alpha]
            cdp_ratio = xf['CD_p'] / xf['CD'] if xf['CD'] > 0 else 0
            print(f"{alpha:5}° | {xf['CL']:8.4f} | {xf['CD']:8.5f} | {xf['CD_f']:8.5f} | {xf['CD_p']:8.5f} | {cdp_ratio:8.1%}")
    
    # Ue at x=0.5 (suction peak region)
    print(f"\nUe at x=0.5 (suction peak indicator):")
    print(f"{'Alpha':>5} | {'Upper Ue':>10} | {'Lower Ue':>10} | {'Ue diff':>10}")
    print("-" * 50)
    for alpha in alphas:
        if alpha in all_xfoil_data:
            xf = all_xfoil_data[alpha]
            upper_s = interpolate_at_x(xf['upper'], 0.5)
            lower_s = interpolate_at_x(xf['lower'], 0.5)
            if upper_s and lower_s:
                diff = upper_s['u'] - lower_s['u']
                print(f"{alpha:5}° | {upper_s['u']:10.4f} | {lower_s['u']:10.4f} | {diff:+10.4f}")
    
    # theta at TE (key for drag)
    print(f"\nθ at TE (key for CD):")
    print(f"{'Alpha':>5} | {'Upper θ':>12} | {'Lower θ':>12} | {'Sum θ':>12}")
    print("-" * 55)
    for alpha in alphas:
        if alpha in all_xfoil_data:
            xf = all_xfoil_data[alpha]
            if xf['upper'] and xf['lower']:
                theta_u = xf['upper'][-1]['theta']
                theta_l = xf['lower'][-1]['theta']
                print(f"{alpha:5}° | {theta_u:12.2e} | {theta_l:12.2e} | {theta_u+theta_l:12.2e}")

if __name__ == "__main__":
    main()
