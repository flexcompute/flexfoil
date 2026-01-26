#!/usr/bin/env python3
"""
Compare BL distributions (Ue, θ, H, Cf) between RustFoil and XFOIL.
Shows WHERE the differences occur, not just integrated results.
"""

import json
import subprocess
import os
from pathlib import Path

XFOIL_BIN = Path(__file__).parent.parent / "Xfoil-instrumented/bin"
RUSTFOIL = Path(__file__).parent.parent / "target/release/rustfoil"
AIRFOIL = Path(__file__).parent.parent / "testdata/naca0012.dat"

def get_xfoil_distributions(alpha):
    """Extract full distributions from XFOIL debug JSON."""
    json_path = XFOIL_BIN / f"xfoil_debug_a{alpha}.json"
    if not json_path.exists():
        return None, None
    
    with open(json_path) as f:
        data = json.load(f)
    
    blvar = [e for e in data['events'] if e.get('subroutine') == 'BLVAR' and e.get('iteration') == 1]
    
    upper = []
    lower = []
    seen_upper = set()
    seen_lower = set()
    
    for e in blvar:
        side = e.get('side', 0)
        x = e['input']['x']
        
        if side == 1 and x not in seen_upper and x <= 1.05:
            upper.append({
                'x': x,
                'ue': e['input']['u'],
                'theta': e['input']['theta'],
                'dstar': e['input']['delta_star'],
                'h': e['output']['H'],
                'hk': e['output']['Hk'],
                'cf': e['output']['Cf'],
            })
            seen_upper.add(x)
        elif side == 2 and x not in seen_lower and x <= 1.05:
            lower.append({
                'x': x,
                'ue': e['input']['u'],
                'theta': e['input']['theta'],
                'dstar': e['input']['delta_star'],
                'h': e['output']['H'],
                'hk': e['output']['Hk'],
                'cf': e['output']['Cf'],
            })
            seen_lower.add(x)
    
    upper.sort(key=lambda s: s['x'])
    lower.sort(key=lambda s: s['x'])
    
    return upper, lower

def get_rustfoil_distributions(alpha):
    """Run RustFoil and parse full BL output."""
    env = os.environ.copy()
    env['RUSTFOIL_BL_DEBUG'] = '1'
    
    cmd = [str(RUSTFOIL), "viscous", str(AIRFOIL), "--alpha", str(alpha), "--re", "1e6"]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    
    output = result.stderr
    
    upper = []
    lower = []
    current_surface = None
    
    for line in output.split('\n'):
        if 'UPPER_SURFACE' in line:
            current_surface = 'upper'
        elif 'LOWER_SURFACE' in line:
            current_surface = 'lower'
        elif line.startswith('[BL_DEBUG]') and 'x=' in line:
            parts = line.replace('[BL_DEBUG]', '').strip().split()
            data = {}
            for part in parts:
                if '=' in part:
                    key, val = part.split('=')
                    try:
                        data[key] = float(val)
                    except:
                        pass
            
            if data and data.get('x', 0) <= 1.05:
                station = {
                    'x': data.get('x', 0),
                    'ue': data.get('Ue', 0),
                    'theta': data.get('theta', 0),
                    'dstar': data.get('dstar', 0),
                    'h': data.get('H', 0),
                    'hk': data.get('Hk', 0),
                    'cf': data.get('Cf', 0),
                }
                if current_surface == 'upper':
                    upper.append(station)
                elif current_surface == 'lower':
                    lower.append(station)
    
    return upper, lower

def interpolate_at_x(stations, x_target):
    """Linear interpolation at specific x."""
    if not stations:
        return None
    
    # Find bracketing stations
    for i in range(len(stations) - 1):
        if stations[i]['x'] <= x_target <= stations[i+1]['x']:
            x1, x2 = stations[i]['x'], stations[i+1]['x']
            t = (x_target - x1) / (x2 - x1) if x2 != x1 else 0.5
            
            result = {}
            for key in ['ue', 'theta', 'dstar', 'h', 'hk', 'cf']:
                v1 = stations[i].get(key, 0)
                v2 = stations[i+1].get(key, 0)
                result[key] = v1 + t * (v2 - v1)
            return result
    
    return None

def print_distribution(name, rf_stations, xf_stations, key, fmt='.4f'):
    """Print distribution comparison at matched x locations."""
    print(f"\n{name} Distribution:")
    print(f"{'x':>8} | {'RF':>10} | {'XF':>10} | {'Diff':>10} | {'Err%':>8}")
    print("-" * 55)
    
    # Use RF x locations and interpolate XF
    for rf_s in rf_stations:
        x = rf_s['x']
        rf_val = rf_s[key]
        
        xf_s = interpolate_at_x(xf_stations, x)
        if xf_s is None:
            continue
        xf_val = xf_s[key]
        
        diff = rf_val - xf_val
        if abs(xf_val) > 1e-10:
            err = diff / xf_val * 100
            err_str = f"{err:+7.1f}%"
        else:
            err_str = "N/A"
        
        print(f"{x:8.4f} | {rf_val:10{fmt}} | {xf_val:10{fmt}} | {diff:+10{fmt}} | {err_str}")

def main():
    alpha = 4  # Use α=4° as representative case
    
    print("=" * 70)
    print(f"DISTRIBUTION COMPARISON: RustFoil vs XFOIL at α={alpha}°")
    print("=" * 70)
    
    xf_upper, xf_lower = get_xfoil_distributions(alpha)
    rf_upper, rf_lower = get_rustfoil_distributions(alpha)
    
    if xf_upper is None:
        print("No XFOIL data available")
        return
    
    print(f"\nStations: RF upper={len(rf_upper)}, lower={len(rf_lower)}")
    print(f"          XF upper={len(xf_upper)}, lower={len(xf_lower)}")
    
    # ===== UPPER SURFACE =====
    print("\n" + "=" * 70)
    print("UPPER SURFACE")
    print("=" * 70)
    
    print_distribution("Ue (edge velocity)", rf_upper, xf_upper, 'ue', '.4f')
    print_distribution("θ (momentum thickness)", rf_upper, xf_upper, 'theta', '.2e')
    print_distribution("H (shape factor)", rf_upper, xf_upper, 'h', '.3f')
    print_distribution("Cf (skin friction)", rf_upper, xf_upper, 'cf', '.2e')
    
    # ===== LOWER SURFACE =====
    print("\n" + "=" * 70)
    print("LOWER SURFACE")
    print("=" * 70)
    
    print_distribution("Ue (edge velocity)", rf_lower, xf_lower, 'ue', '.4f')
    print_distribution("θ (momentum thickness)", rf_lower, xf_lower, 'theta', '.2e')
    print_distribution("H (shape factor)", rf_lower, xf_lower, 'h', '.3f')
    print_distribution("Cf (skin friction)", rf_lower, xf_lower, 'cf', '.2e')
    
    # ===== FIND LARGEST DEVIATIONS =====
    print("\n" + "=" * 70)
    print("LARGEST DEVIATIONS (where differences concentrate)")
    print("=" * 70)
    
    # Collect all deviations
    deviations = []
    for surf_name, rf_st, xf_st in [('Upper', rf_upper, xf_upper), ('Lower', rf_lower, xf_lower)]:
        for rf_s in rf_st:
            x = rf_s['x']
            xf_s = interpolate_at_x(xf_st, x)
            if xf_s is None:
                continue
            
            for key in ['ue', 'theta', 'h', 'cf']:
                rf_val = rf_s[key]
                xf_val = xf_s[key]
                if abs(xf_val) > 1e-10:
                    err_pct = (rf_val - xf_val) / xf_val * 100
                    deviations.append({
                        'surface': surf_name,
                        'x': x,
                        'quantity': key,
                        'rf': rf_val,
                        'xf': xf_val,
                        'err_pct': err_pct
                    })
    
    # Sort by absolute error
    deviations.sort(key=lambda d: abs(d['err_pct']), reverse=True)
    
    print(f"\nTop 15 largest percentage errors:")
    print(f"{'Surface':<8} | {'x':>6} | {'Quantity':>8} | {'RF':>10} | {'XF':>10} | {'Error':>8}")
    print("-" * 65)
    for d in deviations[:15]:
        print(f"{d['surface']:<8} | {d['x']:6.3f} | {d['quantity']:>8} | {d['rf']:10.4g} | {d['xf']:10.4g} | {d['err_pct']:+7.1f}%")

if __name__ == "__main__":
    main()
