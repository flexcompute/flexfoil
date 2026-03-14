#!/usr/bin/env python3
"""
Compare BL quantities between RustFoil and XFOIL at specific x-locations.
"""

import json
import subprocess
import os
from pathlib import Path

XFOIL_BIN = Path(__file__).parent.parent / "Xfoil-instrumented/bin"
RUSTFOIL = Path(__file__).parent.parent / "target/release/rustfoil"
AIRFOIL = Path(__file__).parent.parent / "testdata/naca0012.dat"

def get_xfoil_bl(alpha):
    """Extract BL data from XFOIL debug JSON."""
    json_path = XFOIL_BIN / f"xfoil_debug_a{alpha}.json"
    if not json_path.exists():
        return None, None
    
    with open(json_path) as f:
        data = json.load(f)
    
    events = data.get('events', [])
    
    # Use iteration 1 for BL data
    blvar_events = [e for e in events if e.get('subroutine') == 'BLVAR' 
                    and e.get('iteration') == 1]
    
    upper = {}
    lower = {}
    seen = set()
    
    for e in blvar_events:
        side = e.get('side', 0)
        ibl = e.get('ibl', 0)
        key = (side, ibl)
        if key in seen:
            continue
        seen.add(key)
        
        inp = e.get('input', {})
        out = e.get('output', {})
        
        station = {
            'x': inp.get('x', 0),
            'u': inp.get('u', 0),
            'theta': inp.get('theta', 0),
            'delta_star': inp.get('delta_star', 0),
            'H': out.get('H', 0),
            'Hk': out.get('Hk', 0),
            'Cf': out.get('Cf', 0),
        }
        
        if side == 1:
            upper[station['x']] = station
        elif side == 2 and station['x'] <= 1.1:  # Skip wake
            lower[station['x']] = station
    
    return upper, lower

def interpolate_xfoil(stations, x_target):
    """Interpolate XFOIL data at specific x."""
    x_vals = sorted(stations.keys())
    if not x_vals:
        return None
    
    # Find bracketing
    for i in range(len(x_vals) - 1):
        if x_vals[i] <= x_target <= x_vals[i+1]:
            x1, x2 = x_vals[i], x_vals[i+1]
            t = (x_target - x1) / (x2 - x1) if x2 != x1 else 0.5
            
            s1 = stations[x1]
            s2 = stations[x2]
            result = {}
            for key in ['u', 'theta', 'delta_star', 'H', 'Hk', 'Cf']:
                result[key] = s1[key] + t * (s2[key] - s1[key])
            result['x'] = x_target
            return result
    
    # Return closest if out of range
    closest = min(x_vals, key=lambda x: abs(x - x_target))
    return stations[closest]

def get_rustfoil_bl(alpha):
    """Run RustFoil with BL debug output."""
    env = os.environ.copy()
    env['RUSTFOIL_BL_DEBUG'] = '1'
    
    cmd = [str(RUSTFOIL), "viscous", str(AIRFOIL), "--alpha", str(alpha), "--re", "1e6"]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    
    output = result.stderr + result.stdout
    
    upper = []
    lower = []
    current_surface = None
    cl = cd = None
    
    for line in output.split('\n'):
        if 'UPPER_SURFACE' in line:
            current_surface = 'upper'
        elif 'LOWER_SURFACE' in line:
            current_surface = 'lower'
        elif line.startswith('[BL_DEBUG]') and 'x=' in line:
            # Parse: x=0.1008 Ue=1.4288 theta=1.78e-4 ...
            parts = line.replace('[BL_DEBUG]', '').strip().split()
            data = {}
            for part in parts:
                if '=' in part:
                    key, val = part.split('=')
                    try:
                        data[key] = float(val)
                    except:
                        pass
            if data:
                station = {
                    'x': data.get('x', 0),
                    'u': data.get('Ue', 0),
                    'theta': data.get('theta', 0),
                    'delta_star': data.get('dstar', 0),
                    'H': data.get('H', 0),
                    'Hk': data.get('Hk', 0),
                    'Cf': data.get('Cf', 0),
                }
                if current_surface == 'upper':
                    upper.append(station)
                elif current_surface == 'lower':
                    lower.append(station)
        elif 'CL:' in line:
            try:
                cl = float(line.split()[-1])
            except:
                pass
        elif 'CD:' in line and 'CD_' not in line:
            try:
                cd = float(line.split()[-1])
            except:
                pass
    
    return upper, lower, cl, cd

def pct_diff(rf, xf):
    """Calculate percentage difference."""
    if xf == 0:
        return float('nan') if rf != 0 else 0
    return (rf - xf) / abs(xf) * 100

def main():
    alphas = [0, 2, 4, 6, 8]
    x_locs = [0.1, 0.3, 0.5, 0.7, 0.9]
    
    # Collect all data for trend analysis
    all_data = []
    
    for alpha in alphas:
        print(f"\n{'='*100}")
        print(f"ALPHA = {alpha}°")
        print(f"{'='*100}")
        
        xf_upper, xf_lower = get_xfoil_bl(alpha)
        rf_upper, rf_lower, rf_cl, rf_cd = get_rustfoil_bl(alpha)
        
        if xf_upper is None:
            print("No XFOIL data")
            continue
        
        # Get XFOIL forces from VISCAL_RESULT
        json_path = XFOIL_BIN / f"xfoil_debug_a{alpha}.json"
        with open(json_path) as f:
            data = json.load(f)
        viscres = [e for e in data['events'] if e.get('subroutine') == 'VISCAL_RESULT']
        xf_cl = viscres[-1]['CL'] if viscres else 0
        xf_cd = viscres[-1]['CD'] if viscres else 0
        
        print(f"\nFORCES: RF_CL={rf_cl:.4f}, XF_CL={xf_cl:.4f} ({pct_diff(rf_cl, xf_cl):+.1f}%)")
        print(f"        RF_CD={rf_cd:.5f}, XF_CD={xf_cd:.5f} ({pct_diff(rf_cd, xf_cd):+.1f}%)")
        
        all_data.append({
            'alpha': alpha,
            'rf_cl': rf_cl, 'xf_cl': xf_cl,
            'rf_cd': rf_cd, 'xf_cd': xf_cd,
        })
        
        # Upper surface comparison
        print(f"\nUPPER SURFACE:")
        print(f"{'x':>6} | {'RF_Ue':>7} {'XF_Ue':>7} {'%':>6} | {'RF_θ':>9} {'XF_θ':>9} {'%':>6} | {'RF_H':>6} {'XF_H':>6} {'%':>6} | {'RF_Cf':>9} {'XF_Cf':>9} {'%':>6}")
        print("-" * 110)
        
        for i, rf_s in enumerate(rf_upper):
            x_rf = rf_s['x']
            xf_s = interpolate_xfoil(xf_upper, x_rf)
            if xf_s is None:
                continue
            
            ue_d = pct_diff(rf_s['u'], xf_s['u'])
            th_d = pct_diff(rf_s['theta'], xf_s['theta'])
            h_d = pct_diff(rf_s['H'], xf_s['H'])
            cf_d = pct_diff(rf_s['Cf'], xf_s['Cf'])
            
            print(f"{x_rf:6.3f} | {rf_s['u']:7.4f} {xf_s['u']:7.4f} {ue_d:+5.1f}% | {rf_s['theta']:9.2e} {xf_s['theta']:9.2e} {th_d:+5.1f}% | {rf_s['H']:6.3f} {xf_s['H']:6.3f} {h_d:+5.1f}% | {rf_s['Cf']:9.2e} {xf_s['Cf']:9.2e} {cf_d:+5.1f}%")
        
        # Lower surface comparison  
        print(f"\nLOWER SURFACE:")
        print(f"{'x':>6} | {'RF_Ue':>7} {'XF_Ue':>7} {'%':>6} | {'RF_θ':>9} {'XF_θ':>9} {'%':>6} | {'RF_H':>6} {'XF_H':>6} {'%':>6} | {'RF_Cf':>9} {'XF_Cf':>9} {'%':>6}")
        print("-" * 110)
        
        for i, rf_s in enumerate(rf_lower):
            x_rf = rf_s['x']
            xf_s = interpolate_xfoil(xf_lower, x_rf)
            if xf_s is None:
                continue
            
            ue_d = pct_diff(rf_s['u'], xf_s['u'])
            th_d = pct_diff(rf_s['theta'], xf_s['theta'])
            h_d = pct_diff(rf_s['H'], xf_s['H'])
            cf_d = pct_diff(rf_s['Cf'], xf_s['Cf'])
            
            print(f"{x_rf:6.3f} | {rf_s['u']:7.4f} {xf_s['u']:7.4f} {ue_d:+5.1f}% | {rf_s['theta']:9.2e} {xf_s['theta']:9.2e} {th_d:+5.1f}% | {rf_s['H']:6.3f} {xf_s['H']:6.3f} {h_d:+5.1f}% | {rf_s['Cf']:9.2e} {xf_s['Cf']:9.2e} {cf_d:+5.1f}%")
    
    # Summary trend analysis
    print(f"\n{'='*100}")
    print("ERROR TREND SUMMARY")
    print(f"{'='*100}")
    print(f"{'Alpha':>5} | {'CL_err%':>8} | {'CD_err%':>8} | Notes")
    print("-" * 60)
    for d in all_data:
        cl_err = pct_diff(d['rf_cl'], d['xf_cl'])
        cd_err = pct_diff(d['rf_cd'], d['xf_cd'])
        notes = []
        if abs(cl_err) > 10:
            notes.append("HIGH CL error")
        if cd_err < -30:
            notes.append("CD too low")
        if cd_err > 10:
            notes.append("CD too high")
        print(f"{d['alpha']:5}° | {cl_err:+8.1f} | {cd_err:+8.1f} | {', '.join(notes)}")

if __name__ == "__main__":
    main()
