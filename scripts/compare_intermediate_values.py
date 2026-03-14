#!/usr/bin/env python3
"""
Compare intermediate values between RustFoil and XFOIL at specific angles.
"""
import json
import numpy as np
from pathlib import Path

def extract_bl_stations(debug_file, side=1):
    """Extract BL station data from debug trace."""
    with open(debug_file) as f:
        data = json.load(f)
    
    stations = []
    for e in data['events']:
        if e.get('subroutine') == 'MRCHUE' and e.get('side') == side:
            stations.append({
                'ibl': e['ibl'],
                'x': e['x'],
                'x_coord': e.get('x_coord', e['x']),
                'Ue': e['Ue'],
                'theta': e['theta'],
                'delta_star': e['delta_star'],
                'Hk': e['Hk'],
                'Cf': e['Cf']
            })
    return sorted(stations, key=lambda s: s['ibl'])

def extract_xfoil_bl(trace_file, side=1):
    """Extract XFOIL BL station data."""
    with open(trace_file) as f:
        data = json.load(f)
    
    # XFOIL emits BLVAR for upper at early iterations, lower at final iteration
    # Use iteration 1 for upper surface (side=1), max iteration for lower (side=2)
    iter_values = [e.get('iteration', 0) for e in data['events'] 
                   if e.get('subroutine') == 'BLVAR' and e.get('side') == side]
    if not iter_values:
        return []
    
    # For upper, use iteration 1; for lower, use max iteration
    target_iter = 1 if side == 1 else max(iter_values)
    
    stations = []
    for e in data['events']:
        if e.get('subroutine') == 'BLVAR' and e.get('iteration') == target_iter and e.get('side') == side:
            inp = e.get('input', {})
            out = e.get('output', {})
            stations.append({
                'ibl': e.get('ibl', 0),
                'x': inp.get('x', 0),
                'Ue': inp.get('u', 0),
                'theta': inp.get('theta', 0),
                'delta_star': inp.get('delta_star', 0),
                'Hk': out.get('Hk', out.get('H', 0)),
                'Cf': out.get('Cf', 0)
            })
    return sorted(stations, key=lambda s: s['ibl'])

def compare_at_alpha(alpha_deg):
    """Compare RustFoil vs XFOIL at specific alpha."""
    print(f"\n{'='*100}")
    print(f"INTERMEDIATE VALUE COMPARISON AT α={alpha_deg}°")
    print('='*100)
    
    # Generate RustFoil debug trace
    import subprocess
    alpha_str = f"{alpha_deg:+.1f}".replace('+', '%2B').replace('.', '')
    debug_file = f"/Users/harry/flexfoil-boundary-layer/debug_a{alpha_deg}.json"
    
    result = subprocess.run([
        "cargo", "run", "--release", "--",
        "viscous", "naca0012_xfoil_paneled.dat",
        f"--alpha={alpha_deg}",
        "--re=3000000",
        "--no-repanel",
        f"--debug={debug_file}"
    ], cwd="/Users/harry/flexfoil-boundary-layer", capture_output=True, text=True, timeout=30)
    
    # Extract RustFoil stations
    rf_upper = extract_bl_stations(debug_file, side=1)
    rf_lower = extract_bl_stations(debug_file, side=2)
    
    # Extract XFOIL stations
    # Format: alpha_+004.0.json or alpha_-007.0.json
    alpha_sign = '+' if alpha_deg >= 0 else ''
    xf_trace = f"/Users/harry/flexfoil-boundary-layer/traces/xfoil/naca0012/re3e06/alpha_{alpha_sign}{abs(alpha_deg):03.0f}.0.json"
    if not Path(xf_trace).exists():
        print(f"Warning: XFOIL trace not found: {xf_trace}")
        print(f"Skipping comparison for α={alpha_deg}°")
        return
    
    xf_upper = extract_xfoil_bl(xf_trace, side=1)
    xf_lower = extract_xfoil_bl(xf_trace, side=2)
    
    # Compare key stations
    print(f"\nUPPER SURFACE - Key Stations:")
    print(f"{'Station':<10} {'Source':<10} {'x/c':>8} {'Ue':>10} {'theta':>12} {'H':>8} {'Cf':>10}")
    print('-'*75)
    
    # Match by approximate x location
    key_x = [0.01, 0.1, 0.5, 0.9, 1.0]
    for target_x in key_x:
        # Find closest RF station
        rf_closest = min(rf_upper, key=lambda s: abs(s['x_coord'] - target_x)) if rf_upper else None
        xf_closest = min(xf_upper, key=lambda s: abs(s['x'] - target_x)) if xf_upper else None
        
        if rf_closest:
            print(f"x≈{target_x:<7.2f} {'RustFoil':<10} {rf_closest['x_coord']:>8.6f} {rf_closest['Ue']:>10.6f} "
                  f"{rf_closest['theta']:>12.6e} {rf_closest['Hk']:>8.4f} {rf_closest['Cf']:>10.6f}")
        if xf_closest:
            print(f"{'':10} {'XFOIL':<10} {xf_closest['x']:>8.6f} {xf_closest['Ue']:>10.6f} "
                  f"{xf_closest['theta']:>12.6e} {xf_closest['Hk']:>8.4f} {xf_closest['Cf']:>10.6f}")
            if rf_closest:
                theta_err = (rf_closest['theta'] - xf_closest['theta']) / xf_closest['theta'] * 100
                ue_err = (rf_closest['Ue'] - xf_closest['Ue']) / xf_closest['Ue'] * 100
                print(f"{'':10} {'Δ%':<10} {'':>8} {ue_err:>10.1f}% {theta_err:>12.1f}%")
            print()
    
    # Summary
    print(f"\nSUMMARY FOR α={alpha_deg}°:")
    print(f"  RustFoil: {len(rf_upper)} upper stations, {len(rf_lower)} lower stations")
    print(f"  XFOIL:    {len(xf_upper)} upper stations, {len(xf_lower)} lower stations")
    
    # Extract CL/CD from CLI output
    for line in result.stdout.split('\n'):
        if 'CL:' in line or 'CD:' in line:
            print(f"  RustFoil {line.strip()}")

def main():
    # Compare at key angles
    angles = [4, 8, 12, -7]  # Good case, mid, high, anomaly
    
    for alpha in angles:
        try:
            compare_at_alpha(alpha)
        except Exception as e:
            print(f"Error at α={alpha}°: {e}")

if __name__ == '__main__':
    main()
