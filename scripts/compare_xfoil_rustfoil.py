#!/usr/bin/env python3
"""
Compare XFOIL and RustFoil viscous results across a range of AoA.
"""

import subprocess
import json
import sys
import re
from pathlib import Path

# XFOIL reference data from xfoil_polar_full.txt
# Format: alpha, CL, CD, CDp, CM, xtr_upper, xtr_lower
XFOIL_DATA = [
    (-4.0, -0.4278, 0.007279, 0.001183, -0.006024, 0.9685, 0.2537),
    (-3.0, -0.3199, 0.006390, 0.000869, -0.004852, 0.9285, 0.3642),
    (-2.0, -0.2142, 0.005801, 0.000642, -0.003007, 0.8676, 0.4743),
    (-1.0, -0.1074, 0.005488, 0.000504, -0.001428, 0.7848, 0.5826),
    ( 0.0,  0.0000, 0.005405, 0.000457,  0.000000, 0.6870, 0.6870),
    ( 1.0,  0.1074, 0.005488, 0.000504,  0.001428, 0.5826, 0.7848),
    ( 2.0,  0.2142, 0.005801, 0.000642,  0.003005, 0.4743, 0.8676),
    ( 3.0,  0.3200, 0.006390, 0.000870,  0.004849, 0.3642, 0.9285),
    ( 4.0,  0.4278, 0.007279, 0.001184,  0.006020, 0.2537, 0.9685),
    ( 5.0,  0.5580, 0.008479, 0.001653,  0.001721, 0.1485, 0.9849),
    ( 6.0,  0.6948, 0.009728, 0.002238, -0.004279, 0.0812, 0.9940),
    ( 7.0,  0.8264, 0.010940, 0.002893, -0.009233, 0.0506, 1.0000),
    ( 8.0,  0.9098, 0.012115, 0.003558, -0.003939, 0.0381, 1.0000),
    ( 9.0,  0.9947, 0.013408, 0.004343,  0.000976, 0.0307, 1.0000),
    (10.0,  1.0809, 0.014977, 0.005334,  0.005269, 0.0256, 1.0000),
]

def run_rustfoil(alpha_deg):
    """Run RustFoil at a specific angle of attack."""
    # Run the test and capture output
    cmd = [
        "cargo", "test", 
        "--package", "rustfoil-solver",
        "test_viscous_at_alpha",
        "--",
        "--nocapture",
        "--test-threads=1"
    ]
    
    # Set environment variable for alpha
    env = {"RUSTFOIL_TEST_ALPHA": str(alpha_deg)}
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parent.parent,
        env={**subprocess.os.environ, **env},
        timeout=60
    )
    
    # Parse output for CL, CD values
    output = result.stdout + result.stderr
    
    # Look for output like "CL = 0.2417" etc
    cl_match = re.search(r'CL\s*=\s*([-\d.]+)', output)
    cd_match = re.search(r'CD\s*=\s*([-\d.]+)', output)
    cf_match = re.search(r'Cf\s*=\s*([-\d.]+)', output)
    cp_match = re.search(r'Cp\s*=\s*([-\d.]+)', output)
    xtr_upper_match = re.search(r'x_tr upper\s*=\s*([-\d.]+)', output)
    xtr_lower_match = re.search(r'x_tr lower\s*=\s*([-\d.]+)', output)
    
    if cl_match and cd_match:
        return {
            'cl': float(cl_match.group(1)),
            'cd': float(cd_match.group(1)),
            'cf': float(cf_match.group(1)) if cf_match else None,
            'cp': float(cp_match.group(1)) if cp_match else None,
            'xtr_upper': float(xtr_upper_match.group(1)) if xtr_upper_match else None,
            'xtr_lower': float(xtr_lower_match.group(1)) if xtr_lower_match else None,
        }
    
    return None

def main():
    print("=" * 80)
    print("XFOIL vs RustFoil Comparison - NACA 0012, Re=1e6, Ncrit=9")
    print("=" * 80)
    print()
    
    # Header
    print(f"{'Alpha':>6} | {'XFOIL CL':>10} {'RF CL':>10} {'Err%':>7} | "
          f"{'XFOIL CD':>10} {'RF CD':>10} {'Ratio':>7} | "
          f"{'xtr_u XF':>8} {'xtr_u RF':>8}")
    print("-" * 100)
    
    results = []
    
    for alpha, xf_cl, xf_cd, xf_cdp, xf_cm, xf_xtr_u, xf_xtr_l in XFOIL_DATA:
        rf = run_rustfoil(alpha)
        
        if rf:
            cl_err = ((rf['cl'] - xf_cl) / max(abs(xf_cl), 0.01)) * 100 if xf_cl != 0 else rf['cl'] * 100
            cd_ratio = rf['cd'] / xf_cd
            
            print(f"{alpha:>6.1f} | {xf_cl:>10.4f} {rf['cl']:>10.4f} {cl_err:>6.1f}% | "
                  f"{xf_cd:>10.6f} {rf['cd']:>10.6f} {cd_ratio:>6.2f}x | "
                  f"{xf_xtr_u:>8.4f} {rf['xtr_upper'] or 0:>8.4f}")
            
            results.append({
                'alpha': alpha,
                'xfoil': {'cl': xf_cl, 'cd': xf_cd, 'cdp': xf_cdp, 'xtr_u': xf_xtr_u, 'xtr_l': xf_xtr_l},
                'rustfoil': rf,
                'cl_err_pct': cl_err,
                'cd_ratio': cd_ratio,
            })
        else:
            print(f"{alpha:>6.1f} | {xf_cl:>10.4f} {'FAIL':>10} {'---':>7} | "
                  f"{xf_cd:>10.6f} {'FAIL':>10} {'---':>7}")
    
    print("-" * 100)
    
    # Summary statistics
    if results:
        cd_ratios = [r['cd_ratio'] for r in results]
        cl_errs = [r['cl_err_pct'] for r in results]
        
        print()
        print("Summary:")
        print(f"  CD ratio: min={min(cd_ratios):.2f}x, max={max(cd_ratios):.2f}x, avg={sum(cd_ratios)/len(cd_ratios):.2f}x")
        print(f"  CL error: min={min(cl_errs):.1f}%, max={max(cl_errs):.1f}%, avg={sum(cl_errs)/len(cl_errs):.1f}%")
        
        # Target: CD ratio should be 1.0-1.1x
        if max(cd_ratios) < 1.15:
            print("  STATUS: GOOD - CD ratio within 15% of XFOIL")
        elif max(cd_ratios) < 1.5:
            print("  STATUS: NEEDS WORK - CD ratio within 50% of XFOIL")
        else:
            print("  STATUS: SIGNIFICANT GAP - CD ratio > 50% higher than XFOIL")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
