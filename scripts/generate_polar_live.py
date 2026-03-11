#!/usr/bin/env python3
"""Generate live CL polar comparison between XFOIL and RustFoil."""

import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os
import re as regex

XFOIL_PATH = '/Users/harry/flexfoil-boundary-layer/Xfoil/bin/xfoil'
RUSTFOIL_DIR = '/Users/harry/flexfoil-boundary-layer'

def run_xfoil_polar(naca="0012", alphas=None, re=1e6):
    """Run XFOIL for multiple alpha values and get CL, CD."""
    if alphas is None:
        alphas = list(range(-6, 17, 2))
    
    results = {}
    
    # Generate airfoil and save it for RustFoil
    save_script = f"""NACA {naca}
PSAV /tmp/naca{naca}.dat

QUIT
"""
    subprocess.run([XFOIL_PATH], input=save_script, capture_output=True, text=True, timeout=10)
    
    for alpha in alphas:
        script = f"""NACA {naca}
PANE
OPER
VISC {re:.0f}
ITER 200
ALFA {alpha}

QUIT
"""
        try:
            result = subprocess.run(
                [XFOIL_PATH],
                input=script,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            # Parse output - get the last converged CL/CD
            cl, cd = None, None
            for line in result.stdout.split('\n'):
                if 'CL =' in line:
                    match = regex.search(r'CL\s*=\s*([-\d.]+)', line)
                    if match:
                        cl = float(match.group(1))
                if 'CD =' in line:
                    match = regex.search(r'CD\s*=\s*([-\d.]+)', line)
                    if match:
                        cd = float(match.group(1))
            
            if cl is not None and cd is not None:
                results[alpha] = (cl, cd)
        except Exception as e:
            print(f"XFOIL error at alpha={alpha}: {e}")
    
    return results

def run_rustfoil(airfoil_file, alpha, re=1e6):
    """Run RustFoil and get CL, CD."""
    try:
        result = subprocess.run(
            ['cargo', 'run', '--release', '-q', '-p', 'rustfoil-cli', '--',
             'viscous', airfoil_file, f'--alpha={alpha}', f'--re={re:.0f}', '--format', 'json'],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=RUSTFOIL_DIR
        )
        
        if result.returncode == 0:
            data = json.loads(result.stdout)
            return data.get('cl'), data.get('cd')
    except Exception as e:
        print(f"RustFoil error at alpha={alpha}: {e}")
    return None, None

def main():
    # Run both symmetric and cambered airfoils
    airfoils = ["0012", "4412"]
    alphas = list(range(-25, 26, 2))  # -25 to +25 in 2° increments
    re = 1e6
    
    for naca in airfoils:
        run_polar_for_airfoil(naca, alphas, re)

def run_polar_for_airfoil(naca, alphas, re):
    
    print(f"Running polar sweep for NACA {naca} at Re={re:.0e}")
    print("=" * 80)
    
    # Run XFOIL for all alphas
    print("Running XFOIL...")
    xfoil_results = run_xfoil_polar(naca, alphas, re)
    
    airfoil_file = f"/tmp/naca{naca}.dat"
    
    # Run RustFoil for all alphas  
    print("Running RustFoil...")
    rustfoil_results = {}
    for alpha in alphas:
        cl, cd = run_rustfoil(airfoil_file, alpha, re)
        if cl is not None:
            rustfoil_results[alpha] = (cl, cd)
    
    # Print comparison table
    print()
    print(f"{'Alpha':>6} | {'XFOIL CL':>10} {'XFOIL CD':>10} | {'RustFoil CL':>12} {'RustFoil CD':>12} | {'CL Err%':>8}")
    print("-" * 85)
    
    valid_alphas = []
    xfoil_cl_v, xfoil_cd_v = [], []
    rustfoil_cl_v, rustfoil_cd_v = [], []
    
    for alpha in alphas:
        xcl, xcd = xfoil_results.get(alpha, (None, None))
        rcl, rcd = rustfoil_results.get(alpha, (None, None))
        
        xcl_s = f"{xcl:.4f}" if xcl is not None else "FAIL"
        xcd_s = f"{xcd:.5f}" if xcd is not None else "FAIL"
        rcl_s = f"{rcl:.4f}" if rcl is not None else "FAIL"
        rcd_s = f"{rcd:.5f}" if rcd is not None else "FAIL"
        
        if xcl is not None and rcl is not None:
            if abs(xcl) > 0.01:
                err = (rcl - xcl) / abs(xcl) * 100
            else:
                err = (rcl - xcl) * 100
            err_s = f"{err:+.1f}%"
            
            valid_alphas.append(alpha)
            xfoil_cl_v.append(xcl)
            xfoil_cd_v.append(xcd)
            rustfoil_cl_v.append(rcl)
            rustfoil_cd_v.append(rcd)
        else:
            err_s = "-"
        
        print(f"{alpha:6.1f} | {xcl_s:>10} {xcd_s:>10} | {rcl_s:>12} {rcd_s:>12} | {err_s:>8}")
    
    if not valid_alphas:
        print("\nNo valid data points!")
        return
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # CL vs Alpha
    ax1.plot(valid_alphas, xfoil_cl_v, 'b-o', label='XFOIL', linewidth=2, markersize=6)
    ax1.plot(valid_alphas, rustfoil_cl_v, 'r--s', label='RustFoil', linewidth=2, markersize=6)
    ax1.set_xlabel('Angle of Attack (deg)', fontsize=12)
    ax1.set_ylabel('Lift Coefficient CL', fontsize=12)
    ax1.set_title(f'CL vs Alpha - NACA {naca}, Re={re/1e6:.0f}M', fontsize=14)
    ax1.legend(loc='lower right', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax1.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
    
    # CD vs Alpha
    ax2.plot(valid_alphas, [cd * 10000 for cd in xfoil_cd_v], 'b-o', label='XFOIL', linewidth=2, markersize=6)
    ax2.plot(valid_alphas, [cd * 10000 for cd in rustfoil_cd_v], 'r--s', label='RustFoil', linewidth=2, markersize=6)
    ax2.set_xlabel('Angle of Attack (deg)', fontsize=12)
    ax2.set_ylabel('Drag Coefficient CD (counts)', fontsize=12)
    ax2.set_title(f'CD vs Alpha - NACA {naca}, Re={re/1e6:.0f}M', fontsize=14)
    ax2.legend(loc='upper left', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f'{RUSTFOIL_DIR}/scripts/cl_polar_naca{naca}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved to {output_file}")

if __name__ == "__main__":
    main()
