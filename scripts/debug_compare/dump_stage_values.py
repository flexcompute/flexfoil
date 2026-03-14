#!/usr/bin/env python3
"""
Stage-by-stage comparison of XFOIL vs RustFoil intermediate values.
Finds the FIRST divergence point.
"""

import json
import subprocess
import sys
from pathlib import Path
import numpy as np

WORKSPACE = Path(__file__).parent.parent.parent
XFOIL_BIN = WORKSPACE / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
AIRFOIL = WORKSPACE / "testdata" / "naca0012.dat"

def run_xfoil_dump(alpha=4.0, re=1e6):
    """Run XFOIL and dump intermediate values."""
    script = f"""LOAD {AIRFOIL}
PANE
OPER
VISC {re:.0f}
ITER 50
ALFA {alpha}

QUIT
"""
    script_path = WORKSPACE / "scripts" / "debug_compare" / "xfoil_script.txt"
    with open(script_path, 'w') as f:
        f.write(script)
    
    result = subprocess.run(
        [str(XFOIL_BIN)],
        stdin=open(script_path),
        capture_output=True,
        text=True,
        timeout=60,
        cwd=WORKSPACE
    )
    
    return result.stdout + result.stderr

def load_xfoil_debug():
    """Load XFOIL debug JSON."""
    path = WORKSPACE / "xfoil_debug.json"
    if not path.exists():
        return None
    with open(path) as f:
        data = json.load(f)
    return data.get('events', data)

def compare_value(name, xfoil, rustfoil, tol=0.01):
    """Compare two values with tolerance."""
    if xfoil is None or rustfoil is None:
        return False, "missing"
    
    diff = abs(xfoil - rustfoil)
    if abs(xfoil) > 1e-10:
        rel = diff / abs(xfoil)
    else:
        rel = diff
    
    passed = rel < tol
    return passed, f"xf={xfoil:.6e}, rf={rustfoil:.6e}, diff={rel*100:.2f}%"

def main():
    alpha = 4.0
    re = 1e6
    
    print("=" * 70)
    print(f"XFOIL vs RustFoil Stage-by-Stage Comparison")
    print(f"Alpha = {alpha}°, Re = {re:.0e}")
    print("=" * 70)
    
    # Run XFOIL
    print("\n[1] Running XFOIL...")
    xfoil_out = run_xfoil_dump(alpha, re)
    xfoil_events = load_xfoil_debug()
    
    if not xfoil_events:
        print("ERROR: No XFOIL debug output")
        return 1
    
    print(f"    Loaded {len(xfoil_events)} XFOIL events")
    
    # Extract key XFOIL values
    xfoil_data = {}
    for event in xfoil_events:
        sub = event.get('subroutine', '')
        
        if sub == 'VISCAL':
            xfoil_data['alpha_rad'] = event.get('alpha_rad')
            xfoil_data['reynolds'] = event.get('reynolds')
            
        elif sub == 'BLVAR' and event.get('iteration') == 1:
            side = event.get('side', 0)
            ibl = event.get('ibl', 0)
            key = f"blvar_{side}_{ibl}"
            if key not in xfoil_data:
                xfoil_data[key] = event
                
        elif sub == 'UPDATE':
            side = event.get('side', 0)
            ibl = event.get('ibl', 0)
            key = f"update_{side}_{ibl}"
            if key not in xfoil_data:
                xfoil_data[key] = {
                    'due': event.get('delta_Ue', 0),
                    'rlx': event.get('relaxation', 1),
                }
    
    # Check DIJ values by looking at UPDATE deltas
    print("\n[2] Comparing UPDATE dUe values (DIJ coupling result)...")
    
    # Get first few UPDATE events from each side
    for side in [1, 2]:
        surface = "upper" if side == 1 else "lower"
        print(f"\n    {surface.title()} surface:")
        
        for ibl in [2, 3, 4, 5, 10]:
            key = f"update_{side}_{ibl}"
            if key in xfoil_data:
                xf_due = xfoil_data[key]['due']
                print(f"      IBL {ibl}: dUe = {xf_due:+.6e}")
    
    print("\n[3] Key insight from XFOIL debug:")
    print("-" * 60)
    
    # Look for pattern in dUe
    upper_dues = []
    lower_dues = []
    for key, val in xfoil_data.items():
        if key.startswith('update_1_'):
            upper_dues.append(val['due'])
        elif key.startswith('update_2_'):
            lower_dues.append(val['due'])
    
    if upper_dues:
        print(f"    Upper dUe: min={min(upper_dues):.6e}, max={max(upper_dues):.6e}")
    if lower_dues:
        print(f"    Lower dUe: min={min(lower_dues):.6e}, max={max(lower_dues):.6e}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
