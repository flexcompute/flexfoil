#!/usr/bin/env python3
"""
Compare XFOIL and RustFoil transition at alpha=10°
"""
import subprocess
import json
import tempfile
import os

def run_xfoil(alpha, re):
    """Run XFOIL and extract transition location"""
    xfoil_script = f"""
NACA 0012
OPER
v {re}
iter 200
pacc
polar.txt
alfa {alpha}
quit
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(xfoil_script)
        script_file = f.name
    
    try:
        result = subprocess.run(
            ['./xfoil'],
            input=xfoil_script,
            capture_output=True,
            text=True,
            cwd='/Users/harry/flexfoil-boundary-layer/Xfoil/bin',
            timeout=30
        )
        
        # Parse output for transition
        lines = result.stdout.split('\n')
        xtr_upper = None
        xtr_lower = None
        for line in lines:
            if 'top' in line.lower() and 'xtr' in line.lower():
                parts = line.split()
                try:
                    xtr_upper = float([p for p in parts if '.' in p][0])
                except:
                    pass
            if 'bottom' in line.lower() and 'xtr' in line.lower():
                parts = line.split()
                try:
                    xtr_lower = float([p for p in parts if '.' in p][0])
                except:
                    pass
        
        # Also try to read from polar file
        polar_file = '/Users/harry/flexfoil-boundary-layer/Xfoil/bin/polar.txt'
        if os.path.exists(polar_file):
            with open(polar_file, 'r') as pf:
                for line in pf:
                    if line.strip() and not line.startswith('---'):
                        parts = line.split()
                        if len(parts) >= 7:
                            try:
                                xtr_upper = float(parts[4])
                                xtr_lower = float(parts[5])
                            except:
                                pass
        
        return xtr_upper, xtr_lower
    finally:
        if os.path.exists(script_file):
            os.unlink(script_file)

def run_rustfoil(alpha):
    """Run RustFoil test and extract transition"""
    env = os.environ.copy()
    env['RUSTFOIL_TEST_ALPHA'] = str(int(alpha))
    env['RUSTFOIL_DEBUG_TRACE'] = '0'
    
    result = subprocess.run(
        ['cargo', 'test', '--package', 'rustfoil-solver', '--test', 
         'xfoil_viscous_comparison', 'test_newton_iteration_trace', '--', '--nocapture'],
        capture_output=True,
        text=True,
        cwd='/Users/harry/flexfoil-boundary-layer',
        env=env,
        timeout=60
    )
    
    # Parse output for transition
    lines = result.stdout.split('\n')
    xtr_upper = None
    xtr_lower = None
    for line in lines:
        if 'x_tr upper' in line:
            parts = line.split('=')
            if len(parts) >= 2:
                try:
                    xtr_upper = float(parts[1].strip())
                except:
                    pass
        if 'x_tr lower' in line:
            parts = line.split('=')
            if len(parts) >= 2:
                try:
                    xtr_lower = float(parts[1].strip())
                except:
                    pass
    
    return xtr_upper, xtr_lower

if __name__ == '__main__':
    alpha = 10
    re = 20000000
    
    print(f"Comparing transition at α = {alpha}°, Re = {re:.0e}")
    print("="*60)
    
    print("\nRunning RustFoil...")
    rust_upper, rust_lower = run_rustfoil(alpha)
    print(f"  RustFoil upper: {rust_upper}")
    print(f"  RustFoil lower: {rust_lower}")
    
    print("\nRunning XFOIL...")
    xfoil_upper, xfoil_lower = run_xfoil(alpha, re)
    print(f"  XFOIL upper: {xfoil_upper}")
    print(f"  XFOIL lower: {xfoil_lower}")
    
    if rust_upper and xfoil_upper:
        error = abs(rust_upper - xfoil_upper)
        error_pct = 100 * error / xfoil_upper
        print(f"\nUpper surface transition:")
        print(f"  XFOIL:    x/c = {xfoil_upper:.4f}")
        print(f"  RustFoil: x/c = {rust_upper:.4f}")
        print(f"  Error:    {error:.4f} ({error_pct:.1f}%)")
        
        if error_pct > 5:
            print(f"\n❌ ERROR > 5%! Transition prediction is significantly wrong!")
        elif error_pct > 1:
            print(f"\n⚠️  Error > 1%. Transition prediction needs improvement.")
        else:
            print(f"\n✓ Error < 1%. Transition prediction is acceptable.")
