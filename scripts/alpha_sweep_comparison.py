#!/usr/bin/env python3
"""
Run RustFoil viscous analysis over a range of angles and compare with XFOIL data.
"""
import json
import subprocess
import sys
from pathlib import Path

def run_rustfoil_alpha(alpha_deg, re=3e6, airfoil="naca0012_xfoil_paneled.dat"):
    """Run RustFoil at a specific angle of attack."""
    cmd = [
        "cargo", "run", "--release", "--",
        "viscous", airfoil,
        f"--alpha={alpha_deg}",
        f"--re={int(re)}",
        "--no-repanel"
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            cwd="/Users/harry/flexfoil-boundary-layer"
        )
        
        # Parse output for CL and CD
        cl, cd, cm = None, None, None
        converged = False
        
        for line in result.stdout.split('\n'):
            if line.startswith('CL:'):
                cl = float(line.split()[1])
            elif line.startswith('CD:'):
                cd = float(line.split()[1])
            elif line.startswith('CM:'):
                cm = float(line.split()[1])
            elif 'Converged:' in line:
                converged = 'true' in line.lower()
        
        return {
            'alpha': alpha_deg,
            'cl': cl,
            'cd': cd,
            'cm': cm,
            'converged': converged,
            'success': result.returncode == 0
        }
    except subprocess.TimeoutExpired:
        return {
            'alpha': alpha_deg,
            'cl': None,
            'cd': None,
            'cm': None,
            'converged': False,
            'success': False,
            'error': 'timeout'
        }
    except Exception as e:
        return {
            'alpha': alpha_deg,
            'cl': None,
            'cd': None,
            'cm': None,
            'converged': False,
            'success': False,
            'error': str(e)
        }

def main():
    # Alpha range: -15 to 15 in 1° steps
    alphas = range(-15, 16, 1)
    
    results = []
    
    print(f"Running RustFoil sweep from α={min(alphas)}° to {max(alphas)}° in 1° steps...")
    print(f"{'Alpha':>6} {'CL':>10} {'CD':>10} {'CM':>10} {'Conv':>6}")
    print("-" * 50)
    
    for alpha in alphas:
        result = run_rustfoil_alpha(alpha)
        results.append(result)
        
        if result['success']:
            conv_str = "✓" if result['converged'] else "✗"
            print(f"{alpha:>6.1f} {result['cl']:>10.4f} {result['cd']:>10.6f} "
                  f"{result['cm']:>10.6f} {conv_str:>6}")
        else:
            print(f"{alpha:>6.1f} {'FAILED':>10} {'':>10} {'':>10} {'✗':>6}")
    
    # Save results
    output_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/rustfoil_alpha_sweep.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'airfoil': 'naca0012',
            'reynolds': 3e6,
            'alpha_range': list(alphas),
            'results': results
        }, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    
    # Summary statistics
    successful = [r for r in results if r['success']]
    converged = [r for r in successful if r['converged']]
    
    print(f"\nSummary:")
    print(f"  Total runs: {len(results)}")
    print(f"  Successful: {len(successful)}")
    print(f"  Converged: {len(converged)}")
    
    if successful:
        cl_values = [r['cl'] for r in successful if r['cl'] is not None]
        if cl_values:
            print(f"  CL range: {min(cl_values):.4f} to {max(cl_values):.4f}")

if __name__ == '__main__':
    main()
