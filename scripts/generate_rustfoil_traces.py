#!/usr/bin/env python3
"""
Generate RustFoil debug traces for the full test matrix.

Test Matrix (same as XFOIL):
- NACA 0012: Re = 1e6, 3e6, 6e6
- NACA 2412: Re = 3e6
- NACA 4412: Re = 3e6
- Alpha: -15 to +15 degrees in 1 degree steps

Usage:
    python scripts/generate_rustfoil_traces.py [--output-dir traces/rustfoil] [--dry-run]
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional


def run_rustfoil_viscous(
    rustfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    reynolds: float,
    output_json: Path,
    ncrit: float = 9.0,
    mach: float = 0.0,
    timeout: int = 120,
) -> Optional[dict]:
    """
    Run RustFoil viscous for a single alpha/Re combination.
    
    Returns dict with results or None if failed.
    """
    cmd = [
        str(rustfoil_bin),
        "viscous",
        str(foil_file),
        f"--alpha={alpha}",  # Use = syntax for negative values
        f"--re={int(reynolds)}",
        f"--mach={mach}",
        f"--ncrit={ncrit}",
        f"--debug={output_json}",
        "--no-repanel",  # Use pre-paneled coordinates
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT at alpha={alpha}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    
    if result.returncode != 0:
        # Check if it's a convergence failure vs actual error
        if "Converged: false" in result.stdout:
            converged = False
        else:
            print(f"    ERROR: RustFoil failed with code {result.returncode}")
            return None
    else:
        converged = True
    
    # Parse output for CL/CD
    cl = None
    cd = None
    x_tr_upper = None
    x_tr_lower = None
    
    for line in result.stdout.split('\n'):
        line = line.strip()
        if line.startswith("CL:"):
            try:
                cl = float(line.split()[1])
            except (ValueError, IndexError):
                pass
        elif line.startswith("CD:"):
            try:
                cd = float(line.split()[1])
            except (ValueError, IndexError):
                pass
        elif line.startswith("x_tr (U):"):
            try:
                x_tr_upper = float(line.split()[2])
            except (ValueError, IndexError):
                pass
        elif line.startswith("x_tr (L):"):
            try:
                x_tr_lower = float(line.split()[2])
            except (ValueError, IndexError):
                pass
        elif line.startswith("Converged:"):
            converged = "true" in line.lower()
    
    # Add metadata to debug output
    if output_json.exists():
        try:
            with open(output_json) as f:
                debug_data = json.load(f)
            
            debug_data['metadata'] = {
                'foil': foil_file.stem,
                'alpha_deg': alpha,
                'reynolds': reynolds,
                'ncrit': ncrit,
                'mach': mach,
                'cl': cl,
                'cd': cd,
                'converged': converged,
            }
            
            with open(output_json, 'w') as f:
                json.dump(debug_data, f, indent=2)
        except (json.JSONDecodeError, IOError):
            pass
    
    return {
        'alpha': alpha,
        'reynolds': reynolds,
        'cl': cl,
        'cd': cd,
        'x_tr_upper': x_tr_upper,
        'x_tr_lower': x_tr_lower,
        'converged': converged,
        'output_file': str(output_json),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate RustFoil traces for test matrix"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("traces/rustfoil"),
        help="Output directory for traces (default: traces/rustfoil)"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Print what would be run without executing"
    )
    parser.add_argument(
        "--alpha-start",
        type=float,
        default=-15.0,
        help="Starting alpha (default: -15)"
    )
    parser.add_argument(
        "--alpha-end",
        type=float,
        default=15.0,
        help="Ending alpha (default: 15)"
    )
    parser.add_argument(
        "--alpha-step",
        type=float,
        default=1.0,
        help="Alpha step (default: 1)"
    )
    parser.add_argument(
        "--single-alpha",
        type=float,
        help="Run single alpha only (for debugging)"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    rustfoil_bin = workspace / "target" / "release" / "rustfoil"
    output_dir = workspace / args.output_dir
    
    if not rustfoil_bin.exists():
        print(f"ERROR: RustFoil binary not found at {rustfoil_bin}")
        print("Build it with: cargo build --release")
        return 1
    
    # Test matrix definition (same as XFOIL script)
    test_cases = [
        ("naca0012", workspace / "naca0012_xfoil_paneled.dat", [1e6, 3e6, 6e6]),
        ("naca2412", workspace / "naca2412_xfoil_paneled.dat", [3e6]),
        ("naca4412", workspace / "naca4412_xfoil_paneled.dat", [3e6]),
    ]
    
    # Alpha range
    if args.single_alpha is not None:
        alphas = [args.single_alpha]
    else:
        alphas = []
        alpha = args.alpha_start
        while alpha <= args.alpha_end + 1e-6:
            alphas.append(alpha)
            alpha += args.alpha_step
    
    print("=" * 70)
    print(" RUSTFOIL TRACE GENERATION")
    print("=" * 70)
    print(f"  Output directory: {output_dir}")
    print(f"  Alpha range: {args.alpha_start} to {args.alpha_end} step {args.alpha_step}")
    print(f"  Total alphas: {len(alphas)}")
    print()
    
    if args.dry_run:
        print("DRY RUN - would generate:")
        for foil_name, panel_file, re_list in test_cases:
            for re in re_list:
                print(f"  {foil_name} Re={re:.0e}: {len(alphas)} traces")
        return 0
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Summary results
    all_results = []
    
    for foil_name, panel_file, re_list in test_cases:
        print(f"\n{'='*60}")
        print(f" {foil_name.upper()}")
        print("=" * 60)
        
        if not panel_file.exists():
            print(f"  WARNING: Panel file not found: {panel_file}")
            continue
        
        for re in re_list:
            print(f"\n  Re = {re:.0e}")
            re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
            
            foil_results = []
            
            for alpha in alphas:
                output_file = output_dir / foil_name / re_str / f"alpha_{alpha:+06.1f}.json"
                output_file.parent.mkdir(parents=True, exist_ok=True)
                
                result = run_rustfoil_viscous(
                    rustfoil_bin, panel_file, alpha, re, output_file
                )
                
                if result:
                    foil_results.append(result)
                    status = "OK" if result['converged'] else "NC"
                    cl_str = f"CL={result['cl']:.4f}" if result['cl'] else "CL=?"
                    cd_str = f" CD={result['cd']:.5f}" if result.get('cd') else ""
                    print(f"    α={alpha:+6.1f}° {status} {cl_str}{cd_str}")
                else:
                    print(f"    α={alpha:+6.1f}° FAILED")
            
            all_results.extend(foil_results)
            
            # Save summary for this foil/Re
            summary_file = output_dir / foil_name / re_str / "summary.json"
            with open(summary_file, 'w') as f:
                json.dump({
                    'foil': foil_name,
                    'reynolds': re,
                    'results': foil_results,
                }, f, indent=2)
    
    # Overall summary
    print(f"\n{'='*70}")
    print(" SUMMARY")
    print("=" * 70)
    
    converged = sum(1 for r in all_results if r.get('converged', False))
    
    print(f"  Total traces generated: {len(all_results)}")
    print(f"  Converged: {converged}")
    print(f"  Not converged: {len(all_results) - converged}")
    print(f"  Output directory: {output_dir}")
    
    # Save master summary
    master_summary = {
        'test_matrix': [
            {'foil': f, 'reynolds': r, 'n_alphas': len(alphas)}
            for f, _, r_list in test_cases
            for r in r_list
        ],
        'alphas': alphas,
        'total_traces': len(all_results),
        'converged': converged,
    }
    
    with open(output_dir / "master_summary.json", 'w') as f:
        json.dump(master_summary, f, indent=2)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
