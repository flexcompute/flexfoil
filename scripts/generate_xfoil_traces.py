#!/usr/bin/env python3
"""
Generate instrumented XFOIL traces for the full test matrix.

Test Matrix:
- NACA 0012: Re = 1e6, 3e6, 6e6
- NACA 2412: Re = 3e6
- NACA 4412: Re = 3e6
- Alpha: -15 to +15 degrees in 1 degree steps

Usage:
    python scripts/generate_xfoil_traces.py [--output-dir traces] [--dry-run]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


def run_xfoil_viscous(
    xfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    reynolds: float,
    output_json: Path,
    ncrit: float = 9.0,
    mach: float = 0.0,
    timeout: int = 60,
) -> Optional[dict]:
    """
    Run instrumented XFOIL for a single alpha/Re combination.
    
    Returns dict with results or None if failed.
    """
    # Build XFOIL command script
    xfoil_commands = f"""LOAD {foil_file}
PANE
OPER
VISC {reynolds:.0f}
MACH {mach}
VPAR
N {ncrit}

ALFA {alpha}

QUIT
"""
    
    # Run XFOIL from the bin directory (so xfoil_debug.json is created there)
    work_dir = xfoil_bin.parent
    
    try:
        result = subprocess.run(
            [str(xfoil_bin)],
            input=xfoil_commands,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=work_dir,
        )
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT at alpha={alpha}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    
    # Check for debug output
    debug_file = work_dir / "xfoil_debug.json"
    if not debug_file.exists():
        print(f"    No debug output generated")
        return None
    
    # Read and parse debug output
    try:
        with open(debug_file) as f:
            debug_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"    JSON parse error: {e}")
        return None
    
    # Add metadata
    debug_data['metadata'] = {
        'foil': foil_file.stem,
        'alpha_deg': alpha,
        'reynolds': reynolds,
        'ncrit': ncrit,
        'mach': mach,
    }
    
    # Save to output location
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(debug_data, f, indent=2)
    
    # Extract key results from XFOIL output
    cl = None
    cd = None
    for line in result.stdout.split('\n'):
        if 'CL =' in line:
            try:
                parts = line.split()
                cl_idx = parts.index('CL') + 2
                cd_idx = parts.index('CD') + 2 if 'CD' in parts else None
                cl = float(parts[cl_idx])
                if cd_idx:
                    cd = float(parts[cd_idx])
            except (ValueError, IndexError):
                pass
    
    return {
        'alpha': alpha,
        'reynolds': reynolds,
        'cl': cl,
        'cd': cd,
        'converged': cl is not None,
        'output_file': str(output_json),
    }


def run_xfoil_inviscid(
    xfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    output_json: Path,
    timeout: int = 30,
) -> Optional[dict]:
    """
    Run instrumented XFOIL in inviscid mode for a single alpha.
    
    Returns dict with results or None if failed.
    """
    xfoil_commands = f"""LOAD {foil_file}
PANE
OPER
ALFA {alpha}

QUIT
"""
    
    work_dir = xfoil_bin.parent
    
    try:
        result = subprocess.run(
            [str(xfoil_bin)],
            input=xfoil_commands,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=work_dir,
        )
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT at alpha={alpha}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    
    debug_file = work_dir / "xfoil_debug.json"
    if not debug_file.exists():
        print(f"    No debug output generated")
        return None
    
    try:
        with open(debug_file) as f:
            debug_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"    JSON parse error: {e}")
        return None
    
    debug_data['metadata'] = {
        'foil': foil_file.stem,
        'alpha_deg': alpha,
        'inviscid': True,
    }
    
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(debug_data, f, indent=2)
    
    # Extract CL from output
    cl = None
    for line in result.stdout.split('\n'):
        if 'CL =' in line:
            try:
                parts = line.split()
                cl_idx = parts.index('CL') + 2
                cl = float(parts[cl_idx])
            except (ValueError, IndexError):
                pass
    
    return {
        'alpha': alpha,
        'cl': cl,
        'output_file': str(output_json),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate instrumented XFOIL traces for test matrix"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("traces/xfoil"),
        help="Output directory for traces (default: traces/xfoil)"
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
        "--inviscid-only",
        action="store_true",
        help="Run inviscid analysis only (no viscous)"
    )
    parser.add_argument(
        "--single-alpha",
        type=float,
        help="Run single alpha only (for debugging)"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    xfoil_bin = workspace / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
    output_dir = workspace / args.output_dir
    
    if not xfoil_bin.exists():
        print(f"ERROR: Instrumented XFOIL not found at {xfoil_bin}")
        print("Build it with: cd Xfoil-instrumented/bin && make")
        return 1
    
    # Test matrix definition
    test_cases = [
        # (foil_name, panel_file, reynolds_list)
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
    print(" XFOIL INSTRUMENTED TRACE GENERATION")
    print("=" * 70)
    print(f"  Output directory: {output_dir}")
    print(f"  Alpha range: {args.alpha_start} to {args.alpha_end} step {args.alpha_step}")
    print(f"  Total alphas: {len(alphas)}")
    print(f"  Mode: {'Inviscid only' if args.inviscid_only else 'Viscous'}")
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
        
        if args.inviscid_only:
            re_list = [None]  # Single inviscid run
        
        for re in re_list:
            if re is not None:
                print(f"\n  Re = {re:.0e}")
                re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
            else:
                print(f"\n  Inviscid")
                re_str = "inviscid"
            
            foil_results = []
            
            for alpha in alphas:
                output_file = output_dir / foil_name / re_str / f"alpha_{alpha:+06.1f}.json"
                
                if re is not None:
                    result = run_xfoil_viscous(
                        xfoil_bin, panel_file, alpha, re, output_file
                    )
                else:
                    result = run_xfoil_inviscid(
                        xfoil_bin, panel_file, alpha, output_file
                    )
                
                if result:
                    foil_results.append(result)
                    status = "OK" if result.get('converged', True) else "NC"
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
    
    total_cases = sum(len(r) * len(a) for _, _, r in test_cases for a in [alphas])
    converged = sum(1 for r in all_results if r.get('converged', True))
    
    print(f"  Total traces generated: {len(all_results)}")
    print(f"  Converged: {converged}")
    print(f"  Failed: {len(all_results) - converged}")
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
