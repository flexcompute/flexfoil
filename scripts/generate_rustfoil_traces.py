#!/usr/bin/env python3
"""
Generate RustFoil debug traces for the full test matrix.

Test Matrix (same as XFOIL):
- NACA 0012: Re = 1e6, 3e6, 6e6
- NACA 2412: Re = 3e6
- NACA 4412: Re = 3e6
- Alpha: -15 to +15 degrees in 1 degree steps

Usage:
    python scripts/generate_rustfoil_traces.py [--output-dir traces/rustfoil_faithful] [--dry-run]
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

from sweep_matrix import get_sweep_cases


def run_rustfoil_viscous(
    rustfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    reynolds: float,
    output_json: Path,
    ncrit: float = 9.0,
    mach: float = 0.0,
    timeout: int = 120,
    no_repanel: bool = True,
    summary_only: bool = False,
) -> Optional[dict]:
    """
    Run RustFoil viscous for a single alpha/Re combination.
    
    Returns dict with results or None if failed.
    """
    cmd = [
        str(rustfoil_bin),
        "faithful-viscous",
        str(foil_file),
        f"--alpha={alpha}",  # Use = syntax for negative values
        f"--re={int(reynolds)}",
        f"--mach={mach}",
        f"--ncrit={ncrit}",
        "--format=json",
    ]
    if not summary_only:
        cmd.append(f"--debug={output_json}")
    if no_repanel:
        cmd.append("--no-repanel")
    
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
        print(f"    ERROR: RustFoil failed with code {result.returncode}")
        return None
    converged = True

    cl = None
    cd = None
    x_tr_upper = None
    x_tr_lower = None

    try:
        parsed = json.loads(result.stdout)
        cl = parsed.get("cl")
        cd = parsed.get("cd")
        x_tr_upper = parsed.get("x_tr_upper")
        x_tr_lower = parsed.get("x_tr_lower")
        converged = bool(parsed.get("converged", converged))
    except json.JSONDecodeError:
        print("    ERROR: Failed to parse faithful JSON output")
        return None
    
    metadata = {
        'foil': foil_file.stem,
        'alpha_deg': alpha,
        'reynolds': reynolds,
        'ncrit': ncrit,
        'mach': mach,
        'cl': cl,
        'cd': cd,
        'converged': converged,
        'solver_path': 'faithful-xfoil-side',
        'summary_only': summary_only,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    if summary_only:
        with open(output_json, 'w') as f:
            json.dump({'metadata': metadata}, f, indent=2)
    elif output_json.exists():
        try:
            with open(output_json) as f:
                debug_data = json.load(f)
            debug_data['metadata'] = metadata
            with open(output_json, 'w') as f:
                json.dump(debug_data, f, indent=2)
        except (json.JSONDecodeError, IOError):
            with open(output_json, 'w') as f:
                json.dump({'metadata': metadata}, f, indent=2)
    
    return {
        'alpha': alpha,
        'reynolds': reynolds,
        'cl': cl,
        'cd': cd,
        'x_tr_upper': x_tr_upper,
        'x_tr_lower': x_tr_lower,
        'converged': converged,
        'solver_path': 'faithful-xfoil-side',
        'output_file': str(output_json),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate RustFoil traces for test matrix"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("traces/rustfoil_faithful"),
        help="Output directory for traces (default: traces/rustfoil_faithful)"
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
        default=25.0,
        help="Ending alpha (default: 25)"
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
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Number of parallel RustFoil jobs (default: 8)"
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Store only per-alpha metadata instead of full debug event payloads"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    rustfoil_bin = workspace / "target" / "release" / "rustfoil"
    output_dir = workspace / args.output_dir
    
    if not rustfoil_bin.exists():
        print(f"ERROR: RustFoil binary not found at {rustfoil_bin}")
        print("Build it with: cargo build --release")
        return 1
    
    test_cases = get_sweep_cases(workspace)
    
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
    print(f"  Solver path: faithful-xfoil-side")
    print(f"  Parallel jobs: {args.jobs}")
    print(f"  Trace mode: {'summary-only' if args.summary_only else 'full-debug'}")
    print()
    
    if args.dry_run:
        print("DRY RUN - would generate:")
        for case in test_cases:
            for re in case["reynolds"]:
                print(f"  {case['foil']} Re={re:.0e}: {len(alphas)} traces")
        return 0
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tasks = []
    grouped_results: dict[tuple[str, str], list[dict]] = {}

    for case in test_cases:
        foil_name = case["foil"]
        panel_file = case["file"]
        re_list = case["reynolds"]
        if not panel_file.exists():
            print(f"WARNING: Panel file not found: {panel_file}")
            continue
        for re in re_list:
            re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
            grouped_results[(foil_name, re_str)] = []
            for alpha in alphas:
                output_file = output_dir / foil_name / re_str / f"alpha_{alpha:+06.1f}.json"
                output_file.parent.mkdir(parents=True, exist_ok=True)
                tasks.append((foil_name, panel_file, re, re_str, alpha, output_file, case["no_repanel"]))

    print(f"Submitting {len(tasks)} faithful RustFoil runs...")

    all_results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        future_map = {
            executor.submit(
                run_rustfoil_viscous,
                rustfoil_bin,
                panel_file,
                alpha,
                re,
                output_file,
                no_repanel=no_repanel,
                summary_only=args.summary_only,
            ): (foil_name, re, re_str, alpha)
            for foil_name, panel_file, re, re_str, alpha, output_file, no_repanel in tasks
        }
        for future in as_completed(future_map):
            foil_name, re, re_str, alpha = future_map[future]
            result = future.result()
            if result:
                grouped_results[(foil_name, re_str)].append(result)
                all_results.append(result)
                status = "OK" if result['converged'] else "NC"
                cl_str = f"CL={result['cl']:.4f}" if result['cl'] else "CL=?"
                cd_str = f" CD={result['cd']:.5f}" if result.get('cd') else ""
                print(f"[{foil_name} Re={re:.0e}] α={alpha:+6.1f}° {status} {cl_str}{cd_str}")
            else:
                print(f"[{foil_name} Re={re:.0e}] α={alpha:+6.1f}° FAILED")

    for case in test_cases:
        foil_name = case["foil"]
        re_list = case["reynolds"]
        for re in re_list:
            re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
            foil_results = sorted(grouped_results.get((foil_name, re_str), []), key=lambda row: row["alpha"])
            if not foil_results:
                continue
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
            {'foil': case['foil'], 'reynolds': r, 'n_alphas': len(alphas)}
            for case in test_cases
            for r in case['reynolds']
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
