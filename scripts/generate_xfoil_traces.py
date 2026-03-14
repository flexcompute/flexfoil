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
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

from sweep_matrix import get_sweep_cases


def extract_final_coefficients(debug_data: dict, stdout: str) -> tuple[Optional[float], Optional[float]]:
    cl = None
    cd = None

    metadata = debug_data.get("metadata", {})
    if metadata.get("cl") is not None:
        cl = metadata.get("cl")
    if metadata.get("cd") is not None:
        cd = metadata.get("cd")

    for event in reversed(debug_data.get("events", [])):
        if cl is None:
            cl = event.get("CL", event.get("cl"))
        if cd is None:
            cd = event.get("CD", event.get("cd"))
        if cl is not None and cd is not None:
            break

    if cl is None or cd is None:
        for line in stdout.split('\n'):
            if 'CL =' in line:
                try:
                    parts = line.split()
                    cl_idx = parts.index('CL') + 2
                    cl = float(parts[cl_idx]) if cl is None else cl
                    if 'CD' in parts:
                        cd_idx = parts.index('CD') + 2
                        cd = float(parts[cd_idx]) if cd is None else cd
                except (ValueError, IndexError):
                    pass

    return cl, cd


def run_xfoil_viscous(
    xfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    reynolds: float,
    output_json: Path,
    ncrit: float = 9.0,
    mach: float = 0.0,
    timeout: int = 60,
    summary_only: bool = False,
) -> Optional[dict]:
    """
    Run instrumented XFOIL for a single alpha/Re combination.
    
    Returns dict with results or None if failed.
    """
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
    try:
        with tempfile.TemporaryDirectory(prefix="xfoil-trace-") as tmpdir:
            work_dir = Path(tmpdir)
            result = subprocess.run(
                [str(xfoil_bin)],
                input=xfoil_commands,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=work_dir,
            )
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
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT at alpha={alpha}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    
    metadata = {
        'foil': foil_file.stem,
        'alpha_deg': alpha,
        'reynolds': reynolds,
        'ncrit': ncrit,
        'mach': mach,
    }

    cl, cd = extract_final_coefficients(debug_data, result.stdout)

    metadata.update({
        'cl': cl,
        'cd': cd,
        'converged': cl is not None,
        'solver_path': 'xfoil-instrumented',
        'summary_only': summary_only,
    })

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_payload = {'metadata': metadata} if summary_only else debug_data
    if not summary_only:
        output_payload['metadata'] = metadata
    with open(output_json, 'w') as f:
        json.dump(output_payload, f, indent=2)
    
    return {
        'alpha': alpha,
        'reynolds': reynolds,
        'cl': cl,
        'cd': cd,
        'converged': cl is not None,
        'solver_path': 'xfoil-instrumented',
        'output_file': str(output_json),
    }


def run_xfoil_inviscid(
    xfoil_bin: Path,
    foil_file: Path,
    alpha: float,
    output_json: Path,
    timeout: int = 30,
    summary_only: bool = False,
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
    
    try:
        with tempfile.TemporaryDirectory(prefix="xfoil-inviscid-") as tmpdir:
            work_dir = Path(tmpdir)
            result = subprocess.run(
                [str(xfoil_bin)],
                input=xfoil_commands,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=work_dir,
            )
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
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT at alpha={alpha}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    
    metadata = {
        'foil': foil_file.stem,
        'alpha_deg': alpha,
        'inviscid': True,
        'solver_path': 'xfoil-instrumented',
        'summary_only': summary_only,
    }

    cl, _ = extract_final_coefficients(debug_data, result.stdout)

    metadata['cl'] = cl

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_payload = {'metadata': metadata} if summary_only else debug_data
    if not summary_only:
        output_payload['metadata'] = metadata
    with open(output_json, 'w') as f:
        json.dump(output_payload, f, indent=2)
    
    return {
        'alpha': alpha,
        'cl': cl,
        'solver_path': 'xfoil-instrumented',
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
        "--inviscid-only",
        action="store_true",
        help="Run inviscid analysis only (no viscous)"
    )
    parser.add_argument(
        "--single-alpha",
        type=float,
        help="Run single alpha only (for debugging)"
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Number of parallel XFOIL jobs (default: 4)"
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Store only per-alpha metadata instead of full debug event payloads"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    xfoil_bin = workspace / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
    output_dir = workspace / args.output_dir
    
    if not xfoil_bin.exists():
        print(f"ERROR: Instrumented XFOIL not found at {xfoil_bin}")
        print("Build it with: cd Xfoil-instrumented/bin && make")
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
    print(" XFOIL INSTRUMENTED TRACE GENERATION")
    print("=" * 70)
    print(f"  Output directory: {output_dir}")
    print(f"  Alpha range: {args.alpha_start} to {args.alpha_end} step {args.alpha_step}")
    print(f"  Total alphas: {len(alphas)}")
    print(f"  Mode: {'Inviscid only' if args.inviscid_only else 'Viscous'}")
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
        effective_re_list = [None] if args.inviscid_only else re_list
        for re in effective_re_list:
            re_str = "inviscid" if re is None else f"re{re:.0e}".replace("+", "").replace(".", "")
            grouped_results[(foil_name, re_str)] = []
            for alpha in alphas:
                output_file = output_dir / foil_name / re_str / f"alpha_{alpha:+06.1f}.json"
                tasks.append((foil_name, panel_file, re, re_str, alpha, output_file))

    print(f"Submitting {len(tasks)} XFOIL runs...")

    all_results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        future_map = {}
        for foil_name, panel_file, re, re_str, alpha, output_file in tasks:
            if re is not None:
                future = executor.submit(
                    run_xfoil_viscous,
                    xfoil_bin,
                    panel_file,
                    alpha,
                    re,
                    output_file,
                    summary_only=args.summary_only,
                )
            else:
                future = executor.submit(
                    run_xfoil_inviscid,
                    xfoil_bin,
                    panel_file,
                    alpha,
                    output_file,
                    summary_only=args.summary_only,
                )
            future_map[future] = (foil_name, re, re_str, alpha)

        for future in as_completed(future_map):
            foil_name, re, re_str, alpha = future_map[future]
            result = future.result()
            if result:
                grouped_results[(foil_name, re_str)].append(result)
                all_results.append(result)
                status = "OK" if result.get('converged', True) else "NC"
                cl_str = f"CL={result['cl']:.4f}" if result['cl'] else "CL=?"
                cd_str = f" CD={result['cd']:.5f}" if result.get('cd') else ""
                re_label = "inviscid" if re is None else f"{re:.0e}"
                print(f"[{foil_name} Re={re_label}] α={alpha:+6.1f}° {status} {cl_str}{cd_str}")
            else:
                re_label = "inviscid" if re is None else f"{re:.0e}"
                print(f"[{foil_name} Re={re_label}] α={alpha:+6.1f}° FAILED")

    for case in test_cases:
        foil_name = case["foil"]
        panel_file = case["file"]
        re_list = case["reynolds"]
        if not panel_file.exists():
            continue
        effective_re_list = [None] if args.inviscid_only else re_list
        for re in effective_re_list:
            re_str = "inviscid" if re is None else f"re{re:.0e}".replace("+", "").replace(".", "")
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
    
    converged = sum(1 for r in all_results if r.get('converged', True))
    
    print(f"  Total traces generated: {len(all_results)}")
    print(f"  Converged: {converged}")
    print(f"  Failed: {len(all_results) - converged}")
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
