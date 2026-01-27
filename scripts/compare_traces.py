#!/usr/bin/env python3
"""
Compare XFOIL and RustFoil traces to identify divergence points.

Reads the summary files from traces/xfoil and traces/rustfoil directories
and produces a detailed comparison report.

Usage:
    python scripts/compare_traces.py [--output-dir comparison_results]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional
import numpy as np


def load_summary(summary_file: Path) -> Optional[dict]:
    """Load a summary.json file."""
    if not summary_file.exists():
        return None
    try:
        with open(summary_file) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return None


def compare_foil_re(
    foil: str,
    re: float,
    xfoil_dir: Path,
    rustfoil_dir: Path,
) -> dict:
    """
    Compare XFOIL and RustFoil results for a specific foil/Re combination.
    
    Returns dict with comparison metrics.
    """
    re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
    
    xfoil_summary = load_summary(xfoil_dir / foil / re_str / "summary.json")
    rustfoil_summary = load_summary(rustfoil_dir / foil / re_str / "summary.json")
    
    if not xfoil_summary or not rustfoil_summary:
        return {'error': f'Missing summary for {foil} Re={re}'}
    
    # Build alpha -> result maps
    xfoil_results = {r['alpha']: r for r in xfoil_summary.get('results', [])}
    rustfoil_results = {r['alpha']: r for r in rustfoil_summary.get('results', [])}
    
    # Compare at each alpha
    comparisons = []
    
    for alpha in sorted(set(xfoil_results.keys()) | set(rustfoil_results.keys())):
        xf = xfoil_results.get(alpha)
        rf = rustfoil_results.get(alpha)
        
        if not xf or not rf:
            comparisons.append({
                'alpha': alpha,
                'error': 'Missing data',
            })
            continue
        
        xf_cl = xf.get('cl')
        rf_cl = rf.get('cl')
        xf_cd = xf.get('cd')
        rf_cd = rf.get('cd')
        
        # Handle None values
        if xf_cl is None or rf_cl is None:
            cl_diff = None
            cl_rel = None
        else:
            cl_diff = rf_cl - xf_cl
            cl_rel = abs(cl_diff / xf_cl) * 100 if abs(xf_cl) > 0.01 else None
        
        if xf_cd is None or rf_cd is None:
            cd_diff = None
            cd_rel = None
        else:
            cd_diff = rf_cd - xf_cd
            cd_rel = abs(cd_diff / xf_cd) * 100 if xf_cd > 1e-6 else None
        
        comparisons.append({
            'alpha': alpha,
            'xfoil_cl': xf_cl,
            'rustfoil_cl': rf_cl,
            'cl_diff': cl_diff,
            'cl_rel_pct': cl_rel,
            'xfoil_cd': xf_cd,
            'rustfoil_cd': rf_cd,
            'cd_diff': cd_diff,
            'cd_rel_pct': cd_rel,
            'xfoil_converged': xf.get('converged', True),
            'rustfoil_converged': rf.get('converged', True),
        })
    
    # Summary statistics
    cl_diffs = [c['cl_diff'] for c in comparisons if c.get('cl_diff') is not None]
    cd_diffs = [c['cd_diff'] for c in comparisons if c.get('cd_diff') is not None]
    cl_rels = [c['cl_rel_pct'] for c in comparisons if c.get('cl_rel_pct') is not None]
    cd_rels = [c['cd_rel_pct'] for c in comparisons if c.get('cd_rel_pct') is not None]
    
    return {
        'foil': foil,
        'reynolds': re,
        'comparisons': comparisons,
        'n_points': len(comparisons),
        'cl_max_abs_diff': max(abs(d) for d in cl_diffs) if cl_diffs else None,
        'cl_mean_abs_diff': np.mean([abs(d) for d in cl_diffs]) if cl_diffs else None,
        'cl_max_rel_pct': max(cl_rels) if cl_rels else None,
        'cl_mean_rel_pct': np.mean(cl_rels) if cl_rels else None,
        'cd_max_abs_diff': max(abs(d) for d in cd_diffs) if cd_diffs else None,
        'cd_mean_abs_diff': np.mean([abs(d) for d in cd_diffs]) if cd_diffs else None,
        'cd_max_rel_pct': max(cd_rels) if cd_rels else None,
        'cd_mean_rel_pct': np.mean(cd_rels) if cd_rels else None,
    }


def find_worst_cases(all_results: list[dict], metric: str, top_n: int = 10) -> list[dict]:
    """
    Find the worst N cases based on a given metric.
    
    metric: 'cl_rel_pct', 'cd_rel_pct', 'cl_diff', 'cd_diff'
    """
    cases = []
    for result in all_results:
        if 'error' in result:
            continue
        for comp in result.get('comparisons', []):
            if comp.get(metric) is not None:
                cases.append({
                    'foil': result['foil'],
                    'reynolds': result['reynolds'],
                    'alpha': comp['alpha'],
                    'value': comp[metric],
                    'xfoil_cl': comp.get('xfoil_cl'),
                    'rustfoil_cl': comp.get('rustfoil_cl'),
                    'xfoil_cd': comp.get('xfoil_cd'),
                    'rustfoil_cd': comp.get('rustfoil_cd'),
                })
    
    # Sort by absolute value
    cases.sort(key=lambda x: abs(x['value']), reverse=True)
    return cases[:top_n]


def print_comparison_table(result: dict):
    """Print detailed comparison table for a foil/Re combination."""
    print(f"\n{'='*90}")
    print(f" {result['foil'].upper()} at Re = {result['reynolds']:.0e}")
    print("=" * 90)
    
    if 'error' in result:
        print(f"  ERROR: {result['error']}")
        return
    
    print(f"{'Alpha':>7} | {'XFOIL CL':>10} | {'RF CL':>10} | {'ΔCL':>10} | {'%CL':>7} | "
          f"{'XFOIL CD':>10} | {'RF CD':>10} | {'ΔCD':>10} | {'%CD':>7}")
    print("-" * 90)
    
    for comp in result.get('comparisons', []):
        if 'error' in comp:
            print(f"{comp['alpha']:>7.1f}° | ERROR: {comp['error']}")
            continue
        
        xf_cl = f"{comp['xfoil_cl']:>10.4f}" if comp.get('xfoil_cl') is not None else "         ?"
        rf_cl = f"{comp['rustfoil_cl']:>10.4f}" if comp.get('rustfoil_cl') is not None else "         ?"
        cl_diff = f"{comp['cl_diff']:>+10.4f}" if comp.get('cl_diff') is not None else "         -"
        cl_rel = f"{comp['cl_rel_pct']:>6.1f}%" if comp.get('cl_rel_pct') is not None else "      -"
        
        xf_cd = f"{comp['xfoil_cd']:>10.5f}" if comp.get('xfoil_cd') is not None else "         ?"
        rf_cd = f"{comp['rustfoil_cd']:>10.5f}" if comp.get('rustfoil_cd') is not None else "         ?"
        cd_diff = f"{comp['cd_diff']:>+10.5f}" if comp.get('cd_diff') is not None else "         -"
        cd_rel = f"{comp['cd_rel_pct']:>6.1f}%" if comp.get('cd_rel_pct') is not None else "      -"
        
        # Flag large errors
        flag = ""
        if comp.get('cl_rel_pct') and comp['cl_rel_pct'] > 20:
            flag = " *CL"
        if comp.get('cd_rel_pct') and comp['cd_rel_pct'] > 50:
            flag += " *CD"
        
        print(f"{comp['alpha']:>7.1f}° | {xf_cl} | {rf_cl} | {cl_diff} | {cl_rel} | "
              f"{xf_cd} | {rf_cd} | {cd_diff} | {cd_rel}{flag}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare XFOIL and RustFoil traces"
    )
    parser.add_argument(
        "--xfoil-dir",
        type=Path,
        default=Path("traces/xfoil"),
        help="XFOIL traces directory"
    )
    parser.add_argument(
        "--rustfoil-dir",
        type=Path,
        default=Path("traces/rustfoil"),
        help="RustFoil traces directory"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("comparison_results"),
        help="Output directory for comparison results"
    )
    parser.add_argument(
        "--brief", "-b",
        action="store_true",
        help="Brief output (summary only)"
    )
    args = parser.parse_args()
    
    workspace = Path(__file__).parent.parent
    xfoil_dir = workspace / args.xfoil_dir
    rustfoil_dir = workspace / args.rustfoil_dir
    output_dir = workspace / args.output_dir
    
    if not xfoil_dir.exists():
        print(f"ERROR: XFOIL traces directory not found: {xfoil_dir}")
        return 1
    if not rustfoil_dir.exists():
        print(f"ERROR: RustFoil traces directory not found: {rustfoil_dir}")
        return 1
    
    # Test matrix
    test_cases = [
        ("naca0012", [1e6, 3e6, 6e6]),
        ("naca2412", [3e6]),
        ("naca4412", [3e6]),
    ]
    
    print("=" * 90)
    print(" XFOIL vs RUSTFOIL COMPARISON REPORT")
    print("=" * 90)
    print(f"  XFOIL traces:    {xfoil_dir}")
    print(f"  RustFoil traces: {rustfoil_dir}")
    print()
    
    all_results = []
    
    for foil, re_list in test_cases:
        for re in re_list:
            result = compare_foil_re(foil, re, xfoil_dir, rustfoil_dir)
            all_results.append(result)
            
            if not args.brief:
                print_comparison_table(result)
    
    # Summary
    print(f"\n{'='*90}")
    print(" SUMMARY STATISTICS")
    print("=" * 90)
    
    print(f"\n{'Foil':<12} | {'Re':>10} | {'Max ΔCL':>10} | {'Mean ΔCL':>10} | "
          f"{'Max %CL':>8} | {'Max ΔCD':>10} | {'Mean ΔCD':>10} | {'Max %CD':>8}")
    print("-" * 90)
    
    for result in all_results:
        if 'error' in result:
            print(f"{result['foil']:<12} | {result['reynolds']:>10.0e} | ERROR")
            continue
        
        cl_max = f"{result['cl_max_abs_diff']:>10.4f}" if result.get('cl_max_abs_diff') else "         -"
        cl_mean = f"{result['cl_mean_abs_diff']:>10.4f}" if result.get('cl_mean_abs_diff') else "         -"
        cl_max_rel = f"{result['cl_max_rel_pct']:>7.1f}%" if result.get('cl_max_rel_pct') else "       -"
        cd_max = f"{result['cd_max_abs_diff']:>10.5f}" if result.get('cd_max_abs_diff') else "         -"
        cd_mean = f"{result['cd_mean_abs_diff']:>10.5f}" if result.get('cd_mean_abs_diff') else "         -"
        cd_max_rel = f"{result['cd_max_rel_pct']:>7.1f}%" if result.get('cd_max_rel_pct') else "       -"
        
        print(f"{result['foil']:<12} | {result['reynolds']:>10.0e} | {cl_max} | {cl_mean} | "
              f"{cl_max_rel} | {cd_max} | {cd_mean} | {cd_max_rel}")
    
    # Worst cases
    print(f"\n{'='*90}")
    print(" WORST 10 CASES BY CL ERROR")
    print("=" * 90)
    
    worst_cl = find_worst_cases(all_results, 'cl_rel_pct', 10)
    print(f"{'Foil':<12} | {'Re':>10} | {'Alpha':>7} | {'%CL':>8} | {'XFOIL CL':>10} | {'RF CL':>10}")
    print("-" * 70)
    for case in worst_cl:
        print(f"{case['foil']:<12} | {case['reynolds']:>10.0e} | {case['alpha']:>6.1f}° | "
              f"{case['value']:>7.1f}% | {case['xfoil_cl']:>10.4f} | {case['rustfoil_cl']:>10.4f}")
    
    print(f"\n{'='*90}")
    print(" WORST 10 CASES BY CD ERROR")
    print("=" * 90)
    
    worst_cd = find_worst_cases(all_results, 'cd_rel_pct', 10)
    print(f"{'Foil':<12} | {'Re':>10} | {'Alpha':>7} | {'%CD':>8} | {'XFOIL CD':>10} | {'RF CD':>10}")
    print("-" * 70)
    for case in worst_cd:
        print(f"{case['foil']:<12} | {case['reynolds']:>10.0e} | {case['alpha']:>6.1f}° | "
              f"{case['value']:>7.1f}% | {case['xfoil_cd']:>10.5f} | {case['rustfoil_cd']:>10.5f}")
    
    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(output_dir / "comparison_summary.json", 'w') as f:
        json.dump({
            'results': all_results,
            'worst_cl_cases': worst_cl,
            'worst_cd_cases': worst_cd,
        }, f, indent=2)
    
    print(f"\n\nResults saved to: {output_dir / 'comparison_summary.json'}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
