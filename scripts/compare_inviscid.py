#!/usr/bin/env python3
"""
Compare inviscid solver outputs between XFOIL and RustFoil.

Extracts FULLGAMMA and FULLAIC events from debug JSON files and compares:
- Gamma distribution (circulation)
- Cp distribution (pressure coefficient)
- Qinv distribution (inviscid velocity)
- Inviscid CL
- AIC base solutions (gamu_0, gamu_90)

Usage:
    python scripts/compare_inviscid.py xfoil_debug.json rustfoil_debug.json [--tolerance 0.01] [--plot]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np


def load_json(path: Path) -> dict:
    """Load and parse a debug JSON file."""
    with open(path) as f:
        return json.load(f)


def extract_fullgamma(data: dict, source: str = "unknown") -> Optional[dict]:
    """
    Extract FULLGAMMA event from debug JSON.
    
    Expected format:
    {
        "subroutine": "FULLGAMMA",
        "gamma": [...],      # Full gamma distribution
        "qinv": [...],       # Inviscid velocity at each panel
        "cp": [...],         # Pressure coefficient at each panel
        "cl_inv": float,     # Inviscid lift coefficient
        "x": [...],          # Optional: x coordinates
        "y": [...],          # Optional: y coordinates
        "s": [...],          # Optional: arc length coordinates
    }
    """
    events = data.get('events', [])
    
    # Find FULLGAMMA event (use last one if multiple)
    fullgamma_events = [e for e in events if e.get('subroutine') == 'FULLGAMMA']
    
    if not fullgamma_events:
        print(f"  [{source}] Warning: No FULLGAMMA event found")
        return None
    
    # Use the last FULLGAMMA event (most converged state)
    event = fullgamma_events[-1]
    
    return {
        'gamma': np.array(event.get('gamma', [])),
        'qinv': np.array(event.get('qinv', [])),
        'cp': np.array(event.get('cp', [])),
        'cl_inv': event.get('cl_inv', None),
        'x': np.array(event.get('x', [])) if 'x' in event else None,
        'y': np.array(event.get('y', [])) if 'y' in event else None,
        's': np.array(event.get('s', [])) if 's' in event else None,
        'n_panels': len(event.get('gamma', [])),
        'iteration': event.get('iteration', None),
    }


def extract_fullaic(data: dict, source: str = "unknown") -> Optional[dict]:
    """
    Extract FULLAIC event from debug JSON.
    
    Expected format:
    {
        "subroutine": "FULLAIC",
        "gamu_0": [...],     # Gamma for alpha=0 base solution
        "gamu_90": [...],    # Gamma for alpha=90 base solution
        "qinf_0": [...],     # Optional: qinf for alpha=0
        "qinf_90": [...],    # Optional: qinf for alpha=90
    }
    """
    events = data.get('events', [])
    
    # Find FULLAIC event
    fullaic_events = [e for e in events if e.get('subroutine') == 'FULLAIC']
    
    if not fullaic_events:
        print(f"  [{source}] Warning: No FULLAIC event found")
        return None
    
    event = fullaic_events[-1]
    
    return {
        'gamu_0': np.array(event.get('gamu_0', [])),
        'gamu_90': np.array(event.get('gamu_90', [])),
        'qinf_0': np.array(event.get('qinf_0', [])) if 'qinf_0' in event else None,
        'qinf_90': np.array(event.get('qinf_90', [])) if 'qinf_90' in event else None,
    }


def compute_differences(arr1: np.ndarray, arr2: np.ndarray, name: str) -> dict:
    """
    Compute absolute and relative differences between two arrays.
    
    Returns dict with statistics and per-element differences.
    """
    if len(arr1) == 0 or len(arr2) == 0:
        return {'error': 'Empty array(s)'}
    
    if len(arr1) != len(arr2):
        return {
            'error': f'Array length mismatch: {len(arr1)} vs {len(arr2)}',
            'len_xfoil': len(arr1),
            'len_rustfoil': len(arr2),
        }
    
    abs_diff = np.abs(arr1 - arr2)
    
    # Relative difference (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff = np.where(
            np.abs(arr1) > 1e-12,
            abs_diff / np.abs(arr1),
            np.where(abs_diff > 1e-12, np.inf, 0.0)
        )
    
    # Replace inf with a large but finite number for statistics
    rel_diff_finite = np.where(np.isinf(rel_diff), np.nan, rel_diff)
    
    return {
        'name': name,
        'n_points': len(arr1),
        'abs_diff': abs_diff,
        'rel_diff': rel_diff,
        'max_abs': float(np.max(abs_diff)),
        'mean_abs': float(np.mean(abs_diff)),
        'rms_abs': float(np.sqrt(np.mean(abs_diff**2))),
        'max_rel': float(np.nanmax(rel_diff_finite)) if not np.all(np.isnan(rel_diff_finite)) else np.inf,
        'mean_rel': float(np.nanmean(rel_diff_finite)) if not np.all(np.isnan(rel_diff_finite)) else np.inf,
    }


def find_first_divergence(
    arr1: np.ndarray,
    arr2: np.ndarray,
    threshold: float,
    x_coords: Optional[np.ndarray] = None
) -> Optional[dict]:
    """
    Find the first station where relative difference exceeds threshold.
    
    Returns dict with index, values, and optional x location.
    """
    if len(arr1) != len(arr2) or len(arr1) == 0:
        return None
    
    for i in range(len(arr1)):
        ref = abs(arr1[i])
        diff = abs(arr1[i] - arr2[i])
        
        if ref > 1e-12:
            rel_diff = diff / ref
            if rel_diff > threshold:
                result = {
                    'index': i,
                    'xfoil_value': float(arr1[i]),
                    'rustfoil_value': float(arr2[i]),
                    'abs_diff': float(diff),
                    'rel_diff': float(rel_diff),
                    'rel_diff_pct': float(rel_diff * 100),
                }
                if x_coords is not None and i < len(x_coords):
                    result['x'] = float(x_coords[i])
                return result
    
    return None


def print_comparison_table(diffs: list[dict], title: str):
    """Print a formatted comparison table for multiple quantities."""
    print(f"\n{'='*80}")
    print(f" {title}")
    print(f"{'='*80}")
    print(f"{'Quantity':<12} | {'N':>6} | {'Max Abs':>12} | {'Mean Abs':>12} | {'RMS':>12} | {'Max Rel%':>10}")
    print("-" * 80)
    
    for d in diffs:
        if 'error' in d:
            print(f"{d.get('name', '?'):<12} | ERROR: {d['error']}")
        else:
            max_rel_pct = d['max_rel'] * 100 if d['max_rel'] != np.inf else float('inf')
            print(f"{d['name']:<12} | {d['n_points']:>6} | {d['max_abs']:>12.6g} | "
                  f"{d['mean_abs']:>12.6g} | {d['rms_abs']:>12.6g} | {max_rel_pct:>9.4f}%")


def print_divergence_report(divergences: list[dict], tolerance: float):
    """Print report of first divergences for each quantity."""
    print(f"\n{'='*80}")
    print(f" FIRST DIVERGENCE REPORT (threshold = {tolerance*100:.2f}%)")
    print(f"{'='*80}")
    
    has_any = False
    for div in divergences:
        name = div.get('name', '?')
        info = div.get('info')
        
        if info is None:
            print(f"  {name:<12}: No divergence above threshold")
        else:
            has_any = True
            x_str = f", x={info['x']:.6f}" if 'x' in info else ""
            print(f"  {name:<12}: DIVERGES at index {info['index']}{x_str}")
            print(f"      XFOIL:    {info['xfoil_value']:>15.8g}")
            print(f"      RustFoil: {info['rustfoil_value']:>15.8g}")
            print(f"      Diff:     {info['rel_diff_pct']:>14.4f}%")
    
    return has_any


def print_scalar_comparison(label: str, xf_val, rf_val):
    """Print comparison of a single scalar value."""
    if xf_val is None or rf_val is None:
        print(f"  {label:<20}: XFOIL={xf_val}, RustFoil={rf_val} (missing data)")
        return None
    
    diff = abs(xf_val - rf_val)
    rel = diff / abs(xf_val) if abs(xf_val) > 1e-12 else (float('inf') if diff > 1e-12 else 0.0)
    
    print(f"  {label:<20}: XFOIL={xf_val:>12.8g}  RustFoil={rf_val:>12.8g}  "
          f"Diff={diff:>10.6g}  ({rel*100:>7.4f}%)")
    
    return rel


def create_plots(
    xf_gamma: dict,
    rf_gamma: dict,
    output_prefix: str = "inviscid_comparison"
):
    """Generate matplotlib comparison plots."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nWarning: matplotlib not available, skipping plots")
        return
    
    # Determine x-axis (prefer x coordinates, fall back to index)
    if xf_gamma.get('x') is not None and len(xf_gamma['x']) > 0:
        x_xf = xf_gamma['x']
        x_rf = rf_gamma.get('x', np.arange(len(rf_gamma['gamma'])))
        x_label = "x/c"
    elif xf_gamma.get('s') is not None and len(xf_gamma['s']) > 0:
        x_xf = xf_gamma['s']
        x_rf = rf_gamma.get('s', np.arange(len(rf_gamma['gamma'])))
        x_label = "s (arc length)"
    else:
        x_xf = np.arange(len(xf_gamma['gamma']))
        x_rf = np.arange(len(rf_gamma['gamma']))
        x_label = "Panel index"
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Inviscid Solver Comparison: XFOIL vs RustFoil", fontsize=14)
    
    # 1. Gamma overlay
    ax = axes[0, 0]
    if len(xf_gamma['gamma']) > 0:
        ax.plot(x_xf, xf_gamma['gamma'], 'b-', label='XFOIL', linewidth=1.5)
    if len(rf_gamma['gamma']) > 0:
        ax.plot(x_rf, rf_gamma['gamma'], 'r--', label='RustFoil', linewidth=1.5)
    ax.set_xlabel(x_label)
    ax.set_ylabel("γ (gamma)")
    ax.set_title("Gamma Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Cp overlay
    ax = axes[0, 1]
    if len(xf_gamma['cp']) > 0:
        ax.plot(x_xf, xf_gamma['cp'], 'b-', label='XFOIL', linewidth=1.5)
    if len(rf_gamma['cp']) > 0:
        ax.plot(x_rf, rf_gamma['cp'], 'r--', label='RustFoil', linewidth=1.5)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Cp")
    ax.set_title("Pressure Coefficient")
    ax.invert_yaxis()  # Convention: negative Cp up
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Qinv overlay
    ax = axes[0, 2]
    if len(xf_gamma['qinv']) > 0:
        ax.plot(x_xf, xf_gamma['qinv'], 'b-', label='XFOIL', linewidth=1.5)
    if len(rf_gamma['qinv']) > 0:
        ax.plot(x_rf, rf_gamma['qinv'], 'r--', label='RustFoil', linewidth=1.5)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Qinv")
    ax.set_title("Inviscid Velocity")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Gamma difference
    ax = axes[1, 0]
    if len(xf_gamma['gamma']) == len(rf_gamma['gamma']) and len(xf_gamma['gamma']) > 0:
        diff = xf_gamma['gamma'] - rf_gamma['gamma']
        ax.plot(x_xf, diff, 'g-', linewidth=1.5)
        ax.axhline(y=0, color='k', linestyle=':', alpha=0.5)
        ax.fill_between(x_xf, diff, alpha=0.3)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Δγ (XFOIL - RustFoil)")
    ax.set_title("Gamma Difference")
    ax.grid(True, alpha=0.3)
    
    # 5. Cp difference
    ax = axes[1, 1]
    if len(xf_gamma['cp']) == len(rf_gamma['cp']) and len(xf_gamma['cp']) > 0:
        diff = xf_gamma['cp'] - rf_gamma['cp']
        ax.plot(x_xf, diff, 'g-', linewidth=1.5)
        ax.axhline(y=0, color='k', linestyle=':', alpha=0.5)
        ax.fill_between(x_xf, diff, alpha=0.3)
    ax.set_xlabel(x_label)
    ax.set_ylabel("ΔCp (XFOIL - RustFoil)")
    ax.set_title("Cp Difference")
    ax.grid(True, alpha=0.3)
    
    # 6. Relative differences
    ax = axes[1, 2]
    if len(xf_gamma['gamma']) == len(rf_gamma['gamma']) and len(xf_gamma['gamma']) > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_gamma = np.abs(xf_gamma['gamma'] - rf_gamma['gamma']) / (np.abs(xf_gamma['gamma']) + 1e-12) * 100
            rel_gamma = np.clip(rel_gamma, 0, 100)  # Cap at 100% for visualization
        ax.plot(x_xf, rel_gamma, 'm-', linewidth=1.5, label='γ rel. error')
    if len(xf_gamma['cp']) == len(rf_gamma['cp']) and len(xf_gamma['cp']) > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_cp = np.abs(xf_gamma['cp'] - rf_gamma['cp']) / (np.abs(xf_gamma['cp']) + 1e-12) * 100
            rel_cp = np.clip(rel_cp, 0, 100)
        ax.plot(x_xf, rel_cp, 'c-', linewidth=1.5, label='Cp rel. error')
    ax.set_xlabel(x_label)
    ax.set_ylabel("Relative Error (%)")
    ax.set_title("Relative Differences")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    plot_path = f"{output_prefix}.png"
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlots saved to: {plot_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Compare inviscid solver outputs between XFOIL and RustFoil",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python scripts/compare_inviscid.py xfoil_debug.json rustfoil_debug.json --tolerance 0.01 --plot
        """
    )
    parser.add_argument("xfoil_json", type=Path, help="Path to XFOIL debug JSON file")
    parser.add_argument("rustfoil_json", type=Path, help="Path to RustFoil debug JSON file")
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=0.01,
        help="Relative difference threshold for divergence detection (default: 0.01 = 1%%)"
    )
    parser.add_argument(
        "--plot", "-p",
        action="store_true",
        help="Generate matplotlib comparison plots"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="inviscid_comparison",
        help="Output prefix for plot files (default: inviscid_comparison)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-station comparison"
    )
    
    args = parser.parse_args()
    
    # Validate input files
    if not args.xfoil_json.exists():
        print(f"Error: XFOIL JSON file not found: {args.xfoil_json}")
        return 1
    if not args.rustfoil_json.exists():
        print(f"Error: RustFoil JSON file not found: {args.rustfoil_json}")
        return 1
    
    print("=" * 80)
    print(" INVISCID SOLVER COMPARISON: XFOIL vs RustFoil")
    print("=" * 80)
    print(f"  XFOIL file:    {args.xfoil_json}")
    print(f"  RustFoil file: {args.rustfoil_json}")
    print(f"  Tolerance:     {args.tolerance*100:.2f}%")
    
    # Load data
    print("\nLoading JSON files...")
    xf_data = load_json(args.xfoil_json)
    rf_data = load_json(args.rustfoil_json)
    
    # Extract FULLGAMMA events
    print("\nExtracting FULLGAMMA events...")
    xf_gamma = extract_fullgamma(xf_data, "XFOIL")
    rf_gamma = extract_fullgamma(rf_data, "RustFoil")
    
    # Extract FULLAIC events
    print("Extracting FULLAIC events...")
    xf_aic = extract_fullaic(xf_data, "XFOIL")
    rf_aic = extract_fullaic(rf_data, "RustFoil")
    
    # Check if we have data to compare
    has_gamma_data = xf_gamma is not None and rf_gamma is not None
    has_aic_data = xf_aic is not None and rf_aic is not None
    
    if not has_gamma_data and not has_aic_data:
        print("\nError: No FULLGAMMA or FULLAIC events found in either file.")
        print("Make sure both solvers are configured to output inviscid debug events.")
        return 1
    
    # Summary metrics
    print(f"\n{'='*80}")
    print(" DATA SUMMARY")
    print(f"{'='*80}")
    
    if xf_gamma:
        print(f"  XFOIL FULLGAMMA:    {xf_gamma['n_panels']} panels, CL_inv={xf_gamma['cl_inv']}")
    if rf_gamma:
        print(f"  RustFoil FULLGAMMA: {rf_gamma['n_panels']} panels, CL_inv={rf_gamma['cl_inv']}")
    if xf_aic:
        print(f"  XFOIL FULLAIC:      gamu_0[{len(xf_aic['gamu_0'])}], gamu_90[{len(xf_aic['gamu_90'])}]")
    if rf_aic:
        print(f"  RustFoil FULLAIC:   gamu_0[{len(rf_aic['gamu_0'])}], gamu_90[{len(rf_aic['gamu_90'])}]")
    
    # =========================================================================
    # COMPARISON: FULLGAMMA quantities
    # =========================================================================
    all_diffs = []
    all_divergences = []
    overall_pass = True
    
    if has_gamma_data:
        # Compute differences for each distribution
        gamma_diff = compute_differences(xf_gamma['gamma'], rf_gamma['gamma'], "gamma")
        cp_diff = compute_differences(xf_gamma['cp'], rf_gamma['cp'], "Cp")
        qinv_diff = compute_differences(xf_gamma['qinv'], rf_gamma['qinv'], "Qinv")
        
        all_diffs.extend([gamma_diff, cp_diff, qinv_diff])
        
        # Print comparison table
        print_comparison_table([gamma_diff, cp_diff, qinv_diff], "FULLGAMMA DISTRIBUTION COMPARISON")
        
        # Scalar comparisons
        print(f"\n{'='*80}")
        print(" SCALAR COMPARISONS")
        print(f"{'='*80}")
        
        cl_rel = print_scalar_comparison("CL_inv", xf_gamma['cl_inv'], rf_gamma['cl_inv'])
        if cl_rel is not None and cl_rel > args.tolerance:
            overall_pass = False
        
        # Find first divergences
        x_coords = xf_gamma.get('x')
        
        gamma_div = find_first_divergence(xf_gamma['gamma'], rf_gamma['gamma'], args.tolerance, x_coords)
        cp_div = find_first_divergence(xf_gamma['cp'], rf_gamma['cp'], args.tolerance, x_coords)
        qinv_div = find_first_divergence(xf_gamma['qinv'], rf_gamma['qinv'], args.tolerance, x_coords)
        
        all_divergences.extend([
            {'name': 'gamma', 'info': gamma_div},
            {'name': 'Cp', 'info': cp_div},
            {'name': 'Qinv', 'info': qinv_div},
        ])
        
        # Check overall pass/fail based on max relative error
        for d in [gamma_diff, cp_diff, qinv_diff]:
            if 'error' not in d and d['max_rel'] > args.tolerance:
                overall_pass = False
    
    # =========================================================================
    # COMPARISON: FULLAIC quantities
    # =========================================================================
    if has_aic_data:
        gamu0_diff = compute_differences(xf_aic['gamu_0'], rf_aic['gamu_0'], "gamu_0")
        gamu90_diff = compute_differences(xf_aic['gamu_90'], rf_aic['gamu_90'], "gamu_90")
        
        all_diffs.extend([gamu0_diff, gamu90_diff])
        
        # Print AIC comparison table
        print_comparison_table([gamu0_diff, gamu90_diff], "FULLAIC BASE SOLUTION COMPARISON")
        
        # Find divergences in AIC
        gamu0_div = find_first_divergence(xf_aic['gamu_0'], rf_aic['gamu_0'], args.tolerance)
        gamu90_div = find_first_divergence(xf_aic['gamu_90'], rf_aic['gamu_90'], args.tolerance)
        
        all_divergences.extend([
            {'name': 'gamu_0', 'info': gamu0_div},
            {'name': 'gamu_90', 'info': gamu90_div},
        ])
        
        for d in [gamu0_diff, gamu90_diff]:
            if 'error' not in d and d['max_rel'] > args.tolerance:
                overall_pass = False
    
    # =========================================================================
    # DIVERGENCE REPORT
    # =========================================================================
    has_divergence = print_divergence_report(all_divergences, args.tolerance)
    
    # =========================================================================
    # VERBOSE: Per-station details
    # =========================================================================
    if args.verbose and has_gamma_data:
        print(f"\n{'='*80}")
        print(" DETAILED PER-STATION COMPARISON (first 20 stations)")
        print(f"{'='*80}")
        
        n_show = min(20, xf_gamma['n_panels'], rf_gamma['n_panels'])
        
        print(f"{'Idx':>4} | {'x':>10} | {'γ_XF':>12} | {'γ_RF':>12} | {'Δγ':>12} | {'Err%':>8}")
        print("-" * 75)
        
        for i in range(n_show):
            x_val = xf_gamma['x'][i] if xf_gamma['x'] is not None and i < len(xf_gamma['x']) else float('nan')
            g_xf = xf_gamma['gamma'][i] if i < len(xf_gamma['gamma']) else float('nan')
            g_rf = rf_gamma['gamma'][i] if i < len(rf_gamma['gamma']) else float('nan')
            diff = g_xf - g_rf
            rel = abs(diff) / abs(g_xf) * 100 if abs(g_xf) > 1e-12 else 0.0
            
            flag = " *" if rel > args.tolerance * 100 else ""
            print(f"{i:>4} | {x_val:>10.6f} | {g_xf:>12.6g} | {g_rf:>12.6g} | {diff:>+12.6g} | {rel:>7.3f}%{flag}")
    
    # =========================================================================
    # FINAL STATUS
    # =========================================================================
    print(f"\n{'='*80}")
    print(" OVERALL STATUS")
    print(f"{'='*80}")
    
    if overall_pass:
        print(f"  ✓ PASS - All quantities within {args.tolerance*100:.2f}% tolerance")
        status = 0
    else:
        print(f"  ✗ FAIL - Some quantities exceed {args.tolerance*100:.2f}% tolerance")
        status = 1
    
    # Summary of worst differences
    print("\n  Worst relative differences:")
    for d in all_diffs:
        if 'error' not in d:
            pct = d['max_rel'] * 100
            flag = " ← EXCEEDS" if pct > args.tolerance * 100 else ""
            print(f"    {d['name']:<12}: {pct:>10.4f}%{flag}")
    
    # =========================================================================
    # PLOTS
    # =========================================================================
    if args.plot and has_gamma_data:
        create_plots(xf_gamma, rf_gamma, args.output)
    
    return status


if __name__ == "__main__":
    sys.exit(main())
