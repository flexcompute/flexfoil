#!/usr/bin/env python3
"""
Compare boundary layer station values between XFOIL and RustFoil.

This script extracts BLVAR events from both debug traces and compares
the BL variables (H, Hk, Cf, theta, etc.) at each station.

Usage:
    python scripts/compare_bl_stations.py <xfoil_trace.json> <rustfoil_trace.json>
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional
import numpy as np


def load_trace(filepath: Path) -> dict:
    """Load debug JSON trace."""
    with open(filepath) as f:
        return json.load(f)


def extract_blvar_events(data: dict, iteration: int = 1) -> list[dict]:
    """
    Extract BLVAR events for a specific iteration.
    
    Returns list of {side, ibl, input, output} dicts.
    """
    events = data.get('events', [])
    
    blvar_events = []
    for e in events:
        if e.get('subroutine') == 'BLVAR':
            # Get iteration number (different field names)
            iter_num = e.get('iteration') or e.get('iter') or e.get('newton_iter') or 1
            if iter_num == iteration:
                blvar_events.append({
                    'side': e.get('side', e.get('SIDE', 0)),
                    'ibl': e.get('ibl', e.get('IBL', 0)),
                    'input': e.get('input', {}),
                    'output': e.get('output', {}),
                })
    
    return blvar_events


def extract_final_blvar(data: dict) -> list[dict]:
    """
    Extract final BLVAR events (last iteration).
    
    Returns list sorted by side, then ibl.
    """
    events = data.get('events', [])
    
    # Find max iteration
    max_iter = 0
    for e in events:
        if e.get('subroutine') == 'BLVAR':
            iter_num = e.get('iteration') or e.get('iter') or e.get('newton_iter') or 1
            if iter_num > max_iter:
                max_iter = iter_num
    
    # Get events from last iteration
    blvar_events = extract_blvar_events(data, max_iter)
    
    # Sort by side, then ibl
    blvar_events.sort(key=lambda x: (x['side'], x['ibl']))
    
    return blvar_events, max_iter


def compare_blvar(xf_events: list[dict], rf_events: list[dict]) -> dict:
    """
    Compare BLVAR events from both solvers.
    
    Returns comparison results.
    """
    # Build lookup by (side, ibl)
    xf_lookup = {(e['side'], e['ibl']): e for e in xf_events}
    rf_lookup = {(e['side'], e['ibl']): e for e in rf_events}
    
    comparisons = []
    
    # Compare all stations
    all_keys = sorted(set(xf_lookup.keys()) | set(rf_lookup.keys()))
    
    for key in all_keys:
        side, ibl = key
        xf = xf_lookup.get(key)
        rf = rf_lookup.get(key)
        
        comp = {'side': side, 'ibl': ibl}
        
        if not xf:
            comp['error'] = 'Missing in XFOIL'
            comparisons.append(comp)
            continue
        if not rf:
            comp['error'] = 'Missing in RustFoil'
            comparisons.append(comp)
            continue
        
        # Extract values (handle different field names)
        xf_out = xf.get('output', {})
        rf_out = rf.get('output', {})
        xf_in = xf.get('input', {})
        rf_in = rf.get('input', {})
        
        # Key variables to compare
        for var in ['H', 'Hk', 'Hs', 'Cf', 'Cd', 'Rtheta', 'Us', 'Cq', 'De']:
            xf_val = xf_out.get(var) or xf_out.get(var.lower())
            rf_val = rf_out.get(var) or rf_out.get(var.lower())
            
            if xf_val is not None and rf_val is not None:
                diff = rf_val - xf_val
                rel = abs(diff / xf_val) * 100 if abs(xf_val) > 1e-12 else None
                comp[f'{var}_xf'] = xf_val
                comp[f'{var}_rf'] = rf_val
                comp[f'{var}_diff'] = diff
                comp[f'{var}_rel'] = rel
        
        # Input variables
        for var in ['x', 'u', 'theta', 'delta_star', 'ctau', 'ampl']:
            xf_val = xf_in.get(var) or xf_in.get(var.upper())
            rf_val = rf_in.get(var) or rf_in.get(var.upper())
            
            if xf_val is not None and rf_val is not None:
                diff = rf_val - xf_val
                comp[f'in_{var}_xf'] = xf_val
                comp[f'in_{var}_rf'] = rf_val
                comp[f'in_{var}_diff'] = diff
        
        comparisons.append(comp)
    
    return comparisons


def find_first_divergence(comparisons: list[dict], var: str, threshold: float) -> Optional[dict]:
    """
    Find first station where a variable diverges beyond threshold.
    
    threshold: relative error threshold (percentage)
    """
    rel_key = f'{var}_rel'
    
    for comp in comparisons:
        if 'error' in comp:
            continue
        rel = comp.get(rel_key)
        if rel is not None and rel > threshold:
            return {
                'side': comp['side'],
                'ibl': comp['ibl'],
                'variable': var,
                'xfoil': comp[f'{var}_xf'],
                'rustfoil': comp[f'{var}_rf'],
                'diff': comp[f'{var}_diff'],
                'rel_pct': rel,
            }
    
    return None


def print_station_table(comparisons: list[dict], vars_to_show: list[str]):
    """Print comparison table for specified variables."""
    # Header
    header = f"{'Side':>4} {'IBL':>4}"
    for var in vars_to_show:
        header += f" | {var+'_XF':>10} {var+'_RF':>10} {'Δ'+var:>10} {'%':>6}"
    print(header)
    print("-" * len(header))
    
    for comp in comparisons:
        if 'error' in comp:
            print(f"{comp['side']:>4} {comp['ibl']:>4} | ERROR: {comp['error']}")
            continue
        
        line = f"{comp['side']:>4} {comp['ibl']:>4}"
        for var in vars_to_show:
            xf = comp.get(f'{var}_xf')
            rf = comp.get(f'{var}_rf')
            diff = comp.get(f'{var}_diff')
            rel = comp.get(f'{var}_rel')
            
            if xf is not None:
                xf_str = f"{xf:>10.6f}"
                rf_str = f"{rf:>10.6f}"
                diff_str = f"{diff:>+10.6f}"
                rel_str = f"{rel:>5.1f}%" if rel is not None else "     -"
            else:
                xf_str = rf_str = diff_str = rel_str = "         -"
            
            line += f" | {xf_str} {rf_str} {diff_str} {rel_str}"
        
        print(line)


def main():
    parser = argparse.ArgumentParser(
        description="Compare BL station values between XFOIL and RustFoil"
    )
    parser.add_argument("xfoil_trace", type=Path, help="XFOIL debug JSON")
    parser.add_argument("rustfoil_trace", type=Path, help="RustFoil debug JSON")
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=5.0,
        help="Divergence threshold (percent, default: 5%%)"
    )
    parser.add_argument(
        "--iteration", "-i",
        type=int,
        default=None,
        help="Specific iteration to compare (default: last)"
    )
    parser.add_argument(
        "--side", "-s",
        type=int,
        default=None,
        help="Filter by side (1=upper, 2=lower)"
    )
    args = parser.parse_args()
    
    print("=" * 100)
    print(" BOUNDARY LAYER STATION COMPARISON")
    print("=" * 100)
    print(f"  XFOIL:    {args.xfoil_trace}")
    print(f"  RustFoil: {args.rustfoil_trace}")
    print(f"  Threshold: {args.threshold}%")
    print()
    
    # Load traces
    xf_data = load_trace(args.xfoil_trace)
    rf_data = load_trace(args.rustfoil_trace)
    
    # Extract BLVAR events
    if args.iteration:
        xf_events = extract_blvar_events(xf_data, args.iteration)
        rf_events = extract_blvar_events(rf_data, args.iteration)
        iter_str = f"iteration {args.iteration}"
    else:
        xf_events, xf_iter = extract_final_blvar(xf_data)
        rf_events, rf_iter = extract_final_blvar(rf_data)
        iter_str = f"XFOIL iter {xf_iter}, RustFoil iter {rf_iter}"
    
    print(f"  Comparing {iter_str}")
    print(f"  XFOIL stations: {len(xf_events)}")
    print(f"  RustFoil stations: {len(rf_events)}")
    
    # Filter by side if requested
    if args.side:
        xf_events = [e for e in xf_events if e['side'] == args.side]
        rf_events = [e for e in rf_events if e['side'] == args.side]
        print(f"  Filtered to side {args.side}: XFOIL={len(xf_events)}, RustFoil={len(rf_events)}")
    
    # Compare
    comparisons = compare_blvar(xf_events, rf_events)
    
    # Find first divergences
    print(f"\n{'='*100}")
    print(" FIRST DIVERGENCE POINTS")
    print("=" * 100)
    
    vars_to_check = ['H', 'Hk', 'Cf', 'Cd', 'Rtheta']
    first_divs = []
    
    for var in vars_to_check:
        div = find_first_divergence(comparisons, var, args.threshold)
        if div:
            first_divs.append(div)
            print(f"  {var:<8}: Side={div['side']} IBL={div['ibl']} "
                  f"XF={div['xfoil']:.6f} RF={div['rustfoil']:.6f} "
                  f"Δ={div['diff']:+.6f} ({div['rel_pct']:.1f}%)")
        else:
            print(f"  {var:<8}: No divergence above {args.threshold}%")
    
    # Summary statistics
    print(f"\n{'='*100}")
    print(" SUMMARY STATISTICS")
    print("=" * 100)
    
    for var in ['H', 'Hk', 'Cf']:
        diffs = [abs(c.get(f'{var}_diff', 0)) for c in comparisons if f'{var}_diff' in c]
        rels = [c.get(f'{var}_rel', 0) for c in comparisons if c.get(f'{var}_rel') is not None]
        
        if diffs:
            print(f"  {var}:")
            print(f"    Max abs diff: {max(diffs):.6f}")
            print(f"    Mean abs diff: {np.mean(diffs):.6f}")
            if rels:
                print(f"    Max rel error: {max(rels):.1f}%")
                print(f"    Mean rel error: {np.mean(rels):.1f}%")
    
    # Detailed table for upper surface first 20 stations
    print(f"\n{'='*100}")
    print(" DETAILED COMPARISON (first 20 upper surface stations)")
    print("=" * 100)
    
    upper_comps = [c for c in comparisons if c.get('side') == 1][:20]
    if upper_comps:
        print_station_table(upper_comps, ['Hk', 'Cf'])
    else:
        print("  No upper surface data found")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
