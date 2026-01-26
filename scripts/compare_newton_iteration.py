#!/usr/bin/env python3
"""
Compare Newton system setup and solutions between XFOIL and RustFoil.

This script analyzes the viscous-inviscid coupling Newton iteration:
1. SETBL events - System matrix setup (VA, VB blocks)
2. BLDIF events - Jacobian matrices and residuals
3. Solution deltas (VDEL) - Changes in BL variables

For each Newton iteration (or specified iteration):
    a. Extract SETBL events (VA, VB, VDEL matrices)
    b. Extract BLDIF events (VS1, VS2, VSREZ)
3. Compare:
    - VA blocks: compare each 3x3 block
    - VB blocks: compare each 3x3 block
    - RHS vectors: compare each station
    - Solution deltas: compare [dCtau, dTheta, dMass] per station
4. Find first station/element where divergence exceeds threshold
5. Report convergence pattern comparison (residuals over iterations)

Usage:
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json --iteration 1 --tolerance 0.01
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json --all-iterations
"""

import json
import argparse
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass


@dataclass
class SetblData:
    """Data from a SETBL event."""
    side: int
    ibl: int
    iv: int
    va: np.ndarray  # 3x2 or 3x3 matrix
    vb: np.ndarray  # 3x2 or 3x3 matrix
    vdel: np.ndarray  # 3x2 solution deltas
    vm_sample: Optional[np.ndarray] = None  # VM row samples


@dataclass
class BldifData:
    """Data from a BLDIF event."""
    side: int
    ibl: int
    flow_type: int
    newton_iter: int
    x1: float
    x2: float
    vs1: np.ndarray  # 4x5 Jacobian for station 1
    vs2: np.ndarray  # 4x5 Jacobian for station 2
    vsrez: np.ndarray  # 4-element residual vector


@dataclass
class StationComparison:
    """Comparison result for a single station."""
    station: int
    side: int
    max_va_diff: float
    max_vb_diff: float
    max_vdel_diff: float
    max_vsrez_diff: float
    passed: bool


def load_debug_events(json_path: Path) -> List[Dict[str, Any]]:
    """Load events from debug JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    
    if isinstance(data, dict) and 'events' in data:
        return data['events']
    elif isinstance(data, list):
        return data
    else:
        raise ValueError(f"Unexpected JSON structure in {json_path}")


def extract_setbl_events(events: List[Dict], iteration: Optional[int] = None) -> List[SetblData]:
    """Extract SETBL events from debug output."""
    result = []
    
    for event in events:
        if event.get('subroutine') != 'SETBL':
            continue
        if iteration is not None and event.get('iteration') != iteration:
            continue
        
        va = np.array(event.get('VA', []))
        vb = np.array(event.get('VB', []))
        vdel = np.array(event.get('VDEL', []))
        vm_sample = np.array(event.get('VM_sample', [])) if 'VM_sample' in event else None
        
        result.append(SetblData(
            side=event.get('side', 0),
            ibl=event.get('ibl', 0),
            iv=event.get('iv', 0),
            va=va,
            vb=vb,
            vdel=vdel,
            vm_sample=vm_sample,
        ))
    
    return result


def extract_bldif_events(events: List[Dict], iteration: Optional[int] = None) -> List[BldifData]:
    """Extract BLDIF events from debug output."""
    result = []
    
    for event in events:
        if event.get('subroutine') != 'BLDIF':
            continue
        if iteration is not None and event.get('iteration') != iteration:
            continue
        
        vs1 = np.array(event.get('VS1', []))
        vs2 = np.array(event.get('VS2', []))
        vsrez = np.array(event.get('VSREZ', []))
        
        result.append(BldifData(
            side=event.get('side', 0),
            ibl=event.get('ibl', 0),
            flow_type=event.get('flow_type', 0),
            newton_iter=event.get('newton_iter', 0),
            x1=event.get('X1', 0.0),
            x2=event.get('X2', 0.0),
            vs1=vs1,
            vs2=vs2,
            vsrez=vsrez,
        ))
    
    return result


def extract_viscal_results(events: List[Dict]) -> List[Dict[str, float]]:
    """Extract VISCAL_RESULT events for convergence tracking."""
    results = []
    
    for event in events:
        if event.get('subroutine') in ('VISCAL_RESULT', 'ViscalResult'):
            results.append({
                'iteration': event.get('iteration', len(results) + 1),
                'rms_residual': event.get('rms_residual', 0.0),
                'max_residual': event.get('max_residual', 0.0),
                'CL': event.get('CL', 0.0),
                'CD': event.get('CD', 0.0),
            })
    
    return results


def compare_setbl(
    xf_setbl: List[SetblData],
    rf_setbl: List[SetblData],
    tolerance: float
) -> Tuple[List[StationComparison], Dict[str, Any]]:
    """Compare SETBL events between XFOIL and RustFoil."""
    comparisons = []
    
    # Index by (side, ibl)
    xf_by_key = {(s.side, s.ibl): s for s in xf_setbl}
    rf_by_key = {(s.side, s.ibl): s for s in rf_setbl}
    
    all_keys = sorted(set(xf_by_key.keys()) | set(rf_by_key.keys()))
    
    max_va_global = 0.0
    max_vb_global = 0.0
    max_vdel_global = 0.0
    first_exceed = None
    
    for key in all_keys:
        xf = xf_by_key.get(key)
        rf = rf_by_key.get(key)
        
        if xf is None or rf is None:
            # One side missing
            continue
        
        # Compare VA blocks
        if xf.va.size > 0 and rf.va.size > 0:
            # Align shapes
            min_rows = min(xf.va.shape[0], rf.va.shape[0])
            min_cols = min(xf.va.shape[1], rf.va.shape[1]) if len(xf.va.shape) > 1 else 1
            va_diff = np.abs(xf.va[:min_rows, :min_cols] - rf.va[:min_rows, :min_cols]).max()
        else:
            va_diff = 0.0
        
        # Compare VB blocks
        if xf.vb.size > 0 and rf.vb.size > 0:
            min_rows = min(xf.vb.shape[0], rf.vb.shape[0])
            min_cols = min(xf.vb.shape[1], rf.vb.shape[1]) if len(xf.vb.shape) > 1 else 1
            vb_diff = np.abs(xf.vb[:min_rows, :min_cols] - rf.vb[:min_rows, :min_cols]).max()
        else:
            vb_diff = 0.0
        
        # Compare VDEL (solution deltas)
        if xf.vdel.size > 0 and rf.vdel.size > 0:
            min_rows = min(xf.vdel.shape[0], rf.vdel.shape[0])
            min_cols = min(xf.vdel.shape[1], rf.vdel.shape[1]) if len(xf.vdel.shape) > 1 else 1
            vdel_diff = np.abs(xf.vdel[:min_rows, :min_cols] - rf.vdel[:min_rows, :min_cols]).max()
        else:
            vdel_diff = 0.0
        
        max_diff = max(va_diff, vb_diff, vdel_diff)
        passed = max_diff <= tolerance
        
        max_va_global = max(max_va_global, va_diff)
        max_vb_global = max(max_vb_global, vb_diff)
        max_vdel_global = max(max_vdel_global, vdel_diff)
        
        if not passed and first_exceed is None:
            first_exceed = {
                'side': key[0],
                'ibl': key[1],
                'va_diff': va_diff,
                'vb_diff': vb_diff,
                'vdel_diff': vdel_diff,
            }
        
        comparisons.append(StationComparison(
            station=key[1],
            side=key[0],
            max_va_diff=va_diff,
            max_vb_diff=vb_diff,
            max_vdel_diff=vdel_diff,
            max_vsrez_diff=0.0,  # Not in SETBL
            passed=passed,
        ))
    
    summary = {
        'total_stations': len(comparisons),
        'passed_stations': sum(1 for c in comparisons if c.passed),
        'max_va_diff': max_va_global,
        'max_vb_diff': max_vb_global,
        'max_vdel_diff': max_vdel_global,
        'first_exceed': first_exceed,
    }
    
    return comparisons, summary


def compare_bldif(
    xf_bldif: List[BldifData],
    rf_bldif: List[BldifData],
    tolerance: float
) -> Tuple[List[StationComparison], Dict[str, Any]]:
    """Compare BLDIF events between XFOIL and RustFoil."""
    comparisons = []
    
    # Index by (side, ibl)
    xf_by_key = {(b.side, b.ibl): b for b in xf_bldif}
    rf_by_key = {(b.side, b.ibl): b for b in rf_bldif}
    
    all_keys = sorted(set(xf_by_key.keys()) | set(rf_by_key.keys()))
    
    max_vs1_global = 0.0
    max_vs2_global = 0.0
    max_vsrez_global = 0.0
    first_exceed = None
    
    for key in all_keys:
        xf = xf_by_key.get(key)
        rf = rf_by_key.get(key)
        
        if xf is None or rf is None:
            continue
        
        # Compare VS1 (Jacobian wrt station 1)
        if xf.vs1.size > 0 and rf.vs1.size > 0:
            min_r = min(xf.vs1.shape[0], rf.vs1.shape[0])
            min_c = min(xf.vs1.shape[1], rf.vs1.shape[1]) if len(xf.vs1.shape) > 1 else 1
            vs1_diff = np.abs(xf.vs1[:min_r, :min_c] - rf.vs1[:min_r, :min_c]).max()
        else:
            vs1_diff = 0.0
        
        # Compare VS2 (Jacobian wrt station 2)
        if xf.vs2.size > 0 and rf.vs2.size > 0:
            min_r = min(xf.vs2.shape[0], rf.vs2.shape[0])
            min_c = min(xf.vs2.shape[1], rf.vs2.shape[1]) if len(xf.vs2.shape) > 1 else 1
            vs2_diff = np.abs(xf.vs2[:min_r, :min_c] - rf.vs2[:min_r, :min_c]).max()
        else:
            vs2_diff = 0.0
        
        # Compare VSREZ (residuals)
        if xf.vsrez.size > 0 and rf.vsrez.size > 0:
            min_len = min(len(xf.vsrez), len(rf.vsrez))
            vsrez_diff = np.abs(xf.vsrez[:min_len] - rf.vsrez[:min_len]).max()
        else:
            vsrez_diff = 0.0
        
        max_diff = max(vs1_diff, vs2_diff, vsrez_diff)
        passed = max_diff <= tolerance
        
        max_vs1_global = max(max_vs1_global, vs1_diff)
        max_vs2_global = max(max_vs2_global, vs2_diff)
        max_vsrez_global = max(max_vsrez_global, vsrez_diff)
        
        if not passed and first_exceed is None:
            first_exceed = {
                'side': key[0],
                'ibl': key[1],
                'flow_type': xf.flow_type,
                'vs1_diff': vs1_diff,
                'vs2_diff': vs2_diff,
                'vsrez_diff': vsrez_diff,
                'xf_vsrez': xf.vsrez.tolist() if xf.vsrez.size > 0 else [],
                'rf_vsrez': rf.vsrez.tolist() if rf.vsrez.size > 0 else [],
            }
        
        comparisons.append(StationComparison(
            station=key[1],
            side=key[0],
            max_va_diff=vs1_diff,  # Reusing fields for VS1/VS2
            max_vb_diff=vs2_diff,
            max_vdel_diff=0.0,
            max_vsrez_diff=vsrez_diff,
            passed=passed,
        ))
    
    summary = {
        'total_stations': len(comparisons),
        'passed_stations': sum(1 for c in comparisons if c.passed),
        'max_vs1_diff': max_vs1_global,
        'max_vs2_diff': max_vs2_global,
        'max_vsrez_diff': max_vsrez_global,
        'first_exceed': first_exceed,
    }
    
    return comparisons, summary


def print_setbl_report(
    comparisons: List[StationComparison],
    summary: Dict[str, Any],
    tolerance: float,
    verbose: bool = False
):
    """Print SETBL comparison report."""
    print("\n" + "=" * 70)
    print("SETBL COMPARISON (Newton System Setup)")
    print("=" * 70)
    
    print(f"\nStations compared: {summary['total_stations']}")
    print(f"Stations passed:   {summary['passed_stations']}")
    print(f"Tolerance:         {tolerance}")
    
    print("\n--- Global Statistics ---")
    print(f"  Max VA difference:   {summary['max_va_diff']:.6e}")
    print(f"  Max VB difference:   {summary['max_vb_diff']:.6e}")
    print(f"  Max VDEL difference: {summary['max_vdel_diff']:.6e}")
    
    if summary['first_exceed']:
        fe = summary['first_exceed']
        print(f"\n--- First Exceedance ---")
        print(f"  Side {fe['side']}, Station {fe['ibl']}")
        print(f"  VA diff:   {fe['va_diff']:.6e}")
        print(f"  VB diff:   {fe['vb_diff']:.6e}")
        print(f"  VDEL diff: {fe['vdel_diff']:.6e}")
    
    if verbose:
        print("\n--- Station Details ---")
        print(f"{'Side':>4} {'IBL':>4} | {'VA diff':>12} | {'VB diff':>12} | {'VDEL diff':>12} | Status")
        print("-" * 65)
        for c in sorted(comparisons, key=lambda x: (x.side, x.station)):
            status = "✓" if c.passed else "✗"
            print(f"{c.side:4d} {c.station:4d} | {c.max_va_diff:12.4e} | {c.max_vb_diff:12.4e} | {c.max_vdel_diff:12.4e} | {status}")


def print_bldif_report(
    comparisons: List[StationComparison],
    summary: Dict[str, Any],
    tolerance: float,
    verbose: bool = False
):
    """Print BLDIF comparison report."""
    print("\n" + "=" * 70)
    print("BLDIF COMPARISON (Jacobians and Residuals)")
    print("=" * 70)
    
    print(f"\nStations compared: {summary['total_stations']}")
    print(f"Stations passed:   {summary['passed_stations']}")
    print(f"Tolerance:         {tolerance}")
    
    print("\n--- Global Statistics ---")
    print(f"  Max VS1 difference:   {summary['max_vs1_diff']:.6e}")
    print(f"  Max VS2 difference:   {summary['max_vs2_diff']:.6e}")
    print(f"  Max VSREZ difference: {summary['max_vsrez_diff']:.6e}")
    
    if summary['first_exceed']:
        fe = summary['first_exceed']
        print(f"\n--- First Exceedance ---")
        print(f"  Side {fe['side']}, Station {fe['ibl']} (flow_type={fe['flow_type']})")
        print(f"  VS1 diff:   {fe['vs1_diff']:.6e}")
        print(f"  VS2 diff:   {fe['vs2_diff']:.6e}")
        print(f"  VSREZ diff: {fe['vsrez_diff']:.6e}")
        print(f"  XFOIL VSREZ:    {fe['xf_vsrez']}")
        print(f"  RustFoil VSREZ: {fe['rf_vsrez']}")
    
    if verbose:
        print("\n--- Station Details ---")
        print(f"{'Side':>4} {'IBL':>4} | {'VS1 diff':>12} | {'VS2 diff':>12} | {'VSREZ diff':>12} | Status")
        print("-" * 70)
        for c in sorted(comparisons, key=lambda x: (x.side, x.station)):
            status = "✓" if c.passed else "✗"
            print(f"{c.side:4d} {c.station:4d} | {c.max_va_diff:12.4e} | {c.max_vb_diff:12.4e} | {c.max_vsrez_diff:12.4e} | {status}")


def print_convergence_comparison(xf_results: List[Dict], rf_results: List[Dict]):
    """Print convergence history comparison."""
    print("\n" + "=" * 70)
    print("CONVERGENCE HISTORY")
    print("=" * 70)
    
    max_iter = max(len(xf_results), len(rf_results))
    
    print(f"\n{'Iter':>4} | {'XF RMS':>12} | {'RF RMS':>12} | {'XF Max':>12} | {'RF Max':>12} | {'XF CL':>8} | {'RF CL':>8}")
    print("-" * 90)
    
    for i in range(max_iter):
        xf = xf_results[i] if i < len(xf_results) else {}
        rf = rf_results[i] if i < len(rf_results) else {}
        
        xf_rms = f"{xf.get('rms_residual', 0):.4e}" if xf else "---"
        rf_rms = f"{rf.get('rms_residual', 0):.4e}" if rf else "---"
        xf_max = f"{xf.get('max_residual', 0):.4e}" if xf else "---"
        rf_max = f"{rf.get('max_residual', 0):.4e}" if rf else "---"
        xf_cl = f"{xf.get('CL', 0):.5f}" if xf else "---"
        rf_cl = f"{rf.get('CL', 0):.5f}" if rf else "---"
        
        print(f"{i+1:4d} | {xf_rms:>12} | {rf_rms:>12} | {xf_max:>12} | {rf_max:>12} | {xf_cl:>8} | {rf_cl:>8}")


def print_detailed_station(
    xf_setbl: List[SetblData],
    rf_setbl: List[SetblData],
    xf_bldif: List[BldifData],
    rf_bldif: List[BldifData],
    side: int,
    ibl: int
):
    """Print detailed comparison for a specific station."""
    print("\n" + "=" * 70)
    print(f"DETAILED STATION: Side {side}, IBL {ibl}")
    print("=" * 70)
    
    # Find matching SETBL
    xf_s = next((s for s in xf_setbl if s.side == side and s.ibl == ibl), None)
    rf_s = next((s for s in rf_setbl if s.side == side and s.ibl == ibl), None)
    
    if xf_s and rf_s:
        print("\n--- VA Block ---")
        print("XFOIL:")
        print(xf_s.va)
        print("RustFoil:")
        print(rf_s.va)
        
        print("\n--- VB Block ---")
        print("XFOIL:")
        print(xf_s.vb)
        print("RustFoil:")
        print(rf_s.vb)
        
        print("\n--- VDEL (Solution Deltas) ---")
        print("XFOIL:")
        print(xf_s.vdel)
        print("RustFoil:")
        print(rf_s.vdel)
    
    # Find matching BLDIF
    xf_b = next((b for b in xf_bldif if b.side == side and b.ibl == ibl), None)
    rf_b = next((b for b in rf_bldif if b.side == side and b.ibl == ibl), None)
    
    if xf_b and rf_b:
        print("\n--- VS1 Jacobian (wrt station 1) ---")
        print("XFOIL:")
        print(xf_b.vs1)
        print("RustFoil:")
        print(rf_b.vs1)
        
        print("\n--- VS2 Jacobian (wrt station 2) ---")
        print("XFOIL:")
        print(xf_b.vs2)
        print("RustFoil:")
        print(rf_b.vs2)
        
        print("\n--- VSREZ (Residuals) ---")
        print(f"XFOIL:    {xf_b.vsrez}")
        print(f"RustFoil: {rf_b.vsrez}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare Newton system between XFOIL and RustFoil',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json --iteration 1 --tolerance 0.01
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json --all-iterations
    python scripts/compare_newton_iteration.py xfoil_debug.json rustfoil_debug.json --station 1:38 --verbose
        """
    )
    parser.add_argument('xfoil_json', help='Path to XFOIL debug JSON file')
    parser.add_argument('rustfoil_json', help='Path to RustFoil debug JSON file')
    parser.add_argument('--iteration', '-i', type=int, default=1,
                        help='Newton iteration to compare (default: 1)')
    parser.add_argument('--all-iterations', '-a', action='store_true',
                        help='Compare all iterations')
    parser.add_argument('--tolerance', '-t', type=float, default=0.01,
                        help='Relative tolerance for comparison (default: 0.01)')
    parser.add_argument('--station', '-s', type=str,
                        help='Show detailed comparison for station (format: side:ibl, e.g., 1:38)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show all station comparisons')
    parser.add_argument('--convergence-only', '-c', action='store_true',
                        help='Only show convergence history comparison')
    
    args = parser.parse_args()
    
    # Load debug files
    print(f"Loading XFOIL debug: {args.xfoil_json}")
    xf_events = load_debug_events(Path(args.xfoil_json))
    print(f"  {len(xf_events)} events")
    
    print(f"Loading RustFoil debug: {args.rustfoil_json}")
    rf_events = load_debug_events(Path(args.rustfoil_json))
    print(f"  {len(rf_events)} events")
    
    # Extract convergence history
    xf_results = extract_viscal_results(xf_events)
    rf_results = extract_viscal_results(rf_events)
    
    if args.convergence_only:
        print_convergence_comparison(xf_results, rf_results)
        return 0
    
    # Determine iterations to compare
    if args.all_iterations:
        iterations = sorted(set(
            e.get('iteration') for e in xf_events 
            if e.get('subroutine') in ('SETBL', 'BLDIF') and e.get('iteration')
        ))
    else:
        iterations = [args.iteration]
    
    print(f"\nComparing iterations: {iterations}")
    
    overall_passed = True
    
    for iteration in iterations:
        print(f"\n{'='*70}")
        print(f"ITERATION {iteration}")
        print('='*70)
        
        # Extract events for this iteration
        xf_setbl = extract_setbl_events(xf_events, iteration)
        rf_setbl = extract_setbl_events(rf_events, iteration)
        
        xf_bldif = extract_bldif_events(xf_events, iteration)
        rf_bldif = extract_bldif_events(rf_events, iteration)
        
        print(f"\nSETBL events: XFOIL={len(xf_setbl)}, RustFoil={len(rf_setbl)}")
        print(f"BLDIF events: XFOIL={len(xf_bldif)}, RustFoil={len(rf_bldif)}")
        
        # Compare SETBL
        if xf_setbl and rf_setbl:
            setbl_comparisons, setbl_summary = compare_setbl(xf_setbl, rf_setbl, args.tolerance)
            print_setbl_report(setbl_comparisons, setbl_summary, args.tolerance, args.verbose)
            if setbl_summary['first_exceed']:
                overall_passed = False
        else:
            print("\n  No SETBL events to compare")
        
        # Compare BLDIF
        if xf_bldif and rf_bldif:
            bldif_comparisons, bldif_summary = compare_bldif(xf_bldif, rf_bldif, args.tolerance)
            print_bldif_report(bldif_comparisons, bldif_summary, args.tolerance, args.verbose)
            if bldif_summary['first_exceed']:
                overall_passed = False
        else:
            print("\n  No BLDIF events to compare")
        
        # Show detailed station if requested
        if args.station:
            parts = args.station.split(':')
            if len(parts) == 2:
                side = int(parts[0])
                ibl = int(parts[1])
                print_detailed_station(xf_setbl, rf_setbl, xf_bldif, rf_bldif, side, ibl)
    
    # Show convergence comparison
    if xf_results or rf_results:
        print_convergence_comparison(xf_results, rf_results)
    
    # Summary
    print("\n" + "=" * 70)
    if overall_passed:
        print("✓ All comparisons within tolerance")
        return 0
    else:
        print("✗ Some comparisons exceeded tolerance")
        return 1


if __name__ == '__main__':
    exit(main())
