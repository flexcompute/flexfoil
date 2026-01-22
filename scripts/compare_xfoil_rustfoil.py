#!/usr/bin/env python3
"""
Compare XFOIL instrumented output with RustFoil.

Usage:
    python compare_xfoil_rustfoil.py xfoil_debug.json
    python compare_xfoil_rustfoil.py xfoil_debug.json --phase 1
    python compare_xfoil_rustfoil.py xfoil_debug.json --phase 3 --station 10
"""

import json
import argparse
import numpy as np
from pathlib import Path


def load_xfoil_events(json_path):
    """Load events from xfoil_debug.json."""
    with open(json_path) as f:
        return json.load(f)['events']


def filter_events(events, subroutine, iteration=None, side=None, ibl=None):
    """Filter events by criteria."""
    result = [e for e in events if e['subroutine'] == subroutine]
    if iteration is not None:
        result = [e for e in result if e.get('iteration') == iteration]
    if side is not None:
        result = [e for e in result if e.get('side') == side]
    if ibl is not None:
        result = [e for e in result if e.get('ibl') == ibl]
    return result


def summarize_events(events):
    """Print summary of event types."""
    counts = {}
    for e in events:
        sub = e['subroutine']
        counts[sub] = counts.get(sub, 0) + 1
    
    print("=" * 60)
    print("Event Summary")
    print("=" * 60)
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k:20s}: {v:6d}")
    print(f"  {'TOTAL':20s}: {len(events):6d}")
    print()


def phase1_closures(events):
    """Phase 1: Compare closure functions."""
    print("\n" + "=" * 60)
    print("Phase 1: Closure Function Values")
    print("=" * 60)
    
    # HKIN examples
    print("\n--- HKIN (H -> Hk) ---")
    hkin_events = filter_events(events, 'HKIN')[:10]
    print(f"{'H':>10s} {'Msq':>10s} {'Hk':>10s} {'Hk_H':>10s} {'Hk_Msq':>10s}")
    for ev in hkin_events:
        print(f"{ev['H']:10.4f} {ev['Msq']:10.4f} {ev['Hk']:10.4f} {ev['Hk_H']:10.4f} {ev['Hk_Msq']:10.4f}")
    
    # HSL examples
    print("\n--- HSL (Laminar Hs) ---")
    hsl_events = filter_events(events, 'HSL')[:10]
    print(f"{'Hk':>10s} {'Rtheta':>12s} {'Hs':>10s} {'Hs_Hk':>10s}")
    for ev in hsl_events:
        print(f"{ev['Hk']:10.4f} {ev['Rtheta']:12.2f} {ev['Hs']:10.4f} {ev['Hs_Hk']:10.4f}")
    
    # CFL examples
    print("\n--- CFL (Laminar Cf) ---")
    cfl_events = filter_events(events, 'CFL')[:10]
    print(f"{'Hk':>10s} {'Rtheta':>12s} {'Cf':>12s} {'Cf_Hk':>12s}")
    for ev in cfl_events:
        print(f"{ev['Hk']:10.4f} {ev['Rtheta']:12.2f} {ev['Cf']:12.6f} {ev['Cf_Hk']:12.6f}")
    
    # DAMPL examples
    print("\n--- DAMPL (Amplification Rate) ---")
    dampl_events = filter_events(events, 'DAMPL')[:10]
    print(f"{'Hk':>10s} {'theta':>12s} {'Rtheta':>12s} {'Ax':>12s}")
    for ev in dampl_events:
        print(f"{ev['Hk']:10.4f} {ev['theta']:12.6f} {ev['Rtheta']:12.2f} {ev['Ax']:12.4f}")


def phase2_blvar(events, station=None):
    """Phase 2: BLVAR secondary variables."""
    print("\n" + "=" * 60)
    print("Phase 2: BLVAR Secondary Variables")
    print("=" * 60)
    
    blvar_events = filter_events(events, 'BLVAR', iteration=1)
    if station is not None:
        blvar_events = [e for e in blvar_events if e.get('ibl') == station]
    
    for ev in blvar_events[:20]:
        print(f"\n--- Side {ev.get('side', '?')}, Station {ev.get('ibl', '?')}, Type {ev.get('flow_type', '?')} ---")
        inp = ev.get('input', {})
        out = ev.get('output', {})
        
        print("  Inputs:")
        print(f"    x        = {inp.get('x', 0):12.6e}")
        print(f"    u (Ue)   = {inp.get('u', 0):12.6e}")
        print(f"    theta    = {inp.get('theta', 0):12.6e}")
        print(f"    delta*   = {inp.get('delta_star', 0):12.6e}")
        print(f"    ctau     = {inp.get('ctau', 0):12.6e}")
        
        print("  Outputs:")
        print(f"    H        = {out.get('H', 0):12.6f}")
        print(f"    Hk       = {out.get('Hk', 0):12.6f}")
        print(f"    Hs       = {out.get('Hs', 0):12.6f}")
        print(f"    Cf       = {out.get('Cf', 0):12.6e}")
        print(f"    Cd       = {out.get('Cd', 0):12.6e}")
        print(f"    Rtheta   = {out.get('Rtheta', 0):12.2f}")


def phase3_jacobians(events, station=None):
    """Phase 3: BLDIF Jacobian analysis."""
    print("\n" + "=" * 60)
    print("Phase 3: BLDIF Jacobian Analysis")
    print("=" * 60)
    
    bldif_events = filter_events(events, 'BLDIF', iteration=1)
    if station is not None:
        bldif_events = [e for e in bldif_events if e.get('ibl') == station]
    
    for ev in bldif_events[:10]:
        print(f"\n--- Side {ev.get('side', '?')}, Station {ev.get('ibl', '?')}, Type {ev.get('flow_type', '?')} ---")
        
        vs1 = np.array(ev['VS1'])
        vs2 = np.array(ev['VS2'])
        vsrez = np.array(ev['VSREZ'])
        
        print(f"  VS1 (upstream station derivatives):")
        print(f"    Shape: {vs1.shape}")
        for i in range(min(3, vs1.shape[0])):
            print(f"    Row {i}: [{', '.join(f'{v:12.4e}' for v in vs1[i])}]")
        
        print(f"  VS2 (current station derivatives):")
        for i in range(min(3, vs2.shape[0])):
            print(f"    Row {i}: [{', '.join(f'{v:12.4e}' for v in vs2[i])}]")
        
        print(f"  VSREZ (residuals): [{', '.join(f'{v:12.4e}' for v in vsrez)}]")
        
        # Check for problematic patterns
        warnings = []
        
        # Check column 2 (delta* derivatives) - was previously zero in RustFoil
        if vs1.shape[1] > 2:
            col2_vs1 = vs1[:3, 2] if vs1.shape[0] >= 3 else vs1[:, 2]
            col2_vs2 = vs2[:3, 2] if vs2.shape[0] >= 3 else vs2[:, 2]
            if np.allclose(col2_vs1, 0, atol=1e-20):
                warnings.append("VS1 column 2 (δ* derivs) is all zeros!")
            if np.allclose(col2_vs2, 0, atol=1e-20):
                warnings.append("VS2 column 2 (δ* derivs) is all zeros!")
        
        # Check 3x3 VA block determinant
        if vs1.shape[0] >= 3 and vs1.shape[1] >= 3:
            va = vs1[:3, :3]
            try:
                det = np.linalg.det(va)
                if abs(det) < 1e-10:
                    warnings.append(f"Near-singular VA block! det={det:.2e}")
                else:
                    print(f"  VA block det: {det:.4e} (OK)")
            except:
                warnings.append("Could not compute VA determinant")
        
        for w in warnings:
            print(f"  ⚠️  WARNING: {w}")


def phase4_march(events):
    """Phase 4: Initial march state."""
    print("\n" + "=" * 60)
    print("Phase 4: MRCHUE Initial March")
    print("=" * 60)
    
    mrchue_events = filter_events(events, 'MRCHUE')
    
    # Group by side
    side1 = [e for e in mrchue_events if e.get('side') == 1]
    side2 = [e for e in mrchue_events if e.get('side') == 2]
    
    for side, side_events in [(1, side1), (2, side2)]:
        print(f"\n--- Side {side} ({len(side_events)} stations) ---")
        print(f"{'IBL':>4s} {'x':>12s} {'Ue':>12s} {'theta':>12s} {'delta*':>12s} {'Hk':>8s} {'Cf':>12s}")
        for ev in side_events[:15]:
            print(f"{ev.get('ibl', 0):4d} {ev['x']:12.6f} {ev['Ue']:12.6f} "
                  f"{ev['theta']:12.6e} {ev['delta_star']:12.6e} "
                  f"{ev['Hk']:8.4f} {ev['Cf']:12.6e}")


def phase5_iteration(events):
    """Phase 5: Newton iteration convergence."""
    print("\n" + "=" * 60)
    print("Phase 5: Iteration Convergence")
    print("=" * 60)
    
    viscal_results = filter_events(events, 'VISCAL_RESULT')
    
    print(f"{'Iter':>4s} {'RMS':>12s} {'Max':>12s} {'CL':>10s} {'CD':>12s} {'CM':>10s}")
    for ev in viscal_results:
        print(f"{ev['iteration']:4d} {ev['rms_residual']:12.4e} {ev['max_residual']:12.4e} "
              f"{ev['CL']:10.6f} {ev['CD']:12.6e} {ev['CM']:10.6f}")
    
    # UPDATE statistics per iteration
    print("\n--- UPDATE deltas (first iteration) ---")
    update_events = filter_events(events, 'UPDATE', iteration=1)[:10]
    print(f"{'Side':>4s} {'IBL':>4s} {'Δctau':>12s} {'Δtheta':>12s} {'Δmass':>12s} {'ΔUe':>12s} {'RLX':>8s}")
    for ev in update_events:
        print(f"{ev['side']:4d} {ev['ibl']:4d} {ev['delta_ctau']:12.4e} {ev['delta_theta']:12.4e} "
              f"{ev['delta_mass']:12.4e} {ev['delta_Ue']:12.4e} {ev['relaxation']:8.4f}")


def phase6_dij(events):
    """Phase 6: DIJ influence matrix."""
    print("\n" + "=" * 60)
    print("Phase 6: QDCALC DIJ Matrix")
    print("=" * 60)
    
    qdcalc_events = filter_events(events, 'QDCALC')
    
    for ev in qdcalc_events:
        print(f"  n_airfoil: {ev['n_airfoil']}")
        print(f"  n_wake:    {ev['n_wake']}")
        print(f"  n_total:   {ev['n_total']}")
        
        diag = ev.get('DIJ_diagonal_sample', [])
        row1 = ev.get('DIJ_row1_sample', [])
        
        print(f"  Diagonal (first 10): {[f'{v:.4e}' for v in diag]}")
        print(f"  Row 1 (first 10):    {[f'{v:.4e}' for v in row1]}")


def export_test_vectors(events, output_dir):
    """Export test vectors for RustFoil unit tests."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Export closure function test cases
    hkin_events = filter_events(events, 'HKIN')[:50]
    with open(output_dir / 'hkin_test_vectors.json', 'w') as f:
        json.dump(hkin_events, f, indent=2)
    
    hsl_events = filter_events(events, 'HSL')[:50]
    with open(output_dir / 'hsl_test_vectors.json', 'w') as f:
        json.dump(hsl_events, f, indent=2)
    
    cfl_events = filter_events(events, 'CFL')[:50]
    with open(output_dir / 'cfl_test_vectors.json', 'w') as f:
        json.dump(cfl_events, f, indent=2)
    
    blvar_events = filter_events(events, 'BLVAR', iteration=1)[:50]
    with open(output_dir / 'blvar_test_vectors.json', 'w') as f:
        json.dump(blvar_events, f, indent=2)
    
    # Export DAMPL - mix of zero and non-zero Ax values
    all_dampl = filter_events(events, 'DAMPL')
    dampl_zero = [e for e in all_dampl if e.get('Ax', 0) == 0][:25]
    dampl_nonzero = [e for e in all_dampl if e.get('Ax', 0) != 0][:25]
    dampl_events = dampl_zero + dampl_nonzero
    with open(output_dir / 'dampl_test_vectors.json', 'w') as f:
        json.dump(dampl_events, f, indent=2)
    
    print(f"Test vectors exported to {output_dir}")


def phase_mrchue_stations(xfoil_events, rustfoil_debug_path=None):
    """
    Compare theta/dstar evolution station by station.
    
    This compares XFOIL's MRCHUE_ITER final values with RustFoil's march results.
    """
    print("=" * 70)
    print("MRCHUE Station-by-Station Comparison")
    print("=" * 70)
    
    # Group MRCHUE_ITER events by side and ibl
    mrchue_events = [e for e in xfoil_events if e.get('subroutine') == 'MRCHUE_ITER']
    
    if not mrchue_events:
        print("No MRCHUE_ITER events found in XFOIL data")
        return
    
    # Organize by side -> ibl -> final iteration
    stations = {}
    for e in mrchue_events:
        side = e['side']
        ibl = e['ibl']
        key = (side, ibl)
        
        # Keep only the final (converged) iteration for each station
        if key not in stations or e['newton_iter'] > stations[key]['newton_iter']:
            stations[key] = e
    
    # Load RustFoil debug if available
    rust_stations = {}
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_data = json.load(f)
            # Extract march stations from RustFoil debug
            for e in rust_data.get('events', []):
                if e.get('subroutine') == 'MRCHUE':
                    side = e['side']
                    ibl = e['ibl']
                    rust_stations[(side, ibl)] = e
    
    for side in [1, 2]:
        side_name = "Upper" if side == 1 else "Lower"
        side_stations = {k: v for k, v in stations.items() if k[0] == side}
        
        if not side_stations:
            continue
        
        print(f"\n=== {side_name} Surface (Side {side}) ===")
        
        # Header depends on whether we have RustFoil data
        if rust_stations:
            print(f"{'IBL':>4} {'X':>10} {'XFOIL_θ':>12} {'RUST_θ':>12} {'Err%':>8}")
        else:
            print(f"{'IBL':>4} {'X':>10} {'θ_init':>12} {'θ_final':>12} {'N_iter':>6}")
        
        for (s, ibl), e in sorted(side_stations.items()):
            x = e['x']
            theta_in = e['input']['theta']
            theta_out = e['output']['theta']
            n_iter = e['newton_iter']
            
            if rust_stations and (side, ibl) in rust_stations:
                rust = rust_stations[(side, ibl)]
                rust_theta = rust['theta']
                err = abs(theta_out - rust_theta) / theta_out * 100 if theta_out != 0 else 0
                print(f"{ibl:4d} {x:10.4f} {theta_out:12.4e} {rust_theta:12.4e} {err:8.2f}")
            else:
                print(f"{ibl:4d} {x:10.4f} {theta_in:12.4e} {theta_out:12.4e} {n_iter:6d}")
        
        # Summary stats
        thetas = [e['output']['theta'] for (s, ibl), e in side_stations.items()]
        print(f"\n  Stations: {len(thetas)}")
        print(f"  θ range: [{min(thetas):.4e}, {max(thetas):.4e}]")
        
        # Check for transition (where flow type changes)
        for (s, ibl), e in sorted(side_stations.items()):
            if e.get('converged') and e.get('dmax', 1.0) < 1e-4:
                # Found converged station
                pass


def main():
    parser = argparse.ArgumentParser(description='Compare XFOIL instrumented output with RustFoil')
    parser.add_argument('xfoil_json', help='Path to xfoil_debug.json')
    parser.add_argument('--phase', type=int, default=0,
                        help='Comparison phase (0=summary, 1=closures, 2=blvar, 3=jacobians, 4=march, 5=iteration, 6=dij, 7=mrchue)')
    parser.add_argument('--station', type=int, default=None,
                        help='Filter to specific station number')
    parser.add_argument('--export', type=str, default=None,
                        help='Export test vectors to directory')
    parser.add_argument('--rustfoil', type=str, default=None,
                        help='Path to RustFoil debug JSON for comparison')
    args = parser.parse_args()
    
    print(f"Loading {args.xfoil_json}...")
    events = load_xfoil_events(args.xfoil_json)
    
    summarize_events(events)
    
    if args.export:
        export_test_vectors(events, args.export)
        return
    
    if args.phase == 0:
        print("Use --phase 1-7 for detailed comparison")
        print("  1: Closure functions (HKIN, HSL, CFL, DAMPL)")
        print("  2: BLVAR secondary variables")
        print("  3: BLDIF Jacobian matrices")
        print("  4: MRCHUE initial march")
        print("  5: Newton iteration convergence")
        print("  6: QDCALC DIJ matrix")
        print("  7: MRCHUE station-by-station theta/dstar")
    elif args.phase == 1:
        phase1_closures(events)
    elif args.phase == 2:
        phase2_blvar(events, args.station)
    elif args.phase == 3:
        phase3_jacobians(events, args.station)
    elif args.phase == 4:
        phase4_march(events)
    elif args.phase == 5:
        phase5_iteration(events)
    elif args.phase == 6:
        phase6_dij(events)
    elif args.phase == 7:
        phase_mrchue_stations(events, args.rustfoil)


if __name__ == '__main__':
    main()
