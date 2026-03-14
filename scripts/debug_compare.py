#!/usr/bin/env python3
"""
Side-by-side comparison of XFOIL and RustFoil debug events.

Aligns events by type and station, highlights divergences.
"""

import json
import argparse
from collections import defaultdict
from pathlib import Path


def load_events(path):
    """Load events from JSON file."""
    with open(path) as f:
        data = json.load(f)
    return data.get('events', data)


def group_by_type(events):
    """Group events by subroutine name."""
    groups = defaultdict(list)
    for ev in events:
        groups[ev['subroutine']].append(ev)
    return groups


def compare_mrchue(xfoil_events, rust_events, tolerance=0.01):
    """Compare initial march profiles."""
    print("\n" + "=" * 70)
    print("MRCHUE (Initial March) Comparison")
    print("=" * 70)

    # Group by station
    xf_by_station = {(e.get('side', 1), e['ibl']): e for e in xfoil_events}
    rf_by_station = {(e.get('side', 1), e['ibl']): e for e in rust_events}

    all_stations = sorted(set(xf_by_station.keys()) | set(rf_by_station.keys()))

    print(f"XFOIL has {len(xfoil_events)} MRCHUE events")
    print(f"RustFoil has {len(rust_events)} MRCHUE events")
    print()

    fields = ['theta', 'delta_star', 'Hk', 'Cf', 'Ue', 'x']
    errors = []

    # Print header
    print(f"{'Side':<5} {'IBL':<5} {'Field':<12} {'XFOIL':>14} {'RustFoil':>14} {'Diff%':>10} {'OK':<4}")
    print("-" * 70)

    for station in all_stations[:30]:  # First 30 stations
        xf = xf_by_station.get(station, {})
        rf = rf_by_station.get(station, {})

        if not xf:
            print(f"{station[0]:<5} {station[1]:<5} {'MISSING in XFOIL':<40}")
            continue
        if not rf:
            print(f"{station[0]:<5} {station[1]:<5} {'MISSING in RustFoil':<40}")
            continue

        for field in fields:
            xf_val = xf.get(field, 0)
            rf_val = rf.get(field, 0)

            if abs(xf_val) > 1e-15:
                diff_pct = abs(xf_val - rf_val) / abs(xf_val) * 100
            else:
                diff_pct = 0 if abs(rf_val) < 1e-15 else 100

            status = "OK" if diff_pct < tolerance * 100 else "BAD"
            if diff_pct >= tolerance * 100:
                errors.append((station, field, xf_val, rf_val, diff_pct))

            # Only print first field for each station, or bad ones
            print(f"{station[0]:<5} {station[1]:<5} {field:<12} {xf_val:>14.6e} {rf_val:>14.6e} {diff_pct:>9.2f}% {status:<4}")

    print(f"\nTotal errors (>{tolerance*100}% diff): {len(errors)}")
    return errors


def compare_blvar(xfoil_events, rust_events, iteration=1, tolerance=0.01):
    """Compare BLVAR secondary variables for a specific iteration."""
    print(f"\n" + "=" * 70)
    print(f"BLVAR (Secondary Variables) - Iteration {iteration}")
    print("=" * 70)

    # Filter by iteration
    xf_iter = [e for e in xfoil_events if e.get('iteration') == iteration]
    rf_iter = [e for e in rust_events if e.get('iteration') == iteration]

    print(f"XFOIL has {len(xf_iter)} BLVAR events for iteration {iteration}")
    print(f"RustFoil has {len(rf_iter)} BLVAR events for iteration {iteration}")

    if not xf_iter:
        print("No XFOIL BLVAR events found")
        return []
    if not rf_iter:
        print("No RustFoil BLVAR events found")
        return []

    # Group by station
    xf_by_station = {(e.get('side', 1), e['ibl']): e for e in xf_iter}
    rf_by_station = {(e.get('side', 1), e['ibl']): e for e in rf_iter}

    errors = []
    fields = ['H', 'Hk', 'Hs', 'Cf', 'Cd', 'Rtheta']

    print(f"\n{'Station':<12} {'Field':<10} {'XFOIL':>14} {'RustFoil':>14} {'Diff%':>10} {'OK':<4}")
    print("-" * 70)

    for station in sorted(xf_by_station.keys())[:20]:
        xf = xf_by_station.get(station, {})
        rf = rf_by_station.get(station, {})

        if not rf:
            print(f"{str(station):<12} MISSING in RustFoil")
            continue

        xf_out = xf.get('output', {})
        rf_out = rf.get('output', {})

        for field in fields:
            xf_val = xf_out.get(field, 0)
            rf_val = rf_out.get(field, 0)

            if abs(xf_val) > 1e-15:
                diff_pct = abs(xf_val - rf_val) / abs(xf_val) * 100
            else:
                diff_pct = 0

            status = "OK" if diff_pct < tolerance * 100 else "BAD"
            print(f"{str(station):<12} {field:<10} {xf_val:>14.6e} {rf_val:>14.6e} {diff_pct:>9.2f}% {status:<4}")

            if diff_pct >= tolerance * 100:
                errors.append((station, field, xf_val, rf_val, diff_pct))

    print(f"\nTotal errors (>{tolerance*100}% diff): {len(errors)}")
    return errors


def compare_update(xfoil_events, rust_events, iteration=1, tolerance=0.10):
    """Compare UPDATE deltas - this is where divergence likely starts."""
    print(f"\n" + "=" * 70)
    print(f"UPDATE (Newton Deltas) - Iteration {iteration}")
    print("=" * 70)

    xf_iter = [e for e in xfoil_events if e.get('iteration') == iteration]
    rf_iter = [e for e in rust_events if e.get('iteration') == iteration]

    print(f"XFOIL has {len(xf_iter)} UPDATE events for iteration {iteration}")
    print(f"RustFoil has {len(rf_iter)} UPDATE events for iteration {iteration}")

    if not xf_iter:
        print("No XFOIL UPDATE events found")
        return []
    if not rf_iter:
        print("No RustFoil UPDATE events found - need to instrument update_stations()")
        return []

    fields = ['delta_ctau', 'delta_theta', 'delta_mass', 'delta_Ue', 'relaxation']

    print(f"\n{'Station':<12} {'Field':<14} {'XFOIL':>14} {'RustFoil':>14}")
    print("-" * 60)

    # Compare first few stations
    for i, xf in enumerate(xf_iter[:15]):
        station = (xf.get('side', 1), xf.get('ibl', i))
        rf = rf_iter[i] if i < len(rf_iter) else None

        for field in fields:
            xf_val = xf.get(field, 0)
            rf_val = rf.get(field, 0) if rf else 0

            print(f"{str(station):<12} {field:<14} {xf_val:>14.6e} {rf_val:>14.6e}")

        print()  # Blank line between stations

    return []


def compare_bldif(xfoil_events, rust_events, iteration=1):
    """Compare BLDIF Jacobian matrices."""
    print(f"\n" + "=" * 70)
    print(f"BLDIF (Jacobian Matrices) - Iteration {iteration}")
    print("=" * 70)

    xf_iter = [e for e in xfoil_events if e.get('iteration') == iteration]
    rf_iter = [e for e in rust_events if e.get('iteration') == iteration]

    print(f"XFOIL has {len(xf_iter)} BLDIF events for iteration {iteration}")
    print(f"RustFoil has {len(rf_iter)} BLDIF events for iteration {iteration}")

    if not xf_iter or not rf_iter:
        return []

    # Compare first few stations
    for i, xf in enumerate(xf_iter[:5]):
        station = (xf.get('side', 1), xf.get('ibl', i))
        rf = rf_iter[i] if i < len(rf_iter) else None

        print(f"\n--- Station {station} ---")

        if not rf:
            print("MISSING in RustFoil")
            continue

        # Compare residuals
        xf_rez = xf.get('VSREZ', [0, 0, 0, 0])
        rf_rez = rf.get('VSREZ', [0, 0, 0, 0])

        print("VSREZ (residuals):")
        print(f"  XFOIL:   [{xf_rez[0]:>10.4e}, {xf_rez[1]:>10.4e}, {xf_rez[2]:>10.4e}, {xf_rez[3]:>10.4e}]")
        print(f"  RustFoil:[{rf_rez[0]:>10.4e}, {rf_rez[1]:>10.4e}, {rf_rez[2]:>10.4e}, {rf_rez[3]:>10.4e}]")

        # Compare VS2 diagonal (most important)
        xf_vs2 = xf.get('VS2', [[0]*5]*4)
        rf_vs2 = rf.get('VS2', [[0]*5]*4)

        print("VS2 diagonal (self-derivatives):")
        for row in range(3):
            xf_diag = xf_vs2[row][row] if row < len(xf_vs2) and row < len(xf_vs2[row]) else 0
            rf_diag = rf_vs2[row][row] if row < len(rf_vs2) and row < len(rf_vs2[row]) else 0
            print(f"  [{row},{row}] XFOIL={xf_diag:>12.4e}  Rust={rf_diag:>12.4e}")

    return []


def compare_convergence(xfoil_events, rust_events):
    """Compare iteration convergence history."""
    print(f"\n" + "=" * 70)
    print("VISCAL_RESULT (Convergence History)")
    print("=" * 70)

    print(f"XFOIL has {len(xfoil_events)} VISCAL_RESULT events")
    print(f"RustFoil has {len(rust_events)} VISCAL_RESULT events")
    print()

    print(f"{'Iter':<6} {'XFOIL RMS':>12} {'Rust RMS':>12} {'XFOIL CL':>10} {'Rust CL':>10} {'XFOIL CD':>12} {'Rust CD':>12}")
    print("-" * 80)

    for xf in xfoil_events:
        it = xf.get('iteration', 0)
        rf = next((e for e in rust_events if e.get('iteration') == it), None)

        xf_rms = xf.get('rms_residual', 0)
        xf_cl = xf.get('CL', 0)
        xf_cd = xf.get('CD', 0)

        rf_rms = rf.get('rms_residual', 0) if rf else 0
        rf_cl = rf.get('CL', 0) if rf else 0
        rf_cd = rf.get('CD', 0) if rf else 0

        print(f"{it:<6} {xf_rms:>12.4e} {rf_rms:>12.4e} {xf_cl:>10.4f} {rf_cl:>10.4f} {xf_cd:>12.6e} {rf_cd:>12.6e}")


def compare_qdcalc(xfoil_events, rust_events):
    """Compare QDCALC DIJ matrix samples."""
    print(f"\n" + "=" * 70)
    print("QDCALC (DIJ Influence Matrix)")
    print("=" * 70)

    xf = xfoil_events[0] if xfoil_events else {}
    rf = rust_events[0] if rust_events else {}

    if not xf:
        print("No XFOIL QDCALC events")
        return
    if not rf:
        print("No RustFoil QDCALC events")
        return

    print(f"Dimensions: XFOIL n={xf.get('n_total', 0)}, RustFoil n={rf.get('n_total', 0)}")
    print(f"Airfoil pts: XFOIL={xf.get('n_airfoil', 0)}, RustFoil={rf.get('n_airfoil', 0)}")
    print(f"Wake pts: XFOIL={xf.get('n_wake', 0)}, RustFoil={rf.get('n_wake', 0)}")

    xf_diag = xf.get('DIJ_diagonal_sample', [])
    rf_diag = rf.get('DIJ_diagonal_sample', [])

    print("\nDIJ Diagonal elements (first 10):")
    print(f"{'Index':<6} {'XFOIL':>14} {'RustFoil':>14} {'Diff%':>10}")
    print("-" * 50)

    for i in range(min(10, len(xf_diag), len(rf_diag))):
        xf_val = xf_diag[i]
        rf_val = rf_diag[i]
        diff = abs(xf_val - rf_val) / abs(xf_val) * 100 if xf_val != 0 else 0
        print(f"{i:<6} {xf_val:>14.6e} {rf_val:>14.6e} {diff:>9.2f}%")


def main():
    parser = argparse.ArgumentParser(description="Compare XFOIL and RustFoil debug output")
    parser.add_argument('xfoil_json', help='Path to XFOIL debug JSON')
    parser.add_argument('rust_json', help='Path to RustFoil debug JSON')
    parser.add_argument('--stage', choices=['mrchue', 'blvar', 'bldif', 'update', 'convergence', 'qdcalc', 'all'],
                        default='all', help='Which stage to compare')
    parser.add_argument('--iteration', type=int, default=1, help='Iteration to compare (for blvar, bldif, update)')
    parser.add_argument('--tolerance', type=float, default=0.01, help='Tolerance for "OK" (default 1%%)')
    args = parser.parse_args()

    print(f"Loading XFOIL: {args.xfoil_json}")
    xf_events = load_events(args.xfoil_json)
    print(f"  {len(xf_events)} events")

    print(f"Loading RustFoil: {args.rust_json}")
    rf_events = load_events(args.rust_json)
    print(f"  {len(rf_events)} events")

    xf_groups = group_by_type(xf_events)
    rf_groups = group_by_type(rf_events)

    print("\n" + "=" * 70)
    print("Event Counts")
    print("=" * 70)
    all_types = sorted(set(xf_groups.keys()) | set(rf_groups.keys()))
    for t in all_types:
        xf_count = len(xf_groups.get(t, []))
        rf_count = len(rf_groups.get(t, []))
        status = "OK" if xf_count > 0 and rf_count > 0 else "MISSING"
        print(f"  {t:<18} XFOIL={xf_count:>5}  Rust={rf_count:>5}  {status}")

    if args.stage in ['mrchue', 'all']:
        compare_mrchue(xf_groups.get('MRCHUE', []), rf_groups.get('MRCHUE', []), args.tolerance)

    if args.stage in ['blvar', 'all']:
        compare_blvar(xf_groups.get('BLVAR', []), rf_groups.get('BLVAR', []), args.iteration, args.tolerance)

    if args.stage in ['bldif', 'all']:
        compare_bldif(xf_groups.get('BLDIF', []), rf_groups.get('BLDIF', []), args.iteration)

    if args.stage in ['update', 'all']:
        compare_update(xf_groups.get('UPDATE', []), rf_groups.get('UPDATE', []), args.iteration)

    if args.stage in ['qdcalc', 'all']:
        compare_qdcalc(xf_groups.get('QDCALC', []), rf_groups.get('QDCALC', []))

    if args.stage in ['convergence', 'all']:
        compare_convergence(xf_groups.get('VISCAL_RESULT', []), rf_groups.get('VISCAL_RESULT', []))


if __name__ == '__main__':
    main()
