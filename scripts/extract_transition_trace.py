#!/usr/bin/env python3
"""
Extract transition trace data from XFOIL instrumented debug output.

This script processes the xfoil_debug.json file and extracts detailed
DAMPL, AXSET, and TRCHEK data for transition analysis.

Usage:
    python extract_transition_trace.py xfoil_debug.json
    python extract_transition_trace.py xfoil_debug.json -o testdata/transition_trace.json
"""

import json
import argparse
from pathlib import Path
from collections import defaultdict


def load_xfoil_events(json_path):
    """Load events from xfoil_debug.json."""
    with open(json_path) as f:
        return json.load(f)['events']


def extract_transition_trace(events):
    """
    Extract transition-related events organized by station.
    """
    # Extract metadata from VISCAL event
    viscal_events = [e for e in events if e.get('subroutine') == 'VISCAL']
    metadata = {}
    if viscal_events:
        v = viscal_events[0]
        metadata = {
            "reynolds": v.get('reynolds', 1e6),
            "mach": v.get('mach', 0.0),
            "ncrit": v.get('ncrit', 9.0),
            "alpha_rad": v.get('alpha_rad', 0.0)
        }
    
    # Extract TRCHEK events (final results per station)
    trchek_events = [e for e in events if e.get('subroutine') == 'TRCHEK']
    
    # Extract TRCHEK2_FINAL events (detailed per-station results)
    trchek2_final = [e for e in events if e.get('subroutine') == 'TRCHEK2_FINAL']
    # Extract TRCHEK2_ITER events (intermediate data, includes wf and T1/T2)
    trchek2_iter = [e for e in events if e.get('subroutine') == 'TRCHEK2_ITER']
    
    # Extract AXSET events
    axset_events = [e for e in events if e.get('subroutine') == 'AXSET']
    
    # Extract DAMPL events
    dampl_events = [e for e in events if e.get('subroutine') == 'DAMPL']
    
    # Group by side and station
    stations_by_side = defaultdict(lambda: defaultdict(dict))
    
    # Process TRCHEK events (use the last iteration for each station)
    for e in trchek_events:
        side = e['side']
        ibl = e['ibl']
        # Keep the latest (highest iteration) data
        current = stations_by_side[side][ibl]
        if not current or e.get('iteration', 0) >= current.get('iteration', 0):
            stations_by_side[side][ibl] = {
                'ibl': ibl,
                'side': side,
                'x1': e.get('x1'),
                'x2': e.get('x2'),
                'ampl1': e.get('ampl1'),
                'ampl2': e.get('ampl2'),
                'Hk1': e.get('Hk1'),
                'Hk2': e.get('Hk2'),
                'Rt1': e.get('Rt1'),
                'Rt2': e.get('Rt2'),
                'transition': e.get('transition', False),
                'x_transition': e.get('x_transition'),
                'iteration': e.get('iteration', 0)
            }
    
    # Index TRCHEK2_ITER by (side, ibl) using latest trchek_iter
    trchek2_iter_map = {}
    for e in trchek2_iter:
        key = (e.get('side'), e.get('ibl'))
        if key not in trchek2_iter_map or e.get('trchek_iter', 0) >= trchek2_iter_map[key].get('trchek_iter', 0):
            trchek2_iter_map[key] = e

    # Add TRCHEK2_FINAL data
    for e in trchek2_final:
        side = e['side']
        ibl = e['ibl']
        if ibl in stations_by_side[side]:
            stations_by_side[side][ibl]['ax_final'] = e.get('ax_final')
            stations_by_side[side][ibl]['converged'] = e.get('converged')
            stations_by_side[side][ibl]['n_iterations'] = e.get('n_iterations')
            stations_by_side[side][ibl]['xt_final'] = e.get('xt_final', e.get('xt'))
            stations_by_side[side][ibl]['ampl2_final'] = e.get('ampl2_final', e.get('ampl2'))
            # Optional RustFoil extra fields
            for field in ['tt', 'dt', 'ut', 'Hk_t', 'Rt_t', 'St', 'Cq_t']:
                if field in e:
                    stations_by_side[side][ibl][field] = e.get(field)

    # Add TRCHEK2_ITER interpolation data (WF, TT, DT, UT)
    for side in [1, 2]:
        for ibl, station in stations_by_side[side].items():
            key = (side, ibl)
            iter_ev = trchek2_iter_map.get(key)
            if not iter_ev:
                continue
            station['wf1'] = iter_ev.get('wf1')
            station['wf2'] = iter_ev.get('wf2')
            station['xt_iter'] = iter_ev.get('xt')
            station['T1'] = iter_ev.get('T1')
            station['T2'] = iter_ev.get('T2')
            t1 = iter_ev.get('T1')
            t2 = iter_ev.get('T2')
            wf1 = iter_ev.get('wf1')
            wf2 = iter_ev.get('wf2')
            if t1 is not None and t2 is not None and wf1 is not None and wf2 is not None:
                station['tt_iter'] = t1 * wf1 + t2 * wf2
    
    # Build output structure
    result = {
        "metadata": metadata,
        "source": "xfoil_instrumented",
        "description": "Transition trace data for debugging N-factor evolution",
        "sides": {}
    }
    
    for side in [1, 2]:
        side_data = {"stations": []}
        
        for ibl in sorted(stations_by_side[side].keys()):
            station = stations_by_side[side][ibl]
            side_data["stations"].append(station)
        
        result["sides"][str(side)] = side_data
    
    return result


def print_summary(data):
    """Print a summary of the extracted data."""
    print("=" * 70)
    print("Transition Trace Data Summary")
    print("=" * 70)
    print()
    print(f"Metadata:")
    for k, v in data['metadata'].items():
        print(f"  {k}: {v}")
    print()
    
    for side in ["1", "2"]:
        side_name = "Upper" if side == "1" else "Lower"
        stations = data['sides'][side]['stations']
        print(f"{side_name} Surface (Side {side}): {len(stations)} stations")
        
        # Find transition station
        trans_station = None
        for s in stations:
            if s.get('transition'):
                trans_station = s
                break
        
        if trans_station:
            print(f"  Transition at IBL {trans_station['ibl']}, x={trans_station.get('x_transition', 'N/A')}")
        else:
            print(f"  No transition detected")
        
        # Show N-factor evolution
        print(f"\n  N-factor evolution:")
        print(f"  {'IBL':>4} {'x':>10} {'ampl1':>12} {'ampl2':>12} {'ax':>12} {'Hk2':>8} {'Rt2':>10}")
        for s in stations[:20]:
            ax = s.get('ax_final', 0) or 0
            print(f"  {s['ibl']:4d} {s.get('x2', 0):10.4f} {s.get('ampl1', 0):12.4e} {s.get('ampl2', 0):12.4e} {ax:12.4e} {s.get('Hk2', 0):8.4f} {s.get('Rt2', 0):10.2f}")
        if len(stations) > 20:
            print(f"  ... ({len(stations) - 20} more stations)")
        print()


def main():
    parser = argparse.ArgumentParser(
        description="Extract transition trace data from XFOIL debug output"
    )
    parser.add_argument("input", help="Path to xfoil_debug.json")
    parser.add_argument("-o", "--output", 
                       help="Output JSON file (default: testdata/transition_trace.json)")
    parser.add_argument("--summary", action="store_true",
                       help="Print summary only, don't write output file")
    
    args = parser.parse_args()
    
    # Load events
    print(f"Loading {args.input}...")
    events = load_xfoil_events(args.input)
    print(f"Loaded {len(events)} events")
    
    # Extract transition trace
    data = extract_transition_trace(events)
    
    # Print summary
    print_summary(data)
    
    # Write output
    if not args.summary:
        output_path = args.output or "testdata/transition_trace.json"
        output_path = Path(output_path)
        
        # Ensure directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
