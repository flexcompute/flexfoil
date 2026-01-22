#!/usr/bin/env python3
"""
Extract MRCHUE Newton iteration data from XFOIL instrumented debug output.

This script processes the xfoil_debug.json file and extracts the BL station
data organized by station, showing the Newton iteration progression at each
station.

Usage:
    python extract_mrchue_iterations.py xfoil_debug.json
    python extract_mrchue_iterations.py xfoil_debug.json -o testdata/mrchue_iterations.json
"""

import json
import argparse
from pathlib import Path
from collections import defaultdict


def load_xfoil_events(json_path):
    """Load events from xfoil_debug.json."""
    with open(json_path) as f:
        return json.load(f)['events']


def extract_mrchue_iterations(events):
    """
    Extract MRCHUE_ITER events and organize by station.
    
    Returns a structure like:
    {
        "metadata": {...},
        "sides": {
            1: {"stations": [...]},
            2: {"stations": [...]}
        }
    }
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
    
    # Extract MRCHUE_ITER events
    mrchue_events = [e for e in events if e.get('subroutine') == 'MRCHUE_ITER']
    
    # Group by side and ibl
    stations_by_side = defaultdict(lambda: defaultdict(list))
    for e in mrchue_events:
        side = e['side']
        ibl = e['ibl']
        stations_by_side[side][ibl].append(e)
    
    # Build output structure
    result = {
        "metadata": metadata,
        "source": "xfoil_instrumented",
        "description": "MRCHUE Newton iteration data for BL march debugging",
        "sides": {}
    }
    
    for side in [1, 2]:
        side_data = {"stations": []}
        
        for ibl in sorted(stations_by_side[side].keys()):
            iterations = stations_by_side[side][ibl]
            
            # Sort by newton_iter
            iterations = sorted(iterations, key=lambda e: e['newton_iter'])
            
            # Get first and final state
            first = iterations[0]
            final = iterations[-1]
            
            station = {
                "side": side,
                "ibl": ibl,
                "x": first['x'],
                "Ue": first['Ue'],
                "n_iterations": len(iterations),
                "converged": final['converged'],
                "initial": {
                    "theta": first['input']['theta'],
                    "delta_star": first['input']['delta_star'],
                    "ctau": first['input']['ctau'],
                    "ampl": first['input']['ampl']
                },
                "final": {
                    "theta": final['output']['theta'],
                    "delta_star": final['output']['delta_star'],
                    "ctau": final['output']['ctau'],
                    "ampl": final['output']['ampl'],
                    "dmax": final['dmax']
                },
                "iterations": []
            }
            
            # Add iteration details
            for it in iterations:
                station["iterations"].append({
                    "iter": it['newton_iter'],
                    "theta_in": it['input']['theta'],
                    "delta_star_in": it['input']['delta_star'],
                    "theta_out": it['output']['theta'],
                    "delta_star_out": it['output']['delta_star'],
                    "vsrez": it['vsrez'],
                    "VS2": it['VS2'],
                    "dmax": it['dmax'],
                    "relaxation": it['relaxation'],
                    "converged": it['converged']
                })
            
            side_data["stations"].append(station)
        
        result["sides"][side] = side_data
    
    return result


def print_summary(data):
    """Print a summary of the extracted data."""
    print("=" * 70)
    print("MRCHUE Newton Iteration Data Summary")
    print("=" * 70)
    print()
    print(f"Metadata:")
    for k, v in data['metadata'].items():
        print(f"  {k}: {v}")
    print()
    
    for side in [1, 2]:
        side_name = "Upper" if side == 1 else "Lower"
        stations = data['sides'][side]['stations']
        print(f"{side_name} Surface (Side {side}): {len(stations)} stations")
        
        # Summary stats
        n_iters = [s['n_iterations'] for s in stations]
        converged = sum(1 for s in stations if s['converged'])
        
        print(f"  Newton iterations: min={min(n_iters)}, max={max(n_iters)}, avg={sum(n_iters)/len(n_iters):.1f}")
        print(f"  Converged: {converged}/{len(stations)}")
        
        # Show theta evolution
        print(f"\n  Station theta evolution:")
        print(f"  {'IBL':>4} {'X':>10} {'θ_init':>12} {'θ_final':>12} {'Iters':>6} {'Conv':>5}")
        for s in stations[:10]:  # First 10 stations
            print(f"  {s['ibl']:4d} {s['x']:10.4f} {s['initial']['theta']:12.4e} {s['final']['theta']:12.4e} {s['n_iterations']:6d} {'Y' if s['converged'] else 'N':>5}")
        if len(stations) > 10:
            print(f"  ... ({len(stations) - 10} more stations)")
        print()


def main():
    parser = argparse.ArgumentParser(
        description="Extract MRCHUE Newton iteration data from XFOIL debug output"
    )
    parser.add_argument("input", help="Path to xfoil_debug.json")
    parser.add_argument("-o", "--output", 
                       help="Output JSON file (default: testdata/mrchue_iterations.json)")
    parser.add_argument("--summary", action="store_true",
                       help="Print summary only, don't write output file")
    
    args = parser.parse_args()
    
    # Load events
    print(f"Loading {args.input}...")
    events = load_xfoil_events(args.input)
    print(f"Loaded {len(events)} events")
    
    # Extract MRCHUE iterations
    data = extract_mrchue_iterations(events)
    
    # Print summary
    print_summary(data)
    
    # Write output
    if not args.summary:
        output_path = args.output or "testdata/mrchue_iterations.json"
        output_path = Path(output_path)
        
        # Ensure directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
