#!/usr/bin/env python3
"""
Extract XFOIL CL/CD from actual trace files.
"""
import json
from pathlib import Path

def extract_final_values(trace_file):
    """Extract final CL/CD/CM from XFOIL trace."""
    with open(trace_file) as f:
        data = json.load(f)
    
    events = data.get('events', [])
    
    # Find the last converged iteration
    # Look for CLCALC or CL_DETAIL events at the final iteration
    max_iter = 0
    cl_events = []
    
    for e in events:
        if e.get('subroutine') in ['CL_DETAIL', 'CLCALC']:
            iter_num = e.get('iteration', 0)
            if iter_num > max_iter:
                max_iter = iter_num
            cl_events.append(e)
    
    # Get the last CL_DETAIL event
    if cl_events:
        final_events = [e for e in cl_events if e.get('iteration') == max_iter]
        if not final_events:
            final_events = [cl_events[-1]]
        
        last = final_events[-1]
        return {
            'cl': last.get('cl'),
            'cd': last.get('cdp', last.get('cd')),  # Try cdp (pressure) first
            'cm': last.get('cm', 0),
            'alpha': last.get('alpha_deg', last.get('alpha', 0)),
            'iteration': max_iter
        }
    
    return None

def main():
    traces_dir = Path("/Users/harry/flexfoil-boundary-layer/traces/xfoil/naca0012/re3e06")
    
    results = []
    for trace_file in sorted(traces_dir.glob("alpha_*.json")):
        try:
            data = extract_final_values(trace_file)
            if data:
                results.append(data)
                print(f"{trace_file.name}: α={data['alpha']:>6.1f}°, "
                      f"CL={data['cl']:>8.4f}, CD={data['cd']:>9.6f}, iter={data['iteration']}")
        except Exception as e:
            print(f"Error processing {trace_file.name}: {e}")
    
    # Save results
    output_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/xfoil_actual_data.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'airfoil': 'naca0012',
            'reynolds': 3e6,
            'results': sorted(results, key=lambda x: x['alpha'])
        }, f, indent=2)
    
    print(f"\n✓ Extracted {len(results)} XFOIL results")
    print(f"✓ Saved to {output_file}")

if __name__ == '__main__':
    main()
