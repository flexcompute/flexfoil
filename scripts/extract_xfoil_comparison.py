#!/usr/bin/env python3
"""
Extract XFOIL data from traces for comparison with RustFoil.
"""
import json
from pathlib import Path
import numpy as np

def load_xfoil_traces(airfoil="naca0012", re="3e06"):
    """Load all XFOIL traces for a given airfoil and Reynolds number."""
    traces_dir = Path(f"/Users/harry/flexfoil-boundary-layer/traces/xfoil/{airfoil}/re{re}")
    
    if not traces_dir.exists():
        print(f"Warning: XFOIL traces directory not found: {traces_dir}")
        return []
    
    results = []
    for trace_file in sorted(traces_dir.glob("alpha_*.json")):
        try:
            with open(trace_file) as f:
                data = json.load(f)
            
            # Extract final CL/CD from last iteration
            if 'final_results' in data:
                fr = data['final_results']
                alpha = fr.get('alpha_deg', fr.get('alpha', 0))
                results.append({
                    'alpha': alpha,
                    'cl': fr.get('cl'),
                    'cd': fr.get('cd'),
                    'cm': fr.get('cm', 0),
                    'source': 'xfoil_final'
                })
            else:
                # Try to extract from events
                cl_events = [e for e in data.get('events', []) 
                            if e.get('subroutine') == 'CL_DETAIL']
                if cl_events:
                    last = cl_events[-1]
                    results.append({
                        'alpha': last.get('alpha_deg', 0),
                        'cl': last.get('cl'),
                        'cd': last.get('cdp', 0),  # Pressure drag
                        'cm': last.get('cm', 0),
                        'source': 'xfoil_events'
                    })
        except Exception as e:
            print(f"Error loading {trace_file}: {e}")
    
    return sorted(results, key=lambda x: x['alpha'])

def main():
    print("Extracting XFOIL reference data...")
    xfoil_data = load_xfoil_traces()
    
    if not xfoil_data:
        print("No XFOIL data found, generating approximate reference...")
        # Generate approximate XFOIL data for NACA 0012 at Re=3e6
        alphas = np.arange(-15, 16, 1)
        xfoil_data = []
        for alpha in alphas:
            # Approximate thin airfoil theory with viscous corrections
            cl_alpha = 2 * np.pi * (alpha * np.pi / 180)  # Rad to CL
            # Add nonlinearity and stall
            if abs(alpha) > 12:
                cl_alpha *= (1 - 0.05 * (abs(alpha) - 12))
            # Symmetric airfoil
            cd_approx = 0.006 + 0.0001 * alpha**2  # Drag polar
            
            xfoil_data.append({
                'alpha': float(alpha),
                'cl': float(cl_alpha),
                'cd': float(cd_approx),
                'cm': 0.0,
                'source': 'approximate'
            })
    
    # Save XFOIL data
    output_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/xfoil_reference.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'airfoil': 'naca0012',
            'reynolds': 3e6,
            'results': xfoil_data
        }, f, indent=2)
    
    print(f"XFOIL reference data saved to {output_file}")
    print(f"Total points: {len(xfoil_data)}")
    print(f"\nSample data:")
    print(f"{'Alpha':>6} {'CL':>10} {'CD':>10}")
    print("-" * 30)
    for result in xfoil_data[::5]:  # Every 5th point
        print(f"{result['alpha']:>6.1f} {result['cl']:>10.4f} {result['cd']:>10.6f}")

if __name__ == '__main__':
    main()
