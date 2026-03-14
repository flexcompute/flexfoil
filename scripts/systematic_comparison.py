#!/usr/bin/env python3
"""
Systematic comparison of RustFoil vs XFOIL intermediate quantities.
Compares: transition locations, edge velocities, BL quantities, forces.
"""

import json
import subprocess
import sys
from pathlib import Path

# Paths
XFOIL_BIN = Path(__file__).parent.parent / "Xfoil-instrumented/bin"
RUSTFOIL = Path(__file__).parent.parent / "target/release/rustfoil"
AIRFOIL = Path(__file__).parent.parent / "testdata/naca0012.dat"

def get_xfoil_data(alpha):
    """Extract key quantities from XFOIL debug JSON."""
    json_path = XFOIL_BIN / f"xfoil_debug_a{alpha}.json"
    if not json_path.exists():
        return None
    
    with open(json_path) as f:
        data = json.load(f)
    
    events = data.get('events', [])
    
    # Get final iteration result
    viscres = [e for e in events if e.get('subroutine') == 'VISCAL_RESULT']
    if not viscres:
        return None
    final = viscres[-1]
    
    # Get BLVAR events from final iteration to extract Ue, theta, etc.
    final_iter = final.get('iteration', 1)
    blvar_events = [e for e in events if e.get('subroutine') == 'BLVAR' 
                    and e.get('iteration') == final_iter]
    
    # Organize by surface and IBL
    upper_stations = {}
    lower_stations = {}
    
    for e in blvar_events:
        side = e.get('side', 0)
        ibl = e.get('ibl', 0)
        inp = e.get('input', {})
        out = e.get('output', {})
        
        station_data = {
            'x': inp.get('x', 0),
            'u': inp.get('u', 0),
            'theta': inp.get('theta', 0),
            'delta_star': inp.get('delta_star', 0),
            'H': out.get('H', 0),
            'Hk': out.get('Hk', 0),
            'Cf': out.get('Cf', 0),
        }
        
        if side == 1:
            upper_stations[ibl] = station_data
        elif side == 2:
            lower_stations[ibl] = station_data
    
    # Get transition events
    transitions = [e for e in events if e.get('subroutine') == 'TRANSITION']
    x_tr_upper = None
    x_tr_lower = None
    for t in transitions:
        if t.get('side') == 1:
            x_tr_upper = t.get('XT', None)
        elif t.get('side') == 2:
            x_tr_lower = t.get('XT', None)
    
    # Find stations at specific x locations (approximate)
    def find_station_at_x(stations, target_x, tolerance=0.05):
        for ibl, s in sorted(stations.items()):
            if abs(s['x'] - target_x) < tolerance:
                return s
        return None
    
    # Get TE stations (highest IBL)
    upper_te = upper_stations.get(max(upper_stations.keys())) if upper_stations else None
    lower_te = lower_stations.get(max(lower_stations.keys())) if lower_stations else None
    
    return {
        'alpha': alpha,
        'CL': final.get('CL', 0),
        'CD': final.get('CD', 0),
        'CD_f': final.get('CD_friction', 0),
        'CD_p': final.get('CD_viscous', 0),  # pressure drag component
        'x_tr_upper': x_tr_upper,
        'x_tr_lower': x_tr_lower,
        'upper_te': upper_te,
        'lower_te': lower_te,
        'upper_mid': find_station_at_x(upper_stations, 0.5),
        'lower_mid': find_station_at_x(lower_stations, 0.5),
        'iterations': final.get('iteration', 0),
        'residual': final.get('rms_residual', 0),
    }

def get_rustfoil_data(alpha):
    """Run RustFoil and extract key quantities."""
    cmd = [str(RUSTFOIL), "viscous", str(AIRFOIL), "--alpha", str(alpha), "--re", "1e6"]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        output = result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return None
    
    # Parse output
    data = {'alpha': alpha}
    
    for line in output.split('\n'):
        if 'CL:' in line:
            try:
                data['CL'] = float(line.split()[-1])
            except:
                pass
        elif 'CD:' in line and 'CD_f' not in line and 'CD_p' not in line:
            try:
                data['CD'] = float(line.split()[-1])
            except:
                pass
        elif 'x_tr (U):' in line:
            try:
                data['x_tr_upper'] = float(line.split()[-1])
            except:
                pass
        elif 'x_tr (L):' in line:
            try:
                data['x_tr_lower'] = float(line.split()[-1])
            except:
                pass
        elif 'Iters:' in line:
            try:
                data['iterations'] = int(line.split()[-1])
            except:
                pass
    
    return data

def compare_quantity(name, rf_val, xf_val, unit=""):
    """Compare a single quantity and return formatted string."""
    if rf_val is None or xf_val is None:
        return f"  {name:15s}: RF={rf_val}, XF={xf_val} (N/A)"
    
    if xf_val == 0:
        if rf_val == 0:
            err = 0
        else:
            err = float('inf')
    else:
        err = (rf_val - xf_val) / abs(xf_val) * 100
    
    return f"  {name:15s}: RF={rf_val:8.5f}, XF={xf_val:8.5f}, Err={err:+6.1f}% {unit}"

def main():
    print("=" * 70)
    print("SYSTEMATIC COMPARISON: RustFoil vs XFOIL")
    print("=" * 70)
    
    alphas = [0, 2, 4, 6, 8]
    
    # Collect all data first
    all_data = []
    
    for alpha in alphas:
        xf = get_xfoil_data(alpha)
        rf = get_rustfoil_data(alpha)
        
        if xf is None:
            print(f"\nAlpha {alpha}°: No XFOIL data")
            continue
        if rf is None:
            print(f"\nAlpha {alpha}°: No RustFoil data")
            continue
        
        all_data.append((alpha, rf, xf))
        
        print(f"\n{'='*70}")
        print(f"ALPHA = {alpha}°")
        print(f"{'='*70}")
        
        # Forces
        print("\nFORCES:")
        print(compare_quantity("CL", rf.get('CL'), xf.get('CL')))
        print(compare_quantity("CD", rf.get('CD'), xf.get('CD')))
        if xf.get('CD_f'):
            print(f"  {'CD_f (XFOIL)':15s}: {xf['CD_f']:.6f}")
        if xf.get('CD_p'):
            print(f"  {'CD_p (XFOIL)':15s}: {xf['CD_p']:.6f}")
        
        # Transition
        print("\nTRANSITION:")
        print(compare_quantity("x_tr (upper)", rf.get('x_tr_upper'), xf.get('x_tr_upper')))
        print(compare_quantity("x_tr (lower)", rf.get('x_tr_lower'), xf.get('x_tr_lower')))
        
        # Convergence
        print("\nCONVERGENCE:")
        print(f"  {'Iterations':15s}: RF={rf.get('iterations', 'N/A')}, XF={xf.get('iterations', 'N/A')}")
        if xf.get('residual'):
            print(f"  {'Final residual':15s}: XF={xf['residual']:.2e}")
        
        # TE quantities (if available)
        if xf.get('upper_te'):
            te = xf['upper_te']
            print(f"\nUPPER TE (XFOIL):")
            print(f"  x={te['x']:.4f}, Ue={te['u']:.4f}, θ={te['theta']:.6e}, H={te['H']:.3f}, Cf={te['Cf']:.6e}")
        
        if xf.get('lower_te'):
            te = xf['lower_te']
            print(f"\nLOWER TE (XFOIL):")
            print(f"  x={te['x']:.4f}, Ue={te['u']:.4f}, θ={te['theta']:.6e}, H={te['H']:.3f}, Cf={te['Cf']:.6e}")
    
    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Alpha':>5} | {'CL_RF':>8} | {'CL_XF':>8} | {'CL_Err%':>8} | {'CD_RF':>8} | {'CD_XF':>8} | {'CD_Err%':>8}")
    print("-" * 70)
    
    for alpha, rf, xf in all_data:
        cl_rf = rf.get('CL', 0)
        cl_xf = xf.get('CL', 0)
        cd_rf = rf.get('CD', 0)
        cd_xf = xf.get('CD', 0)
        
        cl_err = (cl_rf - cl_xf) / cl_xf * 100 if cl_xf != 0 else float('nan')
        cd_err = (cd_rf - cd_xf) / cd_xf * 100 if cd_xf != 0 else float('nan')
        
        print(f"{alpha:5}° | {cl_rf:8.4f} | {cl_xf:8.4f} | {cl_err:+7.1f}% | {cd_rf:8.5f} | {cd_xf:8.5f} | {cd_err:+7.1f}%")
    
    # Transition comparison
    print("\n" + "=" * 70)
    print("TRANSITION COMPARISON")
    print("=" * 70)
    print(f"{'Alpha':>5} | {'xtr_U_RF':>8} | {'xtr_U_XF':>8} | {'xtr_L_RF':>8} | {'xtr_L_XF':>8}")
    print("-" * 70)
    
    for alpha, rf, xf in all_data:
        xtr_u_rf = rf.get('x_tr_upper', 0)
        xtr_u_xf = xf.get('x_tr_upper', 0)
        xtr_l_rf = rf.get('x_tr_lower', 0)
        xtr_l_xf = xf.get('x_tr_lower', 0)
        
        print(f"{alpha:5}° | {xtr_u_rf:8.4f} | {xtr_u_xf if xtr_u_xf else 'N/A':>8} | {xtr_l_rf:8.4f} | {xtr_l_xf if xtr_l_xf else 'N/A':>8}")

if __name__ == "__main__":
    main()
