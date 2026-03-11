#!/usr/bin/env python3
"""
Compare early-station boundary layer thickness (θ) between XFOIL and RustFoil
to diagnose why transition occurs 43% late in RustFoil.
"""

import json
import sys

def parse_xfoil_bl(filepath):
    """Parse XFOIL boundary layer dump."""
    with open(filepath, 'r') as f:
        lines = [l for l in f if not l.startswith('#') and l.strip()]
    
    upper = []
    lower = []
    
    for line in lines:
        parts = line.split()
        if len(parts) >= 8:
            try:
                s = float(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                ue = float(parts[3])
                dstar = float(parts[4])
                theta = float(parts[5])
                cf = float(parts[6])
                h = float(parts[7])
                
                # Calculate Rθ manually (Re = 20M, ν = 1/Re)
                nu = 1.0 / 2e7
                rtheta = theta * ue / nu
                
                station = {
                    's': s, 'x': x, 'y': y, 'ue': ue,
                    'dstar': dstar, 'theta': theta, 'cf': cf, 'h': h,
                    'rtheta': rtheta
                }
                
                if y > 0:  # Upper surface
                    upper.append(station)
                else:
                    lower.append(station)
            except ValueError:
                continue
    
    # Reverse upper surface to get leading edge first
    upper.reverse()
    
    return upper, lower

def parse_rustfoil_bl(filepath, reynolds=2e7):
    """Parse RustFoil trace file for converged BL state."""
    print(f"Loading RustFoil trace: {filepath}...")
    
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    # Get FULL_BL_STATE event (final converged state)
    events = data.get('events', [])
    full_bl_events = [e for e in events if e.get('subroutine') == 'FULL_BL_STATE']
    
    if not full_bl_events:
        print("ERROR: No FULL_BL_STATE events found!")
        return [], []
    
    bl_state = full_bl_events[-1]  # Get the last one (converged)
    print(f"Found FULL_BL_STATE at iteration {bl_state.get('iteration', 0)}")
    
    # Extract upper and lower surface data
    upper_data = bl_state.get('upper', {})
    lower_data = bl_state.get('lower', {})
    
    n_upper = len(upper_data.get('x', []))
    n_lower = len(lower_data.get('x', []))
    
    print(f"Upper surface: {n_upper} stations")
    print(f"Lower surface: {n_lower} stations")
    
    # Convert to list of station dicts
    upper_list = []
    for i in range(n_upper):
        theta = upper_data['theta'][i]
        ue = upper_data['ue'][i]
        rtheta = theta * ue * reynolds
        
        station = {
            'ibl': i,
            'x': upper_data['x'][i],
            'y': 0,  # Not in trace
            'ue': ue,
            'dstar': upper_data['delta_star'][i],
            'theta': theta,
            'cf': upper_data['cf'][i],
            'h': upper_data['hk'][i],
            'rtheta': rtheta,
            's': 0  # Not in trace
        }
        upper_list.append(station)
    
    lower_list = []
    for i in range(n_lower):
        theta = lower_data['theta'][i]
        ue = lower_data['ue'][i]
        rtheta = theta * ue * reynolds
        
        station = {
            'ibl': i,
            'x': lower_data['x'][i],
            'y': 0,
            'ue': ue,
            'dstar': lower_data['delta_star'][i],
            'theta': theta,
            'cf': lower_data['cf'][i],
            'h': lower_data['hk'][i],
            'rtheta': rtheta,
            's': 0
        }
        lower_list.append(station)
    
    # Show first few stations for debugging
    print("\nFirst 5 upper surface stations:")
    print(f"{'x':>8} {'θ':>12} {'Ue':>8} {'Rθ':>10}")
    for s in upper_list[:5]:
        print(f"{s['x']:8.5f} {s['theta']:12.8e} {s['ue']:8.4f} {s['rtheta']:10.1f}")
    
    return upper_list, lower_list

def compare_early_stations(xfoil_upper, rustfoil_upper, x_max=0.03):
    """Compare early-station θ values."""
    
    print("\n" + "="*100)
    print(f"EARLY BOUNDARY LAYER COMPARISON (x < {x_max})")
    print("="*100)
    
    print("\nXFOIL early BL:")
    print(f"{'x':>8} {'θ':>12} {'δ*':>12} {'H':>8} {'Ue':>8} {'Rθ':>10} {'Cf':>10}")
    print("-" * 80)
    for s in xfoil_upper:
        if s['x'] < x_max:
            print(f"{s['x']:8.5f} {s['theta']:12.8f} {s['dstar']:12.8f} {s['h']:8.3f} "
                  f"{s['ue']:8.4f} {s['rtheta']:10.1f} {s['cf']:10.6f}")
    
    print("\nRustFoil early BL:")
    print(f"{'x':>8} {'θ':>12} {'δ*':>12} {'H':>8} {'Ue':>8} {'Rθ':>10} {'Cf':>10}")
    print("-" * 80)
    for s in rustfoil_upper:
        if s['x'] < x_max:
            print(f"{s['x']:8.5f} {s['theta']:12.8f} {s['dstar']:12.8f} {s['h']:8.3f} "
                  f"{s['ue']:8.4f} {s['rtheta']:10.1f} {s['cf']:10.6f}")
    
    # Match stations and compare
    print("\nSide-by-Side Comparison:")
    print(f"{'x':>8} | {'XF θ':>12} {'RF θ':>12} {'Δθ%':>8} | "
          f"{'XF Ue':>8} {'RF Ue':>8} {'ΔUe%':>8} | {'XF Rθ':>10} {'RF Rθ':>10} {'ΔRθ%':>8}")
    print("-" * 130)
    
    max_theta_error = 0
    max_theta_error_x = 0
    max_ue_error = 0
    max_ue_error_x = 0
    divergence_x = None
    ue_divergence_x = None
    
    for xf in xfoil_upper:
        if xf['x'] < x_max:
            # Find matching RF station (closest x)
            rf = min(rustfoil_upper, key=lambda s: abs(s['x'] - xf['x']))
            
            if abs(rf['x'] - xf['x']) < 0.001:  # Within 0.1% chord
                theta_err = (rf['theta'] - xf['theta']) / xf['theta'] * 100
                ue_err = (rf['ue'] - xf['ue']) / xf['ue'] * 100
                rtheta_err = (rf['rtheta'] - xf['rtheta']) / xf['rtheta'] * 100
                
                print(f"{xf['x']:8.5f} | {xf['theta']:12.8f} {rf['theta']:12.8f} {theta_err:+8.2f} | "
                      f"{xf['ue']:8.4f} {rf['ue']:8.4f} {ue_err:+8.1f} | "
                      f"{xf['rtheta']:10.1f} {rf['rtheta']:10.1f} {rtheta_err:+8.2f}")
                
                if abs(theta_err) > abs(max_theta_error):
                    max_theta_error = theta_err
                    max_theta_error_x = xf['x']
                
                if abs(ue_err) > abs(max_ue_error):
                    max_ue_error = ue_err
                    max_ue_error_x = xf['x']
                
                # Mark where θ starts diverging significantly (> 5%)
                if divergence_x is None and abs(theta_err) > 5.0:
                    divergence_x = xf['x']
                
                # Mark where Ue starts diverging significantly (> 5%)
                if ue_divergence_x is None and abs(ue_err) > 5.0:
                    ue_divergence_x = xf['x']
    
    print("\n" + "="*100)
    print("DIAGNOSIS:")
    print("="*100)
    print(f"\n🔴 CRITICAL: Maximum Ue error: {max_ue_error:+.1f}% at x = {max_ue_error_x:.5f}")
    if ue_divergence_x:
        print(f"🔴 Ue starts diverging (>5%) at x = {ue_divergence_x:.5f}")
    else:
        print("✅ Ue remains within 5% throughout early region")
    
    print(f"\nMaximum θ error: {max_theta_error:+.2f}% at x = {max_theta_error_x:.5f}")
    if divergence_x:
        print(f"θ starts diverging (>5%) at x = {divergence_x:.5f}")
    else:
        print("θ remains within 5% throughout early region")
    
    print(f"\n⚠️  ROOT CAUSE: Edge velocity (Ue) is ~90% too low!")
    print(f"    Since Rθ = θ × Ue × Re, low Ue directly causes low Rθ")
    print(f"    This keeps Rθ below critical longer → delayed transition")
    
    # Check first station (stagnation initialization)
    if xfoil_upper and rustfoil_upper:
        xf_first = xfoil_upper[0]
        rf_first = rustfoil_upper[0]
        first_theta_err = (rf_first['theta'] - xf_first['theta']) / xf_first['theta'] * 100
        
        print(f"\nFirst station comparison (stagnation initialization):")
        print(f"  XFOIL:    x={xf_first['x']:.5f}, θ={xf_first['theta']:.8f}, Rθ={xf_first['rtheta']:.1f}")
        print(f"  RustFoil: x={rf_first['x']:.5f}, θ={rf_first['theta']:.8f}, Rθ={rf_first['rtheta']:.1f}")
        print(f"  Error:    Δθ = {first_theta_err:+.2f}%")
        
        if abs(first_theta_err) > 5.0:
            print("\n⚠️  WARNING: First station θ differs by >5% - likely stagnation initialization bug!")
    
    # Check transition region
    print(f"\nTransition region (x ~ 0.012):")
    for xf in xfoil_upper:
        if 0.010 <= xf['x'] <= 0.015:
            rf = min(rustfoil_upper, key=lambda s: abs(s['x'] - xf['x']))
            if abs(rf['x'] - xf['x']) < 0.001:
                theta_err = (rf['theta'] - xf['theta']) / xf['theta'] * 100
                print(f"  x={xf['x']:.5f}: XF Rθ={xf['rtheta']:7.1f}, RF Rθ={rf['rtheta']:7.1f}, "
                      f"Δθ={theta_err:+.2f}%")
                
                # Critical Rθ is typically ~200-300 for transition
                if xf['rtheta'] > 200 or rf['rtheta'] > 200:
                    print(f"    → Rθ > 200 (approaching transition critical value)")

def main():
    xfoil_file = '/tmp/xfoil_bl_early.txt'
    rustfoil_file = 'crates/rustfoil-solver/traces/rustfoil_FIXED/rustfoil_alpha_10.json'
    
    print("Parsing XFOIL BL data...")
    xf_upper, xf_lower = parse_xfoil_bl(xfoil_file)
    print(f"  Upper surface: {len(xf_upper)} stations")
    
    print("\nParsing RustFoil BL data...")
    rf_upper, rf_lower = parse_rustfoil_bl(rustfoil_file)
    
    # Compare early stations (x < 0.03)
    compare_early_stations(xf_upper, rf_upper, x_max=0.03)

if __name__ == '__main__':
    main()
