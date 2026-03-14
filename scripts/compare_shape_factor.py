#!/usr/bin/env python3
"""
Phase 1: Shape Factor (H) Distribution Analysis

Compare H distribution along the chord at different angles of attack
to understand separation behavior.

Key thresholds:
- Laminar: H > 2.5 indicates separation imminent
- Turbulent: H > 2.0-2.5 indicates separation
- H > 4.0 typically means fully separated
"""

import json
import subprocess
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def run_rustfoil_with_trace(alpha_deg: float, airfoil: str = "0012") -> str:
    """Run RustFoil and return path to debug JSON."""
    env = os.environ.copy()
    env['RUSTFOIL_TEST_ALPHA'] = str(alpha_deg)
    env['RUSTFOIL_TEST_AIRFOIL'] = airfoil
    env['RUSTFOIL_FULL_TRACE'] = '1'  # Enable full debug trace
    
    # Output file is written relative to the test crate's working directory
    debug_file = f"crates/rustfoil-solver/traces/rustfoil_new/rustfoil_alpha_{int(alpha_deg)}.json"
    
    # Create output directory
    os.makedirs("crates/rustfoil-solver/traces/rustfoil_new", exist_ok=True)
    
    result = subprocess.run(
        ['cargo', 'test', '--package', 'rustfoil-solver',
         '--test', 'xfoil_viscous_comparison',
         'test_newton_iteration_trace', '--', '--nocapture'],
        capture_output=True,
        text=True,
        timeout=120,
        env=env
    )
    
    if result.returncode != 0:
        print(f"Warning: Test returned non-zero exit code")
        print(result.stderr[:500] if result.stderr else "")
    
    return debug_file

def extract_final_h_distribution(debug_file: str):
    """Extract H values from the converged state using MRCHUE events.
    
    MRCHUE events contain the final BL state for each station after marching.
    We use iteration 1 (first full march) to get complete data.
    """
    
    if not os.path.exists(debug_file):
        print(f"Debug file not found: {debug_file}")
        return None, None
        
    with open(debug_file, 'r') as f:
        data = json.load(f)
    
    # Use MRCHUE events which have complete station data
    # Group by iteration and side, take the highest ibl per station to get final state
    iterations = defaultdict(lambda: {'upper': {}, 'lower': {}})
    
    for event in data.get('events', []):
        if event.get('subroutine') == 'MRCHUE':
            iter_num = event.get('iteration', 1)  # MRCHUE doesn't always have iteration
            side = event.get('side', 0)
            ibl = event.get('ibl', 0)
            
            x = event.get('x', 0)
            Hk = event.get('Hk', 0)
            Cf = event.get('Cf', 0)
            theta = event.get('theta', 0)
            delta_star = event.get('delta_star', 0)
            
            # Compute H from delta_star/theta
            H = delta_star / theta if theta > 1e-12 else Hk
            
            key = 'upper' if side == 1 else 'lower'
            # Only keep the station with highest ibl (final march state)
            if ibl not in iterations[iter_num][key] or iterations[iter_num][key][ibl]['x'] < x:
                iterations[iter_num][key][ibl] = {
                    'x': x, 'H': H, 'Hk': Hk, 'Cf': Cf, 'theta': theta, 'dstar': delta_star
                }
    
    # Fall back to BLVAR events if no MRCHUE
    if not any(iterations.values()):
        print("  No MRCHUE events, falling back to BLVAR...")
        for event in data.get('events', []):
            if event.get('subroutine') == 'BLVAR':
                iter_num = event.get('iteration', 0)
                side = event.get('side', 0)
                ibl = event.get('ibl', 0)
                
                x = event.get('input', {}).get('x', 0)
                H = event.get('output', {}).get('H', 0)
                Hk = event.get('output', {}).get('Hk', 0)
                Cf = event.get('output', {}).get('Cf', 0)
                
                key = 'upper' if side == 1 else 'lower'
                iterations[iter_num][key][ibl] = {
                    'x': x, 'H': H, 'Hk': Hk, 'Cf': Cf
                }
    
    if not iterations:
        print(f"No BL events found in {debug_file}")
        return None, None
    
    # Use iteration 1 for complete data (later iterations only update some stations)
    # Or use the highest iteration that has significant data
    best_iter = 1
    max_stations = 0
    for it, data_it in iterations.items():
        total = len(data_it['upper']) + len(data_it['lower'])
        if total > max_stations:
            max_stations = total
            best_iter = it
    
    print(f"  Using iteration {best_iter} with {max_stations} stations")
    final_data = iterations[best_iter]
    
    # Convert to sorted arrays
    upper_data = []
    for ibl in sorted(final_data['upper'].keys()):
        d = final_data['upper'][ibl]
        upper_data.append((d['x'], d['H'], d['Hk'], d['Cf']))
    
    lower_data = []
    for ibl in sorted(final_data['lower'].keys()):
        d = final_data['lower'][ibl]
        lower_data.append((d['x'], d['H'], d['Hk'], d['Cf']))
    
    return upper_data, lower_data

def analyze_separation(upper_data, lower_data, alpha):
    """Analyze separation indicators."""
    print(f"\n{'='*60}")
    print(f"Separation Analysis at α = {alpha}°")
    print(f"{'='*60}")
    
    # Thresholds
    H_SEP_LAM = 2.5  # Laminar separation
    H_SEP_TURB = 2.2  # Turbulent separation onset
    H_SEPARATED = 4.0  # Fully separated
    
    for name, data in [("Upper", upper_data), ("Lower", lower_data)]:
        if not data:
            continue
            
        x_vals = [d[0] for d in data]
        h_vals = [d[1] for d in data]
        cf_vals = [d[3] for d in data]
        
        max_h = max(h_vals)
        max_h_x = x_vals[h_vals.index(max_h)]
        min_cf = min(cf_vals)
        min_cf_x = x_vals[cf_vals.index(min_cf)]
        
        # Find where H exceeds thresholds
        sep_onset_x = None
        for x, h in zip(x_vals, h_vals):
            if h > H_SEP_TURB and sep_onset_x is None:
                sep_onset_x = x
                break
        
        print(f"\n{name} Surface:")
        print(f"  Max H = {max_h:.3f} at x/c = {max_h_x:.4f}")
        print(f"  Min Cf = {min_cf:.6f} at x/c = {min_cf_x:.4f}")
        
        if max_h > H_SEPARATED:
            print(f"  ⚠️  SEPARATION DETECTED: H > {H_SEPARATED}")
        elif max_h > H_SEP_TURB:
            print(f"  ⚠️  Separation onset: H > {H_SEP_TURB} at x/c = {sep_onset_x:.4f}")
        else:
            print(f"  ✓ Attached flow (H < {H_SEP_TURB} everywhere)")
        
        if min_cf < 0:
            print(f"  ⚠️  NEGATIVE Cf detected (reversed flow)")

def main():
    # Alpha values to compare
    alphas = [0, 5, 10, 12, 15]
    airfoil = "0012"
    
    if len(sys.argv) > 1:
        alphas = [float(a) for a in sys.argv[1:]]
    
    print(f"Running shape factor analysis for NACA {airfoil}")
    print(f"Alpha values: {alphas}")
    
    # Collect data
    results = {}
    for alpha in alphas:
        print(f"\n--- Running α = {alpha}° ---")
        debug_file = run_rustfoil_with_trace(alpha, airfoil)
        upper, lower = extract_final_h_distribution(debug_file)
        if upper is not None:
            results[alpha] = {'upper': upper, 'lower': lower}
            analyze_separation(upper, lower, alpha)
        else:
            print(f"  Failed to extract data")
    
    if not results:
        print("No results to plot")
        return
    
    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # === Upper surface H distribution ===
    ax1 = axes[0, 0]
    ax1.set_title('Upper Surface - Shape Factor H', fontweight='bold')
    for alpha, data in sorted(results.items()):
        if data['upper']:
            x = [d[0] for d in data['upper']]
            h = [d[1] for d in data['upper']]
            ax1.plot(x, h, '-o', label=f'α={alpha}°', markersize=3)
    
    ax1.axhline(y=2.2, color='orange', linestyle='--', alpha=0.5, label='Turb sep onset (H=2.2)')
    ax1.axhline(y=4.0, color='red', linestyle='--', alpha=0.5, label='Separated (H=4.0)')
    ax1.set_xlabel('x/c')
    ax1.set_ylabel('Shape Factor H')
    ax1.set_xlim(0, 1.1)
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # === Lower surface H distribution ===
    ax2 = axes[0, 1]
    ax2.set_title('Lower Surface - Shape Factor H', fontweight='bold')
    for alpha, data in sorted(results.items()):
        if data['lower']:
            x = [d[0] for d in data['lower']]
            h = [d[1] for d in data['lower']]
            ax2.plot(x, h, '-o', label=f'α={alpha}°', markersize=3)
    
    ax2.axhline(y=2.2, color='orange', linestyle='--', alpha=0.5, label='Turb sep onset')
    ax2.axhline(y=4.0, color='red', linestyle='--', alpha=0.5, label='Separated')
    ax2.set_xlabel('x/c')
    ax2.set_ylabel('Shape Factor H')
    ax2.set_xlim(0, 1.1)
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    # === Upper surface Cf distribution ===
    ax3 = axes[1, 0]
    ax3.set_title('Upper Surface - Skin Friction Cf', fontweight='bold')
    for alpha, data in sorted(results.items()):
        if data['upper']:
            x = [d[0] for d in data['upper']]
            cf = [d[3] for d in data['upper']]
            ax3.plot(x, cf, '-o', label=f'α={alpha}°', markersize=3)
    
    ax3.axhline(y=0, color='red', linestyle='--', alpha=0.5, label='Cf=0 (separation)')
    ax3.set_xlabel('x/c')
    ax3.set_ylabel('Skin Friction Cf')
    ax3.set_xlim(0, 1.1)
    ax3.legend(loc='upper right', fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    # === Lower surface Cf distribution ===
    ax4 = axes[1, 1]
    ax4.set_title('Lower Surface - Skin Friction Cf', fontweight='bold')
    for alpha, data in sorted(results.items()):
        if data['lower']:
            x = [d[0] for d in data['lower']]
            cf = [d[3] for d in data['lower']]
            ax4.plot(x, cf, '-o', label=f'α={alpha}°', markersize=3)
    
    ax4.axhline(y=0, color='red', linestyle='--', alpha=0.5, label='Cf=0 (separation)')
    ax4.set_xlabel('x/c')
    ax4.set_ylabel('Skin Friction Cf')
    ax4.set_xlim(0, 1.1)
    ax4.legend(loc='upper right', fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle(f'NACA {airfoil} - Shape Factor & Skin Friction Distribution\nRe = 3×10⁶',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('scripts/shape_factor_analysis.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: scripts/shape_factor_analysis.png")
    
    # Summary table
    print("\n" + "="*70)
    print("SUMMARY: Maximum H values")
    print("="*70)
    print(f"{'Alpha':>6} | {'Upper H_max':>12} | {'Lower H_max':>12} | {'Status'}")
    print("-"*70)
    for alpha, data in sorted(results.items()):
        upper_max = max([d[1] for d in data['upper']]) if data['upper'] else 0
        lower_max = max([d[1] for d in data['lower']]) if data['lower'] else 0
        
        status = "✓ Attached"
        if upper_max > 4.0 or lower_max > 4.0:
            status = "❌ SEPARATED"
        elif upper_max > 2.2 or lower_max > 2.2:
            status = "⚠️  Sep onset"
        
        print(f"{alpha:>6}° | {upper_max:>12.3f} | {lower_max:>12.3f} | {status}")

if __name__ == '__main__':
    main()
