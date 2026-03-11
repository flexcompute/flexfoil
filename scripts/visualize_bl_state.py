#!/usr/bin/env python3
"""Visualize boundary layer state comparison: RustFoil vs XFOIL.

Shows δ*, θ, H, Cf along both upper and lower surfaces to diagnose
if surfaces are swapped or BL development is wrong.
"""

import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
import os
import re as regex

XFOIL_PATH = '/Users/harry/flexfoil-boundary-layer/Xfoil/bin/xfoil'
RUSTFOIL_DIR = '/Users/harry/flexfoil-boundary-layer'

def run_xfoil_bl(naca="0012", alpha=14.0, re=1e6):
    """Run XFOIL and get BL data via DUMP command."""
    # Create temp directory for dump files
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        dump_file = os.path.join(tmpdir, 'bl.dat')
        
        script = f"""NACA {naca}
PANE
OPER
VISC {re:.0f}
ITER 200
ALFA {alpha}
DUMP {dump_file}

QUIT
"""
        result = subprocess.run(
            [XFOIL_PATH],
            input=script,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        # Get CL/CD from stdout
        cl, cd = None, None
        for line in result.stdout.split('\n'):
            if 'CL =' in line:
                match = regex.search(r'CL\s*=\s*([-\d.]+)', line)
                if match:
                    cl = float(match.group(1))
            if 'CD =' in line:
                match = regex.search(r'CD\s*=\s*([-\d.]+)', line)
                if match:
                    cd = float(match.group(1))
        
        # Parse dump file
        if os.path.exists(dump_file):
            bl_data = parse_xfoil_dump(dump_file)
            return {'cl': cl, 'cd': cd, 'bl': bl_data}
        
        return {'cl': cl, 'cd': cd, 'bl': None}

def parse_xfoil_dump(filename):
    """Parse XFOIL BL dump file."""
    data = {'upper': {'s': [], 'x': [], 'ue': [], 'dstar': [], 'theta': [], 'cf': [], 'h': []},
            'lower': {'s': [], 'x': [], 'ue': [], 'dstar': [], 'theta': [], 'cf': [], 'h': []}}
    
    with open(filename, 'r') as f:
        content = f.read()
    
    lines = content.strip().split('\n')
    
    # XFOIL DUMP format:
    #     s        x        y       Ue/Vinf    Dstar     Theta      Cf       H
    # Find the header and parse data
    in_upper = False
    in_lower = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
            
        # Detect surface based on context
        if 'Side 1' in line or 'UPPER' in line.upper() or 'Top' in line:
            in_upper = True
            in_lower = False
            continue
        elif 'Side 2' in line or 'LOWER' in line.upper() or 'Bottom' in line:
            in_upper = False
            in_lower = True
            continue
        elif 's' in line.lower() and 'x' in line.lower() and 'Ue' in line:
            # Header line, skip
            continue
        
        # Try to parse data line
        parts = line.split()
        if len(parts) >= 8:
            try:
                s = float(parts[0])
                x = float(parts[1])
                # y = float(parts[2])  # Not used
                ue = float(parts[3])
                dstar = float(parts[4])
                theta = float(parts[5])
                cf = float(parts[6])
                h = float(parts[7])
                
                # Determine surface by Ue sign or position
                # Upper surface typically starts from LE with positive x increasing
                # Use first surface encountered as upper
                if not in_lower and (in_upper or not data['upper']['x']):
                    data['upper']['s'].append(s)
                    data['upper']['x'].append(x)
                    data['upper']['ue'].append(ue)
                    data['upper']['dstar'].append(dstar)
                    data['upper']['theta'].append(theta)
                    data['upper']['cf'].append(cf)
                    data['upper']['h'].append(h)
                    in_upper = True
                elif in_lower:
                    data['lower']['s'].append(s)
                    data['lower']['x'].append(x)
                    data['lower']['ue'].append(ue)
                    data['lower']['dstar'].append(dstar)
                    data['lower']['theta'].append(theta)
                    data['lower']['cf'].append(cf)
                    data['lower']['h'].append(h)
            except (ValueError, IndexError):
                continue
    
    return data

def run_rustfoil_bl(naca="0012", alpha=14.0, re=1e6):
    """Run RustFoil with BL dump."""
    # First generate airfoil file
    script = f"""NACA {naca}
PSAV /tmp/naca{naca}.dat

QUIT
"""
    subprocess.run([XFOIL_PATH], input=script, capture_output=True, text=True, timeout=10)
    
    airfoil_file = f"/tmp/naca{naca}.dat"
    
    result = subprocess.run(
        ['cargo', 'run', '--release', '-q', '-p', 'rustfoil-cli', '--',
         'viscous', airfoil_file, f'--alpha={alpha}', f'--re={re:.0f}', '--format', 'json', '--dump-bl'],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=RUSTFOIL_DIR
    )
    
    if result.returncode == 0:
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError:
            pass
    
    return None

def main():
    naca = "0012"
    alpha = 14.0  # High angle where stall should occur
    re = 1e6
    
    print(f"BL State Comparison: NACA {naca} at α={alpha}°, Re={re:.0e}")
    print("=" * 80)
    
    # Run both solvers
    print("Running XFOIL...")
    xfoil_data = run_xfoil_bl(naca, alpha, re)
    print(f"  XFOIL: CL={xfoil_data.get('cl', 'N/A')}, CD={xfoil_data.get('cd', 'N/A')}")
    
    print("Running RustFoil...")
    rf_data = run_rustfoil_bl(naca, alpha, re)
    if rf_data:
        print(f"  RustFoil: CL={rf_data.get('cl', 'N/A'):.4f}, CD={rf_data.get('cd', 'N/A'):.5f}")
    
    # Check if we have BL data
    xf_bl = xfoil_data.get('bl')
    rf_bl = rf_data.get('bl_state') if rf_data else None
    
    if not xf_bl or not xf_bl['upper']['x']:
        print("\nXFOIL BL data not available from DUMP file.")
        print("Will show RustFoil only.")
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    
    # Line styles: format string or dict of kwargs
    # XFOIL: solid blue/cyan, RustFoil: solid red/orange
    
    # Extract RustFoil data
    if rf_bl:
        rf_upper = rf_bl.get('upper', {})
        rf_lower = rf_bl.get('lower', {})
        
        rf_x_upper = rf_upper.get('x', [])
        rf_x_lower = rf_lower.get('x', [])
        rf_ue_upper = rf_upper.get('ue', [])
        rf_ue_lower = rf_lower.get('ue', [])
        rf_dstar_upper = rf_upper.get('delta_star', [])
        rf_dstar_lower = rf_lower.get('delta_star', [])
        rf_theta_upper = rf_upper.get('theta', [])
        rf_theta_lower = rf_lower.get('theta', [])
        rf_h_upper = rf_upper.get('h', [])
        rf_h_lower = rf_lower.get('h', [])
        rf_cf_upper = rf_upper.get('cf', [])
        rf_cf_lower = rf_lower.get('cf', [])
    
    # Plot 1: Ue distribution
    ax = axes[0, 0]
    if xf_bl and xf_bl['upper']['x']:
        ax.plot(xf_bl['upper']['x'], xf_bl['upper']['ue'], 'b-', label='XFOIL Upper', linewidth=2)
        ax.plot(xf_bl['lower']['x'], xf_bl['lower']['ue'], 'c--', label='XFOIL Lower', linewidth=2)
    if rf_bl:
        ax.plot(rf_x_upper, rf_ue_upper, 'r-', label='RustFoil Upper', linewidth=2, alpha=0.7)
        ax.plot(rf_x_lower, rf_ue_lower, color='orange', linestyle='--', label='RustFoil Lower', linewidth=2, alpha=0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('Ue/V∞')
    ax.set_title('Edge Velocity')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    
    # Plot 2: δ* distribution
    ax = axes[0, 1]
    if xf_bl and xf_bl['upper']['x']:
        ax.plot(xf_bl['upper']['x'], [d*1000 for d in xf_bl['upper']['dstar']], 'b-', label='XFOIL Upper', linewidth=2)
        ax.plot(xf_bl['lower']['x'], [d*1000 for d in xf_bl['lower']['dstar']], 'c--', label='XFOIL Lower', linewidth=2)
    if rf_bl:
        ax.plot(rf_x_upper, [d*1000 for d in rf_dstar_upper], 'r-', label='RustFoil Upper', linewidth=2, alpha=0.7)
        ax.plot(rf_x_lower, [d*1000 for d in rf_dstar_lower], color='orange', linestyle='--', label='RustFoil Lower', linewidth=2, alpha=0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('δ* × 1000')
    ax.set_title('Displacement Thickness')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: θ distribution
    ax = axes[0, 2]
    if xf_bl and xf_bl['upper']['x']:
        ax.plot(xf_bl['upper']['x'], [t*1000 for t in xf_bl['upper']['theta']], 'b-', label='XFOIL Upper', linewidth=2)
        ax.plot(xf_bl['lower']['x'], [t*1000 for t in xf_bl['lower']['theta']], 'c--', label='XFOIL Lower', linewidth=2)
    if rf_bl:
        ax.plot(rf_x_upper, [t*1000 for t in rf_theta_upper], 'r-', label='RustFoil Upper', linewidth=2, alpha=0.7)
        ax.plot(rf_x_lower, [t*1000 for t in rf_theta_lower], color='orange', linestyle='--', label='RustFoil Lower', linewidth=2, alpha=0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('θ × 1000')
    ax.set_title('Momentum Thickness')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: H distribution (CRITICAL)
    ax = axes[1, 0]
    if xf_bl and xf_bl['upper']['x']:
        ax.plot(xf_bl['upper']['x'], xf_bl['upper']['h'], 'b-', label='XFOIL Upper', linewidth=2)
        ax.plot(xf_bl['lower']['x'], xf_bl['lower']['h'], 'c--', label='XFOIL Lower', linewidth=2)
    if rf_bl:
        ax.plot(rf_x_upper, rf_h_upper, 'r-', label='RustFoil Upper', linewidth=2, alpha=0.7)
        ax.plot(rf_x_lower, rf_h_lower, color='orange', linestyle='--', label='RustFoil Lower', linewidth=2, alpha=0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('H = δ*/θ')
    ax.set_title('Shape Factor (H > 4 = SEPARATION)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=4.0, color='k', linestyle=':', linewidth=2, label='Sep threshold')
    ax.set_ylim(0, 10)
    
    # Plot 5: Cf distribution
    ax = axes[1, 1]
    if xf_bl and xf_bl['upper']['x']:
        ax.plot(xf_bl['upper']['x'], [c*1000 for c in xf_bl['upper']['cf']], 'b-', label='XFOIL Upper', linewidth=2)
        ax.plot(xf_bl['lower']['x'], [c*1000 for c in xf_bl['lower']['cf']], 'c--', label='XFOIL Lower', linewidth=2)
    if rf_bl:
        ax.plot(rf_x_upper, [c*1000 for c in rf_cf_upper], 'r-', label='RustFoil Upper', linewidth=2, alpha=0.7)
        ax.plot(rf_x_lower, [c*1000 for c in rf_cf_lower], color='orange', linestyle='--', label='RustFoil Lower', linewidth=2, alpha=0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('Cf × 1000')
    ax.set_title('Skin Friction (Cf < 0 = separated)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linestyle=':', linewidth=2)
    
    # Plot 6: Summary text
    ax = axes[1, 2]
    ax.axis('off')
    
    summary = f"NACA {naca} at α={alpha}°, Re={re/1e6:.0f}M\n\n"
    
    if xfoil_data.get('cl') is not None:
        summary += f"XFOIL:\n"
        summary += f"  CL = {xfoil_data['cl']:.4f}\n"
        summary += f"  CD = {xfoil_data['cd']:.5f}\n"
        if xf_bl and xf_bl['upper']['h']:
            h_max_xf = max(xf_bl['upper']['h'])
            x_h_max_xf = xf_bl['upper']['x'][xf_bl['upper']['h'].index(h_max_xf)]
            summary += f"  Max H (upper) = {h_max_xf:.2f} at x/c={x_h_max_xf:.3f}\n"
        summary += "\n"
    
    if rf_data:
        summary += f"RustFoil:\n"
        summary += f"  CL = {rf_data['cl']:.4f}\n"
        summary += f"  CD = {rf_data['cd']:.5f}\n"
        if rf_h_upper:
            h_max_rf = max(rf_h_upper)
            x_h_max_rf = rf_x_upper[rf_h_upper.index(h_max_rf)]
            summary += f"  Max H (upper) = {h_max_rf:.2f} at x/c={x_h_max_rf:.3f}\n"
        summary += "\n"
    
    if rf_data and xfoil_data.get('cl'):
        cl_err = (rf_data['cl'] - xfoil_data['cl']) / abs(xfoil_data['cl']) * 100
        summary += f"CL Error: {cl_err:+.1f}%\n"
        
        if rf_h_upper and xf_bl and xf_bl['upper']['h']:
            summary += f"\nKEY: XFOIL shows H > 4 (separation)\n"
            summary += f"     RustFoil max H = {max(rf_h_upper):.2f}\n"
            if max(rf_h_upper) < 4:
                summary += f"     --> RustFoil NOT detecting separation!\n"
    
    ax.text(0.1, 0.9, summary, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle(f'Boundary Layer State Comparison: NACA {naca}, α={alpha}°, Re={re/1e6:.0f}M', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_file = f'{RUSTFOIL_DIR}/scripts/bl_comparison_alpha{int(alpha)}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {output_file}")

if __name__ == "__main__":
    main()
