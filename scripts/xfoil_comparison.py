#!/usr/bin/env python3
"""
XFOIL vs RustFoil Polar Comparison

Generates comparison polars for multiple airfoils and Reynolds numbers,
creates plots, and saves data for hardcoding into the app.
"""

import subprocess
import json
import os
import sys
from pathlib import Path
import tempfile
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Configuration
XFOIL_PATH = "/Applications/ESP128/EngSketchPad/bin/xfoil"
PROJECT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_DIR / "testdata" / "xfoil_comparison"

# Test matrix
AIRFOILS = [
    {"name": "NACA 0012", "naca": "0012"},
    {"name": "NACA 2412", "naca": "2412"},
    {"name": "NACA 4412", "naca": "4412"},
]

REYNOLDS = [1e6, 3e6, 6e6]
NCRIT = 9.0
ALPHA_RANGE = (-4, 14, 1)  # start, end, step


def run_xfoil_polar(naca: str, reynolds: float, alpha_start: float, alpha_end: float, alpha_step: float, ncrit: float = 9.0) -> dict:
    """Run XFOIL and extract polar data."""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        polar_filename = "polar.txt"
        polar_file = os.path.join(tmpdir, polar_filename)
        
        # Build XFOIL command sequence
        # Use relative filename and run XFOIL from tmpdir
        commands = f"""NACA {naca}
PANE
OPER
VISC {reynolds:.0f}
VPAR
N {ncrit}

PACC
{polar_filename}

ASEQ {alpha_start} {alpha_end} {alpha_step}

QUIT
"""
        
        # Run XFOIL from the temp directory so relative paths work
        try:
            result = subprocess.run(
                [XFOIL_PATH],
                input=commands,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=tmpdir
            )
        except subprocess.TimeoutExpired:
            print(f"  XFOIL timeout for NACA {naca} at Re={reynolds:.0e}")
            return None
        except Exception as e:
            print(f"  XFOIL error for NACA {naca}: {e}")
            return None
        
        # Check if file exists
        if not os.path.exists(polar_file):
            return None
        
        return parse_xfoil_polar(polar_file)


def parse_xfoil_polar(filepath: str) -> dict:
    """Parse XFOIL polar output file.
    
    XFOIL polar format (version 6.99):
    - Header lines start with spaces or text
    - Data lines start with number (can be scientific notation like -0.200E+01)
    - Columns: alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr, [Top_Itr, Bot_Itr]
    """
    data = {
        "alpha": [],
        "cl": [],
        "cd": [],
        "cdp": [],
        "cm": [],
        "xtr_top": [],
        "xtr_bot": [],
    }
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header line to determine column positions
    in_data = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Skip lines that don't look like data
        # Data lines start with a number (possibly negative, possibly scientific notation)
        if not re.match(r'^\s*-?\d', line):
            continue
        
        # Parse data line
        parts = line.split()
        if len(parts) >= 7:
            try:
                alpha = float(parts[0])
                cl = float(parts[1])
                cd = float(parts[2])
                cdp = float(parts[3])
                cm = float(parts[4])
                xtr_top = float(parts[5])
                xtr_bot = float(parts[6])
                
                # Sanity checks
                if abs(alpha) > 30 or abs(cl) > 5 or cd < 0 or cd > 1:
                    continue
                
                data["alpha"].append(alpha)
                data["cl"].append(cl)
                data["cd"].append(cd)
                data["cdp"].append(cdp)
                data["cm"].append(cm)
                data["xtr_top"].append(xtr_top)
                data["xtr_bot"].append(xtr_bot)
            except (ValueError, IndexError):
                continue
    
    return data if data["alpha"] else None


def run_rustfoil_via_cargo(naca: str, reynolds: float, alpha_start: float, alpha_end: float, alpha_step: float) -> dict:
    """Run RustFoil polar generation via cargo run."""
    
    # Run via cargo
    try:
        result = subprocess.run(
            ["cargo", "run", "--release", "--example", "polar_gen", "--",
             naca, str(reynolds), str(alpha_start), str(alpha_end), str(alpha_step)],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(PROJECT_DIR / "crates" / "rustfoil-solver")
        )
        
        if result.returncode != 0:
            print(f"  RustFoil error: {result.stderr[:200]}")
            return None
        
        # Parse JSON output - find the JSON part (between { and })
        output = result.stdout.strip()
        json_start = output.find('{')
        json_end = output.rfind('}') + 1
        if json_start >= 0 and json_end > json_start:
            json_str = output[json_start:json_end]
            return json.loads(json_str)
        
    except subprocess.TimeoutExpired:
        print(f"  RustFoil timeout for NACA {naca} at Re={reynolds:.0e}")
    except json.JSONDecodeError as e:
        print(f"  JSON parse error: {e}")
    except Exception as e:
        print(f"  RustFoil error: {e}")
    
    return None


def plot_comparison(airfoil_name: str, reynolds: float, xfoil_data: dict, rustfoil_data: dict, output_path: Path):
    """Create comparison plots for a single airfoil/Re combination."""
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)
    
    fig.suptitle(f'{airfoil_name} at Re = {reynolds:.0e}', fontsize=14, fontweight='bold')
    
    # Colors
    xfoil_color = '#1f77b4'  # Blue
    rustfoil_color = '#ff7f0e'  # Orange
    
    # 1. Cl vs Alpha
    ax1 = fig.add_subplot(gs[0, 0])
    if xfoil_data:
        ax1.plot(xfoil_data['alpha'], xfoil_data['cl'], 'o-', color=xfoil_color, 
                 label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_data:
        ax1.plot(rustfoil_data['alpha'], rustfoil_data['cl'], 's--', color=rustfoil_color,
                 label='RustFoil', markersize=4, linewidth=1.5)
    ax1.set_xlabel('α (deg)')
    ax1.set_ylabel('Cl')
    ax1.set_title('Lift Coefficient')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=0, color='k', linewidth=0.5)
    
    # 2. Cd vs Alpha  
    ax2 = fig.add_subplot(gs[0, 1])
    if xfoil_data:
        ax2.plot(xfoil_data['alpha'], [cd * 10000 for cd in xfoil_data['cd']], 'o-', 
                 color=xfoil_color, label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_data:
        ax2.plot(rustfoil_data['alpha'], [cd * 10000 for cd in rustfoil_data['cd']], 's--',
                 color=rustfoil_color, label='RustFoil', markersize=4, linewidth=1.5)
    ax2.set_xlabel('α (deg)')
    ax2.set_ylabel('Cd (counts)')
    ax2.set_title('Drag Coefficient')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Cl vs Cd (Drag Polar)
    ax3 = fig.add_subplot(gs[1, 0])
    if xfoil_data:
        ax3.plot([cd * 10000 for cd in xfoil_data['cd']], xfoil_data['cl'], 'o-',
                 color=xfoil_color, label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_data:
        ax3.plot([cd * 10000 for cd in rustfoil_data['cd']], rustfoil_data['cl'], 's--',
                 color=rustfoil_color, label='RustFoil', markersize=4, linewidth=1.5)
    ax3.set_xlabel('Cd (counts)')
    ax3.set_ylabel('Cl')
    ax3.set_title('Drag Polar')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. L/D vs Cl
    ax4 = fig.add_subplot(gs[1, 1])
    if xfoil_data:
        ld_xfoil = [cl/cd if cd > 1e-6 else 0 for cl, cd in zip(xfoil_data['cl'], xfoil_data['cd'])]
        ax4.plot(xfoil_data['cl'], ld_xfoil, 'o-', color=xfoil_color, 
                 label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_data:
        ld_rustfoil = [cl/cd if cd > 1e-6 else 0 for cl, cd in zip(rustfoil_data['cl'], rustfoil_data['cd'])]
        ax4.plot(rustfoil_data['cl'], ld_rustfoil, 's--', color=rustfoil_color,
                 label='RustFoil', markersize=4, linewidth=1.5)
    ax4.set_xlabel('Cl')
    ax4.set_ylabel('L/D')
    ax4.set_title('Lift-to-Drag Ratio')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_path.name}")


def compute_errors(xfoil_data: dict, rustfoil_data: dict) -> dict:
    """Compute error metrics between XFOIL and RustFoil."""
    if not xfoil_data or not rustfoil_data:
        return None
    
    if not xfoil_data.get('alpha') or not rustfoil_data.get('alpha'):
        return None
    
    errors = {
        "cl_rms": 0.0,
        "cd_rms": 0.0,
        "cl_max_err": 0.0,
        "cd_max_err": 0.0,
        "n_points": 0,
    }
    
    cl_errs = []
    cd_errs = []
    
    n_xfoil = len(xfoil_data['alpha'])
    n_rustfoil = len(rustfoil_data['alpha'])
    
    for i in range(n_xfoil):
        if i >= len(xfoil_data['cd']):
            break
        alpha = xfoil_data['alpha'][i]
        
        # Find matching alpha in rustfoil
        for j in range(n_rustfoil):
            if j >= len(rustfoil_data['cd']):
                break
            rf_alpha = rustfoil_data['alpha'][j]
            
            if abs(alpha - rf_alpha) < 0.1:
                cl_err = rustfoil_data['cl'][j] - xfoil_data['cl'][i]
                cd_err = rustfoil_data['cd'][j] - xfoil_data['cd'][i]
                cl_errs.append(cl_err)
                cd_errs.append(cd_err)
                break
    
    if cl_errs:
        errors["cl_rms"] = np.sqrt(np.mean(np.array(cl_errs)**2))
        errors["cd_rms"] = np.sqrt(np.mean(np.array(cd_errs)**2))
        errors["cl_max_err"] = max(abs(e) for e in cl_errs)
        errors["cd_max_err"] = max(abs(e) for e in cd_errs)
        errors["n_points"] = len(cl_errs)
    
    return errors


def create_summary_plot(all_results: list, output_path: Path):
    """Create a summary comparison plot for all airfoils."""
    
    n_airfoils = len(AIRFOILS)
    fig, axes = plt.subplots(n_airfoils, 2, figsize=(12, 4*n_airfoils))
    if n_airfoils == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle('RustFoil vs XFOIL Comparison Summary', fontsize=14, fontweight='bold')
    
    xfoil_color = '#1f77b4'
    rustfoil_color = '#ff7f0e'
    
    for i, airfoil in enumerate(AIRFOILS):
        ax_cl = axes[i, 0]
        ax_polar = axes[i, 1]
        
        for result in all_results:
            if result['airfoil'] != airfoil['name']:
                continue
            
            re_str = f"Re={result['reynolds']:.0e}"
            linestyle = '-' if result['reynolds'] == 3e6 else '--' if result['reynolds'] == 1e6 else ':'
            
            xf = result.get('xfoil')
            rf = result.get('rustfoil')
            
            if xf:
                ax_cl.plot(xf['alpha'], xf['cl'], linestyle, color=xfoil_color, 
                          alpha=0.7, label=f'XFOIL {re_str}' if result['reynolds'] == 3e6 else None)
                ax_polar.plot([cd*10000 for cd in xf['cd']], xf['cl'], linestyle, 
                             color=xfoil_color, alpha=0.7)
            
            if rf:
                ax_cl.plot(rf['alpha'], rf['cl'], linestyle, color=rustfoil_color,
                          alpha=0.7, label=f'RustFoil {re_str}' if result['reynolds'] == 3e6 else None)
                ax_polar.plot([cd*10000 for cd in rf['cd']], rf['cl'], linestyle,
                             color=rustfoil_color, alpha=0.7)
        
        ax_cl.set_title(f'{airfoil["name"]}')
        ax_cl.set_xlabel('α (deg)')
        ax_cl.set_ylabel('Cl')
        ax_cl.grid(True, alpha=0.3)
        ax_cl.legend(loc='lower right')
        
        ax_polar.set_xlabel('Cd (counts)')
        ax_polar.set_ylabel('Cl')
        ax_polar.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved summary: {output_path.name}")


def main():
    """Main entry point."""
    
    print("=" * 70)
    print("XFOIL vs RustFoil Polar Comparison")
    print("=" * 70)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    all_results = []
    all_data = {}  # For JSON export
    
    alpha_start, alpha_end, alpha_step = ALPHA_RANGE
    
    for airfoil in AIRFOILS:
        print(f"\n{airfoil['name']}")
        print("-" * 40)
        all_data[airfoil['name']] = {}
        
        for reynolds in REYNOLDS:
            print(f"\n  Re = {reynolds:.0e}")
            
            # Run XFOIL
            print("    Running XFOIL...", end=" ", flush=True)
            xfoil_data = run_xfoil_polar(
                airfoil['naca'], reynolds, 
                alpha_start, alpha_end, alpha_step, NCRIT
            )
            if xfoil_data:
                print(f"✓ ({len(xfoil_data['alpha'])} points)")
            else:
                print("✗")
            
            # Run RustFoil
            print("    Running RustFoil...", end=" ", flush=True)
            rustfoil_data = run_rustfoil_via_cargo(
                airfoil['naca'], reynolds,
                alpha_start, alpha_end, alpha_step
            )
            if rustfoil_data:
                print(f"✓ ({len(rustfoil_data['alpha'])} points)")
            else:
                print("✗")
            
            # Store results
            result = {
                'airfoil': airfoil['name'],
                'naca': airfoil['naca'],
                'reynolds': reynolds,
                'xfoil': xfoil_data,
                'rustfoil': rustfoil_data,
            }
            all_results.append(result)
            
            # Compute errors
            if xfoil_data and rustfoil_data:
                errors = compute_errors(xfoil_data, rustfoil_data)
                if errors:
                    print(f"    Cl RMS error: {errors['cl_rms']:.4f}")
                    print(f"    Cd RMS error: {errors['cd_rms']:.6f} ({errors['cd_rms']*10000:.1f} counts)")
                result['errors'] = errors
            
            # Create comparison plot
            plot_path = OUTPUT_DIR / f"{airfoil['naca']}_Re{reynolds:.0e}.png"
            plot_comparison(airfoil['name'], reynolds, xfoil_data, rustfoil_data, plot_path)
            
            # Store for JSON
            re_key = f"Re{reynolds:.0e}"
            all_data[airfoil['name']][re_key] = {
                'xfoil': xfoil_data,
                'rustfoil': rustfoil_data,
                'errors': result.get('errors'),
            }
    
    # Create summary plot
    summary_path = OUTPUT_DIR / "summary_comparison.png"
    create_summary_plot(all_results, summary_path)
    
    # Save JSON data for app
    json_path = OUTPUT_DIR / "xfoil_comparison_data.json"
    with open(json_path, 'w') as f:
        json.dump(all_data, f, indent=2)
    print(f"\nSaved JSON data: {json_path.name}")
    
    # Print error summary
    print("\n" + "=" * 70)
    print("ERROR SUMMARY")
    print("=" * 70)
    print(f"{'Airfoil':<12} {'Re':<10} {'Cl RMS':<10} {'Cd RMS':<12} {'Cd counts':<10}")
    print("-" * 54)
    
    for result in all_results:
        if result.get('errors'):
            e = result['errors']
            print(f"{result['naca']:<12} {result['reynolds']:<10.0e} {e['cl_rms']:<10.4f} {e['cd_rms']:<12.6f} {e['cd_rms']*10000:<10.1f}")
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
