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
PROJECT_DIR = Path(__file__).parent.parent
DEFAULT_XFOIL = PROJECT_DIR / "Xfoil" / "bin" / "xfoil"
XFOIL_PATH = os.environ.get(
    "XFOIL_PATH",
    str(DEFAULT_XFOIL if DEFAULT_XFOIL.exists() else "/Applications/ESP128/EngSketchPad/bin/xfoil")
)
OUTPUT_DIR = PROJECT_DIR / "testdata" / "xfoil_comparison"

# Test matrix
AIRFOILS = [
    {"name": "NACA 0012", "naca": "0012"},
    {"name": "NACA 2412", "naca": "2412"},
    {"name": "NACA 4412", "naca": "4412"},
    {"name": "NACA 6409", "naca": "6409"},
]

REYNOLDS = [0.5e6, 1e6, 3e6]
NCRIT = 9.0
ALPHA_RANGE = (-4, 16, 1)  # start, end, step


def run_xfoil_polar(naca: str, reynolds: float, alpha_start: float, alpha_end: float, alpha_step: float, ncrit: float = 9.0) -> dict:
    """Run XFOIL and extract polar data."""
    
    if not Path(XFOIL_PATH).exists():
        print(f"  XFOIL not found at: {XFOIL_PATH}")
        print("  Set XFOIL_PATH env var if needed.")
        return None

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


def ensure_rustfoil_binary() -> Path:
    """Build or locate the RustFoil polar_gen binary."""
    solver_dir = PROJECT_DIR / "crates" / "rustfoil-solver"
    candidates = [
        PROJECT_DIR / "target" / "release" / "examples" / "polar_gen",
        solver_dir / "target" / "release" / "examples" / "polar_gen",
    ]
    for binary in candidates:
        if binary.exists():
            return binary

    print("  Building RustFoil polar_gen (release)...")
    result = subprocess.run(
        ["cargo", "build", "--release", "--example", "polar_gen"],
        capture_output=True,
        text=True,
        cwd=str(solver_dir),
    )
    if result.returncode != 0:
        print(f"  RustFoil build error: {result.stderr[:200]}")
        return None

    for binary in candidates:
        if binary.exists():
            return binary
    return None


def run_rustfoil_via_cargo(naca: str, reynolds: float, alpha_start: float, alpha_end: float, alpha_step: float) -> dict:
    """Run RustFoil polar generation via compiled binary."""
    binary = ensure_rustfoil_binary()
    if not binary:
        return None

    try:
        result = subprocess.run(
            [str(binary), naca, str(reynolds), str(alpha_start), str(alpha_end), str(alpha_step)],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(binary.parent),
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
    
    # Prepare filtered datasets (drop None/non-finite)
    xfoil_plot = filter_polar_data(xfoil_data)
    rustfoil_plot = filter_polar_data(rustfoil_data)

    # 1. Cl vs Alpha
    ax1 = fig.add_subplot(gs[0, 0])
    if xfoil_plot:
        ax1.plot(xfoil_plot['alpha'], xfoil_plot['cl'], 'o-', color=xfoil_color, 
                 label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_plot:
        ax1.plot(rustfoil_plot['alpha'], rustfoil_plot['cl'], 's--', color=rustfoil_color,
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
    if xfoil_plot:
        ax2.plot(xfoil_plot['alpha'], [cd * 10000 for cd in xfoil_plot['cd']], 'o-', 
                 color=xfoil_color, label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_plot:
        ax2.plot(rustfoil_plot['alpha'], [cd * 10000 for cd in rustfoil_plot['cd']], 's--',
                 color=rustfoil_color, label='RustFoil', markersize=4, linewidth=1.5)
    ax2.set_xlabel('α (deg)')
    ax2.set_ylabel('Cd (counts)')
    ax2.set_title('Drag Coefficient')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Cl vs Cd (Drag Polar)
    ax3 = fig.add_subplot(gs[1, 0])
    if xfoil_plot:
        ax3.plot([cd * 10000 for cd in xfoil_plot['cd']], xfoil_plot['cl'], 'o-',
                 color=xfoil_color, label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_plot:
        ax3.plot([cd * 10000 for cd in rustfoil_plot['cd']], rustfoil_plot['cl'], 's--',
                 color=rustfoil_color, label='RustFoil', markersize=4, linewidth=1.5)
    ax3.set_xlabel('Cd (counts)')
    ax3.set_ylabel('Cl')
    ax3.set_title('Drag Polar')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. L/D vs Cl
    ax4 = fig.add_subplot(gs[1, 1])
    if xfoil_plot:
        ld_xfoil = [cl/cd if cd > 1e-6 else 0 for cl, cd in zip(xfoil_plot['cl'], xfoil_plot['cd'])]
        ax4.plot(xfoil_plot['cl'], ld_xfoil, 'o-', color=xfoil_color, 
                 label='XFOIL', markersize=4, linewidth=1.5)
    if rustfoil_plot:
        ld_rustfoil = [cl/cd if cd > 1e-6 else 0 for cl, cd in zip(rustfoil_plot['cl'], rustfoil_plot['cd'])]
        ax4.plot(rustfoil_plot['cl'], ld_rustfoil, 's--', color=rustfoil_color,
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
        "cl_max_delta": None,
        "stall_alpha_delta": None,
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
                rf_cl = rustfoil_data['cl'][j]
                rf_cd = rustfoil_data['cd'][j]
                if not is_valid_number(rf_cl) or not is_valid_number(rf_cd):
                    continue
                cl_err = rf_cl - xfoil_data['cl'][i]
                cd_err = rf_cd - xfoil_data['cd'][i]
                cl_errs.append(cl_err)
                cd_errs.append(cd_err)
                break
    
    if cl_errs:
        errors["cl_rms"] = np.sqrt(np.mean(np.array(cl_errs)**2))
        errors["cd_rms"] = np.sqrt(np.mean(np.array(cd_errs)**2))
        errors["cl_max_err"] = max(abs(e) for e in cl_errs)
        errors["cd_max_err"] = max(abs(e) for e in cd_errs)
        errors["n_points"] = len(cl_errs)
    
    # Stall metrics
    xfoil_stall = compute_stall_metrics(xfoil_data)
    rustfoil_stall = compute_stall_metrics(rustfoil_data)
    if xfoil_stall and rustfoil_stall:
        errors["cl_max_delta"] = rustfoil_stall["cl_max"] - xfoil_stall["cl_max"]
        errors["stall_alpha_delta"] = rustfoil_stall["alpha_cl_max"] - xfoil_stall["alpha_cl_max"]

    return errors


def compute_stall_metrics(data: dict) -> dict:
    """Return simple stall metrics (peak Cl and its alpha)."""
    if not data or not data.get("alpha") or not data.get("cl"):
        return None
    filtered = filter_polar_data(data)
    if not filtered:
        return None
    cl = np.array(filtered["cl"])
    alpha = np.array(filtered["alpha"])
    if cl.size == 0 or alpha.size != cl.size:
        return None

    max_idx = int(np.argmax(cl))
    cl_max = float(cl[max_idx])
    alpha_max = float(alpha[max_idx])
    cl_drop = 0.0
    if max_idx < len(cl) - 1:
        cl_drop = float(cl_max - np.min(cl[max_idx + 1:]))

    return {
        "cl_max": cl_max,
        "alpha_cl_max": alpha_max,
        "cl_drop": cl_drop,
    }


def is_valid_number(value) -> bool:
    if value is None:
        return False
    try:
        return np.isfinite(value)
    except Exception:
        return False


def filter_polar_data(data: dict) -> dict:
    """Filter out None/non-finite entries from polar data."""
    if not data or not data.get("alpha"):
        return None
    keys = ["alpha", "cl", "cd", "cdp", "cm", "xtr_top", "xtr_bot"]
    length = min(len(data.get(k, [])) for k in keys if k in data)
    if length == 0:
        return None

    filtered = {k: [] for k in keys}
    for i in range(length):
        row = {k: data.get(k, [None] * length)[i] for k in keys}
        if not all(is_valid_number(row[k]) for k in ["alpha", "cl", "cd", "cm"]):
            continue
        for k in keys:
            filtered[k].append(row[k])

    if not filtered["alpha"]:
        return None
    return filtered


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
            if result['reynolds'] == 3e6:
                linestyle = '-'
            elif result['reynolds'] == 1e6:
                linestyle = '--'
            else:
                linestyle = ':'
            
            xf = filter_polar_data(result.get('xfoil'))
            rf = filter_polar_data(result.get('rustfoil'))
            
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
    json_path = OUTPUT_DIR / "xfoil_comparison_data.json"
    if json_path.exists():
        try:
            with open(json_path, "r") as f:
                all_data = json.load(f)
        except Exception:
            all_data = {}
    
    alpha_start, alpha_end, alpha_step = ALPHA_RANGE
    
    for airfoil in AIRFOILS:
        print(f"\n{airfoil['name']}")
        print("-" * 40)
        all_data.setdefault(airfoil["name"], {})
        
        for reynolds in REYNOLDS:
            print(f"\n  Re = {reynolds:.0e}")
            re_key = f"Re{reynolds:.0e}"
            cached = all_data.get(airfoil["name"], {}).get(re_key, {})
            has_cached = cached.get("xfoil") and cached.get("rustfoil")
            
            if has_cached:
                print("    Using cached results.")
                xfoil_data = cached["xfoil"]
                rustfoil_data = cached["rustfoil"]
            else:
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
            all_data[airfoil['name']][re_key] = {
                'xfoil': xfoil_data,
                'rustfoil': rustfoil_data,
                'errors': result.get('errors'),
            }
            with open(json_path, 'w') as f:
                json.dump(all_data, f, indent=2)
    
    # Create summary plot
    summary_path = OUTPUT_DIR / "summary_comparison.png"
    create_summary_plot(all_results, summary_path)
    
    # Save JSON data for app
    print(f"\nSaved JSON data: {json_path.name}")
    
    # Print error summary
    print("\n" + "=" * 70)
    print("ERROR SUMMARY")
    print("=" * 70)
    print(f"{'Airfoil':<12} {'Re':<10} {'Cl RMS':<10} {'Cd RMS':<12} {'Cd counts':<10} {'dCl_max':<10} {'dα_stall':<10}")
    print("-" * 74)
    
    for result in all_results:
        if result.get('errors'):
            e = result['errors']
            dcl_max = e.get('cl_max_delta')
            dstall = e.get('stall_alpha_delta')
            dcl_str = f"{dcl_max:.4f}" if dcl_max is not None else "n/a"
            dstall_str = f"{dstall:.2f}" if dstall is not None else "n/a"
            print(
                f"{result['naca']:<12} {result['reynolds']:<10.0e} "
                f"{e['cl_rms']:<10.4f} {e['cd_rms']:<12.6f} {e['cd_rms']*10000:<10.1f} "
                f"{dcl_str:<10} {dstall_str:<10}"
            )
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
