#!/usr/bin/env python3
"""
Compare H (shape factor) and Cf (skin friction) distributions between XFOIL and RustFoil
at high angle of attack to diagnose stall prediction differences.

Key indicators:
- H > 2.5-3.0 in turbulent region → separation
- Cf < 0 → reversed flow (separation)
- H growing rapidly toward TE → approaching stall
"""

import json
import subprocess
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def run_xfoil_bl_dump(alpha: float, re: float = 3e6, foil: str = "0012") -> dict:
    """Run XFOIL and extract BL parameters."""
    
    xfoil_path = Path("/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin/xfoil_instrumented")
    
    # Create XFOIL commands
    commands = f"""NACA {foil}
OPER
VISC {re:.0f}
ITER 100
ALFA {alpha}

QUIT
"""
    
    try:
        result = subprocess.run(
            [str(xfoil_path)],
            input=commands,
            capture_output=True,
            text=True,
            timeout=60,
            cwd=xfoil_path.parent
        )
    except Exception as e:
        print(f"XFOIL execution error: {e}")
        return None
    
    # Load the debug JSON
    json_path = xfoil_path.parent / "xfoil_debug.json"
    if not json_path.exists():
        print("No XFOIL debug output found")
        return None
    
    with open(json_path) as f:
        return json.load(f)

def extract_bl_params(data: dict, source: str = "xfoil") -> tuple:
    """Extract H, Hk, Cf, x from MRCHUE events.
    
    MRCHUE events contain:
    - side: 1=upper, 2=lower
    - x: arc length position
    - theta, delta_star: BL thicknesses (H = delta_star/theta)
    - Hk: kinematic shape factor
    - Cf: skin friction coefficient
    """
    
    upper_x, upper_h, upper_cf, upper_hk = [], [], [], []
    lower_x, lower_h, lower_cf, lower_hk = [], [], [], []
    
    # Use MRCHUE events (final converged values from march)
    for event in data.get("events", []):
        if event.get("subroutine") == "MRCHUE":
            side = event.get("side", 0)
            x = event.get("x", 0)
            theta = event.get("theta", 0)
            delta_star = event.get("delta_star", 0)
            hk = event.get("Hk", 0)
            cf = event.get("Cf", 0)
            
            # Compute H = delta_star / theta
            h = delta_star / theta if theta > 1e-12 else 0
            
            if x > 0 and h > 0:
                if side == 1:  # Upper surface
                    upper_x.append(x)
                    upper_h.append(h)
                    upper_hk.append(hk)
                    upper_cf.append(cf)
                elif side == 2:  # Lower surface
                    lower_x.append(x)
                    lower_h.append(h)
                    lower_hk.append(hk)
                    lower_cf.append(cf)
    
    # Sort by x and remove duplicates (take last occurrence)
    def dedupe_sort(x_list, h_list, hk_list, cf_list):
        if not x_list:
            return [], [], [], []
        
        # Combine into dict to dedupe by x
        data_dict = {}
        for i, x in enumerate(x_list):
            x_key = round(x, 6)
            data_dict[x_key] = (x, h_list[i], hk_list[i], cf_list[i])
        
        # Sort by x
        sorted_items = sorted(data_dict.values(), key=lambda t: t[0])
        
        return (
            [t[0] for t in sorted_items],
            [t[1] for t in sorted_items],
            [t[2] for t in sorted_items],
            [t[3] for t in sorted_items]
        )
    
    upper_x, upper_h, upper_hk, upper_cf = dedupe_sort(upper_x, upper_h, upper_hk, upper_cf)
    lower_x, lower_h, lower_hk, lower_cf = dedupe_sort(lower_x, lower_h, lower_hk, lower_cf)
    
    return (upper_x, upper_h, upper_hk, upper_cf), (lower_x, lower_h, lower_hk, lower_cf)

def run_rustfoil_test(alpha: float) -> dict:
    """Run RustFoil test at specified alpha."""
    
    env = dict(__import__('os').environ)
    env["RUSTFOIL_DEBUG"] = "1"
    
    result = subprocess.run(
        ["cargo", "test", "--release", "-p", "rustfoil-solver", 
         "--", "--test-threads=1", "test_viscous_naca0012_single", "--nocapture"],
        capture_output=True,
        text=True,
        cwd="/Users/harry/flexfoil-boundary-layer",
        env=env,
        timeout=120
    )
    
    # Load the debug JSON
    json_path = Path("/Users/harry/flexfoil-boundary-layer/debug_test.json")
    if not json_path.exists():
        print("No RustFoil debug output found")
        return None
    
    with open(json_path) as f:
        return json.load(f)

def analyze_separation(x, h, cf, label: str):
    """Analyze separation indicators."""
    
    print(f"\n=== {label} ===")
    
    if not x:
        print("  No data")
        return
    
    # Find max H
    max_h_idx = np.argmax(h)
    print(f"  Max H = {h[max_h_idx]:.3f} at x = {x[max_h_idx]:.4f}")
    
    # Check for H > 2.5 (separation indicator)
    sep_indices = [i for i, hval in enumerate(h) if hval > 2.5]
    if sep_indices:
        print(f"  H > 2.5 at {len(sep_indices)} stations (first at x = {x[sep_indices[0]]:.4f})")
    else:
        print(f"  H never exceeds 2.5 (no separation indicator)")
    
    # Check for H > 3.5 (strong separation / laminar separation bubble)
    strong_sep = [i for i, hval in enumerate(h) if hval > 3.5]
    if strong_sep:
        print(f"  H > 3.5 at {len(strong_sep)} stations (LSB or strong separation)")
    
    # Check for negative Cf
    neg_cf_indices = [i for i, cfval in enumerate(cf) if cfval < 0]
    if neg_cf_indices:
        print(f"  Cf < 0 at {len(neg_cf_indices)} stations (first at x = {x[neg_cf_indices[0]]:.4f})")
        print(f"    → Reversed flow detected!")
    else:
        print(f"  Cf always positive (no flow reversal)")
    
    # Min Cf
    min_cf_idx = np.argmin(cf)
    print(f"  Min Cf = {cf[min_cf_idx]:.6f} at x = {x[min_cf_idx]:.4f}")

def plot_comparison(xfoil_upper, xfoil_lower, rustfoil_upper, rustfoil_lower, alpha: float):
    """Create comparison plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'BL Parameters Comparison at α = {alpha}°\nNACA 0012, Re = 3×10⁶', fontsize=14)
    
    # Upper surface H
    ax = axes[0, 0]
    if xfoil_upper[0]:
        ax.plot(xfoil_upper[0], xfoil_upper[1], 'b-o', markersize=3, label='XFOIL')
    if rustfoil_upper[0]:
        ax.plot(rustfoil_upper[0], rustfoil_upper[1], 'r-s', markersize=3, label='RustFoil')
    ax.axhline(y=2.5, color='orange', linestyle='--', alpha=0.7, label='Separation threshold')
    ax.axhline(y=3.5, color='red', linestyle='--', alpha=0.7, label='Strong separation')
    ax.set_xlabel('x/c')
    ax.set_ylabel('H (Shape Factor)')
    ax.set_title('Upper Surface - Shape Factor H')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1, 6)
    
    # Upper surface Cf
    ax = axes[0, 1]
    if xfoil_upper[0]:
        ax.plot(xfoil_upper[0], xfoil_upper[3], 'b-o', markersize=3, label='XFOIL')
    if rustfoil_upper[0]:
        ax.plot(rustfoil_upper[0], rustfoil_upper[3], 'r-s', markersize=3, label='RustFoil')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Separation (Cf=0)')
    ax.set_xlabel('x/c')
    ax.set_ylabel('Cf (Skin Friction)')
    ax.set_title('Upper Surface - Skin Friction Cf')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Lower surface H
    ax = axes[1, 0]
    if xfoil_lower[0]:
        ax.plot(xfoil_lower[0], xfoil_lower[1], 'b-o', markersize=3, label='XFOIL')
    if rustfoil_lower[0]:
        ax.plot(rustfoil_lower[0], rustfoil_lower[1], 'r-s', markersize=3, label='RustFoil')
    ax.axhline(y=2.5, color='orange', linestyle='--', alpha=0.7, label='Separation threshold')
    ax.set_xlabel('x/c')
    ax.set_ylabel('H (Shape Factor)')
    ax.set_title('Lower Surface - Shape Factor H')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1, 6)
    
    # Lower surface Cf
    ax = axes[1, 1]
    if xfoil_lower[0]:
        ax.plot(xfoil_lower[0], xfoil_lower[3], 'b-o', markersize=3, label='XFOIL')
    if rustfoil_lower[0]:
        ax.plot(rustfoil_lower[0], rustfoil_lower[3], 'r-s', markersize=3, label='RustFoil')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Separation (Cf=0)')
    ax.set_xlabel('x/c')
    ax.set_ylabel('Cf (Skin Friction)')
    ax.set_title('Lower Surface - Skin Friction Cf')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = f'/Users/harry/flexfoil-boundary-layer/scripts/stall_comparison_alpha{int(alpha)}.png'
    plt.savefig(output_path, dpi=150)
    print(f"\nPlot saved to: {output_path}")
    plt.close()

def main():
    alpha = 10.0  # High angle to see stall effects
    if len(sys.argv) > 1:
        alpha = float(sys.argv[1])
    
    print(f"Comparing stall physics at α = {alpha}°")
    print("=" * 60)
    
    # Run XFOIL
    print("\nRunning XFOIL...")
    xfoil_data = run_xfoil_bl_dump(alpha)
    
    if xfoil_data:
        xfoil_upper, xfoil_lower = extract_bl_params(xfoil_data, "xfoil")
        print(f"  XFOIL: {len(xfoil_upper[0])} upper stations, {len(xfoil_lower[0])} lower stations")
        
        analyze_separation(xfoil_upper[0], xfoil_upper[1], xfoil_upper[3], "XFOIL Upper Surface")
        analyze_separation(xfoil_lower[0], xfoil_lower[1], xfoil_lower[3], "XFOIL Lower Surface")
    else:
        xfoil_upper = ([], [], [], [])
        xfoil_lower = ([], [], [], [])
    
    # Run RustFoil
    print("\n" + "=" * 60)
    print("Running RustFoil...")
    rustfoil_data = run_rustfoil_test(alpha)
    
    if rustfoil_data:
        rustfoil_upper, rustfoil_lower = extract_bl_params(rustfoil_data, "rustfoil")
        print(f"  RustFoil: {len(rustfoil_upper[0])} upper stations, {len(rustfoil_lower[0])} lower stations")
        
        analyze_separation(rustfoil_upper[0], rustfoil_upper[1], rustfoil_upper[3], "RustFoil Upper Surface")
        analyze_separation(rustfoil_lower[0], rustfoil_lower[1], rustfoil_lower[3], "RustFoil Lower Surface")
    else:
        rustfoil_upper = ([], [], [], [])
        rustfoil_lower = ([], [], [], [])
    
    # Create plots
    print("\n" + "=" * 60)
    print("Creating comparison plots...")
    plot_comparison(xfoil_upper, xfoil_lower, rustfoil_upper, rustfoil_lower, alpha)
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    if xfoil_upper[0] and rustfoil_upper[0]:
        xfoil_max_h = max(xfoil_upper[1]) if xfoil_upper[1] else 0
        rustfoil_max_h = max(rustfoil_upper[1]) if rustfoil_upper[1] else 0
        
        print(f"\nUpper Surface Max H:")
        print(f"  XFOIL:    {xfoil_max_h:.3f}")
        print(f"  RustFoil: {rustfoil_max_h:.3f}")
        print(f"  Difference: {100*(rustfoil_max_h-xfoil_max_h)/xfoil_max_h:+.1f}%")
        
        if xfoil_max_h > 2.5 and rustfoil_max_h < 2.5:
            print("\n⚠️  XFOIL shows separation (H > 2.5), RustFoil does not!")
            print("   This explains why RustFoil doesn't predict stall.")
        elif rustfoil_max_h > 2.5 and xfoil_max_h > 2.5:
            print("\n✓ Both show H > 2.5, separation should be detected")
        else:
            print("\n  Neither shows separation at this angle")
        
        # Check Cf
        xfoil_min_cf = min(xfoil_upper[3]) if xfoil_upper[3] else 1
        rustfoil_min_cf = min(rustfoil_upper[3]) if rustfoil_upper[3] else 1
        
        print(f"\nUpper Surface Min Cf:")
        print(f"  XFOIL:    {xfoil_min_cf:.6f}")
        print(f"  RustFoil: {rustfoil_min_cf:.6f}")
        
        if xfoil_min_cf < 0 and rustfoil_min_cf >= 0:
            print("\n⚠️  XFOIL shows reversed flow (Cf < 0), RustFoil does not!")

if __name__ == "__main__":
    main()
