#!/usr/bin/env python3
"""
Compare laminar-turbulent transition locations between RustFoil and XFOIL.
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Optional, Dict

# Paths
RUSTFOIL_BIN = Path(__file__).parent.parent / "target" / "release" / "rustfoil"
AIRFOIL_FILE = Path(__file__).parent.parent / "testdata" / "naca0012_repaneled.dat"
TRACES_DIR = Path(__file__).parent.parent / "traces" / "xfoil" / "naca0012" / "re3e06"

# Angles to compare
ANGLES = [0, 2, 4, 6, 8]
REYNOLDS = 3e6


def get_rustfoil_transition(alpha: float) -> Dict[str, Optional[float]]:
    """Run RustFoil and extract transition locations."""
    cmd = [
        str(RUSTFOIL_BIN),
        "viscous",
        str(AIRFOIL_FILE),
        "--alpha", str(alpha),
        "--re", str(int(REYNOLDS)),
        "--format", "json"
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=Path(__file__).parent.parent
        )
        
        if result.returncode != 0:
            print(f"Warning: RustFoil failed for α={alpha}°", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            return {"x_tr_upper": None, "x_tr_lower": None}
        
        # Parse JSON output
        try:
            data = json.loads(result.stdout)
            return {
                "x_tr_upper": data.get("x_tr_upper"),
                "x_tr_lower": data.get("x_tr_lower"),
                "cl": data.get("cl"),
                "cd": data.get("cd"),
            }
        except json.JSONDecodeError:
            print(f"Warning: Failed to parse RustFoil JSON for α={alpha}°", file=sys.stderr)
            return {"x_tr_upper": None, "x_tr_lower": None}
            
    except subprocess.TimeoutExpired:
        print(f"Warning: RustFoil timed out for α={alpha}°", file=sys.stderr)
        return {"x_tr_upper": None, "x_tr_lower": None}


def get_xfoil_transition(alpha: float) -> Dict[str, Optional[float]]:
    """Extract transition locations from XFOIL trace file."""
    # Format: alpha_+000.0.json, alpha_+002.0.json, etc.
    sign = "+" if alpha >= 0 else "-"
    alpha_str = f"{sign}{int(abs(alpha)):03d}.0"
    trace_file = TRACES_DIR / f"alpha_{alpha_str}.json"
    
    if not trace_file.exists():
        print(f"Warning: XFOIL trace file not found: {trace_file}", file=sys.stderr)
        return {"x_tr_upper": None, "x_tr_lower": None}
    
    try:
        with open(trace_file) as f:
            data = json.load(f)
        
        x_tr_upper = None
        x_tr_lower = None
        
        # Find TRANSITION events
        events = data.get("events", [])
        
        # Find maximum iteration number
        max_iteration = 0
        for event in events:
            if event.get("subroutine") == "TRANSITION":
                iteration = event.get("iteration", 0)
                max_iteration = max(max_iteration, iteration)
        
        # Extract transitions from final iteration
        for event in events:
            if (event.get("subroutine") == "TRANSITION" and 
                event.get("iteration") == max_iteration):
                side = event.get("side")
                x_transition = event.get("x_transition")
                
                if side == 1:  # Upper surface
                    x_tr_upper = x_transition
                elif side == 2:  # Lower surface
                    x_tr_lower = x_transition
        
        # Also check final state if available
        final_state = data.get("final_state")
        if final_state:
            if x_tr_upper is None:
                x_tr_upper = final_state.get("x_tr_upper")
            if x_tr_lower is None:
                x_tr_lower = final_state.get("x_tr_lower")
        
        return {
            "x_tr_upper": x_tr_upper,
            "x_tr_lower": x_tr_lower,
        }
        
    except (json.JSONDecodeError, IOError) as e:
        print(f"Warning: Failed to read XFOIL trace for α={alpha}°: {e}", file=sys.stderr)
        return {"x_tr_upper": None, "x_tr_lower": None}


def calculate_error(rf_val: Optional[float], xf_val: Optional[float]) -> Optional[float]:
    """Calculate percentage error."""
    if rf_val is None or xf_val is None:
        return None
    if xf_val == 0:
        return None
    return ((rf_val - xf_val) / xf_val) * 100.0


def format_value(val: Optional[float], fmt: str = ".4f") -> str:
    """Format a value for display."""
    if val is None:
        return "N/A"
    return f"{val:{fmt}}"


def main():
    print("Transition Location Comparison: RustFoil vs XFOIL")
    print("=" * 80)
    print(f"Airfoil: NACA 0012")
    print(f"Reynolds: {REYNOLDS:.0e}")
    print()
    
    results = []
    
    for alpha in ANGLES:
        print(f"Processing α = {alpha}°...", file=sys.stderr)
        
        rf_data = get_rustfoil_transition(alpha)
        xf_data = get_xfoil_transition(alpha)
        
        results.append({
            "alpha": alpha,
            "rf": rf_data,
            "xf": xf_data,
        })
    
    # Print comparison table
    print("\nTransition Location Comparison Table")
    print("-" * 80)
    print(f"{'Alpha':<8} {'RF x_tr_upper':<15} {'XF x_tr_upper':<15} {'Err%':<10} "
          f"{'RF x_tr_lower':<15} {'XF x_tr_lower':<15} {'Err%':<10}")
    print("-" * 80)
    
    for r in results:
        alpha = r["alpha"]
        rf_upper = r["rf"]["x_tr_upper"]
        xf_upper = r["xf"]["x_tr_upper"]
        rf_lower = r["rf"]["x_tr_lower"]
        xf_lower = r["xf"]["x_tr_lower"]
        
        err_upper = calculate_error(rf_upper, xf_upper)
        err_lower = calculate_error(rf_lower, xf_lower)
        
        print(f"{alpha:<8.0f} "
              f"{format_value(rf_upper):<15} {format_value(xf_upper):<15} {format_value(err_upper, '.2f'):<10} "
              f"{format_value(rf_lower):<15} {format_value(xf_lower):<15} {format_value(err_lower, '.2f'):<10}")
    
    # Analysis
    print("\n" + "=" * 80)
    print("Analysis")
    print("=" * 80)
    
    # Check for significant differences
    significant_diffs = []
    for r in results:
        alpha = r["alpha"]
        rf_upper = r["rf"]["x_tr_upper"]
        xf_upper = r["xf"]["x_tr_upper"]
        rf_lower = r["rf"]["x_tr_lower"]
        xf_lower = r["xf"]["x_tr_lower"]
        
        err_upper = calculate_error(rf_upper, xf_upper)
        err_lower = calculate_error(rf_lower, xf_lower)
        
        if err_upper is not None and abs(err_upper) > 10:
            significant_diffs.append({
                "alpha": alpha,
                "surface": "upper",
                "rf": rf_upper,
                "xf": xf_upper,
                "error": err_upper,
            })
        
        if err_lower is not None and abs(err_lower) > 10:
            significant_diffs.append({
                "alpha": alpha,
                "surface": "lower",
                "rf": rf_lower,
                "xf": xf_lower,
                "error": err_lower,
            })
    
    if significant_diffs:
        print("\nSignificant Differences (>10%):")
        for diff in significant_diffs:
            direction = "earlier" if diff["rf"] < diff["xf"] else "later"
            print(f"\n  α = {diff['alpha']}°, {diff['surface']} surface:")
            print(f"    RustFoil: {diff['rf']:.4f}")
            print(f"    XFOIL:    {diff['xf']:.4f}")
            print(f"    Error:    {diff['error']:.2f}%")
            print(f"    RustFoil transitions {direction} than XFOIL")
            
            if diff["rf"] < diff["xf"]:
                print(f"    → RustFoil has longer turbulent region → thicker BL")
            else:
                print(f"    → RustFoil has shorter turbulent region → thinner BL")
    else:
        print("\nNo significant differences found (<10% threshold)")
    
    # Impact on CL
    print("\n" + "-" * 80)
    print("Impact on Lift Coefficient (CL):")
    print("-" * 80)
    
    for r in results:
        alpha = r["alpha"]
        rf_cl = r["rf"].get("cl")
        rf_upper = r["rf"]["x_tr_upper"]
        xf_upper = r["xf"]["x_tr_upper"]
        rf_lower = r["rf"]["x_tr_lower"]
        xf_lower = r["xf"]["x_tr_lower"]
        
        if rf_cl is not None:
            print(f"\nα = {alpha}°:")
            print(f"  CL = {rf_cl:.4f}")
            
            if rf_upper is not None and xf_upper is not None:
                if rf_upper < xf_upper:
                    print(f"  Upper: RF transitions earlier ({rf_upper:.4f} vs {xf_upper:.4f})")
                    print(f"    → More turbulent flow on upper → thicker BL → reduced suction")
                elif rf_upper > xf_upper:
                    print(f"  Upper: RF transitions later ({rf_upper:.4f} vs {xf_upper:.4f})")
                    print(f"    → More laminar flow on upper → thinner BL → increased suction")
            
            if rf_lower is not None and xf_lower is not None:
                if rf_lower < xf_lower:
                    print(f"  Lower: RF transitions earlier ({rf_lower:.4f} vs {xf_lower:.4f})")
                    print(f"    → More turbulent flow on lower → thicker BL")
                elif rf_lower > xf_lower:
                    print(f"  Lower: RF transitions later ({rf_lower:.4f} vs {xf_lower:.4f})")
                    print(f"    → More laminar flow on lower → thinner BL")
    
    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print("\nTransition location differences affect boundary layer thickness:")
    print("  • Earlier transition → longer turbulent region → thicker BL → reduced lift")
    print("  • Later transition → shorter turbulent region → thinner BL → increased lift")
    print("\nIf RustFoil transitions earlier than XFOIL on the upper surface,")
    print("this would explain BL thickness over-prediction and CL under-prediction.")


if __name__ == "__main__":
    main()
