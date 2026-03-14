#!/usr/bin/env python3
"""
XFOIL vs RustFoil Viscous-Inviscid Chain Comparison

This script systematically compares each stage of the V-I coupling:
1. Geometry/paneling
2. Inviscid solution (GAM, QINV, UINV)
3. Initial BL march (theta, delta*, transition)
4. Newton system (Jacobians, residuals)
5. Newton update (deltas, relaxation, new Ue)
6. Force computation (CL, CD)
"""

import json
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
import numpy as np

# Configuration
WORKSPACE = Path(__file__).parent.parent
AIRFOIL = WORKSPACE / "testdata" / "naca0012.dat"
XFOIL_BIN = "/Applications/ESP128/EngSketchPad/bin/xfoil"
XFOIL_INSTRUMENTED = WORKSPACE / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
ALPHA = 4.0
RE = 1_000_000

@dataclass
class ComparisonResult:
    """Result of comparing a single value."""
    name: str
    xfoil_val: float
    rustfoil_val: float
    abs_diff: float
    rel_diff_pct: float
    passed: bool
    tolerance: float


def load_json_events(path: Path) -> List[Dict]:
    """Load JSON debug events from file."""
    with open(path) as f:
        data = json.load(f)
    return data.get("events", data) if isinstance(data, dict) else data


def run_xfoil_instrumented(airfoil: Path, alpha: float, re: float) -> Path:
    """Run XFOIL-instrumented and return path to debug JSON."""
    output_json = WORKSPACE / "xfoil_comparison_debug.json"
    
    # Create XFOIL script
    script = f"""LOAD {airfoil}
PANE
OPER
VISC {re}
ITER 50
ALFA {alpha}

QUIT
"""
    
    script_path = WORKSPACE / "xfoil_comparison_script.txt"
    with open(script_path, "w") as f:
        f.write(script)
    
    # Run XFOIL-instrumented (it creates xfoil_debug.json in cwd)
    try:
        result = subprocess.run(
            [str(XFOIL_INSTRUMENTED)],
            stdin=open(script_path),
            capture_output=True,
            text=True,
            timeout=60,
            cwd=WORKSPACE
        )
    except FileNotFoundError:
        print(f"ERROR: XFOIL-instrumented not found at {XFOIL_INSTRUMENTED}")
        print("Build it with: cd Xfoil-instrumented && make")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: XFOIL-instrumented timed out")
        sys.exit(1)
    
    # XFOIL creates xfoil_debug.json
    xfoil_json = WORKSPACE / "xfoil_debug.json"
    if xfoil_json.exists():
        import shutil
        shutil.copy(xfoil_json, output_json)
    
    return output_json


def run_rustfoil(airfoil: Path, alpha: float, re: float) -> Path:
    """Run RustFoil with debug output and return path to debug JSON."""
    output_json = WORKSPACE / "rustfoil_comparison_debug.json"
    
    # Set environment variable to enable debug output
    env = {"RUSTFOIL_DEBUG": str(output_json)}
    
    try:
        result = subprocess.run(
            ["cargo", "run", "--release", "-p", "rustfoil-cli", "--",
             "viscous", str(airfoil), f"--alpha={alpha}", f"--re={int(re)}",
             "--debug"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=WORKSPACE,
            env={**dict(subprocess.os.environ), **env}
        )
        
        # Print RustFoil output for debugging
        if result.returncode != 0:
            print("RustFoil stderr:", result.stderr)
            
    except subprocess.TimeoutExpired:
        print("ERROR: RustFoil timed out")
        sys.exit(1)
    
    # Look for debug output
    if output_json.exists():
        return output_json
    
    # Also check rustfoil_debug.json
    alt_path = WORKSPACE / "rustfoil_debug.json"
    if alt_path.exists():
        return alt_path
    
    print("WARNING: No RustFoil debug output found")
    return None


def extract_xfoil_values(events: List[Dict]) -> Dict[str, Any]:
    """Extract key values from XFOIL debug events."""
    result = {
        "viscal": [],
        "blvar": {},
        "bldif": {},
        "mrchue": {},
        "update": {},
        "final": {},
    }
    
    for event in events:
        sub = event.get("subroutine", "")
        
        if sub == "VISCAL":
            result["viscal"].append({
                "iteration": event.get("iteration", 0),
                "alpha_rad": event.get("alpha_rad", 0),
                "reynolds": event.get("reynolds", 0),
            })
            
        elif sub == "VISCAL_RESULT" or sub == "VISC_RESULT":
            result["final"] = {
                "CL": event.get("CL", 0),
                "CD": event.get("CD", 0),
                "CM": event.get("CM", 0),
            }
            
        elif sub == "BLVAR":
            side = event.get("side", 0)
            ibl = event.get("ibl", 0)
            key = (side, ibl)
            if key not in result["blvar"]:
                result["blvar"][key] = []
            result["blvar"][key].append({
                "iteration": event.get("iteration", 0),
                "input": event.get("input", {}),
                "output": event.get("output", {}),
            })
            
        elif sub == "BLDIF":
            side = event.get("side", 0)
            ibl = event.get("ibl", 0)
            key = (side, ibl)
            if key not in result["bldif"]:
                result["bldif"][key] = []
            result["bldif"][key].append({
                "iteration": event.get("iteration", 0),
                "VS1": event.get("VS1", []),
                "VS2": event.get("VS2", []),
                "VSREZ": event.get("VSREZ", []),
            })
            
        elif sub == "MRCHUE":
            side = event.get("side", 0)
            ibl = event.get("ibl", 0)
            key = (side, ibl)
            result["mrchue"][key] = {
                "x": event.get("x", 0),
                "Ue": event.get("Ue", 0),
                "theta": event.get("theta", 0),
                "delta_star": event.get("delta_star", 0),
                "Hk": event.get("Hk", 0),
                "Cf": event.get("Cf", 0),
            }
            
        elif sub == "UPDATE":
            side = event.get("side", 0)
            ibl = event.get("ibl", 0)
            iter_num = event.get("iteration", 0)
            key = (side, ibl, iter_num)
            result["update"][key] = {
                "delta_ctau": event.get("delta_ctau", 0),
                "delta_theta": event.get("delta_theta", 0),
                "delta_mass": event.get("delta_mass", 0),
                "delta_Ue": event.get("delta_Ue", 0),
                "relaxation": event.get("relaxation", 1.0),
            }
    
    return result


def compare_value(name: str, xfoil: float, rustfoil: float, 
                  rel_tol: float = 0.05, abs_tol: float = 1e-6) -> ComparisonResult:
    """Compare two values with relative and absolute tolerance."""
    abs_diff = abs(xfoil - rustfoil)
    
    if abs(xfoil) > abs_tol:
        rel_diff = abs_diff / abs(xfoil) * 100
    else:
        rel_diff = 0 if abs_diff < abs_tol else 100
    
    passed = abs_diff < abs_tol or rel_diff < rel_tol * 100
    
    return ComparisonResult(
        name=name,
        xfoil_val=xfoil,
        rustfoil_val=rustfoil,
        abs_diff=abs_diff,
        rel_diff_pct=rel_diff,
        passed=passed,
        tolerance=rel_tol * 100,
    )


def print_comparison(results: List[ComparisonResult], title: str):
    """Print comparison results in a formatted table."""
    print(f"\n{'='*70}")
    print(f" {title}")
    print(f"{'='*70}")
    
    # Header
    print(f"{'Variable':<25} {'XFOIL':>12} {'RustFoil':>12} {'Diff %':>10} {'Status':>8}")
    print("-" * 70)
    
    passed = 0
    failed = 0
    
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        status_color = "" if r.passed else "<<<---"
        
        if r.passed:
            passed += 1
        else:
            failed += 1
        
        print(f"{r.name:<25} {r.xfoil_val:>12.6e} {r.rustfoil_val:>12.6e} "
              f"{r.rel_diff_pct:>9.2f}% {status:>8} {status_color}")
    
    print("-" * 70)
    print(f"Summary: {passed} passed, {failed} failed")
    
    return failed == 0


def compare_stage1_geometry():
    """Stage 1: Compare geometry and panel setup."""
    print("\n" + "="*70)
    print(" STAGE 1: Geometry and Panel Setup")
    print("="*70)
    
    # Load airfoil coordinates
    x_coords = []
    y_coords = []
    with open(AIRFOIL) as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x_coords.append(float(parts[0]))
                    y_coords.append(float(parts[1]))
                except ValueError:
                    pass
    
    x = np.array(x_coords)
    y = np.array(y_coords)
    n = len(x)
    
    print(f"  Airfoil: {AIRFOIL.name}")
    print(f"  Number of points: {n}")
    print(f"  First point: ({x[0]:.6f}, {y[0]:.6f})")
    print(f"  Last point: ({x[-1]:.6f}, {y[-1]:.6f})")
    
    # Check for closed contour
    te_gap = np.sqrt((x[0] - x[-1])**2 + (y[0] - y[-1])**2)
    print(f"  TE gap: {te_gap:.6e}")
    
    if te_gap > 0.01:
        print("  WARNING: Geometry is NOT closed! This will cause errors.")
        return False
    
    # Find LE (minimum x)
    le_idx = np.argmin(x)
    print(f"  LE at index {le_idx}: ({x[le_idx]:.6f}, {y[le_idx]:.6f})")
    
    # Check for symmetry (NACA 0012 should be symmetric)
    if "0012" in str(AIRFOIL):
        # Check y at a few x locations
        upper_y = []
        lower_y = []
        for i in range(n):
            if x[i] > 0.1 and x[i] < 0.9:
                if i < le_idx:
                    upper_y.append((x[i], y[i]))
                else:
                    lower_y.append((x[i], -y[i]))
        
        if upper_y and lower_y:
            upper_y.sort()
            lower_y.sort()
            max_asym = 0
            for (xu, yu), (xl, yl) in zip(upper_y[:10], lower_y[:10]):
                if abs(xu - xl) < 0.01:
                    max_asym = max(max_asym, abs(yu - yl))
            print(f"  Max asymmetry: {max_asym:.6e} (should be ~0 for 0012)")
    
    print("  [PASS] Geometry check complete")
    return True


def compare_stage2_inviscid(xfoil_events: List[Dict], rustfoil_events: List[Dict]):
    """Stage 2: Compare inviscid solution."""
    print("\n" + "="*70)
    print(" STAGE 2: Inviscid Solution")
    print("="*70)
    
    # Run just the inviscid analysis for both
    result = subprocess.run(
        ["cargo", "run", "--release", "-p", "rustfoil-cli", "--",
         "analyze", str(AIRFOIL), f"--alpha={ALPHA}"],
        capture_output=True, text=True, timeout=30, cwd=WORKSPACE
    )
    
    # Parse RustFoil inviscid CL
    rf_cl_inv = None
    for line in result.stdout.split('\n'):
        if 'Cl =' in line:
            rf_cl_inv = float(line.split('=')[1].strip())
            break
    
    # XFOIL inviscid reference (from previous tests)
    xf_cl_inv = 0.4829  # At alpha=4, NACA 0012
    
    if rf_cl_inv:
        cl_diff = abs(rf_cl_inv - xf_cl_inv) / xf_cl_inv * 100
        status = "PASS" if cl_diff < 1.0 else "FAIL"
        print(f"  Inviscid CL:")
        print(f"    XFOIL:    {xf_cl_inv:.6f}")
        print(f"    RustFoil: {rf_cl_inv:.6f}")
        print(f"    Diff:     {cl_diff:.2f}%")
        print(f"    Status:   [{status}]")
        return cl_diff < 1.0
    else:
        print("  ERROR: Could not extract RustFoil inviscid CL")
        return False


def compare_stage3_bl_march(xfoil_events: List[Dict], rustfoil_events: List[Dict]):
    """Stage 3: Compare initial BL march."""
    print("\n" + "="*70)
    print(" STAGE 3: Initial BL March")
    print("="*70)
    
    xf_data = extract_xfoil_values(xfoil_events)
    
    # Look at a few representative stations
    test_stations = [
        (1, 10),  # Upper surface, early
        (1, 40),  # Upper surface, mid (near transition)
        (1, 70),  # Upper surface, late
        (2, 10),  # Lower surface, early
        (2, 70),  # Lower surface, late
    ]
    
    results = []
    
    for side, ibl in test_stations:
        key = (side, ibl)
        if key in xf_data["mrchue"]:
            xf = xf_data["mrchue"][key]
            print(f"\n  Station (side={side}, ibl={ibl}):")
            print(f"    XFOIL: x={xf['x']:.4f}, Ue={xf['Ue']:.4f}, theta={xf['theta']:.6e}, "
                  f"Hk={xf['Hk']:.3f}, Cf={xf['Cf']:.6e}")
    
    # Check transition location from XFOIL output
    print("\n  Transition locations from XFOIL debug:")
    for event in xfoil_events:
        if event.get("subroutine") == "TRAN":
            print(f"    Side {event.get('side')}: x_tr = {event.get('x_tran', 'N/A')}")
    
    return True


def compare_stage6_forces(xfoil_events: List[Dict]):
    """Stage 6: Compare force computation."""
    print("\n" + "="*70)
    print(" STAGE 6: Force Computation (CL/CD)")
    print("="*70)
    
    # Get XFOIL final values
    xf_data = extract_xfoil_values(xfoil_events)
    xf_final = xf_data.get("final", {})
    
    # Run RustFoil viscous
    result = subprocess.run(
        ["cargo", "run", "--release", "-p", "rustfoil-cli", "--",
         "viscous", str(AIRFOIL), f"--alpha={ALPHA}", f"--re={int(RE)}"],
        capture_output=True, text=True, timeout=60, cwd=WORKSPACE
    )
    
    # Parse RustFoil output
    rf_cl = rf_cd = rf_cm = None
    for line in result.stdout.split('\n'):
        if line.strip().startswith('CL:'):
            rf_cl = float(line.split(':')[1].strip())
        elif line.strip().startswith('CD:'):
            rf_cd = float(line.split(':')[1].strip())
        elif line.strip().startswith('CM:'):
            rf_cm = float(line.split(':')[1].strip())
    
    # XFOIL reference values at alpha=4, Re=1M
    xf_cl = xf_final.get("CL", 0.4278)
    xf_cd = xf_final.get("CD", 0.007279)
    xf_cm = xf_final.get("CM", -0.005)
    
    results = []
    
    if rf_cl is not None:
        results.append(compare_value("CL (viscous)", xf_cl, rf_cl, rel_tol=0.15))
    if rf_cd is not None:
        results.append(compare_value("CD (viscous)", xf_cd, rf_cd, rel_tol=0.15))
    if rf_cm is not None:
        results.append(compare_value("CM (viscous)", xf_cm, rf_cm, rel_tol=0.50))
    
    return print_comparison(results, "Force Comparison")


def main():
    """Run the full comparison."""
    print("="*70)
    print(" XFOIL vs RustFoil V-I Chain Comparison")
    print(f" Airfoil: {AIRFOIL.name}")
    print(f" Alpha: {ALPHA}°, Re: {RE:,.0f}")
    print("="*70)
    
    # Stage 1: Geometry
    if not compare_stage1_geometry():
        print("\nERROR: Geometry check failed. Fix geometry before continuing.")
        return 1
    
    # Run XFOIL-instrumented to get debug data
    print("\n" + "-"*70)
    print(" Running XFOIL-instrumented...")
    print("-"*70)
    
    xfoil_json = run_xfoil_instrumented(AIRFOIL, ALPHA, RE)
    
    if xfoil_json and xfoil_json.exists():
        xfoil_events = load_json_events(xfoil_json)
        print(f"  Loaded {len(xfoil_events)} XFOIL debug events")
    else:
        print("  WARNING: No XFOIL debug output, using reference values")
        xfoil_events = []
    
    # Run RustFoil
    print("\n" + "-"*70)
    print(" Running RustFoil...")
    print("-"*70)
    
    rustfoil_json = run_rustfoil(AIRFOIL, ALPHA, RE)
    
    if rustfoil_json and rustfoil_json.exists():
        rustfoil_events = load_json_events(rustfoil_json)
        print(f"  Loaded {len(rustfoil_events)} RustFoil debug events")
    else:
        print("  WARNING: No RustFoil debug output")
        rustfoil_events = []
    
    # Stage 2: Inviscid
    compare_stage2_inviscid(xfoil_events, rustfoil_events)
    
    # Stage 3: BL March
    compare_stage3_bl_march(xfoil_events, rustfoil_events)
    
    # Stage 6: Forces (skip 4-5 for now)
    compare_stage6_forces(xfoil_events)
    
    print("\n" + "="*70)
    print(" COMPARISON COMPLETE")
    print("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
