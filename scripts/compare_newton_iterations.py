#!/usr/bin/env python3
"""
Newton Iteration Brute Force Comparison

Phase 3: Compare XFOIL and RustFoil Newton iterations at each alpha.
Find the exact point of divergence.
"""
import subprocess
import json
import os
import sys
from pathlib import Path

# Configuration - MUST MATCH between XFOIL and RustFoil
RE = 1e6  # Reynolds number
ALPHAS = [0, 2, 4, 8]  # Test angles

def run_xfoil(alpha):
    """Run XFOIL at given alpha and return debug JSON path"""
    print(f"  Running XFOIL at α={alpha}°, Re={RE:.0e}...")
    
    os.chdir("/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin")
    
    # Remove old debug file
    debug_file = "xfoil_debug.json"
    if os.path.exists(debug_file):
        os.remove(debug_file)
    
    script = f"""NACA 0012
PANE
OPER
VISC {RE:.0f}
ITER 100
ALFA {alpha}

QUIT
"""
    proc = subprocess.run(
        ["./xfoil_instrumented"],
        input=script,
        capture_output=True,
        text=True,
        timeout=60
    )
    
    output_path = f"/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin/xfoil_debug_a{alpha}.json"
    if os.path.exists(debug_file):
        os.rename(debug_file, output_path)
        return output_path
    return None

def run_rustfoil(alpha):
    """Run RustFoil at given alpha and return debug JSON path"""
    print(f"  Running RustFoil at α={alpha}°, Re={RE:.0e}...")
    
    os.chdir("/Users/harry/flexfoil-boundary-layer")
    
    env = os.environ.copy()
    env["RUSTFOIL_TEST_ALPHA"] = str(alpha)
    env["RUSTFOIL_FULL_TRACE"] = "1"
    env["RUSTFOIL_RE"] = str(int(RE))  # Pass Re to test
    
    proc = subprocess.run(
        ["cargo", "test", "--release", "-p", "rustfoil-solver",
         "test_newton_iteration_trace", "--", "--nocapture"],
        capture_output=True,
        text=True,
        timeout=120,
        env=env
    )
    
    output_path = f"crates/rustfoil-solver/traces/rustfoil_new/rustfoil_alpha_{alpha}.json"
    if os.path.exists(output_path):
        return output_path
    return None

def load_json(path):
    """Load JSON events from file"""
    print(f"    Loading {path}...")
    with open(path) as f:
        content = f.read()
    print(f"    File size: {len(content)/1e6:.1f} MB")
    data = json.loads(content)
    print(f"    Parsed: {len(data.get('events', []))} events")
    return data

def get_blvar_by_iteration(events):
    """Extract BLVAR events grouped by (iteration, side, ibl)"""
    states = {}
    for e in events:
        if e.get("subroutine") != "BLVAR":
            continue
        key = (e.get("iteration", 0), e.get("side", 0), e.get("ibl", 0))
        if key not in states:
            inp = e.get("input", {})
            out = e.get("output", {})
            states[key] = {
                "x": inp.get("x", e.get("x", 0)),
                "theta": inp.get("theta", e.get("theta", 0)),
                "dstar": inp.get("delta_star", e.get("delta_star", 0)),
                "ue": inp.get("u", e.get("u", e.get("Ue", 0))),
                "hk": out.get("Hk", e.get("Hk", e.get("hk", 0))),
                "cf": out.get("Cf", e.get("Cf", e.get("cf", 0))),
            }
    return states

def compare_states(xf_states, rf_states, threshold_pct=5.0):
    """Compare BLVAR states and find first divergence"""
    divergences = []
    
    for key in sorted(set(xf_states.keys()) & set(rf_states.keys())):
        iter_num, side, ibl = key
        xf = xf_states[key]
        rf = rf_states[key]
        
        for var in ["theta", "ue", "hk"]:
            xf_val = xf.get(var, 0)
            rf_val = rf.get(var, 0)
            
            if xf_val != 0:
                pct_diff = (rf_val - xf_val) / xf_val * 100
                if abs(pct_diff) > threshold_pct:
                    divergences.append({
                        "iter": iter_num,
                        "side": side,
                        "ibl": ibl,
                        "var": var,
                        "xfoil": xf_val,
                        "rustfoil": rf_val,
                        "pct_diff": pct_diff,
                    })
    
    return divergences

def compare_at_alpha(alpha, xf_path, rf_path):
    """Compare XFOIL and RustFoil outputs at given alpha"""
    print(f"\n{'='*70}")
    print(f"ALPHA = {alpha}°")
    print(f"{'='*70}")
    
    if not xf_path or not os.path.exists(xf_path):
        print("  ERROR: No XFOIL data")
        return
    if not rf_path or not os.path.exists(rf_path):
        print("  ERROR: No RustFoil data")
        return
    
    xf_data = load_json(xf_path)
    rf_data = load_json(rf_path)
    
    print(f"  XFOIL events: {len(xf_data['events'])}")
    print(f"  RustFoil events: {len(rf_data['events'])}")
    
    # Check Reynolds numbers match
    xf_re = None
    rf_re = None
    for e in xf_data["events"][:50]:
        if e.get("subroutine") == "VISCAL":
            xf_re = e.get("reynolds")
            break
    for e in rf_data["events"][:50]:
        if e.get("subroutine") == "VISCAL":
            rf_re = e.get("reynolds")
            break
    
    if xf_re and rf_re:
        if abs(xf_re - rf_re) / xf_re > 0.01:
            print(f"  WARNING: Reynolds mismatch! XFOIL={xf_re:.0e}, RustFoil={rf_re:.0e}")
        else:
            print(f"  Reynolds: {xf_re:.0e} (matched)")
    
    # Extract and compare BLVAR states
    xf_states = get_blvar_by_iteration(xf_data["events"])
    rf_states = get_blvar_by_iteration(rf_data["events"])
    
    print(f"  XFOIL BLVAR states: {len(xf_states)}")
    print(f"  RustFoil BLVAR states: {len(rf_states)}")
    
    # Find divergences
    divergences = compare_states(xf_states, rf_states, threshold_pct=5.0)
    
    if divergences:
        # Group by iteration and find first
        by_iter = {}
        for d in divergences:
            iter_num = d["iter"]
            if iter_num not in by_iter:
                by_iter[iter_num] = []
            by_iter[iter_num].append(d)
        
        first_iter = min(by_iter.keys())
        print(f"\n  FIRST DIVERGENCE at iteration {first_iter}:")
        
        # Show first few divergences at that iteration
        for d in by_iter[first_iter][:5]:
            side_name = "upper" if d["side"] == 1 else "lower"
            print(f"    Station {d['ibl']} ({side_name}): {d['var']}")
            print(f"      XFOIL={d['xfoil']:.5e}, RustFoil={d['rustfoil']:.5e} ({d['pct_diff']:+.1f}%)")
        
        if len(by_iter[first_iter]) > 5:
            print(f"    ... and {len(by_iter[first_iter]) - 5} more divergences at this iteration")
    else:
        print("\n  No divergence > 5% found!")

def main():
    print("="*70)
    print("NEWTON ITERATION BRUTE FORCE COMPARISON")
    print(f"Re = {RE:.0e}")
    print("="*70)
    
    results = {}
    
    for alpha in ALPHAS:
        print(f"\n--- Alpha = {alpha}° ---")
        
        # Run both solvers
        xf_path = run_xfoil(alpha)
        rf_path = run_rustfoil(alpha)
        
        results[alpha] = (xf_path, rf_path)
    
    # Compare results
    print("\n" + "="*70)
    print("COMPARISON RESULTS")
    print("="*70)
    
    for alpha in ALPHAS:
        xf_path, rf_path = results[alpha]
        compare_at_alpha(alpha, xf_path, rf_path)

if __name__ == "__main__":
    main()
