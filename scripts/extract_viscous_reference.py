#!/usr/bin/env python3
"""
Extract viscous solver reference data from XFOIL-instrumented.

Runs the FORTRAN XFOIL binary with viscous analysis enabled and parses
the JSON debug output to create reference data for validating the Rust
viscous solver implementation.

Usage:
    python extract_viscous_reference.py [--output testdata/viscous_ref.json]
    python extract_viscous_reference.py --reynolds 3e6 --mach 0.0 --ncrit 9.0
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional


def run_xfoil(xfoil_bin: Path, commands: str, workdir: Path, timeout: int = 60) -> tuple:
    """Run XFOIL with given commands and return stdout, stderr."""
    result = subprocess.run(
        [str(xfoil_bin)],
        input=commands,
        capture_output=True,
        text=True,
        cwd=str(workdir),
        timeout=timeout
    )
    return result.stdout, result.stderr, result.returncode


def extract_viscous_reference(
    panels: int = 160,
    reynolds: float = 3e6,
    mach: float = 0.0,
    ncrit: float = 9.0,
    alphas: list = None,
    airfoil: str = "NACA 0012"
) -> dict:
    """Run XFOIL-instrumented with viscous analysis and extract reference data."""
    
    if alphas is None:
        alphas = [0.0, 4.0]
    
    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    xfoil_bin = project_root / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
    
    if not xfoil_bin.exists():
        raise FileNotFoundError(f"XFOIL binary not found: {xfoil_bin}")
    
    # Build XFOIL commands for viscous analysis
    alpha_commands = "\n".join([f"ALFA {a}" for a in alphas])
    
    commands = f"""{airfoil}
PPAR
N {panels}


OPER
VISC {reynolds}
MACH {mach}
VPAR
N {ncrit}

ITER 50
{alpha_commands}

QUIT
"""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        print(f"Running XFOIL-instrumented (viscous)...", file=sys.stderr)
        print(f"  Re = {reynolds:.2e}, M = {mach}, Ncrit = {ncrit}", file=sys.stderr)
        print(f"  Alphas: {alphas}", file=sys.stderr)
        
        stdout, stderr, rc = run_xfoil(xfoil_bin, commands, tmpdir)
        
        if rc != 0:
            print(f"XFOIL exited with code {rc}", file=sys.stderr)
            print(f"stderr: {stderr[:500]}", file=sys.stderr)
        
        # Look for debug output file
        debug_file = tmpdir / "xfoil_debug.json"
        
        if not debug_file.exists():
            # Check current directory too
            debug_file = Path("xfoil_debug.json")
        
        if debug_file.exists():
            print(f"Found debug output: {debug_file}", file=sys.stderr)
            content = debug_file.read_text()
            
            # Parse JSON
            if content.strip():
                try:
                    data = json.loads(content)
                    events = data.get("events", [])
                    return process_viscous_events(
                        events, panels, reynolds, mach, ncrit, alphas, airfoil
                    )
                except json.JSONDecodeError as e:
                    print(f"JSON parse error: {e}", file=sys.stderr)
                    # Save raw content for debugging
                    with open("/tmp/xfoil_viscous_debug_raw.json", "w") as f:
                        f.write(content)
                    print(f"Raw content saved to /tmp/xfoil_viscous_debug_raw.json", file=sys.stderr)
                    return {"error": str(e), "raw_content": content[:1000]}
        else:
            print("No debug file found", file=sys.stderr)
            print(f"XFOIL stdout:\n{stdout}", file=sys.stderr)
            
        return {"error": "No debug output generated"}


def process_viscous_events(
    events: list,
    panels: int,
    reynolds: float,
    mach: float,
    ncrit: float,
    alphas: list,
    airfoil: str
) -> dict:
    """Process XFOIL debug events into structured viscous reference data."""
    
    reference = {
        "test": "viscous_solver",
        "source": "xfoil_instrumented",
        "n_panels": panels,
        "reynolds": reynolds,
        "mach": mach,
        "ncrit": ncrit,
        "airfoil": airfoil.replace(" ", ""),
        "cases": []
    }
    
    # Group events by alpha
    current_case = None
    
    for event in events:
        if not isinstance(event, dict):
            continue
            
        sub = event.get("subroutine", "")
        
        # Start of a new viscous solve
        if sub == "VISCAL":
            alpha_rad = event.get("alpha_rad", 0.0)
            alpha_deg = alpha_rad * 180.0 / 3.14159265359
            
            # Start new case
            current_case = {
                "alpha_deg": round(alpha_deg, 4),
                "alpha_rad": alpha_rad,
                "inviscid": {},
                "stagnation": {},
                "bl_setup": {},
                "transitions": [],
                "iterations": [],
                "final": {}
            }
            reference["cases"].append(current_case)
            continue
        
        if current_case is None:
            # Store pre-viscous inviscid data
            if sub == "SPECAL":
                reference["inviscid_base"] = {
                    "alpha_rad": event.get("alpha_rad"),
                    "n": event.get("n"),
                    "gam_sample": event.get("GAM", [])[:20],
                    "qinv_sample": event.get("QINV", [])[:20]
                }
            elif sub == "STFIND":
                reference["stagnation_base"] = {
                    "ist": event.get("IST"),
                    "sst": event.get("SST"),
                    "xst": event.get("XST"),
                    "yst": event.get("YST")
                }
            elif sub == "IBLPAN":
                reference["bl_setup_base"] = {
                    "ist": event.get("IST"),
                    "nbl_upper": event.get("NBL_upper"),
                    "nbl_lower": event.get("NBL_lower"),
                    "iblte_upper": event.get("IBLTE_upper"),
                    "iblte_lower": event.get("IBLTE_lower")
                }
            continue
        
        # Within a viscous case
        if sub == "SPECAL":
            current_case["inviscid"]["cl"] = event.get("CL")
            current_case["inviscid"]["gam_sample"] = event.get("GAM", [])[:20]
            current_case["inviscid"]["qinv_sample"] = event.get("QINV", [])[:20]
            
        elif sub == "STFIND":
            current_case["stagnation"]["ist"] = event.get("IST")
            current_case["stagnation"]["sst"] = event.get("SST")
            current_case["stagnation"]["xst"] = event.get("XST")
            current_case["stagnation"]["yst"] = event.get("YST")
            
        elif sub == "IBLPAN":
            current_case["bl_setup"]["ist"] = event.get("IST")
            current_case["bl_setup"]["nbl_upper"] = event.get("NBL_upper")
            current_case["bl_setup"]["nbl_lower"] = event.get("NBL_lower")
            current_case["bl_setup"]["iblte_upper"] = event.get("IBLTE_upper")
            current_case["bl_setup"]["iblte_lower"] = event.get("IBLTE_lower")
            
        elif sub == "TRANSITION":
            current_case["transitions"].append({
                "side": event.get("side"),
                "itran": event.get("ITRAN"),
                "x_transition": event.get("x_transition"),
                "forced": event.get("forced", False)
            })
            
        elif sub == "VISCAL_RESULT":
            current_case["iterations"].append({
                "iteration": event.get("iteration"),
                "rms_residual": event.get("rms_residual"),
                "max_residual": event.get("max_residual"),
                "cl": event.get("CL"),
                "cd": event.get("CD"),
                "cm": event.get("CM")
            })
            
        elif sub == "VISCOUS_FINAL":
            current_case["final"] = {
                "cl": event.get("CL"),
                "cd": event.get("CD"),
                "cm": event.get("CM"),
                "cd_form": event.get("CD_form"),
                "cd_friction": event.get("CD_friction"),
                "x_tr_upper": event.get("x_tr_upper"),
                "x_tr_lower": event.get("x_tr_lower"),
                "x_sep_upper": event.get("x_sep_upper"),
                "x_sep_lower": event.get("x_sep_lower"),
                "iterations": event.get("iterations"),
                "converged": event.get("converged")
            }
            
        elif sub == "MRCHUE":
            # Store first few MRCHUE samples for BL init validation
            if "mrchue_samples" not in current_case:
                current_case["mrchue_samples"] = []
            if len(current_case["mrchue_samples"]) < 20:
                current_case["mrchue_samples"].append({
                    "side": event.get("side"),
                    "ibl": event.get("ibl"),
                    "x": event.get("x"),
                    "ue": event.get("Ue"),
                    "theta": event.get("theta"),
                    "delta_star": event.get("delta_star"),
                    "hk": event.get("Hk"),
                    "cf": event.get("Cf"),
                    "transitional": event.get("transitional")
                })
    
    return reference


def main():
    parser = argparse.ArgumentParser(description="Extract XFOIL viscous reference data")
    parser.add_argument("--panels", type=int, default=160, help="Number of panels")
    parser.add_argument("--reynolds", type=float, default=3e6, help="Reynolds number")
    parser.add_argument("--mach", type=float, default=0.0, help="Mach number")
    parser.add_argument("--ncrit", type=float, default=9.0, help="Ncrit for transition")
    parser.add_argument("--alpha", type=float, nargs="+", default=[0.0, 4.0],
                        help="Angles of attack (deg)")
    parser.add_argument("--airfoil", type=str, default="NACA 0012", help="Airfoil spec")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    
    ref_data = extract_viscous_reference(
        panels=args.panels,
        reynolds=args.reynolds,
        mach=args.mach,
        ncrit=args.ncrit,
        alphas=args.alpha,
        airfoil=args.airfoil
    )
    
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path(__file__).parent.parent / "testdata" / "viscous_ref.json"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(ref_data, f, indent=2)
    
    print(f"\nReference data written to {output_path}", file=sys.stderr)
    
    # Print summary
    if "cases" in ref_data:
        print(f"\nExtracted {len(ref_data['cases'])} viscous cases:", file=sys.stderr)
        for case in ref_data["cases"]:
            final = case.get("final", {})
            alpha = case.get("alpha_deg", "?")
            cl = final.get("cl")
            cd = final.get("cd")
            xtr_u = final.get("x_tr_upper")
            xtr_l = final.get("x_tr_lower")
            conv = "yes" if final.get("converged") else "no"
            
            # Format values, handling None/missing
            def fmt(val, spec):
                if val is None:
                    return "?"
                try:
                    return f"{float(val):{spec}}"
                except (ValueError, TypeError):
                    return str(val)
            
            print(f"  α={fmt(alpha, '6.2f')}°: CL={fmt(cl, '.4f')}, CD={fmt(cd, '.5f')}, "
                  f"xtr_u={fmt(xtr_u, '.4f')}, xtr_l={fmt(xtr_l, '.4f')}, converged={conv}",
                  file=sys.stderr)
    
    if "error" in ref_data:
        print(f"\nError: {ref_data['error']}", file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
