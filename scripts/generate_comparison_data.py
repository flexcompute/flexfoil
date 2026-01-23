#!/usr/bin/env python3
"""
Generate XFOIL reference data for multiple foils at multiple AoA.
Outputs JSON files that can be compared with RustFoil.
"""

import subprocess
import json
import os
import re
from pathlib import Path

# Configuration
FOILS = ["naca0012", "naca2412", "naca4412"]
ALPHAS = [-4, -2, 0, 2, 4, 6, 8, 10]
XFOIL_PATH = "./Xfoil/bin/xfoil"
OUTPUT_DIR = Path("scripts/comparison_data")

def run_xfoil_for_foil(foil_name: str, alphas: list[float]) -> dict:
    """Run XFOIL for a foil at multiple angles, return parsed data."""
    
    # Build XFOIL command sequence
    commands = []
    
    # Check if we have a buffer foil file
    buffer_file = f"{foil_name}_xfoil_paneled.dat"
    if os.path.exists(buffer_file):
        commands.append(f"LOAD {buffer_file}")
        commands.append("")  # Accept default name
        commands.append("PCOP")  # Direct copy to preserve exact geometry from buffer
    else:
        # Generate NACA foil - need to panel it
        naca_digits = foil_name.replace("naca", "")
        commands.append(f"NACA {naca_digits}")
        commands.append("PANE")  # Panel from generated points
    commands.append("OPER")
    
    # Create dump files for each alpha
    dump_files = []
    for alpha in alphas:
        dump_file = OUTPUT_DIR / f"{foil_name}_alpha{alpha:+.0f}.dump"
        commands.append(f"ALFA {alpha}")
        commands.append(f"DUMP {dump_file}")
        dump_files.append((alpha, dump_file))
    
    commands.append("")
    commands.append("QUIT")
    
    # Run XFOIL
    cmd_str = "\n".join(commands)
    try:
        result = subprocess.run(
            [XFOIL_PATH],
            input=cmd_str,
            capture_output=True,
            text=True,
            timeout=30
        )
        stdout = result.stdout
    except subprocess.TimeoutExpired:
        print(f"  XFOIL timeout for {foil_name}")
        return None
    
    # Parse CL values from XFOIL output
    cl_values = {}
    # Look for lines like "  CL =   0.4829" following "a =   4.000"
    alpha_pattern = re.compile(r'a\s*=\s*([-\d.]+)')
    cl_pattern = re.compile(r'CL\s*=\s*([-\d.]+)')
    
    current_alpha = None
    for line in stdout.split('\n'):
        alpha_match = alpha_pattern.search(line)
        if alpha_match:
            current_alpha = float(alpha_match.group(1))
        cl_match = cl_pattern.search(line)
        if cl_match and current_alpha is not None:
            cl_values[int(round(current_alpha))] = float(cl_match.group(1))
    
    # Parse dump files
    foil_data = {
        "foil": foil_name,
        "alphas": {}
    }
    
    for alpha, dump_file in dump_files:
        if dump_file.exists():
            data = parse_dump_file(dump_file)
            if data:
                # Add CL from XFOIL output
                data["cl"] = cl_values.get(alpha, 0.0)
                foil_data["alphas"][alpha] = data
                print(f"  {foil_name} α={alpha:+.0f}°: {len(data['s'])} points, CL={data['cl']:.4f}")
    
    return foil_data

def parse_dump_file(filepath: Path) -> dict:
    """Parse XFOIL DUMP output file."""
    s, x, y, ue_vinf = [], [], [], []
    
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    s.append(float(parts[0]))
                    x.append(float(parts[1]))
                    y.append(float(parts[2]))
                    ue_vinf.append(float(parts[3]))
                except ValueError:
                    continue
    
    if not s:
        return None
    
    # Calculate Cp from gamma (ue/vinf)
    cp = [1.0 - g*g for g in ue_vinf]
    
    return {
        "s": s,
        "x": x, 
        "y": y,
        "gamma": ue_vinf,
        "cp": cp
    }

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    all_data = {}
    
    for foil in FOILS:
        print(f"\nProcessing {foil}...")
        data = run_xfoil_for_foil(foil, ALPHAS)
        if data:
            all_data[foil] = data
    
    # Save combined JSON
    output_file = OUTPUT_DIR / "xfoil_reference.json"
    with open(output_file, "w") as f:
        json.dump(all_data, f, indent=2)
    
    print(f"\nSaved to {output_file}")
    print(f"Foils: {list(all_data.keys())}")
    for foil, data in all_data.items():
        print(f"  {foil}: {list(data['alphas'].keys())} degrees")

if __name__ == "__main__":
    main()
