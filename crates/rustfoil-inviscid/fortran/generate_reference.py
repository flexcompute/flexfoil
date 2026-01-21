#!/usr/bin/env python3
"""
Generate inviscid solver reference data by running Xfoil-instrumented.

This script automates XFOIL to:
1. Load a NACA 0012 airfoil
2. Set up panels
3. Run inviscid analysis at multiple angles
4. Extract the debug JSON output

Usage:
    python generate_reference.py [--panels 160] [--output testdata/inviscid_ref.json]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def run_xfoil_instrumented(panels: int = 160, alphas: list = None) -> dict:
    """Run XFOIL-instrumented and capture debug output."""
    
    if alphas is None:
        alphas = [0.0, 4.0, 8.0]
    
    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent.parent
    xfoil_bin = project_root / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
    
    if not xfoil_bin.exists():
        raise FileNotFoundError(f"XFOIL binary not found: {xfoil_bin}")
    
    # Create input commands for XFOIL
    commands = []
    
    # Load NACA 0012
    commands.append("NACA 0012")
    commands.append("")  # Accept default
    
    # Set panel count (PPAR command)
    commands.append("PPAR")
    commands.append(f"N {panels}")
    commands.append("")
    commands.append("")
    
    # Enter OPER mode
    commands.append("OPER")
    
    # Run at each alpha
    for alpha in alphas:
        commands.append(f"ALFA {alpha}")
    
    # Exit
    commands.append("")
    commands.append("QUIT")
    
    input_text = "\n".join(commands)
    
    # Create temp directory for output
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        debug_file = tmpdir / "xfoil_debug.json"
        
        # Run XFOIL
        result = subprocess.run(
            [str(xfoil_bin)],
            input=input_text,
            capture_output=True,
            text=True,
            cwd=str(tmpdir),
            timeout=30
        )
        
        # Check for debug output
        if debug_file.exists():
            with open(debug_file) as f:
                content = f.read()
                # Parse JSON array
                if content.strip():
                    # Wrap in array if needed
                    if not content.strip().startswith('['):
                        content = '[' + content + ']'
                    # Fix trailing commas
                    content = content.replace(',]', ']').replace(',}', '}')
                    try:
                        return json.loads(content)
                    except json.JSONDecodeError as e:
                        print(f"JSON parse error: {e}", file=sys.stderr)
                        print(f"Content preview: {content[:500]}", file=sys.stderr)
                        return {}
        
        # If no debug file, try to extract data from stdout
        print("No debug file generated, using manual extraction", file=sys.stderr)
        return extract_from_stdout(result.stdout, panels, alphas)


def extract_from_stdout(stdout: str, panels: int, alphas: list) -> dict:
    """Extract reference data manually from XFOIL output."""
    
    # This is a fallback - in practice we'd parse XFOIL's verbose output
    # For now, return empty structure that indicates manual generation needed
    return {
        "error": "Could not extract data automatically",
        "stdout_preview": stdout[:1000] if stdout else "No output"
    }


def generate_reference_data(panels: int = 160) -> dict:
    """Generate complete reference data structure."""
    
    alphas = [0.0, 4.0, 8.0]
    
    # Run XFOIL
    xfoil_data = run_xfoil_instrumented(panels, alphas)
    
    # Structure the output
    reference = {
        "test": "inviscid_solver",
        "n_panels": panels,
        "airfoil": "NACA0012",
        "xfoil_debug": xfoil_data,
        "data": []
    }
    
    # Extract specific data sections from XFOIL debug output
    if isinstance(xfoil_data, list):
        for item in xfoil_data:
            if isinstance(item, dict):
                sub = item.get("subroutine", "")
                if sub in ["NCALC", "APCALC", "TECALC", "PSILIN", "GGCALC", "SPECAL", "CLCALC", "STFIND"]:
                    reference["data"].append(item)
    
    return reference


def main():
    parser = argparse.ArgumentParser(description="Generate inviscid solver reference data")
    parser.add_argument("--panels", type=int, default=160, help="Number of panels")
    parser.add_argument("--output", type=str, default=None, 
                        help="Output file (default: testdata/inviscid_ref.json)")
    args = parser.parse_args()
    
    # Generate reference data
    print(f"Generating reference data with {args.panels} panels...", file=sys.stderr)
    ref_data = generate_reference_data(args.panels)
    
    # Output
    if args.output:
        output_path = Path(args.output)
    else:
        script_dir = Path(__file__).parent
        output_path = script_dir.parent.parent.parent / "testdata" / "inviscid_ref.json"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(ref_data, f, indent=2)
    
    print(f"Reference data written to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
