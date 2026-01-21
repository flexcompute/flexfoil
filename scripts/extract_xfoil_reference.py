#!/usr/bin/env python3
"""
Extract inviscid solver reference data from XFOIL-instrumented.

Runs the actual FORTRAN XFOIL binary and parses its JSON debug output
to create reference data for validating the Rust implementation.

Usage:
    python extract_xfoil_reference.py [--output testdata/inviscid_ref.json]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def run_xfoil(xfoil_bin: Path, commands: str, workdir: Path) -> tuple:
    """Run XFOIL with given commands and return stdout, stderr."""
    result = subprocess.run(
        [str(xfoil_bin)],
        input=commands,
        capture_output=True,
        text=True,
        cwd=str(workdir),
        timeout=30
    )
    return result.stdout, result.stderr, result.returncode


def extract_reference_data(panels: int = 160) -> dict:
    """Run XFOIL-instrumented and extract reference data."""
    
    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    xfoil_bin = project_root / "Xfoil-instrumented" / "bin" / "xfoil_instrumented"
    
    if not xfoil_bin.exists():
        raise FileNotFoundError(f"XFOIL binary not found: {xfoil_bin}")
    
    # XFOIL commands for inviscid analysis
    # Note: XFOIL uses PPAR in GDES menu to change panel count
    commands = f"""NACA 0012
PPAR
N {panels}


OPER
ALFA 0
ALFA 4
ALFA 8

QUIT
"""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        print(f"Running XFOIL-instrumented...", file=sys.stderr)
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
            
            # Parse JSON - handle the comma-separated format
            if content.strip():
                # Wrap in array brackets if needed
                content = content.strip()
                if not content.startswith('['):
                    content = '[' + content + ']'
                # Fix trailing commas
                content = content.replace(',\n]', '\n]')
                content = content.replace(',]', ']')
                
                try:
                    events = json.loads(content)
                    return process_xfoil_events(events, panels)
                except json.JSONDecodeError as e:
                    print(f"JSON parse error: {e}", file=sys.stderr)
                    # Save raw content for debugging
                    with open("/tmp/xfoil_debug_raw.json", "w") as f:
                        f.write(content)
                    print(f"Raw content saved to /tmp/xfoil_debug_raw.json", file=sys.stderr)
                    return {"error": str(e), "raw_content": content[:1000]}
        else:
            print("No debug file found", file=sys.stderr)
            print(f"XFOIL stdout: {stdout[:1000]}", file=sys.stderr)
            
        return {"error": "No debug output generated"}


def process_xfoil_events(events: list, panels: int) -> dict:
    """Process XFOIL debug events into structured reference data."""
    
    reference = {
        "test": "inviscid_solver",
        "source": "xfoil_instrumented",
        "n_panels": panels,
        "airfoil": "NACA0012",
        "data": []
    }
    
    # Extract relevant events by subroutine
    for event in events:
        if not isinstance(event, dict):
            continue
            
        sub = event.get("subroutine", "")
        
        if sub in ["GGCALC", "SPECAL", "CLCALC", "NCALC", "APCALC", "TECALC", "STFIND"]:
            # Add type field for easier parsing
            event["type"] = sub.lower()
            reference["data"].append(event)
    
    return reference


def main():
    parser = argparse.ArgumentParser(description="Extract XFOIL reference data")
    parser.add_argument("--panels", type=int, default=160, help="Number of panels")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    
    ref_data = extract_reference_data(args.panels)
    
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path(__file__).parent.parent / "testdata" / "inviscid_ref.json"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(ref_data, f, indent=2)
    
    print(f"Reference data written to {output_path}", file=sys.stderr)
    
    # Print summary
    if "data" in ref_data:
        print(f"\nExtracted {len(ref_data['data'])} debug events:", file=sys.stderr)
        subs = {}
        for item in ref_data["data"]:
            sub = item.get("subroutine", "unknown")
            subs[sub] = subs.get(sub, 0) + 1
        for sub, count in sorted(subs.items()):
            print(f"  {sub}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
