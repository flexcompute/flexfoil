#!/usr/bin/env python3
"""
Compare Newton system matrices between XFOIL and RustFoil.

Extracts VA, VB, VDEL from XFOIL debug JSON and compares with RustFoil output.
"""

import json
import sys
import subprocess
import re
from pathlib import Path

def extract_xfoil_setbl_data(json_path: Path) -> dict:
    """Extract SETBL entries from XFOIL debug JSON."""
    data = {"events": []}
    
    with open(json_path, 'r') as f:
        content = f.read()
        
    # Find all SETBL entries (they're not proper JSON arrays, need to parse)
    # Look for the pattern: "subroutine": "SETBL" entries
    setbl_pattern = re.compile(
        r'\{\s*"call_id":\s*\d+,\s*"subroutine":\s*"SETBL".*?"VM_sample":\s*\[[^\]]+\]\s*\}',
        re.DOTALL
    )
    
    matches = setbl_pattern.findall(content)
    
    setbl_entries = []
    for match in matches:
        try:
            # Clean up the JSON
            cleaned = match.strip()
            entry = json.loads(cleaned)
            setbl_entries.append(entry)
        except json.JSONDecodeError as e:
            # Try to fix common issues
            pass
    
    return setbl_entries


def parse_xfoil_debug_json_manual(json_path: Path) -> list:
    """Parse XFOIL debug JSON manually, extracting SETBL blocks."""
    with open(json_path, 'r') as f:
        content = f.read()
    
    setbl_entries = []
    lines = content.split('\n')
    
    i = 0
    while i < len(lines):
        line = lines[i]
        if '"subroutine": "SETBL"' in line:
            # Found SETBL entry, go back to find start of JSON object
            start = i
            while start > 0 and '{' not in lines[start]:
                start -= 1
            if '{' in lines[start]:
                # Find call_id line
                for j in range(start, i + 1):
                    if '"call_id"' in lines[j]:
                        start = j - 1  # Include the opening brace
                        break
            
            # Now find the end of the JSON object (look for VM_sample closing bracket)
            end = i + 1
            brace_count = 0
            started = False
            while end < len(lines):
                line_end = lines[end]
                if '{' in line_end:
                    brace_count += line_end.count('{')
                    started = True
                if '}' in line_end:
                    brace_count -= line_end.count('}')
                    if started and brace_count <= 0:
                        break
                end += 1
            
            # Extract and parse the JSON block
            block_lines = lines[start:end+1]
            block = '\n'.join(block_lines)
            
            # Clean up - remove leading/trailing junk
            block = block.strip()
            if block.startswith(','):
                block = block[1:]
            block = block.strip()
            
            # Find the actual JSON object boundaries
            first_brace = block.find('{')
            last_brace = block.rfind('}')
            if first_brace >= 0 and last_brace > first_brace:
                block = block[first_brace:last_brace+1]
            
            try:
                entry = json.loads(block)
                setbl_entries.append(entry)
            except json.JSONDecodeError:
                pass
            
            i = end
        i += 1
    
    return setbl_entries


def format_matrix(matrix, name, width=12):
    """Format a matrix for display."""
    lines = [f"{name}:"]
    for row in matrix:
        row_str = "  [" + ", ".join(f"{v:>{width}.6e}" for v in row) + "]"
        lines.append(row_str)
    return '\n'.join(lines)


def compare_entry(xfoil_entry: dict, station_info: str):
    """Display XFOIL matrix data for a station."""
    print(f"\n{'='*70}")
    print(f"XFOIL SETBL: {station_info}")
    print(f"  side={xfoil_entry.get('side')}, ibl={xfoil_entry.get('ibl')}, iv={xfoil_entry.get('iv')}, iteration={xfoil_entry.get('iteration')}")
    print(f"{'='*70}")
    
    # VA matrix (3x2)
    va = xfoil_entry.get('VA', [])
    print(format_matrix(va, "VA (3x2) - downstream station derivatives"))
    
    # VB matrix (3x2)
    vb = xfoil_entry.get('VB', [])
    print(format_matrix(vb, "\nVB (3x2) - upstream station derivatives"))
    
    # VDEL residuals (3x2, but only first column is residual)
    vdel = xfoil_entry.get('VDEL', [])
    print(format_matrix(vdel, "\nVDEL (residuals with forced changes)"))
    
    # Extract just the residual values
    if vdel:
        print(f"\n  VDEL residuals (first column):")
        print(f"    [0] res_ampl  = {vdel[0][0]:>12.6e}")
        print(f"    [1] res_mom   = {vdel[1][0]:>12.6e}")
        print(f"    [2] res_shape = {vdel[2][0]:>12.6e}")


def main():
    # Path to XFOIL debug JSON
    xfoil_json = Path("/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin/xfoil_debug.json")
    
    if not xfoil_json.exists():
        print(f"XFOIL debug JSON not found at {xfoil_json}")
        print("Run XFOIL first to generate debug output:")
        print("""
cd /Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin
./xfoil_instrumented << 'EOF'
naca 0012
oper
visc 3000000
a 4

quit
EOF
""")
        return 1
    
    print("="*70)
    print("XFOIL Newton System Matrix Comparison")
    print("="*70)
    print(f"\nReading XFOIL debug data from: {xfoil_json}")
    
    # Parse SETBL entries
    setbl_entries = parse_xfoil_debug_json_manual(xfoil_json)
    
    print(f"Found {len(setbl_entries)} SETBL entries")
    
    # Filter for iteration=1, side=1 (upper surface), ibl=2,3,4,5
    upper_iter1 = [e for e in setbl_entries 
                   if e.get('iteration') == 1 and e.get('side') == 1]
    
    print(f"\nUpper surface (side=1), iteration=1 entries: {len(upper_iter1)}")
    
    # Display first few stations
    print("\n" + "="*70)
    print("XFOIL Matrix Values at Key Stations (Upper Surface, Iteration 1)")
    print("="*70)
    
    stations_of_interest = [2, 3, 4, 5, 10, 20]
    for target_ibl in stations_of_interest:
        entries = [e for e in upper_iter1 if e.get('ibl') == target_ibl]
        if entries:
            entry = entries[0]
            compare_entry(entry, f"IBL={target_ibl}, IS=1 (upper)")
    
    # Also show lower surface station 3 for comparison
    lower_iter1 = [e for e in setbl_entries 
                   if e.get('iteration') == 1 and e.get('side') == 2]
    
    print("\n" + "="*70)
    print("XFOIL Matrix Values at Key Stations (Lower Surface, Iteration 1)")
    print("="*70)
    
    for target_ibl in [2, 3]:
        entries = [e for e in lower_iter1 if e.get('ibl') == target_ibl]
        if entries:
            entry = entries[0]
            compare_entry(entry, f"IBL={target_ibl}, IS=2 (lower)")
    
    # Print summary table
    print("\n" + "="*70)
    print("SUMMARY TABLE: XFOIL VDEL Residuals (Iteration 1)")
    print("="*70)
    print(f"{'Surface':<8} {'IBL':>4} {'IV':>4} {'res_ampl':>14} {'res_mom':>14} {'res_shape':>14}")
    print("-"*70)
    
    for entry in sorted(setbl_entries, key=lambda e: (e.get('iteration', 0), e.get('side', 0), e.get('ibl', 0))):
        if entry.get('iteration') == 1 and entry.get('ibl', 0) <= 10:
            side = "upper" if entry.get('side') == 1 else "lower"
            ibl = entry.get('ibl', 0)
            iv = entry.get('iv', 0)
            vdel = entry.get('VDEL', [[0,0], [0,0], [0,0]])
            print(f"{side:<8} {ibl:>4} {iv:>4} {vdel[0][0]:>14.6e} {vdel[1][0]:>14.6e} {vdel[2][0]:>14.6e}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
