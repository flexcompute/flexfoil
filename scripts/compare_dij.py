#!/usr/bin/env python3
"""Direct comparison of DIJ values between XFOIL and RustFoil."""

import subprocess
import sys

# First, get XFOIL DIJ values by parsing the dump file
print("=== Parsing XFOIL DIJ dump ===")

dij_xfoil = {}
with open('/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin/xfoil_dij_dump.txt', 'r') as f:
    lines = f.readlines()
    
    for line in lines:
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit():
            # Diagonal element or single row element
            idx = int(parts[0])
            val = float(parts[1])
            # This is ambiguous - skip for now
        elif len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            # DIJ(i,j) = val
            i = int(parts[0])
            j = int(parts[1])
            val = float(parts[2])
            dij_xfoil[(i, j)] = val

print(f"Parsed {len(dij_xfoil)} DIJ elements from XFOIL dump")

# Key comparisons: look at row 83 (stagnation panel in XFOIL, 1-indexed)
print("\n=== XFOIL DIJ at row 83 (stagnation, 1-indexed) ===")
for j in [82, 83, 84, 85]:
    if (83, j) in dij_xfoil:
        print(f"  DIJ(83, {j}) = {dij_xfoil[(83, j)]:.6e}")

# Row 82 (the panel before stagnation going upper)
print("\n=== XFOIL DIJ at row 82 (1-indexed) ===")
for j in [81, 82, 83, 84]:
    if (82, j) in dij_xfoil:
        print(f"  DIJ(82, {j}) = {dij_xfoil[(82, j)]:.6e}")

# Row 80
print("\n=== XFOIL DIJ at row 80 (1-indexed) ===")
for j in [78, 79, 80, 81, 82]:
    if (80, j) in dij_xfoil:
        print(f"  DIJ(80, {j}) = {dij_xfoil[(80, j)]:.6e}")

print("\n=== RustFoil DIJ (0-indexed) ===")
print("From previous debug output:")
print("  panel 82: dij_diag = -624.5 = DIJ(82,82)")
print("  panel 81: dij_diag = -651.9 = DIJ(81,81)")
print("  panel 80: dij_diag = -674.1 = DIJ(80,80)")
print("  DIJ(82,81) = +202.6 (from trace)")

print("\n=== Comparison (XFOIL 1-indexed to RustFoil 0-indexed) ===")
print("XFOIL row 83 = RustFoil row 82")
print("XFOIL row 82 = RustFoil row 81")

# Get row 83 from XFOIL if we can
row83_elems = {j: dij_xfoil[(83, j)] for j in range(1, 90) if (83, j) in dij_xfoil}
if row83_elems:
    print(f"\nXFOIL row 83 elements: {row83_elems}")
