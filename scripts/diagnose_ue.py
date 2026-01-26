#!/usr/bin/env python3
"""Diagnose the Ue extraction issue."""

import subprocess
import json

# Run RustFoil at single alpha with verbose output
result = subprocess.run(
    ["./target/release/rustfoil", "viscous", "testdata/naca0012.dat", "-a", "4", "-r", "1e6"],
    capture_output=True,
    text=True
)

print("=== RustFoil viscous output ===")
print(result.stdout)
print(result.stderr)

# Compare with inviscid
print("\n=== RustFoil inviscid polar ===")
result2 = subprocess.run(
    ["./target/release/rustfoil", "polar", "testdata/naca0012.dat", 
     "--alpha-start=4", "--alpha-end=4", "--alpha-step=1"],
    capture_output=True,
    text=True
)
print(result2.stdout)
print(result2.stderr)
