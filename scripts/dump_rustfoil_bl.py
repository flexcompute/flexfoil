#!/usr/bin/env python3
"""Dump RustFoil BL values for comparison."""
import subprocess
import re

# Run RustFoil test and capture output
result = subprocess.run(
    ['cargo', 'test', '--release', '-p', 'rustfoil-solver', 
     '--test', 'xfoil_viscous_comparison', 'test_viscous_at_alpha', 
     '--', '--nocapture'],
    capture_output=True, text=True, 
    env={'RUSTFOIL_CL_DEBUG': '1', **__import__('os').environ}
)

output = result.stderr + result.stdout

# Parse Upper/Lower station data
print("=== RustFoil BL State (Initial, before Newton) ===")
print("Station | Surface |    x/c     |    theta     |    dstar     |     Ue      |  mass")
print("-" * 90)

for line in output.split('\n'):
    if 'Upper[' in line or 'Lower[' in line:
        # Parse: Upper[0]: x=2.039e-10, theta=6.77e-8, dstar=1.49e-7, u=0.01, mass=1.49e-9
        m = re.search(r'(Upper|Lower)\[(\d+)\]: x=([0-9.e+-]+), theta=([0-9.e+-]+), dstar=([0-9.e+-]+), u=([0-9.e+-]+), mass=([0-9.e+-]+)', line)
        if m:
            surf = m.group(1)
            idx = int(m.group(2))
            x = float(m.group(3))
            theta = float(m.group(4))
            dstar = float(m.group(5))
            ue = float(m.group(6))
            mass = float(m.group(7))
            print(f"  {idx:3d}   | {surf:6s} | {x:10.6f} | {theta:12.6e} | {dstar:12.6e} | {ue:10.6f} | {mass:12.6e}")

# Also get iteration 5 values
print("\n=== RustFoil BL State (iter 5) ===")
for line in output.split('\n'):
    if 'iter 5' in line and ('upper[' in line or 'lower[' in line):
        print(line)
