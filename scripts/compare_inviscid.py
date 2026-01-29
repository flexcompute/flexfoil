#!/usr/bin/env python3
"""
Compare inviscid panel method between XFOIL and RustFoil.

This script:
1. Runs RustFoil with debug enabled to dump inviscid solution
2. Loads the same geometry and computes gamma from XFOIL-style formula
3. Compares key quantities to find the source of divergence
"""

import json
import math
import subprocess
import sys
from pathlib import Path

def run_rustfoil_debug(airfoil_path: str, alpha_deg: float) -> dict:
    """Run RustFoil with debug enabled and return parsed debug data."""
    # Run rustfoil with debug
    debug_file = Path("/tmp/rustfoil_inviscid_debug.json")
    env = {"RUSTFOIL_DEBUG": str(debug_file), "RUSTFOIL_CL_DEBUG": "1"}
    
    cmd = [
        "cargo", "run", "--release", "-p", "rustfoil-cli", "--",
        "analyze", airfoil_path, "--alpha", str(alpha_deg)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, env={**subprocess.os.environ, **env})
    
    # Parse output for CL
    cl = None
    for line in result.stdout.split('\n'):
        if 'Cl =' in line:
            try:
                cl = float(line.split('=')[1].strip())
            except:
                pass
    
    # Try to load debug file
    debug_data = {}
    if debug_file.exists():
        with open(debug_file) as f:
            debug_data = json.load(f)
    
    return {
        "cl": cl,
        "debug": debug_data,
        "stderr": result.stderr
    }

def compute_theoretical_cl(alpha_deg: float) -> float:
    """Thin airfoil theory: CL = 2*pi*alpha for symmetric airfoil."""
    alpha_rad = math.radians(alpha_deg)
    return 2 * math.pi * alpha_rad

def load_xfoil_reference(alpha_deg: float) -> dict:
    """Load XFOIL reference data."""
    # XFOIL reference values at alpha=4 for NACA 0012 Re=3e6
    # From testdata/naca0012_re3m_xfoil.txt
    xfoil_data = {
        0: {"cl": 0.0000, "cd": 0.00598},
        4: {"cl": 0.4376, "cd": 0.00665},
        8: {"cl": 0.8584, "cd": 0.00986},
    }
    
    return xfoil_data.get(int(alpha_deg), {"cl": None, "cd": None})

def main():
    airfoil = "testdata/naca0012.dat"
    alpha = 4.0
    
    print(f"=== Inviscid Panel Method Comparison ===")
    print(f"Airfoil: {airfoil}")
    print(f"Alpha: {alpha}°")
    print()
    
    # Theoretical
    cl_theory = compute_theoretical_cl(alpha)
    print(f"Thin airfoil theory: CL = {cl_theory:.4f}")
    
    # XFOIL reference
    xfoil_ref = load_xfoil_reference(alpha)
    if xfoil_ref["cl"]:
        print(f"XFOIL (viscous):     CL = {xfoil_ref['cl']:.4f} (error vs theory: {(xfoil_ref['cl']/cl_theory - 1)*100:.1f}%)")
    
    # Run RustFoil
    print()
    print("Running RustFoil...")
    rustfoil = run_rustfoil_debug(airfoil, alpha)
    
    if rustfoil["cl"]:
        print(f"RustFoil (inviscid): CL = {rustfoil['cl']:.4f} (error vs theory: {(rustfoil['cl']/cl_theory - 1)*100:.1f}%)")
        if xfoil_ref["cl"]:
            print(f"RustFoil error vs XFOIL: {(rustfoil['cl']/xfoil_ref['cl'] - 1)*100:.1f}%")
    
    # Check debug data
    if rustfoil["debug"]:
        print()
        print("=== Debug Data ===")
        events = rustfoil["debug"].get("events", [])
        
        # Look for inviscid solution
        for event in events:
            sub = event.get("subroutine", "")
            if "AIC" in sub or "INVISC" in sub:
                print(f"Found: {sub}")
                if "gamu_0" in event:
                    gamu_0 = event["gamu_0"][:5]
                    print(f"  gamu_0 (first 5): {gamu_0}")
                if "gamu_90" in event:
                    gamu_90 = event["gamu_90"][:5]
                    print(f"  gamu_90 (first 5): {gamu_90}")
    
    # Additional analysis
    print()
    print("=== Analysis ===")
    print("The 10% CL over-prediction suggests:")
    print("1. Panel influence coefficients (PSILIN) might differ")
    print("2. Kutta condition implementation might differ")
    print("3. Panel geometry (normals, positions) might differ")
    print()
    print("Next steps:")
    print("1. Compare AIJ matrix elements between XFOIL and RustFoil")
    print("2. Compare GAMU base solutions")
    print("3. Check panel normal calculations")

if __name__ == "__main__":
    main()
