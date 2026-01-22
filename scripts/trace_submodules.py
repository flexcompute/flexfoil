#!/usr/bin/env python3
"""
Trace through all submodules at the divergence point to identify where RustFoil diverges from XFOIL.
"""

import json
import sys

def load_xfoil_data():
    with open("testdata/mrchue_iterations.json") as f:
        return json.load(f)

def main():
    data = load_xfoil_data()
    
    # Get upper surface (side 1) stations
    side1 = data["sides"]["1"]["stations"]
    
    print("=" * 80)
    print("XFOIL Station-by-Station Data (Upper Surface)")
    print("=" * 80)
    print()
    
    # Find stations around the divergence (indices 27-31 based on test output)
    print(f"{'IBL':>4} {'x':>10} {'Ue':>10} {'theta':>12} {'delta*':>12} {'H':>8} {'Hk':>8} {'ampl':>10}")
    print("-" * 90)
    
    for station in side1:
        ibl = station["ibl"]
        x = station["x"]
        ue = station["Ue"]
        final = station["final"]
        theta = final["theta"]
        delta_star = final["delta_star"]
        hk = final.get("Hk", 0)
        ampl = final["ampl"]
        h = delta_star / theta if theta > 0 else 0
        
        # Focus on stations around transition
        if 25 <= ibl <= 35:
            print(f"{ibl:4d} {x:10.6f} {ue:10.6f} {theta:12.6e} {delta_star:12.6e} {h:8.4f} {hk:8.4f} {ampl:10.4f}")
    
    print()
    print("=" * 80)
    print("Detailed Analysis at Station 29 (IBL=31, x~0.1276)")
    print("=" * 80)
    
    # Find station 31 (IBL=31 corresponds to test index 29)
    station_31 = None
    station_30 = None
    for s in side1:
        if s["ibl"] == 31:
            station_31 = s
        if s["ibl"] == 30:
            station_30 = s
    
    if station_31:
        print()
        print("Station 30 (previous, IBL=30):")
        if station_30:
            print(f"  x = {station_30['x']:.6f}")
            print(f"  Ue = {station_30['Ue']:.6f}")
            print(f"  Final theta = {station_30['final']['theta']:.6e}")
            print(f"  Final delta* = {station_30['final']['delta_star']:.6e}")
            print(f"  Final Hk = {station_30['final'].get('Hk', 'N/A')}")
            print(f"  Final ampl = {station_30['final']['ampl']:.4f}")
        
        print()
        print("Station 31 (divergence point, IBL=31):")
        print(f"  x = {station_31['x']:.6f}")
        print(f"  Ue = {station_31['Ue']:.6f}")
        print()
        print("  Initial state (from previous station):")
        print(f"    theta = {station_31['initial']['theta']:.6e}")
        print(f"    delta* = {station_31['initial']['delta_star']:.6e}")
        print(f"    H = {station_31['initial']['delta_star']/station_31['initial']['theta']:.4f}")
        print(f"    ampl = {station_31['initial']['ampl']:.4f}")
        print()
        print("  Final state (after Newton convergence):")
        print(f"    theta = {station_31['final']['theta']:.6e}")
        print(f"    delta* = {station_31['final']['delta_star']:.6e}")
        print(f"    Hk = {station_31['final'].get('Hk', 'N/A')}")
        print(f"    ampl = {station_31['final']['ampl']:.4f}")
        print()
        print(f"  Newton iterations: {station_31['n_iterations']}")
        print(f"  Converged: {station_31['converged']}")
        
        # Show VS2 from first iteration
        if station_31.get("iterations"):
            iter1 = station_31["iterations"][0]
            print()
            print("  First iteration VS2 matrix:")
            vs2 = iter1["VS2"]
            print("    Row 0 (ampl eq):  ", [f"{v:12.4f}" for v in vs2[0]])
            print("    Row 1 (mom eq):   ", [f"{v:12.4f}" for v in vs2[1]])
            print("    Row 2 (shape eq): ", [f"{v:12.4f}" for v in vs2[2]])
            print("    Row 3 (4th eq):   ", [f"{v:12.4f}" for v in vs2[3]])
            print()
            print("  Key Jacobian entries:")
            print(f"    VS2[1][1] (∂mom/∂θ) = {vs2[1][1]:.2f}")
            print(f"    VS2[2][1] (∂shape/∂θ) = {vs2[2][1]:.2f}")
            print(f"    VS2[2][2] (∂shape/∂δ*) = {vs2[2][2]:.2f}")
            print()
            print("  Note: Row 3 has non-zero entries - XFOIL is using INVERSE mode!")
            print("  This means Hk (or separation) is being constrained.")
    
    print()
    print("=" * 80)
    print("KEY OBSERVATION")
    print("=" * 80)
    print()
    print("At station 31 (x=0.1276), XFOIL:")
    print("  - Uses INVERSE mode (VS2 row 3 has non-zero entries)")
    print("  - Hk grows to 5.39 (approaching separation)")
    print("  - Delta* increases 60% in one station")
    print()
    print("RustFoil likely:")
    print("  - Uses DIRECT mode (or wrong inverse mode formulation)")
    print("  - Fails to handle the strong APG properly")
    print("  - Produces collapsed Hk=1.05")
    print()
    print("The issue is likely in the DIRECT vs INVERSE mode logic,")
    print("not in the Jacobian calculation itself.")

if __name__ == "__main__":
    main()
