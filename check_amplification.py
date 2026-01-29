#!/usr/bin/env python3
"""Check amplification rate calculation"""
import math

print("Testing amplification rate critical Rθ calculation:")
print("="*60)

# These are typical values near the leading edge where transition might occur
test_cases = [
    {"hk": 2.6, "theta": 1.84e-5, "u": 1.46, "x": 0.0310, "desc": "RustFoil ibl=10"},
    {"hk": 2.6, "theta": 2.69e-5, "u": 2.68, "x": 0.0511, "desc": "RustFoil ibl=20"},
    {"hk": 2.5, "theta": 2e-5, "u": 1.5, "x": 0.012, "desc": "Typical near LE"},
]

re = 20_000_000

for tc in test_cases:
    hk = tc['hk']
    theta = tc['theta']
    u = tc['u']
    x = tc['x']
    desc = tc['desc']
    rt = re * u * theta
    
    print(f"\n{desc} (x={x:.4f}):")
    print(f"  Hk={hk:.2f}, θ={theta:.2e}, U={u:.2f}")
    print(f"  Rθ = Re*U*θ = {rt:.0f}")
    
    # Manually compute critical Rθ for comparison
    hmi = 1.0 / (hk - 1.0)
    aa = 2.492 * (hmi ** 0.43)
    bb = math.tanh(14.0 * hmi - 9.24)
    grcrit = aa + 0.7 * (bb + 1.0)
    rcrit = 10 ** grcrit
    
    gr = math.log10(rt) if rt > 0 else 0
    DGR = 0.08
    
    print(f"  log10(Rθ)={gr:.3f}, log10(Rcrit)={grcrit:.3f}")
    print(f"  Rcrit={rcrit:.0f}")
    
    if rt < rcrit * (10 ** (-DGR)):  # Below DGR threshold
        print(f"  ❌ SUBCRITICAL: Rθ={rt:.0f} < {rcrit * (10 ** (-DGR)):.0f}, amplification = 0")
    elif rt < rcrit * (10 ** DGR):
        frac = (gr - (grcrit - DGR)) / (2 * DGR)
        rfac = 3 * frac**2 - 2 * frac**3
        print(f"  ⚠️  IN RAMP: Rθ in smooth ramp region, RFAC={rfac:.3f}")
    else:
        print(f"  ✓ SUPERCRITICAL: Rθ={rt:.0f} > {rcrit * (10 ** DGR):.0f}, full amplification")

print("\n" + "="*60)
print("If Rθ is subcritical or barely in the ramp, amplification will be")
print("very low and transition will be delayed or not happen at all.")
