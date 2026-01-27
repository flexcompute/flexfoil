# XFOIL Matching Analysis - Full Test Matrix

## Test Matrix

- **Foils**: NACA 0012, NACA 2412, NACA 4412
- **Reynolds Numbers**: 1e6, 3e6, 6e6 (NACA 0012); 3e6 (others)
- **Alpha Range**: -15° to +15° in 1° steps
- **Total Cases**: 155

## Summary of Results

### CL Error by Case

| Foil | Re | Max ΔCL | Mean ΔCL | Max %CL |
|------|------|---------|----------|---------|
| NACA 0012 | 1e6 | 1.27 | 0.19 | 88.9% |
| NACA 0012 | 3e6 | 0.56 | 0.09 | 36.0% |
| NACA 0012 | 6e6 | 0.26 | 0.05 | 29.8% |
| NACA 2412 | 3e6 | 0.42 | 0.08 | 29.8% |
| NACA 4412 | 3e6 | 0.74 | 0.14 | 180.1% |

### Key Observations

1. **Re dependence**: Lower Re shows larger errors (1e6 worst)
2. **Camber dependence**: Cambered foils show larger errors at low alpha
3. **Stall region**: Extreme alphas (±15°) show largest errors
4. **Best match region**: Mid-range alphas (2-8°) at higher Re

## Root Cause Analysis

### Verified Matching

1. **Panel Geometry**: ✅ Perfect match (< 1e-6 difference)
2. **Inviscid Solution**: ✅ Edge velocity Ue matches exactly
3. **BL Equations**: ✅ BLDIF formulas match XFOIL
4. **Cf Closure**: ✅ cf_laminar() matches CFL exactly

### Identified Divergence Point

**Momentum thickness (θ) evolution is the root cause.**

Station-by-station comparison at α=4°, Re=3e6, NACA 0012:

| Station | θ_XFOIL | θ_RustFoil | Error |
|---------|---------|------------|-------|
| 5 | 1.9e-5 | 1.9e-5 | ~0% |
| 10 | 1.9e-5 | 2.0e-5 | +5% |
| 15 | 2.8e-5 | 3.1e-5 | +11% |
| 20 | 3.9e-5 | 4.3e-5 | +10% |

**Cascade of errors:**
1. θ grows ~10% faster in RustFoil
2. RT (Reynolds θ) = ρUθ/μ → 10% higher
3. Cf ∝ 1/RT → 15% lower (amplified by Hk effect)
4. CD integration → significant underestimation

### Suspected Root Causes

1. **Momentum integral equation**: Different dissipation terms (Di)
2. **Shape parameter (H) evolution**: H grows slightly faster in RustFoil
3. **Transition handling**: May trigger at different locations
4. **Wake evolution**: Different momentum conservation model

## Files Generated

- `traces/xfoil/` - 155 instrumented XFOIL traces
- `traces/rustfoil/` - 155 RustFoil traces with debug output
- `comparison_results/comparison_summary.json` - Full comparison data
- `scripts/verify_paneling_parity.py` - Paneling verification
- `scripts/compare_traces.py` - Trace comparison
- `scripts/compare_bl_stations.py` - Station-by-station BL comparison

## Recommended Next Steps

1. **Debug momentum equation**: Compare Di (dissipation) values at each station
2. **Compare BLDIF Jacobians**: VS1, VS2 matrices at early stations
3. **Check transition location**: Verify transition triggers at same x/c
4. **Compare wake evolution**: θ conservation through wake stations

## Reference Data

Test traces available at:
- XFOIL: `traces/xfoil/{foil}/re{Re}/alpha_{alpha}.json`
- RustFoil: `traces/rustfoil/{foil}/re{Re}/alpha_{alpha}.json`

Each trace contains:
- BLVAR events: BL variables (H, Hk, Cf, θ, δ*, etc.)
- BLDIF events: Jacobian matrices (VS1, VS2, VSREZ)
- MRCHUE events: Marching update details
- TRCHEK events: Transition check data
