# ctau Transition Issue Analysis

## Summary

The boundary layer ctau (shear stress coefficient) values are near-zero at turbulent stations immediately after transition, when they should be ~0.03-0.06. This causes:
- Slow H (shape factor) evolution after transition (takes ~0.15c instead of ~0.05c)
- Incorrect boundary layer profiles in the transition region
- Potential CD underprediction at higher angles of attack

## Investigation Timeline

### What Works Correctly

1. **TRDIF initialization:** The `trdif` and `trdif_full` functions correctly compute initial ctau at transition using:
   ```
   CTR = 1.8 * exp(-3.3/(Hk-1))
   ctau_init = CTR * CQ
   ```
   At Hk≈3.2, this gives CTR≈0.40 and ctau_init≈0.03-0.04.

2. **Single-station Newton march:** The `march_surface` function correctly evolves ctau through the transition:
   - Station 38 (transition): ctau converges to 0.016
   - Station 39: ctau → 0.040
   - Station 40-45: ctau stabilizes around 0.045-0.050

3. **is_laminar flag:** Correctly set to false after transition, ensuring subsequent stations use turbulent equations.

### What Goes Wrong

The **global Newton iterations** (30 iterations in `solve_viscous_two_surfaces`) overwrite the correctly-computed ctau values:

1. After `march_surface` completes, ctau values are reasonable (~0.04-0.05)
2. Global Newton updates stations via `apply_global_updates`:
   ```rust
   station.ctau += rlx * delta[0];
   station.ctau = station.ctau.clamp(1e-6, 0.25);
   ```
3. The `delta[0]` values are large negative, driving ctau to the minimum (1e-6)
4. Final BL_DEBUG output shows ctau ≈ 0.0000 at transition stations

## Root Cause Hypothesis

The shear-lag equation's residual/Jacobian may have an issue in the global Newton context:

1. **Residual sign:** `res.res_third = -REZC` (negative of XFOIL convention)
2. **Jacobian:** Uses `d(REZC)/d(ctau)` but should use `d(-REZC)/d(ctau) = -d(REZC)/d(ctau)`

If the Jacobian sign is wrong relative to the residual, Newton will drive ctau in the wrong direction when `CQ > ctau` (which is the case near transition).

## Key Code Locations

1. **TRDIF equations:** `crates/rustfoil-bl/src/equations.rs:1880-2100`
2. **Single-station Newton:** `crates/rustfoil-coupling/src/march.rs:770-1500`
3. **Global Newton updates:** `crates/rustfoil-coupling/src/global_newton.rs:950-1060`
4. **Shear-lag equation:** `crates/rustfoil-bl/src/equations.rs:1158-1370`

## XFOIL Comparison

At α=4°, upper surface near transition (x≈0.24-0.34):

| x | XFOIL ctau | RustFoil ctau | XFOIL Hk | RustFoil Hk |
|---|------------|---------------|----------|-------------|
| 0.24 | ~0.05 | 0.000 | 2.7 | 3.0 |
| 0.27 | ~0.06 | 0.000 | 1.9 | 3.1 |
| 0.30 | ~0.05 | 0.004 | 1.5 | 3.1 |
| 0.35 | ~0.04 | 0.045 | 1.5 | 1.9 |

The Hk mismatch directly correlates with ctau mismatch - both quantities fail to evolve properly after transition.

## Next Steps

1. **Verify Jacobian signs** in the shear-lag equation for both single-station and global Newton
2. **Compare XFOIL debug output** for Z_S2 (∂res_ctau/∂ctau) with RustFoil values
3. **Consider alternative:** Skip global Newton updates to ctau at transition stations, or use smaller relaxation factor
