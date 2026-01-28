# Boundary Layer Re-March Fix for Stall Prediction

## Problem
RustFoil was unable to predict stall because it marched the boundary layer ONCE before the Newton loop, then used only Jacobian updates for 20 global iterations. XFOIL re-marches the BL at EVERY iteration, allowing H (shape parameter) to exceed htmax nonlinearly.

## Root Cause
The previous implementation in `viscal.rs`:
1. Marched BL once at line 627 (before Newton loop)
2. Copied results to stations at lines 664-708
3. Entered Newton loop at line 809
4. Used only Jacobian updates (blvar + build_global_system + solve)

This prevented nonlinear behavior like H growing beyond htmax at high angles of attack, which is essential for stall prediction.

## Solution
Added boundary layer re-marching INSIDE the Newton loop at line 827 (right after loop start):

```rust
for iter in 0..max_newton_iter {
    // === RE-MARCH BOUNDARY LAYERS WITH CURRENT Ue ===
    // Extract current Ue from stations (updated by previous iteration)
    let current_upper_ue: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
    let current_lower_ue: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
    
    // Re-march with current Ue
    let upper_result_iter = march_surface(&upper_arc, &current_upper_ue, re, msq, &march_config, 1);
    let lower_result_iter = march_surface(&lower_arc, &current_lower_ue, re, msq, &march_config, 2);
    
    // Copy results back to stations (theta, delta_star, h, hk, cf, ctau, etc.)
    // NOTE: Station.u is NOT overwritten - it comes from VI coupling
    
    // Update transition locations
    if let Some(xtr) = upper_result_iter.x_transition {
        x_tr_upper = xtr;
    }
    // ... rest of Newton loop (blvar, Jacobian, solve)
}
```

## Key Implementation Details

1. **Extract Current Ue**: Use `upper_stations[i].u` and `lower_stations[i].u` which were updated by the previous Newton iteration's VI coupling

2. **Re-march Both Surfaces**: Call `march_surface()` with the CURRENT Ue values, not the initial inviscid Ue

3. **Copy BL Variables**: Update theta, delta_star, h, hk, cf, ctau, ampl, flow type flags, mass_defect, r_theta

4. **Preserve Edge Velocity**: Do NOT overwrite `station.u` - it must come from VI coupling (DIJ matrix + Newton solve)

5. **Update Transition**: Copy transition locations from march results

6. **Offset Handling**: Reuse `upper_offset` and `lower_offset` to handle stagnation station properly

## What This Enables

At high angles of attack (α > 15°):
- **Before Fix**: H stays bounded by Jacobian linearization around htmax
- **After Fix**: H can grow nonlinearly (e.g., from 2.4 → 3.0 at TE) over iterations
- **Result**: CL stops increasing, CD rises sharply → STALL PREDICTION

## Testing

Run with debug output to verify:
```bash
RUSTFOIL_CL_DEBUG=1 cargo run --release --bin rustfoil-cli -- \
    --airfoil testdata/naca2412.dat \
    --reynolds 1e6 \
    --alpha-start 0 \
    --alpha-end 18 \
    --alpha-step 1
```

Expected output:
- `[DEBUG Newton] Re-marching BL at iteration N` appears for each iteration
- MRCHUE events appear MULTIPLE times per station (once per iteration)
- At α=15°, ibl=93 (TE) shows H growing from ~2.4 → ~3.0 over iterations
- CL reaches peak around α=15-16°, then decreases (stall)

## Files Modified
- `crates/rustfoil-solver/src/viscous/viscal.rs` (lines 827-910)

## References
- XFOIL's VISCAL: `xoper.f` line 2886-3163
- XFOIL's MRCHUE: Re-marches BL inside Newton loop at every iteration
- Issue: RustFoil CD overprediction at low alpha, inability to predict stall
