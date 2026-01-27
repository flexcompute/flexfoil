# Newton Iteration Explosion Fix

## Problem

RustFoil's Newton iteration was exploding at iteration 7, with residual jumping from ~0.5 to ~2.5. This indicates that the Newton updates were becoming unstable and causing large changes in the boundary layer variables.

## Root Cause

The relaxation safeguards in `apply_global_updates()` were **only checking normalized changes in theta and ctau/ampl**, but **NOT checking changes in Ue or delta_star**. 

When Newton produced large `delta_mass` values, this would cause:
1. Large `due` (edge velocity change) via DIJ coupling: `due = new_u - current_u`
2. Large `d_dstar` (displacement thickness change): `d_dstar = (delta_mass - delta_star * due) / ue`
3. These large changes would be applied without relaxation limiting, causing:
   - Ue to change dramatically
   - delta_star to blow up
   - Mass defect to become inconsistent
   - Residual to explode

## The Bug Location

In `crates/rustfoil-coupling/src/global_newton.rs`, the `apply_global_updates()` function computed relaxation based only on:
- `dn1`: Normalized ctau/ampl change
- `dn2`: Normalized theta change

But it did NOT check:
- `dn3`: Normalized delta_star change  
- `dn4`: Normalized Ue change

## The Fix

Added relaxation checks for Ue and delta_star changes **after** computing the new edge velocities via DIJ coupling:

1. **Compute new Ue** from proposed mass changes (as before)
2. **Compute normalized Ue change**: `dn4 = due / current_u`
3. **Compute normalized delta_star change**: `dn3 = d_dstar / current_delta_star`
4. **Apply relaxation limits** to both `dn3` and `dn4` using the same `dhi=1.5` and `dlo=-0.5` limits

This ensures that if Newton produces large mass changes that would cause explosive Ue or delta_star updates, the relaxation factor is reduced to keep changes within the ±50% to +150% bounds.

## Code Changes

### `crates/rustfoil-coupling/src/global_newton.rs`

- Added Step 1.5: Refine relaxation based on Ue and delta_star changes
- Checks normalized Ue changes (`dn4`) and delta_star changes (`dn3`) for both upper and lower surfaces
- Reduces relaxation if these changes would exceed `dhi=1.5` or `dlo=-0.5` limits
- Added debug output to trace when relaxation is reduced due to Ue/delta_star limits

### `crates/rustfoil-solver/src/viscous/viscal.rs`

- Added debug output for iterations 6-8 to trace explosion
- Logs station values (theta, delta_star, Ue, ctau, mass) at sample stations
- Logs maximum delta magnitudes (theta, mass, ctau)
- Logs relaxation factor used

## Testing

To verify the fix:

1. Run a case that previously exploded at iteration 7
2. Check debug output (if `RUSTFOIL_CL_DEBUG` is set) to see:
   - If relaxation is being reduced due to Ue/delta_star changes
   - Station values at iterations 6-8
   - Maximum delta magnitudes
3. Verify that residual decreases smoothly instead of exploding

## Expected Behavior

- Relaxation should be automatically reduced when large Ue or delta_star changes are detected
- Residual should decrease smoothly without explosions
- Station values should remain physically reasonable (theta < 0.1, Ue > 0.01, etc.)

## Related Issues

This fix addresses the same class of issues as:
- Large theta values causing CD overprediction
- Transition region instabilities
- Near-singular DIJ matrix conditions causing large mass changes

The relaxation safeguards now properly limit ALL variable changes, not just theta and ctau.
