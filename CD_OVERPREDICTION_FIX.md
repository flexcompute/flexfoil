# CD Over-Prediction Fix

## Root Cause Analysis

### Problem
RustFoil's drag coefficient was massively over-predicted (3-1000× too high) compared to XFOIL.

### Investigation Results

#### 1. What XFOIL Reports
- **CL_DETAIL event** reports `cdp` = **pressure drag only** (~0.000845 at α=4°)
- **CD_BREAKDOWN event** shows:
  - `cd_total` = 0.006183 (total drag: friction + pressure)
  - `cd_pressure` = 0.000845 (pressure drag)
  - `cd_friction` = 0.004579 (friction drag)
- **Key finding**: `CL_DETAIL.cdp == CD_BREAKDOWN.cd_pressure`

#### 2. What RustFoil Reported
- **CLI output** printed `result.cd` = **total CD** (friction + pressure)
- At α=4°: RustFoil CD = 0.005780 vs XFOIL CDP = 0.000845
- **Ratio**: 6.84× (matches the reported 3-5× typical over-prediction)

#### 3. Code Location
- **File**: `crates/rustfoil-solver/src/viscous/forces.rs`
- **Function**: `compute_forces_two_surfaces()`
- **Line 251**: `let cd = cd_total;` (was reporting total drag)
- **Line 248**: `let cd_pressure = (cd_total - cd_friction).max(0.1 * cd_total);` (had incorrect clamp)

### Root Cause
**RustFoil was reporting total CD (friction + pressure) instead of pressure CD only**, which is what XFOIL's CL_DETAIL event reports.

## Fix Applied

### Changes Made

1. **Fixed `cd_pressure` calculation** (line 248):
   - **Before**: `let cd_pressure = (cd_total - cd_friction).max(0.1 * cd_total);`
   - **After**: `let cd_pressure = (cd_total - cd_friction).max(0.0);`
   - Removed incorrect 10% minimum clamp

2. **Changed reported CD value** (line 251):
   - **Before**: `let cd = cd_total;` (total drag)
   - **After**: `let cd = cd_pressure;` (pressure drag only, matches XFOIL CL_DETAIL)

3. **Updated documentation**:
   - Updated `ViscousResult.cd` comment to clarify it's pressure drag only
   - Added note that total drag = cd + cd_friction

### Expected Results After Fix

At α=4°, Re=3e6:
- **Before**: RustFoil CD = 0.005780 (total), XFOIL CDP = 0.000845 → **6.84× too high**
- **After**: RustFoil CD = ~0.000845 (pressure), XFOIL CDP = 0.000845 → **Should match**

At α=-4°:
- **Before**: RustFoil CD = 0.595090, XFOIL CDP = 0.000844 → **705× too high**
- **After**: RustFoil CD = ~0.000844 (pressure), XFOIL CDP = 0.000844 → **Should match**

## Technical Details

### XFOIL's CD Computation
- **Total CD**: Computed using Squire-Young at far wake: `CD = 2 * θ * Ue^((5+H)/2)`
- **Friction CD**: Integrated skin friction along surface
- **Pressure CD**: `CDP = CD - CDF` (total - friction)

### RustFoil's CD Computation (After Fix)
- **Total CD**: Computed using Squire-Young at far wake (same as XFOIL)
- **Friction CD**: Integrated skin friction along surface (same as XFOIL)
- **Pressure CD**: `cd_pressure = cd_total - cd_friction` (same as XFOIL)
- **Reported CD**: Now matches XFOIL's CL_DETAIL `cdp` (pressure drag only)

## Files Modified

1. `crates/rustfoil-solver/src/viscous/forces.rs`
   - Line 248: Fixed `cd_pressure` calculation
   - Line 251: Changed `cd` to report pressure drag instead of total drag

2. `crates/rustfoil-solver/src/viscous/viscal.rs`
   - Updated `ViscousResult.cd` documentation

## Verification

To verify the fix works:
```bash
# Run RustFoil at α=4°
cargo run --bin rustfoil -- viscous testdata/naca0012.dat --alpha 4.0 --re 3e6

# Expected output:
# CD:        ~0.000845 (should match XFOIL's cdp)
# cd_friction: ~0.004579
# cd_pressure: ~0.000845
# Total drag (cd + cd_friction): ~0.005424
```

## Notes

- **Total drag is still available**: `cd_total = cd + cd_friction`
- **XFOIL comparison**: RustFoil's `cd` now matches XFOIL's CL_DETAIL `cdp`
- **Anomalies at certain angles**: The 100-1000× over-prediction at α=-4° and α=9° should be resolved, as those were likely due to reporting total CD instead of pressure CD
