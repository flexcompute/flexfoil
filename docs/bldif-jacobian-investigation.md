# BLDIF Shape Equation Jacobian Investigation

## Summary

RustFoil's boundary layer solver has a ~2x error in the shape equation Jacobian compared to XFOIL. This causes the "wrong" Hk derivatives to accidentally work better than the mathematically correct ones.

## Current Status

- **Flat plate tests pass** with "wrong" Hk derivatives
- **Airfoil transition is 58% off** with "wrong" derivatives
- **Correct Hk derivatives break flat plate** (H collapses to 1.02)

## Root Cause Identified

The shape equation Jacobian entry `VS2[2][2]` (∂shape/∂δ*) has ~2x error:

| Metric | XFOIL | RustFoil | Match? |
|--------|-------|----------|--------|
| VS2[2][2] (∂shape/∂δ*) | -17548 | -9309 | ❌ (~53%) |
| VS2[2][1] (∂shape/∂θ) | 18504 | 18500 | ✅ |
| VS2[1][1] (∂mom/∂θ) | 20180 | 20180 | ✅ |

## What Has Been Verified

All individual components match XFOIL:

1. **Closure derivatives** (from `blvar`):
   - `hs_hk`, `cf_hk`, `cd_hk` ✓
   - `h_theta = -H/θ`, `h_delta_star = 1/θ` ✓

2. **Z coefficients** (from `bldif`):
   - `z_hs2 = 1/Hs2` ✓
   - `z_cf2 = upw * xlog/2 * (x2/θ2)` ✓
   - `z_di2 = -upw * xlog * (x2/θ2)` ✓

3. **Derivative chain**:
   - `hs2_d = hs_hk * hk_h * h_d` ✓
   - `cf2_d = cf_hk * hk_h * h_d` ✓
   - `di2_d = cd_hk * hk_h * h_d` ✓

## Known Issues

### 1. UPW Derivatives Hardcoded to Zero

**File**: `crates/rustfoil-bl/src/equations.rs` lines 1051-1056

```rust
let upw_t1 = 0.0;
let upw_d1 = 0.0;
let upw_u1 = 0.0;
let upw_t2 = 0.0;  // Should be: upw_hk2 * hk2_t
let upw_d2 = 0.0;  // Should be: upw_hk2 * hk2_d
let upw_u2 = 0.0;
```

XFOIL computes these dynamically (xblsys.f lines 1637-1642). However, this only contributes ~5 to VS2[2][2], not the ~8000 needed.

### 2. Missing ~8000 Contribution

The sum of all computed terms gives -9309, but XFOIL gives -17548. There's approximately -8000 missing somewhere.

**Possible sources to investigate:**
- Additional terms in XFOIL's BLDIF not replicated in RustFoil
- Different upwinding scheme at certain stations
- Compressibility corrections even at M=0
- Different handling of the `DDLOG` factor

## Files to Examine

### XFOIL Reference
- `Xfoil/src/xblsys.f` lines 1901-1976 (BLDIF shape equation)
- `Xfoil/src/xblsys.f` lines 1605-1644 (UPW upwinding)
- `Xfoil/src/xblsys.f` lines 725-780 (BLKIN - HK derivatives)

### RustFoil Implementation
- `crates/rustfoil-bl/src/equations.rs` - `bldif()` function (~lines 700-1087)
- `crates/rustfoil-coupling/src/march.rs` - `newton_solve_station()` (~lines 350-500)
- `crates/rustfoil-coupling/src/solve.rs` - `build_4x4_system()` (~lines 556-600)

### Test Data
- `testdata/mrchue_iterations.json` - XFOIL instrumented output with VS2 matrices
- `testdata/bldif_test_vectors.json` - BLDIF test cases

## Test to Run

```bash
cargo test --package rustfoil-coupling --lib -- solve::tests::test_bldif_vs2_shape_equation_xfoil_comparison --nocapture
```

This compares RustFoil's VS2 against XFOIL at IBL=11 and prints detailed breakdown of each term.

## What Needs to Be Done

### Step 1: Find the Missing Term

Compare XFOIL's BLDIF line-by-line against RustFoil's `bldif()` for the shape equation (row 3 / index 2). Look for:
- Terms that XFOIL adds but RustFoil doesn't
- Different coefficient values
- Sign errors in multiplication

### Step 2: Fix UPW Derivatives

Implement proper UPW derivative computation in `equations.rs`:
```rust
// Replace hardcoded zeros with:
let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;
let upw_t2 = upw_hk2 * hk2_t;
let upw_d2 = upw_hk2 * hk2_d;
let upw_u2 = upw_hk2 * hk2_u;
```

### Step 3: Verify Fix

Once VS2[2][2] matches XFOIL:
1. Switch to correct Hk derivatives in `march.rs`:
   ```rust
   let hk2_t = station.derivs.hk_h * station.derivs.h_theta;      // -Hk/θ
   let hk2_d = station.derivs.hk_h * station.derivs.h_delta_star; // 1/θ
   ```
2. Run flat plate test - should now pass
3. Run airfoil transition test - should improve from 58% error

### Step 4: Clean Up

- Remove "wrong" derivative workaround comments
- Update transition_validation.rs expected tolerances
- Document the fix

## Key Insight

The "wrong" Hk derivatives (`+Hk/θ`, `-Hk/δ*`) accidentally compensate for the ~2x error in VS2[2][2]. When the Jacobian is fixed, the correct derivatives (`-Hk/θ`, `+1/θ`) should work properly.

## Related Commits

- `e2d5839` - test: Add VS2 shape Jacobian comparison test against XFOIL
- Previous commits documenting the Hk derivative investigation

## Contact

This investigation was done comparing against XFOIL instrumented at `Xfoil-instrumented/`. The instrumentation outputs VS2 matrices at each Newton iteration to `testdata/mrchue_iterations.json`.
