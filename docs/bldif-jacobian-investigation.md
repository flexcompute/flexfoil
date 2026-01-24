# BLDIF Shape Equation Jacobian Investigation

## Summary

RustFoil's boundary layer solver has a ~2x error in the shape equation Jacobian compared to XFOIL. This causes the "wrong" Hk derivatives to accidentally work better than the mathematically correct ones.

## Current Status

- **Flat plate tests pass** with "wrong" Hk derivatives
- **Airfoil transition is 58% off** with "wrong" derivatives
- **Correct Hk derivatives break flat plate** (H collapses to 1.02)
- **BLDIF term instrumentation added** to capture intermediate momentum/shape terms

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

## Implementation Progress (January 2026)

### Completed Fixes

1. **UPW Derivative Framework** (`equations.rs` lines 720-740):
   - Added `upw_hd`, `hd_hk1`, `hd_hk2` calculations matching XFOIL
   - Computed `upw_hk1`, `upw_hk2` via chain rule
   - **Result**: Mathematically correct, but UPW derivatives are disabled (set to 0) because enabling them causes `test_march_adverse_pressure_gradient` to fail. The derivatives only contribute ~5 to VS2[2][2], not enough to explain the ~8000 gap.

2. **Hc Derivative Storage** (`state.rs`, `equations.rs`):
   - Added `hc_hk` and `hc_msq` fields to `BlDerivatives` struct
   - Stored values in `blvar()` after computing from `density_shape()`
   - Updated `bldif()` to use chain rule: `hc2_d = hc_hk * hk2_d`
   - **Result**: At M=0, `hc_hk=0` so no effect on incompressible results. Needed for compressible correctness.

3. **Hk Derivative Documentation** (`march.rs` lines 390-410):
   - Updated comments explaining the "wrong" vs "correct" derivative issue
   - Added TODO for switching to correct derivatives once VS2[2][2] is fixed
   - **Result**: Keeping "wrong" derivatives (`+Hk/θ`, `-Hk/δ*`) as workaround

### Findings

1. **VS2[2][2] still at ~53% of XFOIL** despite UPW and Hc fixes
   - RustFoil: -9319 (was -9309 before fixes)
   - XFOIL: -17548
   - The fixes added only ~12 to the value

2. **UPW derivatives cause solver instability** when enabled
   - Test `test_march_adverse_pressure_gradient` fails
   - Suggests interaction with other Jacobian terms

3. **Transition error remains at 58%** (XFOIL x_tr=0.149, RustFoil=0.235)
   - Same as before because underlying ~2x Jacobian error not resolved

### 2026-01-23 Update: BLDIF Terms vs XFOIL

Added `BLDIF_TERMS` debug events to log the intermediate quantities used in the momentum
and shape equations (xlog/ulog/tlog/hlog, UPW, HA, BTMP, CFX, DIX, etc.). These are
emitted for each station during `march_surface` and written to `rustfoil_debug.json`.

Key observation at lower surface IBL 64–67:
- **HA/BTMP** are ~1.0 lower in RustFoil, driven by smaller `delta_star` (lower H).
- **CFX/CFX_SHAPE** are shifted in sign/magnitude versus XFOIL when recomputed from XFOIL
  states, confirming that the *state* (θ, δ*) mismatch dominates the term differences.
- **UPW/xlog/ulog** differences are small compared to the HA/CFX shifts.

This supports the current hypothesis: the laminar Newton update is producing a
lower `delta_star` state (and thus lower H/Hk), rather than the residual formulas
being wrong.

### 2026-01-23 Update: MRCHUE_ITER Debug

Added `MRCHUE_ITER` debug events during the Newton loop to capture iteration-by-iteration
residuals and updates. Each event logs:
- residuals (`res_third`, `res_mom`, `res_shape`)
- update vector (`delta_s`, `delta_theta`, `delta_delta_star`, `delta_ue`)
- relaxation and `dmax`
- current station state (`theta`, `delta_star`, `Ue`, `H`, `Hk`, `Cf`, `Cd`, `ampl`, `ctau`)

This should let us pinpoint which Newton step drives `delta_star` low at the problematic
lower-surface stations (IBL 64–67).

Initial MRCHUE_ITER observations (lower IBL 65, current run):
- Iter 1 computed in direct mode, then switches to inverse mode without applying updates.
- Converges after 3 iterations with `dmax ≈ 1.3e-2` (tolerance is 0.1).
- Residuals at convergence are still non-trivial (`res_mom ~ 1.8e-3`, `res_shape ~ -7.3e-3`).
- `delta_star` climbs to ~2.87e-3, still below XFOIL’s ~3.61e-3 for this case.

### 2026-01-23 Update: Transition Theta Alignment

After moving TRCHEK2 evaluation into the Newton loop and aligning the transition
system inputs, the TRDIF Newton system now matches XFOIL at ~1e-5 relative.
To validate **theta at transition**, use the TRCHEK2_ITER event (wf1/wf2 + T1/T2)
because XFOIL does not emit `tt` in TRCHEK2_FINAL. The TRCHEK2_ITER-based comparison
matches at ~1e-10 for the transition station in the current case.

This suggests the convergence criterion may be letting the laminar solve stop early
before `delta_star` reaches the XFOIL state; the iteration logs now make that visible.

### 2026-01-23 Update: Convergence Criterion Alignment

XFOIL’s `MRCHUE` loop converges on `DMAX <= 1e-5` and only accepts `DMAX <= 0.1`
after the iteration limit is reached (see `Xfoil/src/xbl.f` lines ~662–786).
RustFoil now uses `tolerance = 1e-5` and only falls back when `dmax > 0.1`.

With the tighter tolerance:
- Rust MRCHUE iteration counts now match XFOIL’s scale (5–7 iterations).
- Lower-surface `delta_star` at IBL 64–67 increased and moved closer to XFOIL’s
  state (Hk now tracks XFOIL in that region).
- Transition location shifted slightly on the upper surface; lower surface stays
  near the prior transition index.

### 2026-01-23 Update: Machine-Level Comparison (Post-Tolerance)

Re-ran `compare_trchek2_full.py` and BLVAR comparisons with the new tolerance.
Results still diverge from XFOIL at machine precision:

- **TRCHEK2_FINAL** AX and XT still high: AX errors ~40–70% and XT shifts
  (`x_tr(U)` ~0.2369 vs 0.2403, `x_tr(L)` ~0.9206 vs 0.9286).
- **Lower BLVAR H/Hk** still offset by ~0.7 at IBL 64–67 even after convergence
  tightening.
- **Lower BLVAR Us/Cf/Cq** remain different (sign and magnitude changes).

Conclusion: convergence tightening alone is insufficient to reach machine-level
agreement; there is still a modeling or state mismatch in the laminar march inputs.

### Remaining Investigation

The ~8000 missing contribution to VS2[2][2] is not from:
- UPW derivatives (contribute ~5)
- Hc derivatives (zero at M=0)
- Closure derivatives (verified to match XFOIL)
- Z coefficients (verified to match XFOIL)

Possible unexplored sources:
- Numerical differences in upwinding calculation at specific flow conditions
- Interaction between equations in the Newton system
- Subtle differences in how XFOIL handles the BLDIF residual vs Jacobian

## Related Commits

- `e2d5839` - test: Add VS2 shape Jacobian comparison test against XFOIL
- Previous commits documenting the Hk derivative investigation

## Contact

This investigation was done comparing against XFOIL instrumented at `Xfoil-instrumented/`. The instrumentation outputs VS2 matrices at each Newton iteration to `testdata/mrchue_iterations.json`.
