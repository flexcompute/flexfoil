# BLVAR Debug Plan

## Problem Statement

Transition location in RustFoil (x_tr ≈ 0.23) differs significantly from XFOIL (x_tr ≈ 0.15). Root cause analysis revealed that the **Hk (kinematic shape factor) values are ~25% lower** in RustFoil than XFOIL:

| Region | XFOIL Hk | RustFoil Hk |
|--------|----------|-------------|
| Pre-transition | 3.3-3.6 | 2.6-2.8 |

Lower Hk → Lower amplification rate → Delayed transition

## CRITICAL FINDING (Phase 2 Complete)

**The `blvar` closure functions are CORRECT.** Testing with XFOIL's exact inputs shows:

| Output | Error vs XFOIL |
|--------|----------------|
| H | 0.00% |
| Hk | 0.00% |
| Hs | 0.00% |
| Cf | 0.00% |
| Cd | 0.00% |
| Rtheta | 0.00% |

**Root cause is UPSTREAM**: The BL march computes different θ and δ* values than XFOIL, which then produce different Hk when passed to blvar.

## Updated Goal

Debug the BL march (Newton iteration) to achieve matching θ and δ* values at each station.

---

## NEW Phase: Debug BL March (Newton Iteration)

Since closures are correct, focus on why θ and δ* differ during marching.

### Hypothesis 1: Initial Conditions
- [ ] Compare stagnation point initialization
- [ ] Check first station values (θ₀, δ*₀)
- [ ] Verify BLDIF residual at first station

### Hypothesis 2: Newton Iteration Convergence
- [ ] Compare iteration count per station
- [ ] Check residual magnitude at convergence
- [ ] Verify Jacobian matrix matches XFOIL

### Hypothesis 3: Step Accumulation Error
- [ ] Compare θ evolution station-by-station
- [ ] Check if error grows with arc length
- [ ] Identify first station with significant divergence

### Key Files to Debug
- `crates/rustfoil-coupling/src/march.rs` - BL marching logic
- `crates/rustfoil-bl/src/equations.rs` - BLDIF residuals and Jacobian
- `Xfoil-instrumented/src/xbl.f` - XFOIL's SETBL/MRCHUE

### Test Plan
1. Extract MRCHUE_ITER events from XFOIL debug output
2. Compare θ_in, δ*_in, θ_out, δ*_out at each Newton iteration
3. Identify where RustFoil diverges from XFOIL

---

## Phase 1: Instrument XFOIL BLVAR (COMPLETE)

### 1.1 Verify existing BLVAR instrumentation
- [ ] Check `xfoil_debug.f` for existing `DBGBLVAR` subroutine
- [ ] Verify it captures all key outputs: H, Hk, Hs, Hc, Rtheta, Cf, Cd, Us, Cq, De
- [ ] Verify it captures inputs: theta, delta_star, ctau, Ue, x, side, ibl

### 1.2 Add missing instrumentation if needed
- [ ] Ensure BLVAR is called at consistent points (after Newton convergence)
- [ ] Capture Mach-related variables (Msq, compressibility corrections)
- [ ] Add flow type indicator (laminar vs turbulent)

### 1.3 Generate fresh test vectors
- [ ] Run instrumented XFOIL at Re=1M, alpha=5°, Ncrit=9
- [ ] Extract BLVAR events from `xfoil_debug.json`
- [ ] Create `testdata/blvar_test_vectors.json` with ~50 representative cases

---

## Phase 2: Create BLVAR Comparison Test

### 2.1 Test structure
```rust
// crates/rustfoil-bl/tests/blvar_comparison.rs

#[test]
fn test_blvar_vs_xfoil() {
    // Load test vectors
    // For each vector:
    //   1. Set up BLStation with XFOIL's inputs
    //   2. Call blvar()
    //   3. Compare all outputs vs XFOIL
}
```

### 2.2 Key comparisons
| Output | Tolerance | Priority |
|--------|-----------|----------|
| Hk | < 1% | **CRITICAL** |
| Rtheta | < 1% | HIGH |
| Cf | < 5% | HIGH |
| Hs | < 2% | MEDIUM |
| Us | < 2% | MEDIUM |
| Cd | < 5% | MEDIUM |
| Cq | < 5% | LOW |
| De | < 5% | LOW |

### 2.3 Separate laminar vs turbulent
- Test laminar closures independently (stations before transition)
- Test turbulent closures independently (stations after transition)
- Identify which closure set has larger errors

---

## Phase 3: Trace Hk Computation

### 3.1 Understand the Hk computation chain
```
H = delta_star / theta
    ↓
Hk = f(H, Msq)  [compressibility correction]
    ↓
Used by: amplification_rate(), Cf correlations, etc.
```

### 3.2 Files to examine
- `crates/rustfoil-bl/src/closures/shape_factors.rs` - Hk computation
- `crates/rustfoil-bl/src/equations.rs` - blvar main logic
- `Xfoil-instrumented/src/xbl.f` - XFOIL's BLVAR subroutine

### 3.3 Specific checks
- [ ] Compare H computation (should be trivial: δ*/θ)
- [ ] Compare Hk compressibility correction formula
- [ ] Check Hk limits (HKLIM in XFOIL)
- [ ] Verify Msq (Mach squared) is computed correctly

---

## Phase 4: Compare Individual Closures

### 4.1 Shape factor closures
- [ ] `hkin()` - Kinematic shape factor Hk from H
- [ ] `hsl()` - Density shape factor Hs from Hk
- [ ] `hct()` - H-star correlation

### 4.2 Friction/dissipation closures
- [ ] `cfl()` - Laminar skin friction Cf
- [ ] `cft()` - Turbulent skin friction Cf
- [ ] `cdl()` - Laminar dissipation Cd
- [ ] `cdt()` - Turbulent dissipation Cd

### 4.3 Velocity defect closures
- [ ] `usl()` - Laminar Us
- [ ] `ust()` - Turbulent Us

### 4.4 For each closure:
1. Extract XFOIL reference values from debug output
2. Create focused unit test with same inputs
3. Compare output and derivatives
4. Fix any discrepancies found

---

## Phase 5: Integration Testing

### 5.1 Station-by-station march comparison
- [ ] Run march with XFOIL's exact Ue distribution
- [ ] Compare BL state at each station
- [ ] Identify where errors accumulate

### 5.2 Transition location validation
- [ ] After fixes, re-run transition validation test
- [ ] Target: x_tr error < 10% on both surfaces

### 5.3 Force coefficient validation
- [ ] Compare Cd (drag coefficient)
- [ ] Target: < 10% error vs XFOIL

---

## Implementation Order (Updated)

1. ~~**Phase 1**: Verify/generate test vectors~~ ✅ COMPLETE
2. ~~**Phase 2**: Create comparison test framework~~ ✅ COMPLETE - closures match 100%
3. ~~**Phase 3**: Check Hk computation~~ ✅ NOT NEEDED - closures are correct
4. ~~**Phase 4**: Fix shape factor closures~~ ✅ NOT NEEDED - closures are correct
5. **NEW**: Debug BL march θ/δ* computation
6. **Phase 5.2**: Validate transition improvement

---

## Success Criteria

- [x] All BLVAR outputs within tolerance vs XFOIL ✅ (0.00% error)
- [x] Hk values match within 1% ✅ (0.00% error)
- [ ] θ and δ* values match XFOIL during march (< 5% error)
- [ ] Transition location within 10% of XFOIL
- [ ] No regression in existing tests

---

## Reference Files

### RustFoil
- `crates/rustfoil-bl/src/equations.rs` - blvar function
- `crates/rustfoil-bl/src/closures/` - individual closure functions
- `crates/rustfoil-bl/src/station.rs` - BLStation struct

### XFOIL
- `Xfoil-instrumented/src/xbl.f` - BLVAR subroutine (lines 2100-2500)
- `Xfoil-instrumented/src/xblsys.f` - Supporting routines
- `Xfoil-instrumented/src/XBL.INC` - Common block definitions
