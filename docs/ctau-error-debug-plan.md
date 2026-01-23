# CTAU Error Debugging Plan

## Status: INVESTIGATION COMPLETE

### Summary of Findings

| Metric | Before | After | Target | Status |
|--------|--------|-------|--------|--------|
| res_third sign | 727% (sign flip) | Correct | - | ✅ Fixed |
| Transition location | - | 0.1% | - | ✅ Excellent |
| Amplification | - | 0.01% | - | ✅ Excellent |
| ctau at transition | 1.9% | 1.9% | < 0.5% | ⚠️ Needs work |
| theta at transition | 13.5% | 13.5% | < 5% | ⚠️ Needs work |
| theta downstream | - | 3-4% | acceptable | ✅ Good |

### Key Verified Items

1. **All Z-coefficients match XFOIL exactly**:
   - `Z_SA = -SCC*DXI*ALD` ✅
   - `Z_SL = -DEA*2.0` ✅
   - `Z_S2 = UPW*Z_SA + Z_SL/S2` ✅
   - `Z_UPW` terms added to Jacobian ✅

2. **UPW calculation is correct**: Both use `ABS()` for Hk ratio ✅

3. **Sign convention is correct**: `VSREZ(1) = -REZC` ✅

4. **Closure relations**: 100% pass rate (BLVAR, HSL, CFL, CFT, HKIN, DAMPL) ✅

### Remaining Issue

The 13.5% theta error at the first turbulent station is localized and diminishes
rapidly downstream (to 3-4% within 1-2 stations). This suggests accurate turbulent
marching but potential differences in:
- Transition initialization (initial theta/delta_star values)
- Inverse mode handling at Hk = 2.5 constraint

---

## Original Problem Statement

The first turbulent station after transition shows:
- **ctau error: 1.9%** (RustFoil: 0.055636, XFOIL: 0.054610)
- **theta error: 13.5%** (RustFoil: 2.769e-4, XFOIL: 2.440e-4)
- **delta_star error: 13.5%** (RustFoil: 6.923e-4, XFOIL: 6.099e-4)

Target: ctau error < 0.5%

## Original Root Cause Analysis (RESOLVED)

From the ctau_debug test output, the discrepancies trace back to `bldif`:

| Component | RustFoil | XFOIL | Error |
|-----------|----------|-------|-------|
| VS2[0][0] (∂shear/∂ctau) | -0.1848 | -0.1053 | 75.5% |
| VS2[0][1] (∂shear/∂θ) | -10.83 | -9.32 | 16.3% |
| res_third (shear-lag) | -7.265e-3 | +1.158e-3 | 727% (sign flip!) |
| Hk (target) | 2.5000 | 2.5000 | 0% ✓ |

**Resolution**: The apparent "sign flip" was due to test setup using wrong input values
that didn't match XFOIL's converged state. When using correct inputs, the formulas match.

---

## Phase 1: Verify Shear-Lag Residual Sign

### Task 1.1: Trace REZC calculation

The shear-lag residual formula (from XFOIL `xblsys.f:1684-1710`):

```fortran
REZC = SCC*(CQA - SA*ALD)*DXI 
     - DEA*2.0*SLOG 
     + DEA*2.0*(UQ*DXI - ULOG)*DUXCON
```

Where:
- `SCC = SCCON * 1.333 / (1 + USA)` — lag constant
- `CQA` — equilibrium shear coefficient (averaged)
- `SA` — actual shear stress τ (averaged)
- `ALD = 1.0` for turbulent flow
- `DXI = X2 - X1` — spatial step
- `SLOG = ln(S2/S1)` — shear ratio log
- `ULOG = ln(U2/U1)` — velocity ratio log
- `UQ` — velocity gradient effect
- `DUXCON = 1.0` — constant

**Hypothesis**: RustFoil may have a sign error in one of these terms.

**Test**: Add debug prints to RustFoil's `bldif` for:
1. Each term in REZC before summing
2. The final REZC value before negation
3. Compare term-by-term against XFOIL instrumented output

### Task 1.2: Check VSREZ sign convention

XFOIL returns `VSREZ(1) = -REZC` (line 1874). Verify RustFoil does the same.

### Task 1.3: Instrument XFOIL for comparison

Modify `Xfoil-instrumented/src/xblsys.f` to print intermediate values at the first turbulent station:

```fortran
IF(IBL.EQ.30 .AND. SIDE.EQ.1) THEN
  WRITE(*,*) 'REZC terms:', SCC*(CQA-SA*ALD)*DXI, -DEA*2.0*SLOG, DEA*2.0*(UQ*DXI-ULOG)*DUXCON
  WRITE(*,*) 'REZC=', REZC
ENDIF
```

---

## Phase 2: Compare Jacobian Term-by-Term

### Task 2.1: Z-coefficient comparison

The Jacobian uses intermediate coefficients (xblsys.f:1796-1841):

```fortran
Z_CFA = DEA*2.0*UQ_CFA*DXI * DUXCON
Z_HKA = DEA*2.0*UQ_HKA*DXI * DUXCON
Z_DE1 = (   -SLOG + (UQ*DXI - ULOG)*DUXCON) * 0.5*2.0
Z_US1 = SCC*(CQA - SA*ALD)*DXI * USA1_USA
Z_S1  = -SCC * ALD * DXI * (1.0-UPW)
```

**Test**: Print all Z-coefficients in RustFoil and compare against XFOIL.

### Task 2.2: Check VS2[0][0] calculation

VS2[0][0] should be `Z_S2` which is:
```fortran
Z_S2 = -SCC * ALD * DXI * UPW
```

RustFoil: -0.1848, XFOIL: -0.1053 (75% error)

This suggests either:
- `SCC` differs
- `DXI` differs  
- `UPW` differs

### Task 2.3: Check VS2[0][1] calculation

VS2[0][1] (∂shear/∂θ) is built from:
```fortran
VS2(1,2) = Z_UPW*UPW_T2 + Z_DE2*DE2_T2 + Z_US2*US2_T2
         + Z_CQ2*CQ2_T2 + Z_CF2*CF2_T2 + Z_HK2*HK2_T2
```

Each term needs verification.

---

## Phase 3: Check Input State Consistency

### Task 3.1: Verify station inputs match

The test uses hardcoded values. Verify:
- `s1` and `s2` x, u, θ, δ*, ctau match XFOIL's state at call_id=2916
- The flow state (laminar→turbulent transition) is correct

### Task 3.2: Check blvar outputs at transition

Even though blvar matches for mid-turbulent stations, verify it matches at the **transition point** where Hk is forced to 2.5:

```rust
// Expected from XFOIL at transition:
station.hk = 2.5;  // Forced
station.cf = ?;    // Should match XFOIL
station.cq = ?;    // Should match XFOIL
```

### Task 3.3: Check UPW (upwinding) calculation

UPW depends on Hk ratio. At transition, Hk changes abruptly from laminar (~5.4) to turbulent (2.5), which affects upwinding.

Verify:
```rust
let arg = (hk2 - 1.0) / (hk1 - 1.0);
let hl = arg.ln();
```

If `arg < 0` (shouldn't happen), this would produce NaN.

---

## Phase 4: Create Minimal Reproduction

### Task 4.1: Extract exact XFOIL state

From the instrumented test output, extract the **exact** station state that XFOIL uses when computing bldif at call_id=2916.

### Task 4.2: Create standalone comparison test

```rust
#[test]
fn test_bldif_transition_exact() {
    // Use exact XFOIL input values
    let s1 = BlStation { /* exact values from XFOIL */ };
    let s2 = BlStation { /* exact values from XFOIL */ };
    
    let (res, jac) = bldif(&s1, &s2, FlowType::Turbulent, msq, re);
    
    // Compare each intermediate value
    // ...
}
```

---

## Phase 5: Fix and Validate

### Task 5.1: Apply fixes based on findings

Common issues to check:
- Sign errors in formulas
- Missing/extra negations
- Different variable interpretations
- Constant mismatches (SCCON, GACON, etc.)

### Task 5.2: Re-run validation tests

After fixes:
```bash
cargo test --package rustfoil-coupling --test ctau_debug -- --nocapture
cargo test --package rustfoil-coupling test_transition_location_upper -- --nocapture
```

### Task 5.3: Update test vectors if needed

The bldif_test_vectors.json may need regeneration with correct XFOIL state.

---

## Success Criteria

1. **res_third sign matches XFOIL** (same sign, not inverted)
2. **VS2[0][0] error < 10%** (currently 75%)
3. **VS2[0][1] error < 10%** (currently 16%)
4. **ctau error < 0.5%** (currently 1.9%)
5. **theta error < 5%** (currently 13.5%)

---

## Test Commands

```bash
# Run ctau debug test
cargo test --package rustfoil-coupling --test ctau_debug -- --nocapture

# Run transition test
cargo test --package rustfoil-coupling test_transition_location_upper -- --nocapture

# Run bldif comparison
cargo test --package rustfoil-bl test_bldif_matches_xfoil -- --nocapture

# Run blvar comparison  
cargo test --package rustfoil-bl test_blvar_matches_xfoil -- --nocapture
```

---

## Files to Modify

| File | Purpose |
|------|---------|
| `crates/rustfoil-bl/src/equations.rs` | BLDIF implementation |
| `crates/rustfoil-coupling/tests/ctau_debug.rs` | Debug test |
| `Xfoil-instrumented/src/xblsys.f` | XFOIL instrumentation |
| `testdata/bldif_test_vectors.json` | Test data (may need regen) |

---

## XFOIL Reference

- **BLDIF**: `Xfoil-instrumented/src/xblsys.f` lines 1584-1874
- **Shear-lag equation**: lines 1684-1710
- **Z-coefficients**: lines 1796-1841
- **Jacobian assembly**: lines 1848-1873
- **Constants**: `Xfoil/src/xbl.f` line ~1587 (SCCON=5.6, GACON=6.70, etc.)
