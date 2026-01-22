# Turbulent BL Solver Fix Plan

## Problem Summary

After transition, the turbulent BL solver diverges:
- First turbulent station: Hk=2.5 ✅ (matches XFOIL)
- Station 2+: Hk collapses to 1.05, ctau hits 0.30 max clamp

## Root Cause: Shear-Lag Equation Bugs

Comparing XFOIL's `xblsys.f` (lines 1686-1841) with RustFoil's `equations.rs` (lines 810-864), several critical bugs were identified:

### Bug 1: Wrong USA Calculation

**XFOIL** (lines 822-825):
```fortran
US2 = 0.5*HS2*(1.0 - (HK2-1.0)/(GBCON*H2))
USA = 0.5*(US1 + US2)
```

**RustFoil** (line 817):
```rust
let usa = 0.5 * (0.5 * s1.hs + 0.5 * s2.hs); // WRONG: 0.25*(hs1+hs2)
```

US is the edge velocity ratio, computed from Hs, Hk, and H. RustFoil is just averaging Hs.

### Bug 2: Wrong CQA Variable

**XFOIL** (line 1689):
```fortran
CQA = (1.0-UPW)*CQ1 + UPW*CQ2  ! CQ is equilibrium shear stress coefficient
```

**RustFoil** (line 813):
```rust
let cqa = (1.0 - upw) * s1.hk + upw * s2.hk; // WRONG: should be CQ, not Hk
```

CQ is the equilibrium maximum shear stress coefficient (√Cτ_eq), which is computed from Hs, Hk, Cf, etc. It is NOT the same as Hk.

### Bug 3: Missing CQ (Equilibrium Ctau) Computation

XFOIL computes CQ in BLVAR (lines 866-884):
```fortran
CQ2 = SQRT(CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2))
```

RustFoil doesn't compute CQ at all. This is the equilibrium shear stress that the actual ctau (S) is trying to approach.

### Bug 4: Wrong DEA Computation

**XFOIL**:
```fortran
DEA = 0.5*(DE1 + DE2)  ! DE is energy thickness from BLVAR
```

**RustFoil** (lines 819-821):
```rust
let dea = 0.5 * ((3.15 + 1.72 / (s1.hk - 1.0).max(0.01)) * s1.theta
              + (3.15 + 1.72 / (s2.hk - 1.0).max(0.01)) * s2.theta);
```

This formula is close but should use the stored DE value from blvar, not recompute it.

## Implementation Plan

### Phase 1: Add Missing CQ Computation (1 day)

1. **Update `blvar()` in equations.rs**:
   - Add CQ (equilibrium shear stress) computation from XFOIL's formula
   - Store in `BlStation` (add `cq: f64` field to state.rs)
   - Include derivatives: `cq_hs`, `cq_us`, `cq_hk`, `cq_rt`, `cq_h`

2. **Add US computation**:
   - US = 0.5 * Hs * (1.0 - (Hk-1.0)/(GBCON*H))
   - Store in `BlStation` (add `us: f64` field)
   - Include derivatives: `us_hs`, `us_hk`, `us_h`

### Phase 2: Fix Shear-Lag Equation (1 day)

1. **Update `bldif()` turbulent section**:
   - Replace `cqa = Hk average` with `cqa = CQ average`
   - Replace `usa = Hs average` with `usa = US average`
   - Use stored `DE` for `dea` instead of recomputing

2. **Update Jacobian entries**:
   - `VS2[0][0]` (∂/∂ctau): currently correct
   - `VS2[0][1]` (∂/∂θ): needs CQ derivatives
   - `VS2[0][2]` (∂/∂δ*): needs CQ and US derivatives
   - `VS2[0][4]` (∂/∂x): currently correct

### Phase 3: Add Turbulent Closure Verification (0.5 day)

1. **Create turbulent closure comparison test**:
   - Compare CQ computation against XFOIL
   - Compare US computation against XFOIL
   - Verify at stations 30-35

2. **Create shear-lag residual comparison test**:
   - Trace REZC term-by-term
   - Compare against XFOIL `mrchue_iterations.json`

### Phase 4: Validate Full Turbulent March (0.5 day)

1. **Update `test_turbulent_station_comparison`**:
   - Expect Hk ~1.48 at stations 31+
   - Expect ctau ~0.04-0.05
   - Expect <10% theta error

2. **Test on multiple cases**:
   - NACA 0012 at Re=1M
   - Flat plate with forced transition

## Key XFOIL Constants

```fortran
SCCON = 5.6     ! Shear lag constant
GACON = 6.70    ! G-beta A constant
GBCON = 0.75    ! G-beta B constant
GCCON = 18.0    ! G-beta wall term
CTCON = 0.5/(GACON**2 * GBCON)  ! = 0.0148
DLCON = 0.9     ! Wake dissipation length factor
DUXCON = 1.0    ! Shear lag UxEQ weight
```

## XFOIL Shear-Lag Residual Formula

```
REZC = SCC*(CQA - SA*ALD)*DXI 
     - DEA*2.0*SLOG 
     + DEA*2.0*(UQ*DXI - ULOG)*DUXCON

where:
  SCC = SCCON * 1.333 / (1.0 + USA)
  CQA = (1-UPW)*CQ1 + UPW*CQ2  (equilibrium shear)
  SA  = (1-UPW)*S1  + UPW*S2   (actual shear = ctau)
  DEA = 0.5*(DE1 + DE2)        (energy thickness)
  SLOG = ln(S2/S1)
  ULOG = ln(U2/U1)
  UQ = equilibrium dUe/dx
```

## Success Criteria

1. Turbulent stations 31-40: Hk stable at ~1.47-1.48
2. ctau evolution: 0.03 → 0.055 → 0.04 (not hitting 0.30 clamp)
3. θ error < 10% at turbulent stations
4. All existing tests continue to pass

## Files to Modify

| File | Changes |
|------|---------|
| `crates/rustfoil-bl/src/state.rs` | Add `cq`, `us` fields to BlStation |
| `crates/rustfoil-bl/src/equations.rs` | Fix CQ, US computation in blvar; Fix shear-lag in bldif |
| `crates/rustfoil-coupling/tests/transition_validation.rs` | Update success criteria |

## Reference

- XFOIL `xblsys.f` lines 1686-1841 (shear-lag equation)
- XFOIL `xblsys.f` lines 822-884 (US and CQ computation)
- XFOIL `xbl.f` lines 1581-1592 (constants)
