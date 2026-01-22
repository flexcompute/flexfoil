# XFOIL vs RustFoil I/O Comparison

This document traces the I/O flow through the boundary layer Newton solver in both XFOIL and RustFoil, identifying the key differences that may cause divergence.

## Overview: Newton System Structure

Both XFOIL and RustFoil solve a 4-variable Newton system at each station:

```
[A] * [dS, dθ, dδ*, dUe]ᵀ = [residuals]ᵀ
```

Where:
- `S` = Amplification (laminar) or Ctau (turbulent)
- `θ` = Momentum thickness
- `δ*` = Displacement thickness  
- `Ue` = Edge velocity

---

## Module Flow Comparison

### XFOIL Flow (xbl.f MRCHUE → xblsys.f BLDIF)

```
MRCHUE (xbl.f:542)
    │
    ├── BLPRV: Set primary variables (θ, δ*, Ue, x, ctau, ampl)
    │
    ├── BLSYS (xblsys.f:583)
    │     │
    │     ├── BLVAR(flow_type): Compute secondary variables and derivatives
    │     │     │
    │     │     ├── BLKIN: Set H, M, V, Hk, Rθ from primaries
    │     │     │         H2_T2 = -H/θ     (NEGATIVE)
    │     │     │         H2_D2 = 1/θ      (POSITIVE)
    │     │     │         HK2_T2 = HK2_H2 * H2_T2  (≈ -Hk/θ at low Mach)
    │     │     │         HK2_D2 = HK2_H2 * H2_D2  (≈ 1/θ at low Mach)
    │     │     │
    │     │     └── Cf, Hs, Cd correlations with derivatives
    │     │
    │     └── BLDIF(flow_type): Build 3x5 Jacobian (VS1, VS2) and residuals
    │           │
    │           ├── Row 0: Amplification/Shear-lag equation
    │           ├── Row 1: Momentum equation  
    │           └── Row 2: Shape parameter equation
    │
    ├── Set Row 3 (4th equation):
    │     DIRECT MODE:  VS2(4,:) = [0, 0, 0, 1], VSREZ(4) = 0
    │     INVERSE MODE: VS2(4,:) = [0, HK2_T2, HK2_D2, HK2_U2]
    │                   VSREZ(4) = HTARG - HK2
    │
    └── GAUSS: Solve 4x4 system
```

### RustFoil Flow (march.rs → equations.rs → solve.rs)

```
newton_solve_station (march.rs:346)
    │
    ├── Initialize station from previous
    │
    ├── blvar(station, flow_type): Compute secondary variables
    │     │
    │     └── Computes: H, Hk, Hs, Cf, Cd, Rθ and derivatives
    │
    ├── bldif(prev, curr, flow_type): Build 3x5 Jacobian and residuals
    │     │
    │     ├── Row 0: Amplification/Shear-lag residual + Jacobian
    │     ├── Row 1: Momentum residual + Jacobian
    │     └── Row 2: Shape parameter residual + Jacobian
    │
    ├── build_4x4_system(vs2, res, direct, hk2_t, hk2_d, ...):
    │     │
    │     ├── Extract columns [0:4] from vs2 (drop x column)
    │     ├── Row 3 (4th equation):
    │     │     DIRECT:  [0, 0, 0, 1] * x = 0
    │     │     INVERSE: [0, hk2_t, hk2_d, hk2_u] * x = hk_target - hk_current
    │     └── NOTE: hk2_t and hk2_d are PASSED IN, not computed from BLVAR!
    │
    └── solve_4x4: Gaussian elimination with pivoting
```

---

## Critical Difference: Hk Derivatives

### XFOIL Computation (xblsys.f:765-769)

```fortran
C---- set shape parameter
      H2    =  D2/T2
      H2_D2 = 1.0/T2          ! ∂H/∂δ* = 1/θ  (POSITIVE)
      H2_T2 = -H2/T2          ! ∂H/∂θ = -H/θ   (NEGATIVE)

C---- set kinematic shape parameter
      CALL HKIN( H2, M2, HK2, HK2_H2, HK2_M2 )

      HK2_T2 = HK2_H2*H2_T2   ! ≈ -Hk/θ at low Mach (NEGATIVE)
      HK2_D2 = HK2_H2*H2_D2   ! ≈ 1/θ at low Mach   (POSITIVE)
```

### RustFoil Computation (march.rs:393-402)

```rust
// In newton_solve_station:
let (mut a, mut b) = build_4x4_system(
    &jac.vs2,
    &[res.res_third, res.res_mom, res.res_shape],
    direct,
    station.hk / station.theta,       // hk2_t: +Hk/θ  (WRONG SIGN!)
    -station.hk / station.delta_star, // hk2_d: -Hk/δ* (WRONG FORMULA!)
    0.0,  // hk2_u
    htarg,
    station.hk,
);
```

### The Problem

| Derivative | XFOIL | RustFoil | Issue |
|------------|-------|----------|-------|
| HK2_T2 (∂Hk/∂θ) | `-Hk/θ` | `+Hk/θ` | **Wrong sign** |
| HK2_D2 (∂Hk/∂δ*) | `+1/θ` | `-Hk/δ*` | **Wrong formula AND sign** |

### Mathematical Derivation

For incompressible flow, `Hk ≈ H = δ*/θ`. Therefore:

```
∂H/∂θ  = ∂(δ*/θ)/∂θ  = -δ*/θ² = -H/θ     (NEGATIVE)
∂H/∂δ* = ∂(δ*/θ)/∂δ* = 1/θ               (POSITIVE)
```

RustFoil's formula `-Hk/δ*` is mathematically wrong:
- `-Hk/δ* = -(δ*/θ)/(δ*) = -1/θ` (negative instead of positive!)

---

## Residual Sign Convention

Both XFOIL and RustFoil use `VSREZ = -residual` convention in the 3 main equations:

| XFOIL | RustFoil |
|-------|----------|
| `VSREZ(1) = -REZC` (amplification) | `res.res_third = -(...)` |
| `VSREZ(2) = -REZT` (momentum) | `res.res_mom = -(...)` |
| `VSREZ(3) = -REZH` (shape) | `res.res_shape = -(...)` |

For the 4th row (inverse mode constraint), XFOIL uses:
```fortran
VSREZ(4) = HTARG - HK2  ! NOT negated
```

RustFoil:
```rust
b[3] = hk_target - hk_current;  // Same convention
```

This appears consistent.

---

## Variable Order in Jacobian

Both use the same variable order:

| Index | Variable | Description |
|-------|----------|-------------|
| 0 | S/Ampl | Ctau (turbulent) or N-factor (laminar) |
| 1 | θ | Momentum thickness |
| 2 | δ* | Displacement thickness |
| 3 | Ue | Edge velocity |
| 4 | x | Arc length (dropped in 4x4 system) |

---

## Observed Behavior

The Hk derivative issue causes a conflict:

1. **With XFOIL-correct derivatives** (`hk2_t = -Hk/θ`, `hk2_d = 1/θ`):
   - Airfoil transition prediction: **9.7% error** (good)
   - Flat plate Blasius test: **FAILS** (H collapses to ~0.0001)

2. **With current RustFoil derivatives** (`hk2_t = +Hk/θ`, `hk2_d = -Hk/δ*`):
   - Airfoil transition prediction: **58% error** (bad)
   - Flat plate Blasius test: **PASSES** (H stays reasonable)

This suggests:
- The flat plate test may have a compensating error elsewhere
- OR the inverse mode is rarely triggered for flat plates (always direct mode)
- OR there's a sign convention difference in how the 4th row interacts with the solver

---

## Recommended Fix

1. **Correct the Hk derivatives in `newton_solve_station`**:
   ```rust
   let hk2_t = -station.hk / station.theta;  // Correct: negative
   let hk2_d = 1.0 / station.theta;          // Correct: 1/θ, positive
   ```

2. **Investigate why flat plate fails with correct derivatives**:
   - Check if inverse mode is ever triggered for flat plates
   - Compare the minimum Hk clamp behavior
   - Verify the similarity station initialization

3. **Consider using `BlStation.derivs` consistently**:
   - The `blvar` function already computes `hk_h * h_theta` and `hk_h * h_delta_star`
   - Use these stored derivatives instead of recomputing in `newton_solve_station`
