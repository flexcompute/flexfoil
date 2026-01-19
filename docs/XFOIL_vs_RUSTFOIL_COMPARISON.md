# Comprehensive Comparison: XFOIL vs RustFoil

## Executive Summary

RustFoil underpredicts stall (Cl_max too low, stall angle too early) due to several fundamental differences from XFOIL's approach. This document identifies the key differences and recommends fixes.

---

## 1. Viscous-Inviscid Coupling Architecture

### XFOIL Approach
- **Simultaneous Global Newton System**: XFOIL solves the inviscid and boundary layer equations **simultaneously** in one large Newton system
- **Variables**: For each BL station: [Ctau/Ampl, θ, m] where m = Ue·δ* (mass defect)
- **Coupling**: The DIJ matrix relates mass defect at station j to edge velocity at station i:
  ```
  Ue_i = Ue_inviscid_i + Σ_j DIJ(i,j) * m_j
  ```
- **Key Insight**: Ue is NOT an independent unknown - it's computed from the mass defect and inviscid solution through the DIJ matrix

### RustFoil Approach  
- **Sequential BL Marching with Iteration**: RustFoil marches the BL with fixed Ue, then updates Ue through transpiration, then re-marches
- **Local Newton**: Each station is solved independently with a local 3x3 Newton iteration
- **Coupling**: Transpiration sources are added to modify the inviscid solution, but not through a true DIJ coupling

### **Impact on Stall**
XFOIL's simultaneous system allows the BL and inviscid solutions to adjust together. When separation occurs, the increased δ* immediately feeds back to reduce Ue through the DIJ matrix, creating a self-consistent solution. RustFoil's sequential approach lags behind.

---

## 2. Inverse Mode Implementation

### XFOIL Approach (xbl.f lines 680-756)
```fortran
C---------- set prescribed Hk for inverse calculation at the current station
IF(IBL.LT.ITRAN(IS)) THEN
C----------- laminar case: relatively slow increase in Hk downstream
  HTARG = HK1 + 0.03*(X2-X1)/T1
ELSE IF(WAKE) THEN
C----------- wake: asymptotic approach to H=1
  HTARG = (Newton iteration for H + const*(H-1)^3 = H_prev)
ELSE
C----------- turbulent case: relatively fast decrease in Hk downstream
  HTARG = HK1 - 0.15*(X2-X1)/T1
ENDIF

C---------- inverse mode: solve 4x4 system for [Ctau, θ, δ*, Ue]
VS2(4,1) = 0.
VS2(4,2) = HK2_T2   ! ∂Hk/∂θ
VS2(4,3) = HK2_D2   ! ∂Hk/∂δ*
VS2(4,4) = HK2_U2   ! ∂Hk/∂Ue
VSREZ(4) = HTARG - HK2
```

**Key Points:**
1. In inverse mode, XFOIL adds Ue as a **4th unknown** to the local Newton system
2. The 4th equation is `Hk - Hk_target = 0`
3. The Jacobian includes derivatives ∂Hk/∂Ue, enabling the solver to find Ue
4. This Ue update then feeds back to the global system through the DIJ matrix

### RustFoil Approach (blsys.rs)
```rust
// In inverse mode, shape equation becomes: Hk - Hk_target = 0
let res_shape = if state2.inverse_mode {
    closures2.hk - state2.hk_target
} else {
    // Normal shape equation
};
```

**Problems:**
1. RustFoil still uses a **3x3 system** [Ctau, θ, mass] - doesn't properly add Ue as unknown
2. The Ue update in inverse mode is done **outside** the Newton system
3. No proper Jacobian coupling between Hk and Ue
4. The local Ue change doesn't feed back to a global system

---

## 3. Cl Computation

### XFOIL Approach
```fortran
CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF,...)
```
- Cl is computed from the **circulation (GAM)** around the airfoil
- The circulation is modified by the viscous solution through:
  - Wake sources (modeling the viscous wake)
  - Transpiration velocities from displacement thickness
- **Key**: The inviscid panel strengths (γ) are directly modified to account for viscous effects

### RustFoil Approach
```rust
let cl_reduction_attached = cl_correction_factor * re_factor * ds_total_te / chord;
let cl_viscous = inviscid.cl * (1.0 - cl_reduction_attached);
```
- Cl is computed as **inviscid Cl × correction factor**
- The correction is based on trailing edge displacement thickness
- No modification of panel strengths or circulation

### **Impact on Stall**
XFOIL's approach naturally captures the lift loss from separation because:
1. Separated flow increases δ* dramatically
2. This increases transpiration, reducing circulation
3. Lower circulation → lower Cl

RustFoil's empirical correction doesn't capture the physics of separation.

---

## 4. Specific Differences in Inverse Mode Triggering

### XFOIL (xbl.f line 682)
```fortran
DIRECT = HKTEST.LT.HMAX
```
The switch is based on `HKTEST` - the shape factor that **would result** from the direct mode solution.

### RustFoil
```rust
let should_inverse = prev_state.h >= hk_max || in_inverse_mode;
```
The switch is based on the **previous station's** shape factor.

### **Impact**
XFOIL tests if the current iteration would push Hk over the limit, then immediately switches. RustFoil only notices after the fact, leading to:
- Late detection of separation
- Overshooting the threshold before switching
- Less smooth transition between modes

---

## 5. Target Hk Evolution

### XFOIL Target Hk Rules
| Regime | Target Hk | Rate |
|--------|-----------|------|
| Laminar separated | HK1 + 0.03·Δs/θ | Slow increase |
| Turbulent separated | HK1 - 0.15·Δs/θ | Decrease toward reattachment |
| Wake | Newton solve for H→1 | Asymptotic |

### RustFoil
Similar formulas, but the target Hk is **clamped** to HK_MAX:
```rust
let hk_new = hk_prev + config.laminar_hk_rate * ds_over_theta;
hk_new.min(HK_MAX_LAMINAR + 1.0) // Allow slight overshoot
```

This clamping prevents Hk from growing sufficiently in separated regions.

---

## 6. Recommended Fixes (Priority Order)

### Fix 1: Implement True 4x4 Newton System in Inverse Mode
**Location**: `blsys.rs`, `solve_station_newton`

Add Ue as true 4th unknown with proper Jacobian:
```rust
if inverse_mode {
    // Build 4x4 system: [Ctau, θ, δ*, Ue]
    // 4th equation: Hk - Hk_target = 0
    // 4th row of Jacobian: [0, ∂Hk/∂θ, ∂Hk/∂δ*, ∂Hk/∂Ue]
}
```

### Fix 2: Use Predictive Inverse Mode Detection
**Location**: `march_bl_surface`

Test if the **proposed solution** would exceed threshold:
```rust
// Solve in direct mode first
let direct_result = solve_direct(...);
let hk_proposed = compute_hk(direct_result);

// Check if inverse mode needed
if hk_proposed >= hk_max {
    // Switch to inverse mode and resolve
    solve_inverse(...)
}
```

### Fix 3: Implement DIJ-Based Ue Update
**Location**: `viscous/coupling.rs`

Instead of transpiration iteration, use:
```rust
// After each BL march, update Ue from mass defect:
for i in 0..n_panels {
    ue[i] = ue_inviscid[i];
    for j in 0..n_panels {
        ue[i] += dij[i][j] * mass_defect[j];
    }
}
```

### Fix 4: Compute Cl from Modified Circulation
**Location**: `viscous/coupling.rs`

Instead of empirical correction:
```rust
// Add transpiration sources to panel system
// Re-solve for modified γ
// Compute Cl from modified γ
```

### Fix 5: Remove Hk Clamping in Inverse Mode
**Location**: `inverse.rs`

Allow Hk to grow naturally in separated regions:
```rust
// Don't clamp to HK_MAX in inverse mode
// Let physics determine equilibrium
```

---

## 7. Why Stall is Underpredicted

1. **Late Inverse Mode Activation**: RustFoil doesn't detect separation until it's already happened
2. **No True Ue Feedback**: The modified Ue in inverse mode doesn't properly feed back to reduce lift
3. **Empirical Cl Correction**: Doesn't capture the physics of separation
4. **Missing DIJ Coupling**: No mechanism for δ* to reduce Ue globally
5. **Aggressive Correction Factor**: The displacement correction reduces Cl too much at low angles, but not enough at high angles

---

## 8. Quick Experiment to Validate

To test if inverse mode is working, add logging:
```rust
if inverse_mode {
    eprintln!("INVERSE: station={}, Hk_target={:.3}, Ue_old={:.4}, Ue_new={:.4}",
        j, hk_target, ue_old, state2.ue);
}
```

If Ue is not changing significantly in inverse mode, the coupling is broken.

---

## References
- XFOIL source: `Xfoil/src/xbl.f` (MRCHUE subroutine, lines 542-840)
- XFOIL source: `Xfoil/src/xsolve.f` (BLSOLV subroutine, global Newton solver)
- XFOIL source: `Xfoil/src/xfoil.f` (CLCALC subroutine, Cl computation)
