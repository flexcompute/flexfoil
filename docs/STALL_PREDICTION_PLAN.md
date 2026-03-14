# Stall Prediction Plan: Matching XFOIL Cl-alpha Behavior

## Problem Statement

RustFoil's Cl-alpha curves don't show stall - lift keeps increasing linearly while XFOIL shows Cl_max followed by lift reduction. This is because:

1. **Current behavior**: RustFoil uses inviscid Cl with a small δ*-based correction
2. **Required behavior**: Stall should emerge naturally from viscous-inviscid coupling when separation occurs

## Root Cause Analysis

XFOIL's stall prediction works through **physics-based coupling**, not an explicit stall model:

```
Separation → Large δ* → Large transpiration Vn → Reduced Ue → Reduced suction peak → Lower Cl
```

### What XFOIL Does (from xbl.f)

1. **Separation thresholds** (lines 556-557):
   ```fortran
   HLMAX = 3.8   ! laminar Hk threshold
   HTMAX = 2.5   ! turbulent Hk threshold
   ```

2. **Mode switching** (lines 680-682):
   ```fortran
   IF(IBL.LT.ITRAN(IS)) HMAX = HLMAX
   IF(IBL.GE.ITRAN(IS)) HMAX = HTMAX
   DIRECT = HKTEST.LT.HMAX
   ```

3. **Inverse mode** (lines 691-728): When Hk exceeds threshold:
   - Prescribe target Hk (HTARG) that changes slowly downstream
   - Solve for Ue that gives this Hk (instead of solving for δ* with fixed Ue)
   - The modified Ue feeds back to the global Newton system

4. **Ue feedback to lift**: The modified Ue distribution from the BL solver is used in the next iteration of the panel method, naturally reducing the suction peak and lift.

### What RustFoil Is Missing

1. **No inverse mode**: Always runs in direct mode (fixed Ue, solve for δ*)
2. **Clamped Hk**: Shape factor is clamped to prevent blow-up, but this prevents proper δ* growth
3. **Weak Ue feedback**: Transpiration may not be strongly enough coupled to modify Ue
4. **Cl from wrong Ue**: May be using inviscid Ue for Cl calculation instead of viscous-modified Ue

---

## Implementation Plan

### Phase 1: Diagnostics (Understand Current State)

**Goal**: Identify exactly where the physics breaks down

#### 1.1 Add BL state logging at high alpha
- Log Hk, δ*, θ, Ue along upper surface at α = 10°, 12°, 14°
- Compare with XFOIL values at same conditions
- Identify: Is Hk reaching threshold? Is δ* growing?

#### 1.2 Compare transpiration velocities
- Extract Vn = d(Ue·δ*)/ds from both solvers
- Check if RustFoil's transpiration is comparable to XFOIL

#### 1.3 Compare Ue distributions
- At high alpha, XFOIL's Ue should be significantly modified from inviscid
- Check if RustFoil's Ue is being modified at all

**Deliverable**: Diagnostic comparison showing where divergence occurs

---

### Phase 2: Implement Inverse Mode Switching

**Goal**: Match XFOIL's separation handling

#### 2.1 Add separation detection
```rust
// In boundary layer marching
const HK_MAX_LAMINAR: f64 = 3.8;
const HK_MAX_TURBULENT: f64 = 2.5;

fn should_use_inverse_mode(hk: f64, is_turbulent: bool) -> bool {
    let hk_max = if is_turbulent { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
    hk >= hk_max
}
```

#### 2.2 Implement inverse mode BL equations
In inverse mode:
- **Given**: Target Hk (prescribed, limited)
- **Solve for**: Ue (edge velocity) instead of δ*

The target Hk evolution (from XFOIL):
```rust
fn compute_target_hk(
    hk_prev: f64,
    x_prev: f64,
    x_curr: f64,
    theta_prev: f64,
    is_turbulent: bool,
    is_wake: bool,
) -> f64 {
    let dx = x_curr - x_prev;
    
    if !is_turbulent {
        // Laminar: slow increase
        (hk_prev + 0.03 * dx / theta_prev).max(HK_MAX_LAMINAR)
    } else if is_wake {
        // Wake: asymptotic approach to H=1
        // Newton iterations for: H + const*(H-1)^3 = H_prev
        let const_val = 0.03 * dx / theta_prev;
        let mut h = hk_prev;
        for _ in 0..3 {
            h = h - (h + const_val * (h - 1.0).powi(3) - hk_prev)
                  / (1.0 + 3.0 * const_val * (h - 1.0).powi(2));
        }
        h.max(1.01)
    } else {
        // Turbulent: faster decrease toward attached
        (hk_prev - 0.15 * dx / theta_prev).max(HK_MAX_TURBULENT)
    }
}
```

#### 2.3 Add inverse mode to Newton system
Modify the BL residual equations:
- Direct mode: R_shape = f(θ, δ*, Ue) - standard shape equation
- Inverse mode: R_shape = Hk - Hk_target

**Deliverable**: BL solver that switches between direct/inverse modes

---

### Phase 3: Strengthen Ue Feedback

**Goal**: Ensure BL solution properly modifies edge velocity

#### 3.1 Review transpiration velocity calculation
Current implementation in `coupling.rs`:
```rust
// Vn = d(Ue·δ*)/ds
let vn = compute_transpiration_velocity(ue, delta_star, s_coords);
```
Verify this matches XFOIL's formulation.

#### 3.2 Ensure Ue is updated from BL solution
In inverse mode, the BL solver outputs a modified Ue. This must:
1. Be stored in the solution state
2. Be used in the next V-I iteration
3. Be used for lift calculation

#### 3.3 Implement proper Ue relaxation
From XFOIL (line 756):
```fortran
UEI = UEI + RLX*VSREZ(4)  ! Update Ue with relaxation
```
Ensure we're doing similar updates.

**Deliverable**: Proper Ue modification from BL in separated regions

---

### Phase 4: Fix Lift Calculation

**Goal**: Use viscous-modified Ue for lift

#### 4.1 Current problem
Lift is calculated from inviscid panel method Cp, not from the viscous-modified Ue distribution.

#### 4.2 Solution: Recompute Cp from modified Ue
After V-I convergence:
```rust
// Use the final viscous Ue to compute Cp
let cp_viscous: Vec<f64> = ue_final.iter()
    .map(|&ue| 1.0 - ue.powi(2))
    .collect();

// Integrate for lift
let cl_viscous = integrate_cp_for_lift(&cp_viscous, &panels);
```

#### 4.3 Alternative: Proper coupled lift
In the Newton iteration, the lift should come from the modified circulation distribution that's consistent with the BL solution.

**Deliverable**: Cl computed from viscous-modified velocity distribution

---

### Phase 5: Validation & Tuning

**Goal**: Match XFOIL Cl-alpha curves

#### 5.1 Test cases
- NACA 0012 at Re = 3M: α = -4° to 16°
- NACA 2412 at Re = 3M: α = -4° to 16°  
- NACA 4412 at Re = 3M: α = -4° to 16°

#### 5.2 Metrics
- Cl-alpha slope in linear region (should match to ~2%)
- Cl_max location (alpha at max lift, within 0.5°)
- Cl_max value (within 5%)
- Post-stall Cl reduction (qualitatively correct shape)

#### 5.3 Tuning parameters
- Hk thresholds (HTMAX, HLMAX)
- Target Hk evolution rates (0.03, 0.15 factors)
- Relaxation factors for inverse mode

**Deliverable**: Validated Cl-alpha curves matching XFOIL

---

## File Changes Summary

| File | Changes |
|------|---------|
| `boundary_layer/mod.rs` | Add inverse mode detection and target Hk calculation |
| `boundary_layer/state.rs` | Add `inverse_mode: bool` flag to state |
| `viscous/blsys.rs` | Add inverse mode residual equations |
| `viscous/coupling.rs` | Implement mode switching, Ue feedback, modified Cl calculation |
| `viscous/newton.rs` | Support inverse mode in Newton system |
| New: `boundary_layer/inverse.rs` | Inverse mode BL equations |

---

## Implementation Order

1. **Phase 1** (Diagnostics) - 1 session
   - Essential to understand current state before making changes
   
2. **Phase 2** (Inverse Mode) - 2-3 sessions
   - Core algorithm change, most complex
   
3. **Phase 3** (Ue Feedback) - 1 session
   - May be partially working, needs verification
   
4. **Phase 4** (Lift Calculation) - 1 session
   - Straightforward once phases 2-3 work
   
5. **Phase 5** (Validation) - 1 session
   - Tuning and verification

---

## Success Criteria

### Minimum Viable
- [ ] Cl_max exists (curve stops increasing at some alpha)
- [ ] Post-stall Cl decreases

### Target
- [ ] Cl_max within 5% of XFOIL
- [ ] Alpha at Cl_max within 1°
- [ ] Cl-alpha slope matches to 3%
- [ ] Post-stall shape qualitatively correct

### Stretch
- [ ] Cl within 2% of XFOIL across full alpha range
- [ ] Works for cambered airfoils (2412, 4412)
- [ ] Works across Re range (1M to 6M)

---

## References

- XFOIL source: `Xfoil/src/xbl.f` (lines 542-840, MRCHUE subroutine)
- XFOIL source: `Xfoil/src/xblsys.f` (BL equations)
- XFOIL source: `Xfoil/src/xsolve.f` (Newton system)
- Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
