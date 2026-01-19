# Viscous Solver Accuracy Improvements Plan

## Current Status

After the transpiration sign fix, results have improved significantly:

| Alpha | Cl Error | Cd Error (before) | Cd Error (after) |
|-------|----------|-------------------|------------------|
| 0°    | 0.2%     | 42%               | 23%              |
| 2°    | 12%      | 69%               | 6%               |
| 4°    | 12%      | 128%              | 12%              |
| 6°    | 3%       | 177%              | 6%               |
| 8°    | 5%       | 302%              | 7%               |

**Remaining Issues:**
1. Cl error ~10-12% at low angles (likely inviscid Cl slope issue)
2. Cd error ~20-30% at α=0° (wake/drag integration issue)
3. No proper wake coupling with inviscid solver

---

## Phase 1: Proper Wake-Based Drag Calculation

**Goal:** Use the existing `wake.rs` infrastructure properly for drag computation.

### 1.1 Use `compute_squire_young_drag` from wake.rs

Currently `coupling.rs` has inline Squire-Young calculation that doesn't use `ue_te` properly.

**File:** `crates/rustfoil-solver/src/viscous/coupling.rs`

**Changes:**
```rust
// Replace inline Squire-Young with wake module
use crate::boundary_layer::wake::compute_squire_young_drag;
use crate::boundary_layer::state::BLState;

// In build_bl_solution or equivalent:
let upper_te = BLState {
    theta: theta_upper.last().copied().unwrap_or(0.001),
    delta_star: delta_star_upper.last().copied().unwrap_or(0.002),
    ue: ue_upper_te,  // Actual edge velocity at TE
    h: h_upper.last().copied().unwrap_or(1.4),
    ..Default::default()
};
let lower_te = BLState { /* similar */ };

let (cd_total, cd_friction, cd_pressure) = compute_squire_young_drag(&upper_te, &lower_te, chord);
```

### 1.2 Extract Proper TE Edge Velocities

The Squire-Young formula depends critically on `(Ue_TE / U∞)^exponent`. Currently using `1.0` which is wrong for cambered airfoils or high angles.

**Changes:**
- Extract `ue_te` from the last BL station on each surface
- Ensure it's properly normalized to freestream

---

## Phase 2: Friction Drag Integration

**Goal:** Compute friction drag by proper surface integration, not approximation.

### 2.1 XFOIL-Style Friction Drag Integration

From `VISCOUS_MODELS.md`:
```
τ = 0.5 × ρ × Uₑ² × Cf
Cd_f = Σ 0.5×(τᵢ + τᵢ₋₁) × Δx × 2/U∞²
```

where Δx is projected onto freestream direction:
```
Δx = (xᵢ - xᵢ₋₁)×cos(α) + (yᵢ - yᵢ₋₁)×sin(α)
```

**File:** `crates/rustfoil-solver/src/viscous/coupling.rs` → `compute_friction_drag`

**Changes:**
```rust
fn compute_friction_drag_xfoil(
    cf: &[f64],
    ue: &[f64],
    x: &[f64],
    y: &[f64],
    alpha: f64,
    chord: f64,
) -> f64 {
    let cos_a = alpha.cos();
    let sin_a = alpha.sin();
    let mut cd_f = 0.0;
    
    for i in 1..cf.len() {
        // Shear stress (normalized)
        let tau_avg = 0.5 * (cf[i] * ue[i].powi(2) + cf[i-1] * ue[i-1].powi(2));
        
        // Projected distance in freestream direction
        let dx = (x[i] - x[i-1]) * cos_a + (y[i] - y[i-1]) * sin_a;
        
        cd_f += tau_avg * dx.abs();
    }
    
    cd_f / chord
}
```

---

## Phase 3: Global Newton-Raphson Coupling (Optional - Higher Complexity)

**Goal:** Replace semi-direct iteration with true simultaneous solution.

### 3.1 Current State

- Block-tridiagonal structures exist in `newton.rs`
- `blsys.rs` has local Newton solver and Jacobian computation
- These are not fully integrated into the VII loop

### 3.2 Full Newton System

The XFOIL approach solves the complete system:
```
[dR_BL/dBL    dR_BL/dUe  ] [ΔBL] = -[R_BL]
[dR_inv/dBL  dR_inv/dUe ] [ΔUe]    [R_inv]
```

**Required Components:**
1. **Surface BL Jacobian** - Already have via `compute_interval_jacobian`
2. **BL-to-Inviscid coupling** - dR_inv/dδ* (transpiration influence)
3. **Inviscid-to-BL coupling** - dR_BL/dUe (edge velocity influence)

### 3.3 Simplified Approach: Lagged Newton

Instead of full simultaneous solve:
1. March BL with current Ue → get residuals and Jacobian
2. Solve block-tridiagonal system for BL corrections
3. Update BL variables
4. Recompute transpiration
5. Re-solve inviscid for new Ue
6. Check convergence

This is what XFOIL actually does (not truly simultaneous).

---

## Phase 4: Cl Accuracy Improvement

**Goal:** Reduce Cl error from ~12% to <5%.

### 4.1 Root Cause Analysis

The Cl error is likely from:
1. Panel count (160 may be insufficient)
2. LE/TE clustering
3. Transpiration velocity magnitude

### 4.2 Panel Refinement Study

Test with 200, 240, 320 panels to see if Cl converges to XFOIL.

### 4.3 Displacement Effect

Verify that the transpiration-based Cl reduction matches expected:
- Cl_viscous ≈ Cl_inviscid × (1 - 2δ*/c) for thin airfoils

---

## Implementation Priority

| Task | Impact | Complexity | Priority |
|------|--------|------------|----------|
| 1.1 Use wake.rs Squire-Young | High | Low | **1** |
| 1.2 Extract proper Ue_TE | High | Low | **2** |
| 2.1 XFOIL friction integration | Medium | Low | **3** |
| 4.2 Panel count study | Medium | Low | **4** |
| 3.3 Lagged Newton | Medium | High | 5 |

---

## Verification

After each change, run:
```bash
cargo test --package rustfoil-solver --test xfoil_comparison -- --nocapture
```

Target metrics:
- Cd error < 15% across all angles
- Cl error < 8% across all angles
- Transition location within ±5% of XFOIL

---

## Files to Modify

1. `crates/rustfoil-solver/src/viscous/coupling.rs`
   - Use `wake::compute_squire_young_drag`
   - Extract proper Ue_TE values
   - Update friction drag calculation

2. `crates/rustfoil-solver/src/boundary_layer/wake.rs`
   - Already complete, just needs to be used

3. `crates/rustfoil-solver/tests/xfoil_comparison.rs`
   - Add panel count sensitivity test
