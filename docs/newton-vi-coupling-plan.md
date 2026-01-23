# Newton Viscous-Inviscid Coupling Implementation Plan

## Status: PLANNING

## Objective
Implement XFOIL's full Newton viscous-inviscid coupling algorithm to achieve converged viscous solutions that match XFOIL's results.

## XFOIL Algorithm Overview

XFOIL's viscous-inviscid coupling uses a global Newton method where:
- **BL equations** provide 3 residuals per station (amplification/shear-lag, momentum, shape)
- **Mass defect coupling** links stations through the DIJ influence matrix
- **Edge velocity updates** propagate mass defect changes to affect all downstream stations

### Key Subroutine Roles

| Subroutine | File | Role |
|------------|------|------|
| VISCAL | xoper.f:2886 | Main iteration loop |
| SETBL | xbl.f:21 | Build Newton system with VM matrix |
| BLSYS | xblsys.f:583 | Build local 3x5 Jacobian |
| BLSOLV | xsolve.f:283 | Solve block system with VM coupling |
| UPDATE | xbl.f:1253 | Apply solution with limiting |
| UESET | xpanel.f:1758 | Update Ue from mass defect via DIJ |
| MRCHDU | xbl.f:851 | Direct march with current Ue and δ* |

### Algorithm Flow

```
VISCAL:
  1. Initialize: XYWAKE, QWCALC, QISET, IBLPAN, IBLSYS
  2. Set initial Ue: UEDG = UINV
  3. Build DIJ matrix if needed (QDCALC)
  4. For iteration = 1 to max:
      a. SETBL: Build Newton system
      b. BLSOLV: Solve block-tridiagonal + VM
      c. UPDATE: Apply solution with limiting
      d. QVFUE, GAMQV: Update Qvis, GAM
      e. STMOVE: Relocate stagnation point
      f. CLCALC, CDCALC: Update forces
      g. Check convergence (RMSBL < 1e-4)
```

---

## Gap Analysis: Current RustFoil vs XFOIL

### What We Have ✓

| Component | Status | Location |
|-----------|--------|----------|
| `solve_viscous()` | ✓ Exists | `viscal.rs:88` |
| `BlNewtonSystem` | ✓ Exists | `newton.rs:48` |
| `build_with_vm()` | ✓ Partial | `newton.rs:196` |
| `solve_coupled_system()` | ✓ Exists | `solve.rs:131` |
| `update_stations()` | ✓ Exists | `update.rs:95` |
| `set_edge_velocities()` | ✓ Exists | `update.rs:367` |
| `march_fixed_ue()` | ✓ Exists | `march.rs` |
| DIJ matrix computation | ✓ Exists | `rustfoil-coupling/src/dij.rs` |

### What's Missing ✗

| Component | Issue | XFOIL Reference |
|-----------|-------|-----------------|
| **VTI Sign Handling** | No surface direction signs | xbl.f:206 |
| **Forced Changes** | No DUE, DDS computation | xbl.f:216-217 |
| **VDEL Column 2** | No α/Re sensitivity | xbl.f:354-357 |
| **XI Sensitivities** | No stagnation point motion | xbl.f:317-322 |
| **Two-Surface VM** | Single surface only | xbl.f:123-132 |
| **VZ Block** | No wake merging | xbl.f:234-268 |
| **MRCHDU** | No direct march preserving transition | xbl.f:851 |
| **Full BLSOLV** | Simplified VM back-substitution | xsolve.f:283-485 |

---

## Implementation Plan

### Phase 1: Single-Surface Newton Coupling (Core Fix)

**Goal**: Get `solve_viscous()` working correctly for a single surface.

#### 1.1 Add VTI Sign Handling
```rust
// In newton.rs
pub struct BlNewtonSystem {
    // ... existing fields ...
    /// Sign array for surface direction (+1 upper, -1 lower)
    pub vti: Vec<f64>,
}
```

The VTI array tracks the sign convention for each station. XFOIL uses:
- VTI = +1 for upper surface (flow goes from LE toward TE)
- VTI = -1 for lower surface (flow goes from LE toward TE, but velocity sign flipped)

**Files to modify**:
- `crates/rustfoil-coupling/src/newton.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`

#### 1.2 Implement Forced Changes (DUE, DDS)

XFOIL computes "forced changes" from the mismatch between current UEDG and what UESET would give:

```fortran
C---- "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
      DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
      DDS2 = D2_U2*DUE2
```

This represents the non-linear update needed because UEDG has been modified outside the Newton linearization.

**Implementation**:
```rust
// In newton.rs build_with_vm()
pub fn build_with_vm_and_forced(
    &mut self,
    stations: &[BlStation],
    flow_types: &[FlowType],
    msq: f64,
    re: f64,
    dij: &DMatrix<f64>,
    ue_inviscid: &[f64],  // NEW: UINV
    vti: &[f64],          // NEW: VTI
) {
    // Save current Ue (USAV)
    let ue_save: Vec<f64> = stations.iter().map(|s| s.u).collect();
    
    // Compute what UESET would give
    let ue_computed = compute_ueset(stations, ue_inviscid, dij, vti);
    
    // Compute forced changes
    let due: Vec<f64> = ue_save.iter().zip(&ue_computed)
        .map(|(save, comp)| save - comp).collect();
    
    // ... use DUE in residual construction ...
}
```

**Files to modify**:
- `crates/rustfoil-coupling/src/newton.rs`

#### 1.3 Fix VM Matrix Construction

Current implementation in `build_with_vm()` is missing:

1. **U_M computation with VTI**:
   ```rust
   // XFOIL: U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
   let u2_m_j = -vti[i] * vti[j] * dij[(i, j)];
   ```

2. **VS1 contributions** (upstream station derivatives):
   ```rust
   // Current only uses VS2, need VS1 too
   vm[i][j][k] = vs2_d[k] * d2_m_j + vs2_u[k] * u2_m_j
               + vs1_d[k] * d1_m_j + vs1_u[k] * u1_m_j;
   ```

**Files to modify**:
- `crates/rustfoil-coupling/src/newton.rs`

#### 1.4 Implement Full BLSOLV Algorithm

Current `solve_coupled_system()` has a simplified back-substitution. XFOIL's BLSOLV:

1. **Forward sweep**: Eliminates VB blocks AND propagates VM through the factorization
2. **Back substitution**: Uses the third equation's delta (mass change) to update upstream stations

Key insight: XFOIL solves for `[dCtau, dθ, d(mass)]` not `[dCtau, dθ, dδ*]`. The mass defect is the coupling variable.

```rust
// In solve.rs - Key fix: solution variable is mass defect change
pub fn solve_blsolv_full(system: &BlNewtonSystem) -> Vec<[f64; 3]> {
    // ... forward sweep ...
    
    // Back substitution with VM coupling
    for iv in (1..n).rev() {
        // Third variable (mass) affects all upstream equations through VM
        let mass_change = solution[iv][2];
        for kv in (1..iv).rev() {
            for k in 0..3 {
                solution[kv][k] -= system.vm[kv][iv][k] * mass_change;
            }
        }
    }
    
    // ... 
}
```

**Files to modify**:
- `crates/rustfoil-coupling/src/solve.rs`

#### 1.5 Fix UPDATE to Match XFOIL

Current `update_stations()` applies deltas directly. XFOIL's UPDATE:

1. Computes **new Ue first** from mass defect changes
2. Computes **new CL** from new Ue to get α/CL sensitivity
3. Determines **global relaxation** from all changes
4. **Applies all updates** with common relaxation factor

```rust
// In update.rs
pub fn update_xfoil_style(
    stations: &mut [BlStation],
    deltas: &[[f64; 3]],
    ue_inviscid: &[f64],
    dij: &DMatrix<f64>,
    vti: &[f64],
    config: &UpdateConfig,
) -> UpdateResult {
    // Step 1: Compute new Ue from mass deltas BEFORE updating stations
    let new_masses: Vec<f64> = stations.iter().zip(deltas.iter())
        .map(|(s, d)| s.mass_defect + d[2]).collect();
    
    let new_ue = compute_ue_from_mass(ue_inviscid, &new_masses, dij, vti);
    
    // Step 2: Compute relaxation factor
    let rlx = compute_relaxation(&stations, &deltas, &new_ue, config);
    
    // Step 3: Apply all updates with same relaxation
    for (i, (station, delta)) in stations.iter_mut().zip(deltas.iter()).enumerate() {
        station.ctau += rlx * delta[0];
        station.theta += rlx * delta[1];
        station.mass_defect += rlx * delta[2];
        station.u = ue_inviscid[i] + (new_ue[i] - ue_inviscid[i]) * rlx;
        station.delta_star = station.mass_defect / station.u;
    }
    // ...
}
```

**Files to modify**:
- `crates/rustfoil-coupling/src/update.rs`

---

### Phase 2: Two-Surface Coupling

**Goal**: Handle upper + lower surfaces with proper wake merging.

#### 2.1 Two-Surface System Structure

XFOIL treats both surfaces in one Newton system:
- Stations 1..N1 = Upper surface (LE to TE)
- Stations N1+1..N1+N2 = Lower surface (LE to TE)
- Stations N1+N2+1..NSYS = Wake

The VM matrix couples ALL stations:
```
       Upper   Lower   Wake
Upper  [VM11]  [VM12]  [VM13]
Lower  [VM21]  [VM22]  [VM23]
Wake   [VM31]  [VM32]  [VM33]
```

#### 2.2 Implement VZ Block for Wake

At the wake start, XFOIL uses a special VZ block to couple TE stations from both surfaces:

```fortran
C----- at wake start, D1 depends on both TE δ* values
       D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
       D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
```

**Implementation**:
```rust
// In newton.rs
pub struct BlNewtonSystem {
    // ... existing fields ...
    /// VZ block for wake coupling (3x2 matrix)
    pub vz: Option<[[f64; 2]; 3]>,
    /// Index of upper surface TE station
    pub ivte1: usize,
    /// Index of lower surface TE station  
    pub ivte2: usize,
}
```

#### 2.3 BLSOLV Two-Surface Extension

BLSOLV has special handling at the wake start:
```fortran
        IF(IV.EQ.IVTE1) THEN
C------- eliminate VZ block
         IVZ = ISYS(IBLTE(2)+1,2)
         DO K=1, 3
           VDEL(K,1,IVZ) = VDEL(K,1,IVZ)
     &         - VZ(K,1)*VDEL(1,1,IV) - VZ(K,2)*VDEL(2,1,IV)
         ENDDO
        ENDIF
```

**Files to modify**:
- `crates/rustfoil-coupling/src/newton.rs`
- `crates/rustfoil-coupling/src/solve.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`

---

### Phase 3: Stagnation Point Motion

**Goal**: Handle XI (arc length) sensitivities for stagnation point motion.

#### 3.1 Stagnation Point Tracking

XFOIL tracks how the stagnation point location depends on LE velocities:
```fortran
C---- set XI sensitivities wrt LE Ue changes
      IF(IS.EQ.1) THEN
       XI_ULE1 =  SST_GO
       XI_ULE2 = -SST_GP
      ELSE
       XI_ULE1 = -SST_GO
       XI_ULE2 =  SST_GP
      ENDIF
```

This affects the VM matrix through VSX terms:
```fortran
     &              + (VS1(1,5) + VS2(1,5) + VSX(1))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
```

#### 3.2 Implementation

```rust
// In newton.rs
pub struct BlNewtonSystem {
    // ... existing fields ...
    /// Arc length sensitivity coefficients
    pub xi_ule1: f64,  // ∂ξ/∂Ue_LE1
    pub xi_ule2: f64,  // ∂ξ/∂Ue_LE2
    /// LE velocity sensitivities to mass defect
    pub ule1_m: Vec<f64>,  // ∂Ue_LE1/∂mass_j
    pub ule2_m: Vec<f64>,  // ∂Ue_LE2/∂mass_j
}
```

**Files to modify**:
- `crates/rustfoil-coupling/src/newton.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`

---

### Phase 4: Alpha/Reynolds Sensitivity

**Goal**: Support fixed-CL mode where α varies.

#### 4.1 VDEL Column 2

XFOIL stores α/Re sensitivity in VDEL column 2:
```fortran
      IF(LALFA) THEN
       VDEL(1,2,IV) = VSR(1)*RE_CLMR + VSM(1)*MSQ_CLMR
      ELSE
       VDEL(1,2,IV) = VS1(1,4)*U1_A + VS2(1,4)*U2_A + ...
      ENDIF
```

This allows solving for Δα along with BL variables when operating in fixed-CL mode.

**Implementation** (lower priority - fixed-α mode works first):
```rust
// In newton.rs
pub struct BlNewtonSystem {
    // ... existing fields ...
    /// Sensitivity to alpha (VDEL column 2)
    pub vdel_alpha: Vec<[f64; 3]>,
}
```

---

## Testing Strategy

### Unit Tests

1. **VM Matrix Construction**
   - Test against XFOIL SETBL output for NACA 0012 at α=0°
   - Verify VM[i][j][k] values at key stations

2. **BLSOLV**
   - Test with known solution (manufactured problem)
   - Verify back-substitution correctness

3. **UPDATE**
   - Test relaxation factor computation
   - Test mass defect → Ue → δ* chain

### Integration Tests

1. **Single Surface Flat Plate**
   - Match XFOIL θ, δ*, Cf at Re=10⁶
   - Verify transition location

2. **NACA 0012 Symmetric**
   - Match XFOIL CL, CD at α=0°, Re=10⁶
   - Verify both surfaces converge

3. **NACA 0012 α=5°**
   - Match CL, CD, Cm
   - Verify asymmetric transition

### Validation Data

Use existing `testdata/mrchue_iterations.json` and create:
- `testdata/setbl_vm_matrix.json` - VM entries at each station
- `testdata/blsolv_solution.json` - BLSOLV solution vectors
- `testdata/update_states.json` - State after each UPDATE call

---

## Success Criteria

| Metric | Target |
|--------|--------|
| CL error vs XFOIL | < 1% |
| CD error vs XFOIL | < 5% |
| Transition location | ±1 station |
| Convergence | < 20 iterations |
| Stability | No divergence for |α| < 10° |

---

## File Change Summary

| File | Changes |
|------|---------|
| `newton.rs` | Add VTI, forced changes, XI sensitivities, VZ block |
| `solve.rs` | Full BLSOLV with two-surface VM handling |
| `update.rs` | XFOIL-style UPDATE with mass-first Ue computation |
| `viscal.rs` | Two-surface iteration loop, wake handling |
| `setup.rs` | VTI initialization, stagnation point sensitivities |

---

## References

- XFOIL Source: `Xfoil/src/xbl.f`, `xblsys.f`, `xsolve.f`, `xoper.f`, `xpanel.f`
- Drela, M., "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
- Drela, M., "Two-Dimensional Transonic Aerodynamic Design and Analysis Using the Euler Equations"
