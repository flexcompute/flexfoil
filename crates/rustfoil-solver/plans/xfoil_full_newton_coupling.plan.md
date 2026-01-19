# XFOIL Full Newton Coupling Implementation Plan

## Current Status

**Problem**: Transpiration-as-sources gives only ~4% Cl reduction vs XFOIL's ~11%

**Root Cause**: Source panels create weak velocity perturbation, not equivalent to true geometric displacement

## XFOIL's Approach

XFOIL solves the coupled V-I system simultaneously using a global Newton method:

```
[dR_BL/d(BL)   dR_BL/dUe  ] [Δ(BL)] = -[R_BL]
[dUe/dm        I          ] [ΔUe  ]    [R_Ue ]
```

Where:
- **BL variables**: Ctau (or N for laminar), θ, m = Ue·δ* (mass defect)
- **R_BL**: BL equation residuals (shear lag, momentum, shape)
- **dUe/dm**: Inviscid influence - how mass defect affects edge velocity via sources

### Key Data Structures (from XFOIL xbl.f, xsolve.f)

```fortran
VA(3,3,IV)    ! "A" block - BL Jacobian at station IV
VB(3,3,IV)    ! "B" block - BL coupling to upstream station
VZ(3,2)       ! "Z" block - wake/TE coupling
VM(3,IV,JV)   ! Mass influence - dRes(IV)/dm(JV) via dUe/dm
VDEL(3,2,IV)  ! RHS: residual + alpha sensitivity
DIJ(I,J)      ! Source influence: dGamma(I)/dSigma(J)
```

### The Magic: DIJ Matrix

XFOIL precomputes `DIJ(I,J)` = sensitivity of surface velocity at node I to source strength at panel J.

Then the V-I coupling is:
```
dUe(I)/dm(J) = -VTI(I) * VTI(J) * DIJ(panel_I, panel_J)
```

Where `VTI` converts between BL stations and panel indices.

## Implementation Phases

### Phase 1: Source Influence Matrix (DIJ)

**Goal**: Compute dGamma/dSigma for all node pairs

**Location**: `crates/rustfoil-solver/src/inviscid/influence.rs`

**Implementation**:
```rust
/// Source influence matrix: dGamma(i)/dSigma(j)
/// 
/// For each panel j with source strength sigma_j, compute the 
/// tangential velocity perturbation at each node i.
pub fn compute_source_influence_matrix(nodes: &[Point]) -> Array2<f64> {
    let n = nodes.len();
    let mut dij = Array2::zeros((n, n));
    
    for i in 0..n {
        for j in 0..n {
            // Source panel influence from XFOIL PSILIN
            // dU_tangent / dSigma = (1/2π) * ln(r2/r1)
            dij[[i, j]] = source_tangent_influence(nodes, i, j);
        }
    }
    dij
}
```

This is similar to `compute_source_velocity_correction` but precomputed as a matrix.

### Phase 2: Mass Influence in Newton System

**Goal**: Add dUe/dm coupling to the block-tridiagonal solver

**Location**: `crates/rustfoil-solver/src/viscous/newton.rs`

**Changes to `BlockTridiagJacobian`**:
```rust
pub struct BlockTridiagJacobian {
    pub a_blocks: Vec<BLBlock>,     // Diagonal blocks
    pub b_blocks: Vec<BLBlock>,     // Sub-diagonal blocks  
    pub z_block: Option<BLBlock>,   // TE/wake coupling
    
    // NEW: Mass influence matrix
    pub mass_influence: Array2<f64>,  // dUe(i)/dm(j) = -VTI(i)*VTI(j)*DIJ
    
    pub rhs: Vec<[f64; 3]>,
    pub solution: Vec<[f64; 3]>,
}
```

**Solver modification**: Instead of tridiagonal solve, use XFOIL's BLSOLV algorithm:
1. Forward elimination of A, B blocks
2. Accumulate mass influence terms into modified RHS
3. Back-substitute with mass coupling

### Phase 3: Update BL Equations to Use Mass Defect

**Goal**: Change third equation variable from δ* to m = Ue·δ*

**Location**: `crates/rustfoil-solver/src/viscous/blsys.rs`

**Current**: Solve for (Ctau, θ, δ*)  
**New**: Solve for (Ctau, θ, m) where m = Ue·δ*

Benefits:
- Mass defect m is conserved in the wake
- Derivatives dR/dm are more stable than dR/dδ*
- Directly couples to inviscid through source distribution

### Phase 4: Newton Iteration with Full Coupling

**Goal**: Global Newton iteration solving BL + inviscid simultaneously

**Location**: `crates/rustfoil-solver/src/viscous/coupling.rs`

**Algorithm**:
```rust
fn full_newton_coupling(&self, ...) -> ViscousSolution {
    // 1. Initial inviscid solve
    let mut ue = solve_inviscid();
    let dij = compute_source_influence_matrix(&nodes);
    
    for iter in 0..max_iter {
        // 2. Compute BL state with current Ue
        let bl_state = march_bl(&ue);
        
        // 3. Build Newton system
        let mut jacobian = build_block_jacobian(&bl_state, &ue, &dij);
        
        // 4. Add mass influence coupling
        jacobian.set_mass_influence(&dij, &vti_mapping);
        
        // 5. Solve Newton system (XFOIL BLSOLV style)
        jacobian.solve_with_mass_coupling();
        
        // 6. Update BL variables
        update_bl_from_solution(&mut bl_state, &jacobian.solution);
        
        // 7. Update Ue from mass defect changes
        let dm = extract_mass_defect_changes(&jacobian.solution);
        ue = update_ue_from_mass(&ue, &dm, &dij);
        
        // 8. Check convergence
        if jacobian.residual_norm() < tol {
            break;
        }
    }
    
    build_solution(...)
}
```

## Expected Accuracy Improvement

| Metric | Current | After Full Newton |
|--------|---------|-------------------|
| Cl reduction from viscous | 1% | ~11% |
| Cl error vs XFOIL | 11% | <3% |
| Cd error at moderate α | 6-13% | <5% |

## Implementation Order

1. **Phase 1**: Source influence matrix (~2-3 hours)
   - Add `compute_source_influence_matrix` to inviscid module
   - Test against existing `compute_source_velocity_correction`

2. **Phase 2**: Mass influence in Newton system (~3-4 hours)
   - Extend `BlockTridiagJacobian` with mass influence
   - Implement XFOIL BLSOLV algorithm

3. **Phase 3**: Mass defect equations (~2-3 hours)
   - Change third variable from δ* to m
   - Update residual and Jacobian calculations

4. **Phase 4**: Full Newton coupling (~3-4 hours)
   - New coupling method `CouplingMethod::FullNewton`
   - Integrate all components

## Files to Modify

1. `crates/rustfoil-solver/src/inviscid/influence.rs` (new file)
2. `crates/rustfoil-solver/src/viscous/newton.rs`
3. `crates/rustfoil-solver/src/viscous/blsys.rs`
4. `crates/rustfoil-solver/src/viscous/coupling.rs`
5. `crates/rustfoil-solver/src/viscous/mod.rs`

## References

- XFOIL source: `xbl.f` (SETBL), `xsolve.f` (BLSOLV), `xpanel.f` (DIJ)
- Drela thesis: "Two-Dimensional Transonic Aerodynamic Design and Analysis Using the Euler Equations" (1987)
