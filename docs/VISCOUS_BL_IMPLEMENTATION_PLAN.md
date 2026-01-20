# Viscous Boundary Layer Solver Implementation Plan

## Quick Start: Copy-Paste Prompts

Use these prompts in new Cursor chats to execute each wave of tasks.

### Wave 1: Test Infrastructure

```
@docs/tasks/TASK_01_TESTKIT.md

Implement this task completely. Create the rustfoil-testkit crate with:
1. Cargo.toml with serde, serde_json, tempfile dependencies
2. src/lib.rs exporting modules
3. src/fortran_runner.rs for compiling/running FORTRAN
4. src/approx.rs for float comparison
5. fortran/Makefile for building test harnesses
6. fortran/test_closures.f with HKIN test generation

Update the workspace Cargo.toml to include the new crate.
Verify with: cd crates/rustfoil-testkit/fortran && make
```

### Wave 2: Constants

```
@docs/tasks/TASK_02_BL_CONSTANTS.md

Implement this task completely. Create the rustfoil-bl crate with:
1. Cargo.toml depending on rustfoil-core
2. src/lib.rs
3. src/constants.rs with all BLPAR.INC constants (SCCON, GACON, GBCON, etc.)
4. src/closures/mod.rs (empty, ready for closure functions)

Update workspace Cargo.toml. Verify with: cargo build -p rustfoil-bl
```

### Wave 3: Closures (Run these 7 prompts in PARALLEL chats)

**Chat 3A - HKIN:**
```
@docs/tasks/TASK_03_HKIN.md

Implement the HKIN closure function in rustfoil-bl/src/closures/hkin.rs.
Port XFOIL's xblsys.f line 2276 exactly. Include:
1. HkinResult struct with hk, hk_h, hk_msq
2. hkin(h, msq) function
3. Unit tests with numerical derivative verification
4. Update closures/mod.rs to export

Verify with: cargo test -p rustfoil-bl hkin
```

**Chat 3B - HS:**
```
@docs/tasks/TASK_04_HS.md

Implement HSL and HST closures in rustfoil-bl/src/closures/hs.rs.
Port XFOIL's xblsys.f lines 2327 (HSL) and 2388 (HST). Include:
1. HsResult struct with hs, hs_hk, hs_rt, hs_msq
2. hs_laminar() and hs_turbulent() functions with full derivatives
3. Unit tests
4. Update closures/mod.rs

Verify with: cargo test -p rustfoil-bl hs
```

**Chat 3C - CF:**
```
@docs/tasks/TASK_05_CF.md

Implement CFL and CFT closures in rustfoil-bl/src/closures/cf.rs.
Port XFOIL's xblsys.f lines 2354 (CFL) and 2483 (CFT). Include:
1. CfResult struct with cf, cf_hk, cf_rt, cf_msq
2. cf_laminar() and cf_turbulent() functions
3. Unit tests with derivative verification
4. Update closures/mod.rs

Verify with: cargo test -p rustfoil-bl cf
```

**Chat 3D - Dissipation:**
```
@docs/tasks/TASK_06_DISSIPATION.md

Implement DIL, DIT, DILW in rustfoil-bl/src/closures/dissipation.rs.
Port XFOIL's xblsys.f lines 2290, 2375, 2308. Include:
1. DissipationResult structs
2. dissipation_laminar(), dissipation_turbulent(), dissipation_wake()
3. Unit tests
4. Update closures/mod.rs

Verify with: cargo test -p rustfoil-bl dissipation
```

**Chat 3E - Density:**
```
@docs/tasks/TASK_07_DENSITY.md

Implement HCT in rustfoil-bl/src/closures/density.rs.
Port XFOIL's xblsys.f line 2514. Include:
1. HctResult struct with hc, hc_hk, hc_msq
2. density_shape() function (Whitfield correlation)
3. Unit tests with numerical derivative checks
4. Update closures/mod.rs

Verify with: cargo test -p rustfoil-bl density
```

**Chat 3F - Transition:**
```
@docs/tasks/TASK_08_TRANSITION.md

Implement DAMPL in rustfoil-bl/src/transition.rs.
Port XFOIL's xblsys.f line 1981 (Drela-Giles correlation). Include:
1. AmplificationResult struct with ax, ax_hk, ax_th, ax_rt
2. amplification_rate() function
3. check_transition() helper
4. Unit tests
5. Update lib.rs to export

Verify with: cargo test -p rustfoil-bl transition
```

**Chat 3G - State:**
```
@docs/tasks/TASK_09_STATE.md

Implement BlStation in rustfoil-bl/src/state.rs.
Match XFOIL's XBL.INC common block. Include:
1. BlStation struct with all primary variables (x, u, theta, delta_star, ctau, ampl)
2. All secondary variables (h, hk, hs, hc, r_theta, cf, cd)
3. BlDerivatives struct for Jacobian partials
4. BlStation::new() and BlStation::stagnation() constructors
5. Update lib.rs to export

Verify with: cargo test -p rustfoil-bl state
```

### Wave 4: Equations

```
@docs/tasks/TASK_10_EQUATIONS.md

Implement BLVAR and BLDIF in rustfoil-bl/src/equations.rs.
Port XFOIL's xblsys.f lines 784 and 1552. Include:
1. blvar() - compute all secondary variables from primary
2. BlResiduals struct for equation residuals
3. BlJacobian struct for Newton blocks
4. bldif() - compute residuals and Jacobian between stations
5. Unit tests
6. Update lib.rs

This requires all closures from Wave 3 to be complete.
Verify with: cargo test -p rustfoil-bl equations
```

### Wave 5: Coupling Crate

```
@docs/tasks/TASK_11_COUPLING_DIJ.md

Create rustfoil-coupling crate and implement QDCALC.
1. Create crate structure with Cargo.toml (depends on rustfoil-core, rustfoil-bl, nalgebra)
2. src/lib.rs exporting all modules
3. src/dij.rs with build_dij_matrix() - mass defect influence
4. Update workspace Cargo.toml

Verify with: cargo build -p rustfoil-coupling
```

### Wave 6: Newton System

```
@docs/tasks/TASK_12_COUPLING_NEWTON.md

Implement BLSYS in rustfoil-coupling/src/newton.rs.
Port XFOIL's xblsys.f line 583. Include:
1. BlNewtonSystem struct with VA, VB blocks and RHS
2. CoupledNewtonSystem for full viscous-inviscid coupling
3. build() method to construct system from stations
4. max_residual() for convergence checking

Verify with: cargo test -p rustfoil-coupling newton
```

### Wave 7: Block Solver

```
@docs/tasks/TASK_13_COUPLING_SOLVE.md

Implement BLSOLV in rustfoil-coupling/src/solve.rs.
Port XFOIL's xsolve.f line 283. Include:
1. solve_bl_system() - block Gaussian elimination
2. Forward sweep to eliminate lower diagonal
3. Back substitution
4. Helper functions: invert_3x3, multiply_3x3, multiply_3x3_vec

Verify with: cargo test -p rustfoil-coupling solve
```

### Wave 8: BL Marching

```
@docs/tasks/TASK_14_COUPLING_MARCH.md

Implement MRCHUE and MRCHDU in rustfoil-coupling/src/march.rs.
Port XFOIL's xbl.f lines 542 and 875. Include:
1. MarchResult struct with stations, x_transition, x_separation
2. MarchConfig for ncrit, tolerance, max_iter
3. march_fixed_ue() - BL march with prescribed Ue
4. march_coupled() - BL march with Ue updates

Verify with: cargo test -p rustfoil-coupling march
```

### Wave 9: Update Functions

```
@docs/tasks/TASK_15_COUPLING_UPDATE.md

Implement UPDATE and UESET in rustfoil-coupling/src/update.rs.
Port XFOIL's xbl.f line 1253 and xpanel.f line 1758. Include:
1. UpdateConfig for limiting parameters
2. update_stations() - apply Newton updates with limiting
3. set_edge_velocities() - compute Ue from inviscid + mass defect
4. limit_change() helper for stability

Verify with: cargo test -p rustfoil-coupling update
```

### Wave 10: Main Solver

```
@docs/tasks/TASK_16_VISCAL.md

Implement VISCAL in rustfoil-solver/src/viscous/.
Port XFOIL's xoper.f line 2886. Create:
1. src/viscous/mod.rs
2. src/viscous/config.rs - ViscousSolverConfig
3. src/viscous/viscal.rs - solve_viscous(), solve_viscous_polar_parallel()
4. src/viscous/forces.rs - compute_forces(), AeroForces
5. Update rustfoil-solver Cargo.toml to depend on rustfoil-coupling

Verify with: cargo build -p rustfoil-solver
```

### Wave 11: CLI Commands

```
@docs/tasks/TASK_17_CLI.md

Add viscous commands to rustfoil-cli/src/main.rs. Include:
1. ViscousCmd struct with alpha, re, mach, ncrit options
2. ViscousPolarCmd struct with alpha range, parallel flag
3. run_viscous() and run_viscous_polar() implementations
4. Add to Commands enum and main match

Test with: cargo run -p rustfoil-cli -- viscous --help
```

### Wave 12: Integration Tests

```
@docs/tasks/TASK_18_INTEGRATION_TESTS.md

Create end-to-end tests in rustfoil-solver/tests/. Include:
1. xfoil_comparison.rs - compare with XFOIL binary
2. XfoilRunner helper struct
3. test_naca0012_vs_xfoil() (ignored, run manually)
4. test_blasius_flat_plate() - analytical validation
5. test_transition_location()
6. test_separation_detection()

Run with: cargo test -p rustfoil-solver --test xfoil_comparison
```

---

## Overview

Implement XFOIL's viscous-inviscid coupling in Rust, building on the existing inviscid solver in `/Users/harry/flexfoil/crates`. Every new Rust function will be unit tested against its FORTRAN counterpart using gfortran-compiled test harnesses.

## Architecture

```mermaid
graph TD
    subgraph existing [Existing Code]
        CLI[rustfoil-cli]
        SOLVER[rustfoil-solver]
        INVISCID[inviscid module]
        CORE[rustfoil-core]
    end
    
    subgraph new [New Code]
        VISCOUS[viscous module]
        COUPLING[rustfoil-coupling]
        BL[rustfoil-bl]
        TESTKIT[rustfoil-testkit]
    end
    
    CLI --> SOLVER
    SOLVER --> INVISCID
    SOLVER --> VISCOUS
    VISCOUS --> COUPLING
    COUPLING --> BL
    COUPLING --> INVISCID
    BL --> CORE
    INVISCID --> CORE
    
    TESTKIT -.-> BL
    TESTKIT -.-> COUPLING
```

## FORTRAN Source Files (Reference)

| File | Lines | Key Functions |
|------|-------|---------------|
| Xfoil/src/xblsys.f | 2525 | HKIN, HSL, HST, CFL, CFT, DIL, DIT, HCT, DAMPL, BLVAR, BLDIF, BLSYS |
| Xfoil/src/xbl.f | 1598 | SETBL, MRCHUE, MRCHDU, UPDATE |
| Xfoil/src/xsolve.f | 488 | BLSOLV |
| Xfoil/src/xoper.f | 3161 | VISCAL, CDCALC |
| Xfoil/src/xpanel.f | 1796 | QDCALC, UESET |
| Xfoil/src/BLPAR.INC | 13 | Constants (SCCON, GACON, etc.) |
| Xfoil/src/XBL.INC | 73 | State variables and derivatives |

---

## Phase 1: Test Infrastructure and Closure Functions

### 1.1 Create Test Infrastructure

Create `rustfoil-testkit` crate with FORTRAN test harness compilation:

- **Location**: `/Users/harry/flexfoil/crates/rustfoil-testkit/`
- **Purpose**: Compile FORTRAN test programs, generate reference data, comparison utilities
- **Key files**:
  - `src/lib.rs` - Test utilities and JSON loading
  - `src/fortran_runner.rs` - Compile and run FORTRAN test programs
  - `src/approx.rs` - Floating-point comparison with tolerance
  - `fortran/test_closures.f` - Test harness for all closure functions
  - `fortran/Makefile` - Compile against XFOIL source

### 1.2 Create rustfoil-bl Crate

**Location**: `/Users/harry/flexfoil/crates/rustfoil-bl/`

**Structure**:
```
rustfoil-bl/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── constants.rs          # BLPAR constants
    ├── state.rs              # BlStation struct
    ├── closures/
    │   ├── mod.rs
    │   ├── hkin.rs           # HKIN (line 2276 of xblsys.f)
    │   ├── hs.rs             # HSL (2327), HST (2388)
    │   ├── cf.rs             # CFL (2354), CFT (2483)
    │   ├── dissipation.rs    # DIL (2290), DIT (2375), DILW (2308)
    │   └── density.rs        # HCT (2514)
    ├── transition.rs         # DAMPL (1981), DAMPL2 (2099), TRCHEK2
    └── equations.rs          # BLVAR (784), BLDIF (1552)
```

### 1.3 Closure Function Implementations

Each function returns value + all partial derivatives (critical for Newton Jacobian):

| Rust Function | FORTRAN | Inputs | Outputs | Test Points |
|---------------|---------|--------|---------|-------------|
| `hkin` | HKIN | H, M² | Hk, ∂Hk/∂H, ∂Hk/∂M² | 400 |
| `hs_laminar` | HSL | Hk, Rθ, M² | Hs, 3 derivatives | 500 |
| `hs_turbulent` | HST | Hk, Rθ, M² | Hs, 3 derivatives | 2000 |
| `cf_laminar` | CFL | Hk, Rθ, M² | Cf, 3 derivatives | 500 |
| `cf_turbulent` | CFT | Hk, Rθ, M² | Cf, 3 derivatives | 2000 |
| `dissipation_laminar` | DIL | Hk, Rθ | Di, 2 derivatives | 500 |
| `dissipation_turbulent` | DIT | Hs, Us, Cf, St | Di, 4 derivatives | 1000 |
| `density_shape` | HCT | Hk, M² | Hc, 2 derivatives | 400 |
| `amplification_rate` | DAMPL | Hk, θ, Rθ | Ax, 3 derivatives | 500 |

### 1.4 FORTRAN Test Harness

Create `fortran/test_closures.f`:
```fortran
      PROGRAM TEST_CLOSURES
      INCLUDE '../../Xfoil/src/BLPAR.INC'
C     Initialize BLPAR constants (from xfoil.f INIT)
      SCCON  = 5.6
      GACON  = 6.70
      GBCON  = 0.75
      GCCON  = 18.0
      DLCON  = 0.9
      CTRCON = 1.8
      CTRCEX = 3.3
      DUXCON = 1.0
      CTCON  = 0.5/(GACON**2 * GBCON)
      CFFAC  = 1.0
      
C     Generate test vectors for each function
      CALL TEST_HKIN()
      CALL TEST_HSL()
      ...
      END
```

---

## Phase 2: Transition Prediction

### 2.1 Transition Module

**File**: `rustfoil-bl/src/transition.rs`

Implement:
- `amplification_rate()` - DAMPL (line 1981)
- `amplification_rate_v2()` - DAMPL2 (line 2099) - modified e^n
- `check_transition()` - TRCHEK2 logic
- `axset()` - AXSET (line 35) - average amplification

### 2.2 FORTRAN Test Harness

Create `fortran/test_transition.f`:
- Test DAMPL at various (Hk, θ, Rθ) combinations
- Test TRCHEK2 with known transition scenarios
- Verify Rθ_crit calculation

---

## Phase 3: BL Equations and State

### 3.1 State Structure

**File**: `rustfoil-bl/src/state.rs`

Match XBL.INC common block structure:
```rust
pub struct BlStation {
    // Primary variables
    pub x: f64,      // Arc length
    pub u: f64,      // Edge velocity
    pub t: f64,      // Momentum thickness θ
    pub d: f64,      // Displacement thickness δ*
    pub s: f64,      // Shear stress coefficient (Cτ)
    pub ampl: f64,   // Amplification factor N
    
    // Derived quantities (computed by blvar)
    pub h: f64,      // Shape factor H = δ*/θ
    pub hk: f64,     // Kinematic Hk
    pub hs: f64,     // Energy Hs
    pub hc: f64,     // Density Hc
    pub rt: f64,     // Rθ
    pub cf: f64,     // Skin friction
    pub di: f64,     // Dissipation
    
    // All partial derivatives (73 values per XBL.INC)
    pub derivs: BlDerivatives,
    
    // Mode flags
    pub is_laminar: bool,
    pub is_wake: bool,
}
```

### 3.2 Equations Module

**File**: `rustfoil-bl/src/equations.rs`

- `blvar()` - BLVAR (line 784) - compute all secondary variables
- `bldif()` - BLDIF (line 1552) - compute residuals and Jacobian blocks

---

## Phase 4: Viscous-Inviscid Coupling

### 4.1 Create rustfoil-coupling Crate

**Location**: `/Users/harry/flexfoil/crates/rustfoil-coupling/`

**Structure**:
```
rustfoil-coupling/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── dij.rs           # QDCALC - mass defect influence
    ├── newton.rs        # BLSYS (xblsys.f:583), SETBL (xbl.f:21)
    ├── solve.rs         # BLSOLV (xsolve.f:283)
    ├── march.rs         # MRCHUE (xbl.f:542), MRCHDU (xbl.f:875)
    ├── update.rs        # UPDATE (xbl.f:1253), UESET, DSLIM
    └── stagnation.rs    # STFIND, STMOVE
```

### 4.2 DIJ Matrix

**File**: `rustfoil-coupling/src/dij.rs`

QDCALC builds the N×N source panel influence matrix relating mass defect to edge velocity:
```
ΔUe_i = Σ_j DIJ_ij * Δ(Ue*δ*)_j
```

### 4.3 Newton System

**File**: `rustfoil-coupling/src/newton.rs`

Block-tridiagonal structure from BLSYS:
- VA[i]: 3×2 diagonal block at station i
- VB[i]: 3×2 off-diagonal block
- VDEL[i]: 3×1 residual vector
- VS1, VS2: 4×5 sensitivity arrays per XBL.INC

### 4.4 Block Solver

**File**: `rustfoil-coupling/src/solve.rs`

BLSOLV implements block Gaussian elimination with:
- Forward sweep: eliminate sub-diagonal
- Back substitution: solve for updates
- Handle wake coupling to TE

---

## Phase 5: Top-Level Viscous Solver

### 5.1 Extend rustfoil-solver

**Location**: `/Users/harry/flexfoil/crates/rustfoil-solver/src/viscous/`

**Structure**:
```
viscous/
├── mod.rs
├── viscal.rs      # Main iteration loop (VISCAL)
├── forces.rs      # CDCALC - drag calculation
└── config.rs      # ViscousSolverConfig
```

### 5.2 VISCAL Main Loop

Port VISCAL (xoper.f:2886):
1. Initialize BL from stagnation point (MRCHUE)
2. Build Newton system (SETBL)
3. Solve (BLSOLV)
4. Update variables (UPDATE)
5. Check convergence
6. Repeat until converged or max iterations

### 5.3 Public API

```rust
pub fn solve_viscous(
    body: &Body,
    alpha: f64,
    config: &ViscousSolverConfig,
) -> Result<ViscousResult, SolverError>;

pub fn solve_viscous_polar_parallel(
    body: &Body,
    alphas: &[f64],
    config: &ViscousSolverConfig,
) -> Vec<Result<ViscousResult, SolverError>>;
```

---

## Phase 6: CLI Extension

### 6.1 New Commands

Add to `rustfoil-cli/src/main.rs`:

```
rustfoil viscous <FILE> -a <ALPHA> -r <RE> [-m <MACH>] [-n <NCRIT>]
rustfoil viscous-polar <FILE> --alpha <START:END:STEP> -r <RE> [--parallel]
```

### 6.2 Output Formats

- JSON: Full solution with BL state
- CSV: Polar data (alpha, CL, CD, CM, x_tr)
- Human-readable: Summary table

---

## Testing Strategy

### Unit Tests (Per Function)

Every Rust function tested against FORTRAN:

1. **Generate reference data**: Run `fortran/test_*.f` compiled with gfortran
2. **Save as JSON**: `testdata/reference/*.json`
3. **Rust tests**: Load JSON, compare with tolerance

```rust
#[test]
fn test_hkin_matches_fortran() {
    let tests = load_reference("closures/hkin.json");
    for t in tests {
        let result = hkin(t.h, t.msq);
        assert_relative_eq!(result.value, t.hk, epsilon = 1e-10);
        assert_relative_eq!(result.d_h, t.hk_h, epsilon = 1e-10);
        assert_relative_eq!(result.d_msq, t.hk_msq, epsilon = 1e-10);
    }
}
```

### Integration Tests

| Test Case | Inputs | Expected | Tolerance |
|-----------|--------|----------|-----------|
| Flat plate Blasius | Re=1e6, α=0° | H≈2.59 | 5% |
| NACA 0012 low-α | Re=3e6, α=4° | CD≈0.007 | 0.001 |
| NACA 0012 high-α | Re=3e6, α=8° | x_tr≈0.03 | 0.02 |
| Stall detection | Re=3e6, α=14° | CL_max≈1.2 | 0.1 |

### End-to-End XFOIL Comparison

```rust
#[test]
fn compare_polar_with_xfoil() {
    let xfoil = XfoilRunner::new("Xfoil/bin/xfoil");
    let body = Body::from_naca4(12, 160);
    
    for alpha in [-4.0, 0.0, 4.0, 8.0, 12.0] {
        let xfoil_result = xfoil.analyze(&body, alpha, 3e6);
        let rust_result = solve_viscous(&body, alpha, &config).unwrap();
        
        assert_relative_eq!(rust_result.cl, xfoil_result.cl, epsilon = 0.01);
        assert_relative_eq!(rust_result.cd, xfoil_result.cd, epsilon = 0.001);
    }
}
```

---

## File Summary

### New Files to Create

| Crate | File | FORTRAN Reference | Est. Lines |
|-------|------|-------------------|------------|
| rustfoil-testkit | src/lib.rs | - | 100 |
| rustfoil-testkit | src/fortran_runner.rs | - | 150 |
| rustfoil-testkit | src/approx.rs | - | 50 |
| rustfoil-testkit | fortran/test_closures.f | xblsys.f | 300 |
| rustfoil-testkit | fortran/test_transition.f | xblsys.f | 150 |
| rustfoil-testkit | fortran/Makefile | - | 30 |
| rustfoil-bl | src/lib.rs | - | 50 |
| rustfoil-bl | src/constants.rs | BLPAR.INC | 30 |
| rustfoil-bl | src/state.rs | XBL.INC | 150 |
| rustfoil-bl | src/closures/mod.rs | - | 30 |
| rustfoil-bl | src/closures/hkin.rs | HKIN | 50 |
| rustfoil-bl | src/closures/hs.rs | HSL, HST | 200 |
| rustfoil-bl | src/closures/cf.rs | CFL, CFT | 200 |
| rustfoil-bl | src/closures/dissipation.rs | DIL, DIT, DILW | 150 |
| rustfoil-bl | src/closures/density.rs | HCT | 50 |
| rustfoil-bl | src/transition.rs | DAMPL, TRCHEK2 | 300 |
| rustfoil-bl | src/equations.rs | BLVAR, BLDIF | 400 |
| rustfoil-coupling | src/lib.rs | - | 50 |
| rustfoil-coupling | src/dij.rs | QDCALC | 200 |
| rustfoil-coupling | src/newton.rs | BLSYS, SETBL | 400 |
| rustfoil-coupling | src/solve.rs | BLSOLV | 250 |
| rustfoil-coupling | src/march.rs | MRCHUE, MRCHDU | 350 |
| rustfoil-coupling | src/update.rs | UPDATE, UESET | 200 |
| rustfoil-coupling | src/stagnation.rs | STFIND, STMOVE | 100 |
| rustfoil-solver | src/viscous/mod.rs | - | 30 |
| rustfoil-solver | src/viscous/viscal.rs | VISCAL | 300 |
| rustfoil-solver | src/viscous/forces.rs | CDCALC | 100 |
| rustfoil-solver | src/viscous/config.rs | - | 50 |

**Total**: ~4,000 lines of Rust + 500 lines of FORTRAN test harnesses

### Files to Modify

| File | Changes |
|------|---------|
| `/Users/harry/flexfoil/Cargo.toml` | Add rustfoil-bl, rustfoil-coupling, rustfoil-testkit to workspace |
| `/Users/harry/flexfoil/crates/rustfoil-solver/Cargo.toml` | Add dependency on rustfoil-coupling |
| `/Users/harry/flexfoil/crates/rustfoil-cli/src/main.rs` | Add viscous and viscous-polar commands |

---

## Implementation Order (TODO List)

1. **rustfoil-testkit** + FORTRAN harness for closures
2. **rustfoil-bl::constants** (BLPAR)
3. **rustfoil-bl::closures::hkin** + tests
4. **rustfoil-bl::closures::hs** (HSL, HST) + tests
5. **rustfoil-bl::closures::cf** (CFL, CFT) + tests
6. **rustfoil-bl::closures::dissipation** + tests
7. **rustfoil-bl::closures::density** + tests
8. **rustfoil-bl::transition** (DAMPL, TRCHEK2) + tests
9. **rustfoil-bl::state** (BlStation)
10. **rustfoil-bl::equations** (BLVAR, BLDIF) + tests
11. **rustfoil-coupling::dij** (QDCALC) + tests
12. **rustfoil-coupling::newton** (BLSYS, SETBL) + tests
13. **rustfoil-coupling::solve** (BLSOLV) + tests
14. **rustfoil-coupling::march** (MRCHUE, MRCHDU)
15. **rustfoil-coupling::update** (UPDATE, UESET)
16. **rustfoil-solver::viscous** (VISCAL)
17. **rustfoil-cli** viscous commands
18. Integration tests + XFOIL comparison

---

## Key FORTRAN Function Locations

For reference when implementing each function:

### Closures (xblsys.f)
- HKIN: line 2276
- HSL: line 2327
- HST: line 2388
- CFL: line 2354
- CFT: line 2483
- DIL: line 2290
- DIT: line 2375
- DILW: line 2308
- HCT: line 2514
- DAMPL: line 1981
- DAMPL2: line 2099
- BLVAR: line 784
- BLDIF: line 1552
- BLSYS: line 583

### BL Marching (xbl.f)
- SETBL: line 21
- MRCHUE: line 542
- MRCHDU: line 875
- UPDATE: line 1253

### Solver (xsolve.f)
- BLSOLV: line 283

### Main Loop (xoper.f)
- VISCAL: line 2886

### Panel Influence (xpanel.f)
- QDCALC: line 1149
- UESET: line 1758

---

## Prerequisites

- gfortran available at `/opt/homebrew/bin/gfortran`
- XFOIL source at `/Users/harry/flexfoil-boundary-layer/Xfoil/src/`
- Existing inviscid solver at `/Users/harry/flexfoil/crates/`
