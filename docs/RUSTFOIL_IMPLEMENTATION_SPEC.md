# RustFoil Implementation Specification

> **Document Status**: Implementation Plan (Updated)  
> **Target**: Rust port of XFOIL 6.99 with parallelization & GPU acceleration  
> **Last Updated**: 2026-01-20  
> **Existing Code**: `/Users/harry/flexfoil/crates`

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Existing Implementation Status](#2-existing-implementation-status)
3. [Remaining Work: Viscous Solver](#3-remaining-work-viscous-solver)
4. [CLI Specification](#4-cli-specification)
5. [Module Specifications](#5-module-specifications)
6. [FORTRAN → Rust Function Mapping](#6-fortran--rust-function-mapping)
7. [Parallelization Strategy](#7-parallelization-strategy)
8. [GPU Acceleration Roadmap](#8-gpu-acceleration-roadmap)
9. [Validation & Test Specification](#9-validation--test-specification)
10. [Implementation Phases](#10-implementation-phases)

---

## 1. Executive Summary

### Goals

1. **Correctness First**: Bit-accurate reproduction of XFOIL results for identical inputs
2. **Performance**: 10-100× speedup for polar sweeps via parallelization
3. **Extensibility**: Clean module boundaries enabling future GPU acceleration
4. **Usability**: Modern CLI with structured output (JSON/CSV) for automation

### Key Finding: Inviscid Solver Already Exists

The `/Users/harry/flexfoil/crates` workspace contains a **complete inviscid solver** with:
- XFOIL-exact linear vorticity panel method
- Two-solution factorization for efficient α sweeps
- Stream function and velocity field evaluation
- XFOIL-compatible NACA generation and paneling (PANGEN algorithm)

**The remaining work is focused on the viscous boundary layer solver.**

---

## 2. Existing Implementation Status

### 2.1 Workspace Structure (Already Built)

```
/Users/harry/flexfoil/crates/
├── rustfoil-core/           ✅ COMPLETE
│   ├── body.rs              ✅ Multi-body airfoil representation
│   ├── panel.rs             ✅ Panel discretization with cached geometry
│   ├── point.rs             ✅ 2D points/vectors with nalgebra
│   ├── spline.rs            ✅ CubicSpline with XFOIL PANGEN algorithm
│   ├── xfoil_spline.rs      ✅ Hermite splines (XFOIL's spline.f exact)
│   ├── naca.rs              ✅ NACA 4-digit generator (XFOIL exact)
│   └── error.rs             ✅ Error types
│
├── rustfoil-solver/         ⚠️  INVISCID COMPLETE, VISCOUS NEEDED
│   └── inviscid/
│       ├── mod.rs           ✅ Panel method (PSILIN, Kutta, factorization)
│       ├── velocity.rs      ✅ velocity_at, psi_at, streamlines
│       ├── smoke.rs         ✅ Particle visualization
│       └── error.rs         ✅ Solver errors
│
├── rustfoil-cli/            ⚠️  BASIC CLI, NEEDS VISCOUS COMMANDS
│   └── main.rs              ✅ analyze, polar, repanel, info
│
└── rustfoil-wasm/           ⚠️  EXISTS, NEEDS VISCOUS BINDINGS
    └── lib.rs               ✅ WASM interface
```

### 2.2 What We DON'T Need to Build

| Component | Status | Notes |
|-----------|--------|-------|
| Spline interpolation | ✅ Done | Both CubicSpline and XfoilSpline |
| Arc-length parameterization | ✅ Done | In spline.rs |
| XFOIL PANGEN paneling | ✅ Done | Curvature-based, exact match |
| NACA 4-digit generator | ✅ Done | Matches XFOIL naca.f |
| Panel geometry (midpoints, normals) | ✅ Done | Cached in Panel struct |
| Influence matrix (AIJ) | ✅ Done | Stream function formulation |
| Kutta condition | ✅ Done | γ₀ + γₙ₋₁ = 0 |
| Two-solution factorization | ✅ Done | α=0° and α=90° base solutions |
| CL, CM calculation | ✅ Done | Pressure integration |
| Velocity field | ✅ Done | velocity_at() with K&P VOR2DL |
| Stream function | ✅ Done | psi_at() with PSILIN |
| Point-in-polygon test | ✅ Done | is_inside_airfoil() |
| Streamline integration | ✅ Done | RK4 with adaptive stepping |
| .dat file I/O | ✅ Done | Selig format parsing |

### 2.3 Existing API (Reference)

```rust
// rustfoil-core
pub struct Body { ... }
pub struct Panel { p1, p2, midpoint, normal, tangent, length, ... }
pub struct CubicSpline { ... }
pub struct XfoilSpline { ... }
pub fn naca4(designation: u32, nside: Option<usize>) -> Vec<Point>;

// rustfoil-solver::inviscid
pub struct InviscidSolver { ... }
pub struct FlowConditions { alpha: f64, v_inf: f64 }
pub struct InviscidSolution { gamma, cp, cl, cm, psi_0, ... }
pub struct FactorizedSolution { ... }  // For efficient α sweeps

impl InviscidSolver {
    pub fn solve(&self, bodies: &[Body], flow: &FlowConditions) -> Result<InviscidSolution>;
    pub fn factorize(&self, bodies: &[Body]) -> Result<FactorizedSolution>;
}

impl FactorizedSolution {
    pub fn solve_alpha(&self, flow: &FlowConditions) -> InviscidSolution;
}

// Velocity field
pub fn velocity_at(x, y, nodes, gamma, alpha, v_inf) -> (f64, f64);
pub fn psi_at(x, y, nodes, gamma, alpha, v_inf) -> f64;
pub fn is_inside_airfoil(x, y, nodes) -> bool;
pub fn build_streamlines(...) -> Vec<Vec<(f64, f64)>>;
```

---

## 3. Remaining Work: Viscous Solver

### 3.1 New Crate: `rustfoil-bl`

Boundary layer physics - the core of what needs to be built:

```
rustfoil-bl/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── state.rs              # BlStation, BlSide, BlState structs
    ├── constants.rs          # BLPAR constants (SCCON, GACON, etc.)
    ├── closures/
    │   ├── mod.rs
    │   ├── hkin.rs           # Kinematic shape parameter
    │   ├── hs.rs             # Energy shape (HSL, HST)
    │   ├── cf.rs             # Skin friction (CFL, CFT)
    │   ├── dissipation.rs    # DIL, DIT, DILW
    │   └── density.rs        # HCT density shape
    ├── transition.rs         # e^N method (DAMPL, TRCHEK)
    ├── equations.rs          # BLDIF momentum/shape equations
    └── shear_lag.rs          # Turbulent shear lag (G-β locus)
```

### 3.2 New Crate: `rustfoil-coupling`

Viscous-inviscid coupling - the Newton solver:

```
rustfoil-coupling/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── dij.rs                # DIJ mass defect influence matrix
    ├── march.rs              # MRCHUE, MRCHDU BL marching
    ├── inverse.rs            # Hk-prescribed inverse mode
    ├── newton.rs             # Newton system build (SETBL, BLSYS)
    ├── solve.rs              # Block elimination (BLSOLV)
    ├── update.rs             # UPDATE, UESET, DSSET, DSLIM
    └── stagnation.rs         # STFIND, STMOVE
```

### 3.3 Extend: `rustfoil-solver`

Add viscous solver to existing crate:

```
rustfoil-solver/
└── src/
    ├── inviscid/             # ✅ Already exists
    └── viscous/              # NEW
        ├── mod.rs
        ├── viscal.rs         # VISCAL main loop
        ├── forces.rs         # CDCALC with friction drag
        └── config.rs         # ViscousSolverConfig
```

### 3.4 Updated Dependency Graph

```
rustfoil-cli (updated)
    └── rustfoil-solver (updated)
            ├── inviscid/     ✅ exists
            ├── viscous/      NEW
            │       └── rustfoil-coupling (NEW)
            │               ├── rustfoil-bl (NEW)
            │               └── rustfoil-core ✅
            └── rustfoil-core ✅
```

---

## 4. CLI Specification

### 4.1 Existing Commands (Keep As-Is)

```bash
rustfoil analyze <FILE> --alpha <DEG>           # ✅ Inviscid
rustfoil polar <FILE> --alpha-start/end/step    # ✅ Inviscid
rustfoil repanel <FILE> --panels <N>            # ✅ Geometry
rustfoil info <FILE>                            # ✅ Geometry
```

### 4.2 New Commands for Viscous

```bash
rustfoil viscous <FILE> [OPTIONS]
    -a, --alpha <DEG>       Angle of attack [required unless --cl]
    -c, --cl <CL>           Target lift coefficient
    -r, --re <RE>           Reynolds number [required]
    -m, --mach <M>          Mach number [default: 0.0]
    -n, --ncrit <N>         Transition Ncrit [default: 9.0]
    --max-iter <N>          Max Newton iterations [default: 100]
    --tolerance <TOL>       Convergence tolerance [default: 1e-4]
    --json                  Output as JSON
    -v, --verbose           Show iteration progress

rustfoil viscous-polar <FILE> [OPTIONS]
    --alpha <START:END:STEP>   Alpha range
    --cl <START:END:STEP>      CL range (alternative)
    -r, --re <RE>              Reynolds number [required]
    -m, --mach <M>             Mach number [default: 0.0]
    -n, --ncrit <N>            Transition Ncrit [default: 9.0]
    --parallel                 Parallel alpha sweep [default: true]
    --continue-on-fail         Continue on convergence failure
    -o, --output <FILE>        Output file (CSV or JSON)

rustfoil validate [OPTIONS]
    --suite <NAME>          closures | transition | inviscid | viscous | all
    --xfoil <PATH>          Path to XFOIL binary for comparison
    --tolerance <T>         Acceptable deviation
    --generate              Generate new reference data
```

### 4.3 Output Format (Viscous)

```json
{
  "airfoil": "naca0012",
  "alpha": 8.0,
  "re": 3000000,
  "mach": 0.0,
  "ncrit": 9.0,
  "results": {
    "cl": 0.9012,
    "cd": 0.00982,
    "cm": -0.0234,
    "cdp": 0.00312,
    "cdf": 0.00670,
    "x_tr_upper": 0.0312,
    "x_tr_lower": 0.4521
  },
  "converged": true,
  "iterations": 8,
  "rms_final": 4.2e-5
}
```

---

## 5. Module Specifications

### 5.1 `rustfoil-bl::state`

```rust
/// State at a single BL station
#[derive(Debug, Clone)]
pub struct BlStation {
    // Position
    pub s: f64,               // Arc length from stagnation
    pub panel_idx: usize,     // Corresponding panel index
    
    // Primary variables (Newton unknowns)
    pub theta: f64,           // Momentum thickness θ
    pub dstar: f64,           // Displacement thickness δ*
    pub ctau: f64,            // Shear stress coeff (turb) or N (lam)
    pub ue: f64,              // Edge velocity
    
    // Derived quantities
    pub h: f64,               // Shape factor H = δ*/θ
    pub hk: f64,              // Kinematic shape parameter
    pub hs: f64,              // Energy shape parameter
    pub hss: f64,             // Density shape parameter H**
    pub cf: f64,              // Skin friction coefficient
    pub cd: f64,              // Dissipation coefficient
    pub rtheta: f64,          // Momentum thickness Reynolds number
    
    // Mode flags
    pub is_laminar: bool,
    pub is_wake: bool,
    pub is_inverse: bool,     // Hk prescribed (separated region)
}

/// BL state for one surface side
pub struct BlSide {
    pub stations: Vec<BlStation>,
    pub transition_idx: Option<usize>,
    pub transition_x: Option<f64>,
    pub te_idx: usize,
}

/// Full BL state (both sides + wake)
pub struct BlState {
    pub upper: BlSide,
    pub lower: BlSide,
    pub wake: Vec<BlStation>,
}
```

### 5.2 `rustfoil-bl::closures`

All closure functions return value + derivatives for Newton Jacobian:

```rust
pub struct ClosureOutput {
    pub value: f64,
    pub d_hk: f64,      // ∂/∂Hk
    pub d_rtheta: f64,  // ∂/∂Rθ
    pub d_msq: f64,     // ∂/∂M²
}

/// HKIN: Kinematic shape parameter
/// Hk = (H - 0.29·M²) / (1 + 0.113·M²)
pub fn hkin(h: f64, msq: f64) -> ClosureOutput;

/// HSL: Laminar energy shape parameter
pub fn hs_laminar(hk: f64) -> ClosureOutput;

/// HST: Turbulent energy shape parameter
pub fn hs_turbulent(hk: f64, rtheta: f64, msq: f64) -> ClosureOutput;

/// CFL: Laminar skin friction (Falkner-Skan)
pub fn cf_laminar(hk: f64, rtheta: f64) -> ClosureOutput;

/// CFT: Turbulent skin friction (Coles wall-law)
pub fn cf_turbulent(hk: f64, rtheta: f64, msq: f64) -> ClosureOutput;

/// DIL: Laminar dissipation 2CD/H*
pub fn dissipation_laminar(hk: f64, rtheta: f64) -> ClosureOutput;

/// DIT: Turbulent dissipation
pub fn dissipation_turbulent(hk: f64, hs: f64, cf: f64, rtheta: f64, msq: f64) -> ClosureOutput;

/// HCT: Density shape parameter H**
pub fn density_shape(hk: f64, msq: f64) -> ClosureOutput;
```

### 5.3 `rustfoil-bl::transition`

```rust
/// Amplification rate (DAMPL equivalent)
pub struct AmplificationResult {
    pub ax: f64,              // dN/dx
    pub ax_hk: f64,           // ∂(dN/dx)/∂Hk
    pub ax_theta: f64,        // ∂(dN/dx)/∂θ
    pub rtheta_crit: f64,     // Critical Rθ
}

pub fn amplification_rate(hk: f64, theta: f64, rtheta: f64) -> AmplificationResult;

/// Transition check between two stations (TRCHEK2)
pub struct TransitionCheck {
    pub n2: f64,              // Amplification at station 2
    pub transitioned: bool,   // Did transition occur?
    pub x_transition: f64,    // Interpolated location
}

pub fn check_transition(
    station1: &BlStation,
    station2: &BlStation,
    n1: f64,
    ncrit: f64,
) -> TransitionCheck;
```

### 5.4 `rustfoil-coupling::newton`

```rust
/// Block-tridiagonal Newton system
pub struct NewtonSystem {
    /// VA: Diagonal blocks (3×2 per station)
    pub va: Vec<[[f64; 2]; 3]>,
    /// VB: Off-diagonal blocks
    pub vb: Vec<[[f64; 2]; 3]>,
    /// VM: Mass defect coupling (sparse or dense)
    pub vm: MassInfluence,
    /// VDEL: RHS residuals
    pub vdel: Vec<[f64; 3]>,
    /// VZ: TE coupling
    pub vz: [[f64; 2]; 3],
}

impl NewtonSystem {
    /// Build from current BL state (SETBL equivalent)
    pub fn build(
        bl_state: &BlState,
        dij: &DijMatrix,
        params: &SolverParams,
    ) -> Self;
    
    /// Solve via block elimination (BLSOLV)
    pub fn solve(&mut self) -> Vec<[f64; 3]>;
}
```

### 5.5 `rustfoil-solver::viscous`

```rust
pub struct ViscousSolverConfig {
    pub re: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub max_iterations: usize,
    pub tolerance: f64,
    pub vaccel: f64,
}

pub struct ViscousResult {
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cdp: f64,              // Pressure drag
    pub cdf: f64,              // Friction drag
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
    pub converged: bool,
    pub iterations: usize,
    pub rms_history: Vec<f64>,
    pub bl_state: BlState,
    pub cp: Vec<f64>,
}

/// Main viscous solver (uses existing inviscid solver internally)
pub fn solve_viscous(
    body: &Body,
    alpha: f64,
    config: &ViscousSolverConfig,
) -> Result<ViscousResult, SolverError>;

/// Solve for alpha given CL target
pub fn solve_for_cl(
    body: &Body,
    cl_target: f64,
    config: &ViscousSolverConfig,
) -> Result<ViscousResult, SolverError>;
```

---

## 6. FORTRAN → Rust Function Mapping

### 6.1 Already Implemented (in flexfoil/crates)

| FORTRAN | Location | Rust | Status |
|---------|----------|------|--------|
| `SPLINE` | spline.f | `CubicSpline::from_points` | ✅ |
| `SPLINT` | spline.f | `CubicSpline::evaluate` | ✅ |
| `SEGSPL` | spline.f | `XfoilSpline::new` | ✅ |
| `SEVAL` | spline.f | `XfoilSpline::eval` | ✅ |
| `DEVAL` | spline.f | `XfoilSpline::deriv` | ✅ |
| `CURV` | spline.f | `XfoilSpline::curvature` | ✅ |
| `LEFIND` | spline.f | `XfoilSpline::lefind` | ✅ |
| `NACA4` | naca.f | `naca::naca4` | ✅ |
| `PANGEN` | xfoil.f | `CubicSpline::resample_xfoil` | ✅ |
| `PANCOP` | xpanel.f | `Body::from_points` | ✅ |
| `PSILIN` | xpanel.f | `InviscidSolver::factorize` | ✅ |
| `GAMSOLV` | xpanel.f | (inline in factorize) | ✅ |
| `CLCALC` | xoper.f | `FactorizedSolution::compute_forces` | ✅ |
| Velocity field | - | `velocity_at` | ✅ |
| Stream function | - | `psi_at` | ✅ |

### 6.2 To Be Implemented (Boundary Layer)

| FORTRAN | File | Rust Module | Rust Function | Priority |
|---------|------|-------------|---------------|----------|
| `HKIN` | xblsys.f | rustfoil-bl::closures | `hkin` | P0 |
| `HSL` | xblsys.f | rustfoil-bl::closures | `hs_laminar` | P0 |
| `HST` | xblsys.f | rustfoil-bl::closures | `hs_turbulent` | P0 |
| `CFL` | xblsys.f | rustfoil-bl::closures | `cf_laminar` | P0 |
| `CFT` | xblsys.f | rustfoil-bl::closures | `cf_turbulent` | P0 |
| `DIL` | xblsys.f | rustfoil-bl::closures | `dissipation_laminar` | P0 |
| `DIT` | xblsys.f | rustfoil-bl::closures | `dissipation_turbulent` | P0 |
| `HCT` | xblsys.f | rustfoil-bl::closures | `density_shape` | P0 |
| `DAMPL` | xblsys.f | rustfoil-bl::transition | `amplification_rate` | P1 |
| `TRCHEK2` | xblsys.f | rustfoil-bl::transition | `check_transition` | P1 |
| `BLDIF` | xblsys.f | rustfoil-bl::equations | `bl_residuals` | P1 |
| `BLVAR` | xblsys.f | rustfoil-bl::equations | `compute_all_closures` | P1 |

### 6.3 To Be Implemented (Coupling/Newton)

| FORTRAN | File | Rust Module | Rust Function | Priority |
|---------|------|-------------|---------------|----------|
| `QDCALC` | xpanel.f | rustfoil-coupling::dij | `DijMatrix::build` | P2 |
| `SETBL` | xbl.f | rustfoil-coupling::newton | `NewtonSystem::build` | P2 |
| `BLSYS` | xblsys.f | rustfoil-coupling::newton | `build_station_system` | P2 |
| `MRCHUE` | xbl.f | rustfoil-coupling::march | `initialize_bl` | P2 |
| `MRCHDU` | xbl.f | rustfoil-coupling::march | `establish_transition` | P2 |
| `BLSOLV` | xsolve.f | rustfoil-coupling::solve | `NewtonSystem::solve` | P2 |
| `UPDATE` | xbl.f | rustfoil-coupling::update | `apply_updates` | P2 |
| `UESET` | xpanel.f | rustfoil-coupling::update | `edge_velocity_from_mass` | P2 |
| `DSSET` | xbl.f | rustfoil-coupling::update | `dstar_from_mass` | P2 |
| `DSLIM` | xbl.f | rustfoil-coupling::update | `clamp_shape_factor` | P2 |
| `VISCAL` | xoper.f | rustfoil-solver::viscous | `solve_viscous` | P3 |
| `STFIND` | xoper.f | rustfoil-coupling::stagnation | `find_stagnation` | P2 |
| `STMOVE` | xoper.f | rustfoil-coupling::stagnation | `relocate_stagnation` | P2 |
| `CDCALC` | xoper.f | rustfoil-solver::viscous | `compute_drag` | P3 |

---

## 7. Parallelization Strategy

### 7.1 Already Parallelizable (Existing Code)

| Operation | Current Status | Parallelization |
|-----------|---------------|-----------------|
| Influence matrix (AIJ) | ✅ Sequential | Row-parallel possible |
| α sweep (inviscid) | ✅ Sequential | α-point parallel (trivial) |
| Streamline integration | ✅ Sequential | Per-streamline parallel |

### 7.2 New Parallelization Opportunities

| Operation | Data Size | Strategy | Expected Speedup |
|-----------|-----------|----------|------------------|
| DIJ matrix build | N² | Row-parallel | ~8× |
| Closure evaluation | N stations | Station-parallel | ~2× |
| **Viscous polar sweep** | M points | **α-point parallel** | **~M×** |
| Upper/lower BL march | 2 sides | Side-parallel | ~2× |

### 7.3 Key Implementation: Parallel Polar

```rust
use rayon::prelude::*;

pub fn generate_viscous_polar_parallel(
    body: &Body,
    alpha_points: &[f64],
    config: &ViscousSolverConfig,
) -> Vec<Result<ViscousResult, SolverError>> {
    // Pre-compute shared geometry (inviscid factorization)
    let inviscid = InviscidSolver::new();
    let factorized = inviscid.factorize(&[body.clone()]).unwrap();
    
    // Each alpha point is independent
    alpha_points.par_iter()
        .map(|&alpha| {
            // Each thread gets its own BL state
            solve_viscous_with_inviscid(body, &factorized, alpha, config)
        })
        .collect()
}
```

---

## 8. GPU Acceleration Roadmap

### Phase 1: CPU Baseline (Current Focus)
- Complete viscous solver on CPU
- Use `rayon` for polar parallelism
- Validate against XFOIL

### Phase 2: GPU-Ready Refactor
- Abstract linear algebra operations
- Prepare DIJ matrix for GPU upload

### Phase 3: GPU Implementation (Future)
- DIJ matrix construction on GPU (wgpu compute shaders)
- Batched closure evaluation
- Target: 50× speedup for 100-point polars

---

## 9. Validation & Test Specification

### 9.1 Closure Function Tests (Priority: Highest)

Generate FORTRAN reference data for every closure function:

```fortran
      PROGRAM TEST_CLOSURES
      IMPLICIT NONE
      REAL*8 H, MSQ, HK, HK_H, HK_MSQ
C     Generate test vectors
      DO I = 1, 50
        H = 1.2 + 0.1*I
        MSQ = 0.0
        CALL HKIN(H, MSQ, HK, HK_H, HK_MSQ)
        WRITE(*,'(5E20.12)') H, MSQ, HK, HK_H, HK_MSQ
      END DO
      END
```

**Test file**: `testdata/reference/closures_test_vectors.json`

```rust
#[test]
fn test_hkin_against_fortran() {
    let tests = load_closure_tests().hkin;
    for t in tests {
        let result = hkin(t.h, t.msq);
        assert_relative_eq!(result.value, t.expected_hk, epsilon = 1e-10);
        assert_relative_eq!(result.d_h, t.expected_dhk_dh, epsilon = 1e-10);
    }
}
```

### 9.2 Integration Tests

| Test | Inputs | Expected | Tolerance |
|------|--------|----------|-----------|
| Flat plate (Blasius) | Re=1e6, α=0° | H≈2.59, Cf=0.664/√Rex | 5% |
| NACA 0012 inviscid | α=4° | CL≈0.46 | 0.01 |
| NACA 0012 viscous | Re=3e6, α=0° | CD≈0.006 | 0.001 |
| NACA 0012 viscous | Re=3e6, α=8° | CL≈0.9, x_tr≈0.03 | 0.05, 0.02 |
| Stall (NACA 0012) | Re=3e6, α=14-16° | CL_max≈1.2 | 0.1 |

### 9.3 XFOIL Live Comparison

```rust
#[cfg(feature = "xfoil-compare")]
#[test]
fn compare_with_xfoil() {
    let xfoil = XfoilRunner::new("/path/to/xfoil");
    let body = Body::from_naca4(12, 160);
    
    let xfoil_result = xfoil.analyze(&body, 8.0, 3e6, 0.0, 9.0);
    let rust_result = solve_viscous(&body, 8.0, &config).unwrap();
    
    assert_relative_eq!(rust_result.cl, xfoil_result.cl, epsilon = 0.01);
    assert_relative_eq!(rust_result.cd, xfoil_result.cd, epsilon = 0.001);
}
```

---

## 10. Implementation Phases

### Phase 1: Closure Relations (Week 1) ⬅️ START HERE

**Goal**: All 8 closure functions with 100% test coverage

- [ ] Create `rustfoil-bl` crate
- [ ] Implement `hkin` with derivatives
- [ ] Implement `hs_laminar`, `hs_turbulent`
- [ ] Implement `cf_laminar`, `cf_turbulent`
- [ ] Implement `dissipation_laminar`, `dissipation_turbulent`
- [ ] Implement `density_shape`
- [ ] Generate FORTRAN reference data
- [ ] 100% test coverage

**Deliverable**: `cargo test -p rustfoil-bl` passes all closure tests

### Phase 2: Transition Prediction (Week 2)

- [ ] Implement `amplification_rate` (DAMPL)
- [ ] Implement `check_transition` (TRCHEK2)
- [ ] Validate transition location vs XFOIL

**Deliverable**: Transition location matches XFOIL within 2% chord

### Phase 3: BL Equations & Marching (Weeks 3-4)

- [ ] Create `rustfoil-coupling` crate
- [ ] Implement `bl_residuals` (BLDIF)
- [ ] Implement `march::initialize_bl` (MRCHUE)
- [ ] Implement `march::establish_transition` (MRCHDU)
- [ ] Direct mode only (no separation)

**Deliverable**: Low-α viscous solutions match XFOIL

### Phase 4: Newton System & Coupling (Weeks 5-6)

- [ ] Implement `DijMatrix::build` (QDCALC)
- [ ] Implement `NewtonSystem::build` (SETBL)
- [ ] Implement `NewtonSystem::solve` (BLSOLV)
- [ ] Implement `update` module (UPDATE, UESET, DSSET, DSLIM)

**Deliverable**: Full Newton iteration converges

### Phase 5: Inverse Mode & Stall (Week 7)

- [ ] Implement inverse mode (Hk prescribed)
- [ ] Mode switching logic
- [ ] Test CL_max prediction

**Deliverable**: Post-stall polars match XFOIL

### Phase 6: Parallelization & CLI (Week 8)

- [ ] Add `rustfoil viscous` command
- [ ] Add `rustfoil viscous-polar` command
- [ ] Parallel polar with `rayon`
- [ ] JSON/CSV output

**Deliverable**: `rustfoil viscous-polar naca0012 --alpha="-5:20:0.5" -r 3e6 --parallel`

### Phase 7: Validation & Release (Week 9)

- [ ] Full validation suite
- [ ] Documentation
- [ ] v0.2.0 release (viscous-capable)

---

## Appendix A: BLPAR Constants

These constants must be exactly matched from XFOIL's `BLPAR.INC`:

```rust
// rustfoil-bl/src/constants.rs
pub const SCCON: f64 = 5.6;      // Shear coefficient lag constant
pub const GACON: f64 = 6.70;     // G-β locus constant A
pub const GBCON: f64 = 0.75;     // G-β locus constant B
pub const GCCON: f64 = 18.0;     // G-β wall term constant
pub const DLCON: f64 = 0.9;      // Wake dissipation length ratio
pub const CTRCON: f64 = 1.8;     // Cτ transition constant
pub const CTRCEX: f64 = 3.3;     // Cτ transition exponent
pub const DUXCON: f64 = 1.0;     // dUe/dx weighting
pub const CTCON: f64 = 0.5 / (GACON * GACON * GBCON);  // Derived
pub const CFFAC: f64 = 1.0;      // Cf correction factor

// Mode switching thresholds
pub const HLMAX: f64 = 3.8;      // Laminar Hk threshold
pub const HTMAX: f64 = 2.5;      // Turbulent Hk threshold

// Clamping limits
pub const HK_MIN_AIRFOIL: f64 = 1.02;
pub const HK_MIN_WAKE: f64 = 1.00005;
pub const CTAU_MIN: f64 = 1e-7;
pub const CTAU_MAX: f64 = 0.30;
pub const RTHETA_MIN: f64 = 200.0;
```

## Appendix B: Reference Formulas

### HKIN (Whitfield compressibility correction)
```
Hk = (H - 0.29·M²) / (1 + 0.113·M²)
```

### HSL (Laminar energy shape)
```
If Hk < 4.35:
    H* = 0.0111·(Hk-4.35)²/(Hk+1) - 0.0278·(Hk-4.35)³/(Hk+1) + 1.528 - 0.0002·(Hk·(Hk-4.35))²
Else:
    H* = 0.015·(Hk-4.35)²/Hk + 1.528
```

### HST (Turbulent energy shape)
```
Ho = 3 + 400/Rθ  (Rθ > 400)
Ho = 4           (Rθ ≤ 400)

If Hk < Ho:
    Hr = (Ho - Hk)/(Ho - 1)
    H* = (2 - 1.5 - 4/Rθ)·Hr²·1.5/(Hk+0.5) + 1.5 + 4/Rθ
Else:
    H* = (Hk - Ho)²·[0.007·ln(Rθ)/(Hk-Ho+4/ln(Rθ))² + 0.015/Hk] + 1.5 + 4/Rθ

Compressibility: H* = (H* + 0.028·M²) / (1 + 0.014·M²)
```

### DAMPL (Amplification rate)
```
log₁₀(Rθ_crit) = 2.492/(Hk-1)^0.43 + 0.7·(tanh(14/(Hk-1) - 9.24) + 1)

If Rθ > Rθ_crit:
    DADR = 0.028·(Hk-1) - 0.0345·exp(-(3.87/(Hk-1) - 2.52)²)
    AF = -0.05 + 2.7/(Hk-1) - 5.5/(Hk-1)² + 3/(Hk-1)³
    AX = AF · DADR / θ
```

---

*Document version: 2.0 | Updated to reflect existing codebase | 2026-01-20*
