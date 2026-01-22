# XFOIL-Matching Inviscid Solver Rewrite Plan

## Problem Statement

RustFoil's inviscid solver produces a fundamentally different velocity distribution than XFOIL near the leading edge:

| Metric | XFOIL | RustFoil | Impact |
|--------|-------|----------|--------|
| Ue_max (upper) | 1.59 | 2.85 | 1.8x higher |
| Peak location | 3.41% chord | 0.88% chord | 3.9x closer to stagnation |
| Max dUe/ds | 86.7 | 757.5 | **8.7x steeper gradient** |

This causes the boundary layer solver to fail regardless of how correctly it's implemented.

## Goal

Recreate the inviscid solver to produce **identical** results to XFOIL's panel method.

---

## XFOIL Inviscid Solver Architecture

### Core Data Flow

```
Geometry (X,Y,S) ──► Panel Setup ──► GGCALC ──► SPECAL ──► Solution
                    ├─ NCALC (normals)     └─ GAMU(α=0,90)   ├─ GAM (combined)
                    ├─ APCALC (angles)                        ├─ QINV (velocity)
                    └─ TECALC (TE geometry)                   └─ CL, CM, Cp
```

### Key Insight: gamma = velocity

In XFOIL, the vortex strength γ at each node **equals the surface tangential velocity**:
```fortran
QINV(I) = COSA*GAMU(I,1) + SINA*GAMU(I,2)  ! Tangential velocity
```

---

## Module-by-Module Specification

### Module 1: Geometry Processing (`geometry.rs`)

**XFOIL Reference:** `xpanel.f` (NCALC, APCALC), `XFOIL.INC`

#### Input
- Raw airfoil coordinates: `Vec<(f64, f64)>` (N points, CCW from upper TE)

#### Output - `AirfoilGeometry` struct
```rust
pub struct AirfoilGeometry {
    // Node arrays (indexed 0..N-1)
    pub x: Vec<f64>,      // X coordinates
    pub y: Vec<f64>,      // Y coordinates
    pub s: Vec<f64>,      // Arc length from first point
    pub xp: Vec<f64>,     // dX/dS (spline derivative)
    pub yp: Vec<f64>,     // dY/dS (spline derivative)
    pub nx: Vec<f64>,     // Unit normal X component (outward)
    pub ny: Vec<f64>,     // Unit normal Y component (outward)
    pub apanel: Vec<f64>, // Panel angle (atan2 formulation)
    
    // Trailing edge info
    pub xte: f64,         // TE midpoint X
    pub yte: f64,         // TE midpoint Y
    pub dste: f64,        // TE gap length
    pub ante: f64,        // TE normal projected gap
    pub aste: f64,        // TE tangent projected gap
    pub sharp: bool,      // True if DSTE < 0.0001 * chord
    
    // Leading edge info
    pub xle: f64,         // LE X coordinate
    pub yle: f64,         // LE Y coordinate
    pub sle: f64,         // Arc length at LE
    
    // Dimensions
    pub n: usize,         // Number of nodes
    pub chord: f64,       // Chord length
}
```

#### Implementation Notes

**NCALC (Normal Calculation):**
```fortran
C     Calculates normal unit vector components at airfoil panel nodes
      CALL SEGSPL(X,XN,S,N)   ! XN = dX/dS via cubic spline
      CALL SEGSPL(Y,YN,S,N)   ! YN = dY/dS via cubic spline
      DO I=1, N
        SX =  YN(I)           ! Rotate 90° CCW
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        XN(I) = SX/SMOD       ! Outward normal
        YN(I) = SY/SMOD
      ENDDO
```

**APCALC (Panel Angle):**
```fortran
C     Panel i goes from node i to node i+1
      DO I=1, N-1
        SX = X(I+1) - X(I)
        SY = Y(I+1) - Y(I)
        APANEL(I) = ATAN2(SX, -SY)  ! Angle of panel normal
      ENDDO
C     TE panel (N to 1)
      APANEL(N) = ATAN2(-SX, SY) + PI
```

**State Validation:**
- `N >= 10` (minimum panels)
- `chord > 0`
- Coordinates form closed contour
- No duplicate sequential points

---

### Module 2: Influence Coefficients (`influence.rs`)

**XFOIL Reference:** `xpanel.f` (PSILIN subroutine, lines 99-800)

#### Core Function: `psilin`

**Purpose:** Calculate stream function ψ at point (xi, yi) due to all vortex panels, and the sensitivity ∂ψ/∂γ for each node.

**Input:**
```rust
fn psilin(
    geom: &AirfoilGeometry,
    i: usize,              // Node index (0..N) where we're evaluating
    xi: f64, yi: f64,      // Field point coordinates
    nxi: f64, nyi: f64,    // Unit normal at field point
    gam: &[f64],           // Vorticity distribution (N values)
    siglin: bool,          // Include source (viscous) effects
    sig: Option<&[f64]>,   // Source distribution (if siglin)
) -> PsilinResult
```

**Output:**
```rust
pub struct PsilinResult {
    pub psi: f64,          // Stream function value
    pub psi_ni: f64,       // ∂ψ/∂n (tangential velocity)
    pub dzdg: Vec<f64>,    // ∂ψ/∂γⱼ for each node j
    pub dzdm: Vec<f64>,    // ∂ψ/∂σⱼ for each node j (if siglin)
    pub qtan1: f64,        // Tangential velocity at α=0°
    pub qtan2: f64,        // Tangential velocity at α=90°
    pub z_alfa: f64,       // ∂ψ/∂α (alpha sensitivity)
}
```

**Algorithm (XFOIL's PSILIN):**

The key formula for linear vorticity panel contribution:

```
For each panel from node JO to JP with linearly varying γ:

1. Transform to panel-local coordinates:
   x1 = (P - P_JO) · tangent
   x2 = (P - P_JP) · tangent  
   yy = (P - P_JO) · normal

2. Compute potential integrals:
   g1 = ln(r1²), t1 = atan2(x1, yy)
   g2 = ln(r2²), t2 = atan2(x2, yy)
   
3. Sum/Difference formulation:
   PSIS = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2)  // Sum term
   PSID = ((x1+x2)*PSIS + 0.5*(r2²*g2 - r1²*g1 + x1² - x2²)) / (x1-x2)  // Diff term
   
4. Influence on ψ:
   ψ += QOPI * (PSIS*(γ_JO+γ_JP) + PSID*(γ_JP-γ_JO))
   
5. Derivative ∂ψ/∂γ:
   dzdg[JO] += QOPI * (PSIS - PSID)
   dzdg[JP] += QOPI * (PSIS + PSID)
```

where `QOPI = 1/(4π)`.

**TE Panel Special Treatment:**

For the TE panel (node N-1 to node 0), XFOIL uses a source+vortex formulation:
```fortran
SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))   ! Source strength
GAMTE = -0.5*SDS*(GAM(JP) - GAM(JO))  ! Vortex strength
PSI += HOPI*(PSIG*SIGTE + PGAM*GAMTE)  ! HOPI = 1/(2π)
```

where `SCS = ANTE/DSTE`, `SDS = ASTE/DSTE` for blunt TE (both = 1, 0 for sharp).

**State Validation:**
- Panel lengths > `SEPS = arc_length * 1e-5`
- No field point coinciding with panel endpoint (handle singularity)

---

### Module 3: Matrix Assembly & Solution (`system.rs`)

**XFOIL Reference:** `xpanel.f` (GGCALC), `xsolve.f` (LUDCMP, BAKSUB)

#### Function: `build_system`

**Purpose:** Build and solve the (N+1)×(N+1) linear system for vorticity.

**System Structure:**
```
| A₀₀  A₀₁  ...  A₀,N-1   -1 | | γ₀    |   | -ψ∞(0)   |
| A₁₀  A₁₁  ...  A₁,N-1   -1 | | γ₁    |   | -ψ∞(1)   |
|  :    :   ...    :      :  | |  :    | = |    :     |
| AN-1,0 ...   AN-1,N-1  -1  | | γN-1  |   | -ψ∞(N-1) |
|  1    0   ...   1       0  | | ψ₀    |   |    0     |
```

**Row i (i = 0..N-1):** Boundary condition ψ_induced + ψ_freestream = ψ₀
```rust
for j in 0..n {
    a[i][j] = dzdg[j];  // From psilin
}
a[i][n] = -1.0;  // Coefficient of ψ₀

// RHS for α=0°:  -V∞*y_i
// RHS for α=90°: +V∞*x_i
rhs_0[i] = -y[i];
rhs_90[i] = x[i];
```

**Row N:** Kutta condition γ₀ + γ_{N-1} = 0
```rust
a[n][0] = 1.0;
a[n][n-1] = 1.0;
rhs_0[n] = 0.0;
rhs_90[n] = 0.0;
```

**Sharp TE Enhancement (XFOIL):**
For sharp TE, XFOIL replaces row N-1 with internal velocity = 0:
```fortran
C     Control point on TE bisector just ahead of TE
      XBIS = XTE - BWT*DSMIN*CBIS   ! BWT = 0.1
      YBIS = YTE - BWT*DSMIN*SBIS
C     Set velocity component along bisector line = 0
      CALL PSILIN(0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,...)
      ! Row N-1: dQdG_j * γ_j = -V∞*(COSA*CBIS + SINA*SBIS)
```

**Output:**
```rust
pub struct FactorizedSystem {
    pub gamu_0: Vec<f64>,   // γ for α=0° (N values)
    pub gamu_90: Vec<f64>,  // γ for α=90° (N values)
    pub psi0_0: f64,        // ψ₀ for α=0°
    pub psi0_90: f64,       // ψ₀ for α=90°
}
```

**State Validation:**
- Matrix condition number reasonable
- LU factorization succeeded

---

### Module 4: Solution Computation (`solution.rs`)

**XFOIL Reference:** `xoper.f` (SPECAL, SPECCL, QISET, CLCALC, CPCALC)

#### Function: `compute_solution`

**Purpose:** Combine α=0° and α=90° base solutions for any angle of attack.

**Algorithm:**
```rust
// From SPECAL
let cosa = alpha.cos();
let sina = alpha.sin();

// Vorticity (= tangential velocity)
let gamma: Vec<f64> = (0..n)
    .map(|i| cosa * gamu_0[i] + sina * gamu_90[i])
    .collect();

// Gamma derivative w.r.t. alpha
let gamma_a: Vec<f64> = (0..n)
    .map(|i| -sina * gamu_0[i] + cosa * gamu_90[i])
    .collect();
```

#### Function: `qiset` - Set inviscid velocity QINV

**XFOIL:**
```fortran
DO I=1, N+NW
    QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
    QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
ENDDO
```

**Note:** For airfoil surface, `QINV(I) = GAM(I)` (vorticity = velocity).

#### Function: `clcalc` - Lift and moment calculation

**XFOIL Algorithm:**
```rust
// Integrate around closed contour (including TE panel from N-1 to 0)
let mut cl = 0.0;
let mut cm = 0.0;

for i in 0..n {
    let ip = (i + 1) % n;
    
    // Panel in wind axes
    let dx = (x[ip] - x[i]) * cosa + (y[ip] - y[i]) * sina;
    let dy = (y[ip] - y[i]) * cosa - (x[ip] - x[i]) * sina;
    
    // Cp = 1 - (γ/V∞)² with compressibility correction (skip for M=0)
    let cp_avg = 0.5 * (cp[i] + cp[ip]);
    let dcp = cp[ip] - cp[i];
    
    // Lift: ∫ Cp dx (in wind axes)
    cl += dx * cp_avg;
    
    // Moment about (xref, yref)
    let ax = (0.5*(x[ip]+x[i]) - xref) * cosa + (0.5*(y[ip]+y[i]) - yref) * sina;
    let ay = (0.5*(y[ip]+y[i]) - yref) * cosa - (0.5*(x[ip]+x[i]) - xref) * sina;
    
    cm -= dx * (cp_avg * ax + dcp * dx / 12.0);
    cm -= dy * (cp_avg * ay + dcp * dy / 12.0);
}
```

**Output:**
```rust
pub struct InviscidSolution {
    pub gamma: Vec<f64>,    // Surface vorticity (= velocity for Qinf=1)
    pub gamma_a: Vec<f64>,  // ∂γ/∂α
    pub qinv: Vec<f64>,     // Tangential velocity (= gamma for airfoil)
    pub cp: Vec<f64>,       // Pressure coefficient
    pub cl: f64,            // Lift coefficient
    pub cm: f64,            // Moment coefficient (about 0.25c)
    pub psi0: f64,          // Internal stream function
}
```

---

### Module 5: Stagnation Point Detection (`stagnation.rs`)

**XFOIL Reference:** `xpanel.f` (STFIND, IBLPAN, XICALC)

Already implemented. Verify against XFOIL:

#### Function: `stfind`

**XFOIL Algorithm:**
```fortran
DO I=1, N-1
    IF(GAM(I).GE.0.0 .AND. GAM(I+1).LT.0.0) GO TO 11
ENDDO
11 CONTINUE
IST = I

DGAM = GAM(I+1) - GAM(I)
DS = S(I+1) - S(I)

C---- minimize roundoff for very small GAM
IF(GAM(I) .LT. -GAM(I+1)) THEN
    SST = S(I)   - DS*(GAM(I)  /DGAM)
ELSE
    SST = S(I+1) - DS*(GAM(I+1)/DGAM)
ENDIF

C---- tweak if on node
IF(SST .LE. S(I)  ) SST = S(I)   + 1.0E-7
IF(SST .GE. S(I+1)) SST = S(I+1) - 1.0E-7
```

**Output:**
```rust
pub struct StagnationPoint {
    pub ist: usize,  // Panel index (stag between IST and IST+1)
    pub sst: f64,    // Exact arc length at stagnation
}
```

---

## Implementation Phases

### Phase 1: Geometry Module
- [ ] Port SEGSPL spline derivative calculation
- [ ] Port NCALC normal calculation  
- [ ] Port APCALC panel angle calculation
- [ ] Port TE geometry (ANTE, ASTE, DSTE, SHARP)
- [ ] Unit tests against XFOIL debug output

### Phase 2: Influence Coefficients
- [ ] Port PSILIN exactly (most critical!)
- [ ] Handle all edge cases (panel on/near point)
- [ ] TE panel special treatment (source + vortex)
- [ ] Validate against XFOIL AIJ matrix values

### Phase 3: System Assembly
- [ ] Build system matrix exactly per GGCALC
- [ ] Implement sharp TE bisector condition
- [ ] LU factorization matching XFOIL
- [ ] Solve for α=0° and α=90° base solutions

### Phase 4: Solution & Forces
- [ ] Combine base solutions per SPECAL/QISET
- [ ] Pressure coefficient calculation (CPCALC)
- [ ] Force integration (CLCALC)
- [ ] Validate CL, CM against XFOIL

### Phase 5: Validation
- [ ] Compare gamma at each node for NACA 0012, α=0°, 4°, 8°
- [ ] Compare velocity gradient dUe/ds near LE
- [ ] Verify stagnation point location matches
- [ ] Full comparison script with instrumented XFOIL

---

## Critical Differences to Check

1. **Indexing:** XFOIL uses 1-based (Fortran), Rust uses 0-based
2. **Panel direction:** XFOIL panel i goes from node i to node i+1; "TE panel" is from node N to node 1
3. **TE treatment:** Blunt TE uses source+vortex, sharp TE uses bisector velocity condition
4. **Tangent sign:** XFOIL's VTI array handles sign convention for upper/lower surface

---

## Test Cases

1. **NACA 0012, 160 panels (XFOIL-style), α=0°**
   - CL should be exactly 0
   - Gamma symmetric about LE
   
2. **NACA 0012, 160 panels, α=4°**
   - CL ≈ 0.458 (inviscid)
   - Peak Ue ≈ 1.59 at ≈3.4% chord
   - Max dUe/ds ≈ 87 (not 750!)

3. **NACA 0012, 160 panels, α=-4° to +4° sweep**
   - Verify polar symmetry: CL(-α) = -CL(α)
   - Lift curve slope ≈ 2π/rad

---

## Deliverables

1. `crates/rustfoil-inviscid/` - New clean inviscid solver crate
2. `docs/INVISCID_VALIDATION.md` - Validation report with XFOIL comparisons
3. Updated `rustfoil-cli` using new solver
4. Passing boundary layer analysis at α=4°

---

## References

- XFOIL 6.99 source: `Xfoil/src/xpanel.f`, `xoper.f`, `xsolve.f`
- XFOIL.INC common block definitions
- Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
- Katz & Plotkin, "Low-Speed Aerodynamics", Chapter 11
