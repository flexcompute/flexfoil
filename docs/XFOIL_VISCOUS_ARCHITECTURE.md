# XFOIL Viscous Solver Architecture: Complete Function Reference

## Purpose
This document provides a comprehensive mapping of every XFOIL function related to viscous flow solving. It is intended to serve as the specification for **FlexFoil** - a clean Rust reimplementation of XFOIL's viscous-inviscid coupling.

---

## Table of Contents
1. [Architecture Overview](#1-architecture-overview)
2. [Include Files (Data Structures)](#2-include-files-data-structures)
3. [Top-Level Solution Orchestration](#3-top-level-solution-orchestration)
4. [Boundary Layer System Setup](#4-boundary-layer-system-setup)
5. [Boundary Layer Equations](#5-boundary-layer-equations)
6. [Closure Relations](#6-closure-relations)
7. [Transition Prediction](#7-transition-prediction)
8. [Newton System Solver](#8-newton-system-solver)
9. [Panel Method Interface](#9-panel-method-interface)
10. [Utility Functions](#10-utility-functions)
11. [FlexFoil Implementation Mapping](#11-flexfoil-implementation-mapping)

---

## 1. Architecture Overview

### Data Flow Diagram
```
┌──────────────────────────────────────────────────────────────────────────────┐
│                           XFOIL VISCOUS SOLUTION                              │
└──────────────────────────────────────────────────────────────────────────────┘
                                     │
                              VISCAL (xoper.f)
                           Main iteration driver
                                     │
        ┌────────────────────────────┼────────────────────────────┐
        │                            │                            │
        ▼                            ▼                            ▼
   SETBL (xbl.f)              BLSOLV (xsolve.f)           UPDATE (xbl.f)
   Build Newton system        Solve block system         Apply changes
        │                                                     │
        │                                                     │
        ▼                                                     ▼
   BLSYS (xblsys.f)                                    STMOVE (xoper.f)
   Assemble equations                                  Move stag point
        │
        ├──── BLVAR (xblsys.f) ──── Closure relations
        │
        ├──── BLMID (xblsys.f) ──── Midpoint Cf
        │
        └──── BLDIF (xblsys.f) ──── Finite differences
```

### Key Files
| File | Purpose |
|------|---------|
| `XFOIL.INC` | Global state: geometry, panel data, BL arrays |
| `XBL.INC` | BL-specific: station variables, Newton system coefficients |
| `BLPAR.INC` | BL model parameters (constants) |
| `xoper.f` | High-level operations, `VISCAL` |
| `xbl.f` | BL setup: `SETBL`, `MRCHUE`, `MRCHDU`, `UPDATE` |
| `xblsys.f` | BL equations: `BLSYS`, `BLVAR`, closures, transition |
| `xsolve.f` | Linear solvers: `BLSOLV`, `GAUSS` |
| `xpanel.f` | Panel method: influence coefficients |

---

## 2. Include Files (Data Structures)

### 2.1 `XFOIL.INC` - Global State

#### Geometry Arrays
```fortran
X(IZX), Y(IZX)        ! Airfoil + wake coordinates (i=1..N, N+1..N+NW)
XP(IZX), YP(IZX)      ! Spline derivatives dX/dS, dY/dS
S(IZX)                ! Arc length parameter
NX(IZX), NY(IZX)      ! Unit normal vectors
APANEL(IZX)           ! Panel angles
```

#### Panel Method Arrays
```fortran
GAM(IQX)              ! Vortex panel strengths (γ)
SIG(IZX)              ! Source/sink panel strengths (σ = mass defect)
QINV(IZX)             ! Inviscid tangential velocity
QVIS(IZX)             ! Viscous tangential velocity
CPI(IZX), CPV(IZX)    ! Inviscid/viscous pressure coefficients
```

#### Influence Matrices
```fortran
AIJ(IQX,IQX)          ! dPsi/dGam influence matrix
DIJ(IZX,IZX)          ! dQtan/dSig influence matrix
```

#### Boundary Layer Arrays (per side, per station)
```fortran
XSSI(IVX,ISX)         ! BL arc length coordinate
UEDG(IVX,ISX)         ! BL edge velocity
UINV(IVX,ISX)         ! Inviscid edge velocity (no mass defect)
MASS(IVX,ISX)         ! Mass defect = Ue * δ*
THET(IVX,ISX)         ! Momentum thickness θ
DSTR(IVX,ISX)         ! Displacement thickness δ*
CTAU(IVX,ISX)         ! Shear coefficient √(τmax/ρUe²) OR amplification N
```

#### Newton System Arrays
```fortran
VA(3,2,IZX)           ! Diagonal block coefficients
VB(3,2,IZX)           ! Off-diagonal block coefficients
VM(3,IZX,IZX)         ! Mass influence vectors
VZ(3,2)               ! TE coupling block
VDEL(3,2,IZX)         ! Residuals (col 1) and sensitivities (col 2)
```

#### Key Scalars
```fortran
N                     ! Number of airfoil panels
NW                    ! Number of wake panels
IST                   ! Stagnation point panel index
ALFA                  ! Angle of attack (radians)
CL, CD, CM            ! Lift, drag, moment coefficients
MINF                  ! Freestream Mach number
REINF                 ! Reynolds number
ACRIT(ISX)            ! Critical amplification ratio (e^N method)
ITRAN(ISX)            ! Transition station index (per side)
```

**FlexFoil Mapping:**
```rust
pub struct XfoilState {
    // Geometry
    pub x: Vec<f64>,           // Panel x-coordinates
    pub y: Vec<f64>,           // Panel y-coordinates
    pub s: Vec<f64>,           // Arc length
    pub nx: Vec<f64>,          // Normal x-component
    pub ny: Vec<f64>,          // Normal y-component
    
    // Panel method
    pub gam: Vec<f64>,         // Vortex strengths
    pub sig: Vec<f64>,         // Source strengths (mass defect)
    pub qinv: Vec<f64>,        // Inviscid velocity
    pub qvis: Vec<f64>,        // Viscous velocity
    
    // BL state (indexed by [side][station])
    pub bl: [BlSide; 2],
    
    // Solution state
    pub alfa: f64,
    pub cl: f64,
    pub cd: f64,
    pub minf: f64,
    pub reinf: f64,
}

pub struct BlSide {
    pub xi: Vec<f64>,          // Arc length from stagnation
    pub ue: Vec<f64>,          // Edge velocity
    pub theta: Vec<f64>,       // Momentum thickness
    pub dstar: Vec<f64>,       // Displacement thickness
    pub ctau: Vec<f64>,        // Shear coeff (turb) or amp (lam)
    pub itran: usize,          // Transition index
}
```

---

### 2.2 `XBL.INC` - Station Variables

```fortran
! Station "1" (upstream) and "2" (downstream) variables
! Suffix meanings: _T = ∂/∂θ, _D = ∂/∂δ*, _U = ∂/∂Ue, _MS = ∂/∂M², _RE = ∂/∂Re

! Primary variables
X1, X2                ! Arc length position
U1, U2                ! Edge velocity (compressible)
T1, T2                ! Momentum thickness θ
D1, D2                ! Displacement thickness δ*
S1, S2                ! Ctau (turbulent) or amplification N (laminar)

! Shape parameters
H1, H2                ! H = δ*/θ
HK1, HK2              ! Kinematic shape parameter Hk
HS1, HS2              ! Energy shape parameter H*
HC1, HC2              ! Density shape parameter H**

! Reynolds numbers
RT1, RT2              ! Momentum thickness Re_θ = ρUeθ/μ
M1, M2                ! Local Mach number

! Skin friction and dissipation
CF1, CF2              ! Skin friction coefficient Cf
DI1, DI2              ! Dissipation coefficient 2CD/H*
CQ1, CQ2              ! Equilibrium shear coefficient

! Newton system coefficients (4 equations × 5 unknowns)
VS1(4,5)              ! Jacobian wrt station 1: [dCt/dTh/dDs][dA/dT/dD/dU/dX]
VS2(4,5)              ! Jacobian wrt station 2
VSREZ(4)              ! Residual vector
VSM(4)                ! Mach sensitivity
VSR(4)                ! Reynolds sensitivity
VSX(4)                ! Xi (forced transition) sensitivity
```

**FlexFoil Mapping:**
```rust
/// BL station state with all derivatives for Newton system
#[derive(Clone, Default)]
pub struct BlStation {
    // Primary variables
    pub x: f64,                // Arc length
    pub ue: f64,               // Edge velocity
    pub theta: f64,            // Momentum thickness
    pub dstar: f64,            // Displacement thickness
    pub ctau: f64,             // Shear coeff or amplification
    
    // Derived quantities with derivatives
    pub hk: Derivative,        // Kinematic shape parameter
    pub hs: Derivative,        // Energy shape parameter
    pub ret: Derivative,       // Re_theta
    pub cf: Derivative,        // Skin friction
    pub cd: Derivative,        // Dissipation
    pub cq: Derivative,        // Equilibrium Ctau
}

/// Value with partial derivatives
#[derive(Clone, Default)]
pub struct Derivative {
    pub val: f64,
    pub theta: f64,            // ∂/∂θ
    pub dstar: f64,            // ∂/∂δ*
    pub ue: f64,               // ∂/∂Ue
    pub ctau: f64,             // ∂/∂Ctau (if applicable)
}
```

---

### 2.3 `BLPAR.INC` - Model Constants

```fortran
COMMON /BLPAR/
     & SCCON,    ! Shear lag constant = 5.6
     & GACON,    ! G-beta constant = 6.70
     & GBCON,    ! G-beta constant = 0.75
     & GCCON,    ! G-beta wall term = 18.0
     & DLCON,    ! Wake dissipation ratio = 0.9
     & CTRCON,   ! Ctau transition constant = 1.8
     & CTRCEX,   ! Ctau transition exponent = 3.3
     & DUXCON,   ! dUe/dx weighting = 1.0
     & CTCON,    ! Ctau equilibrium constant = 0.5/(GACON² × GBCON)
     & CFFAC     ! Cf correction factor = 1.0
```

**FlexFoil Mapping:**
```rust
pub struct BlParameters {
    pub sccon: f64,            // 5.6 - shear lag
    pub gacon: f64,            // 6.70 - G-beta
    pub gbcon: f64,            // 0.75 - G-beta
    pub gccon: f64,            // 18.0 - wall term
    pub dlcon: f64,            // 0.9 - wake dissipation
    pub ctrcon: f64,           // 1.8 - transition Ctau
    pub ctrcex: f64,           // 3.3 - transition exponent
    pub ctcon: f64,            // derived
    pub cffac: f64,            // 1.0 - Cf factor
}

impl Default for BlParameters {
    fn default() -> Self {
        let gacon = 6.70;
        let gbcon = 0.75;
        Self {
            sccon: 5.6,
            gacon,
            gbcon,
            gccon: 18.0,
            dlcon: 0.9,
            ctrcon: 1.8,
            ctrcex: 3.3,
            ctcon: 0.5 / (gacon * gacon * gbcon),
            cffac: 1.0,
        }
    }
}
```

---

## 3. Top-Level Solution Orchestration

### 3.1 `VISCAL` (xoper.f:2886)

**Purpose:** Main viscous solution iteration driver

**Algorithm:**
```
1. Initialize wake geometry if needed (XYWAKE)
2. Set wake velocities (QWCALC)
3. Set initial velocity distribution (QISET)
4. Set up BL-panel pointers if needed (STFIND, IBLPAN, XICALC, IBLSYS)
5. Calculate inviscid edge velocities (UICALC)
6. Newton iteration loop:
   a. Build BL Newton system (SETBL)
   b. Solve block system (BLSOLV)
   c. Apply updates with relaxation (UPDATE)
   d. Update Mach/Re or alpha (MRCL/QISET)
   e. Calculate viscous velocities (QVFUE)
   f. Update vorticity from velocities (GAMQV)
   g. Move stagnation point (STMOVE)
   h. Calculate CL, CD (CLCALC, CDCALC)
   i. Check convergence (RMS < 1e-4)
```

**Inputs:**
- `NITER1`: Maximum iterations
- Global state in `XFOIL.INC`

**Outputs:**
- Updated `UEDG`, `THET`, `DSTR`, `CTAU` arrays
- `CL`, `CD`, `CM` coefficients
- `LVCONV`: Convergence flag

**FlexFoil Equivalent:**
```rust
impl FlexFoilSolver {
    pub fn solve_viscous(&mut self, max_iter: usize) -> Result<Solution, SolverError> {
        // 1. Initialize wake and velocities
        self.init_wake()?;
        self.compute_inviscid_ue()?;
        
        // 2. Newton iteration
        for iter in 0..max_iter {
            // Build system
            let system = self.build_bl_system()?;
            
            // Solve
            let deltas = self.solve_bl_system(&system)?;
            
            // Update with relaxation
            let rlx = self.compute_relaxation(&deltas);
            self.apply_update(&deltas, rlx)?;
            
            // Update coupling
            self.update_viscous_velocities()?;
            self.move_stagnation_point()?;
            
            // Check convergence
            let rms = self.compute_rms(&deltas);
            if rms < 1e-4 {
                return Ok(self.extract_solution());
            }
        }
        Err(SolverError::NotConverged)
    }
}
```

---

### 3.2 `SETBL` (xbl.f:21)

**Purpose:** Builds the complete BL Newton system by marching along both surfaces

**Algorithm:**
```
1. Set compressibility parameters (MRCL, COMSET)
2. Initialize BL if first call (MRCHUE)
3. March BL with current Ue and δ* (MRCHDU)
4. For each side IS = 1,2:
   For each station IBL = 2..NBL(IS):
     a. Set station flags (SIMI, WAKE, TRAN, TURB)
     b. Load primary variables (θ, δ*, Ctau, Ue)
     c. Compute secondary variables (BLPRV, BLKIN)
     d. Check for transition (TRCHEK)
     e. Assemble local system:
        - First wake station: TESYS
        - Other stations: BLSYS
     f. Insert into global Jacobian (VS1, VS2 → VA, VB, VM)
```

**Key Logic - Transition handling:**
```fortran
IF(TRAN) THEN
    CALL TRCHEK       ! Check transition, set XT
    AMI = AMPL2       ! Update amplification
ENDIF
```

**Key Logic - Wake starting point:**
```fortran
IF(IBL.EQ.IBLTE(IS)+1) THEN
    TTE = THET(top) + THET(bottom)  ! Combined TE theta
    DTE = DSTR(top) + DSTR(bottom) + ANTE
    CALL TESYS(CTE,TTE,DTE)
ENDIF
```

**FlexFoil Equivalent:**
```rust
fn build_bl_system(&self) -> BlNewtonSystem {
    let mut system = BlNewtonSystem::new(self.nsys);
    
    for side in 0..2 {
        let mut prev_station = BlStation::default();
        
        for ibl in 1..self.nbl[side] {
            // Set flags
            let is_similarity = ibl == 1;
            let is_wake = ibl > self.iblte[side];
            let is_transition = ibl == self.itran[side];
            let is_turbulent = ibl >= self.itran[side];
            
            // Load and compute station
            let station = self.compute_station(side, ibl);
            
            // Check transition
            if !is_turbulent && !is_similarity {
                self.check_transition(&prev_station, &station);
            }
            
            // Assemble local equations
            let (residual, jacobian) = if is_wake && ibl == self.iblte[side] + 1 {
                self.te_wake_equations(&prev_station, &station)
            } else {
                self.bl_equations(&prev_station, &station, is_turbulent, is_wake)
            };
            
            // Insert into global system
            system.insert(ibl, side, residual, jacobian);
            
            prev_station = station;
        }
    }
    system
}
```

---

### 3.3 `MRCHUE` (xbl.f:542)

**Purpose:** Initialize BL by marching in direct mode, switching to inverse when separation detected

**Critical Constants:**
```fortran
HLMAX = 3.8    ! Laminar separation Hk threshold
HTMAX = 2.5    ! Turbulent separation Hk threshold
```

**Algorithm:**
```
For each side:
  1. Initialize at similarity station with Thwaites:
     θ = sqrt(0.45·x/(Ue·(5β+1)·Re))  where β = d(ln Ue)/d(ln x) ≈ 1
     δ* = 2.2·θ
     
  2. For each downstream station:
     a. Try DIRECT mode: prescribe Ue, solve for (θ, δ*)
     b. Compute resulting Hk
     c. IF Hk > threshold:
        - Switch to INVERSE mode
        - Prescribe HTARG (extrapolated from upstream)
        - Solve for Ue (in addition to θ, δ*)
     d. Store results
```

**CRITICAL: Inverse Mode Logic**
```fortran
IF(DIRECT) THEN
    ! Try direct mode (Ue prescribed)
    VS2(4,1) = 0.
    VS2(4,2) = 0.
    VS2(4,3) = 0.
    VS2(4,4) = 1.0        ! dUe = 0
    VSREZ(4) = 0.
    CALL GAUSS(4,4,VS2,VSREZ,1)
    
    ! Check if Hk exceeded threshold
    HKTEST = compute_hk(DSI + RLX*VSREZ(3), THI + RLX*VSREZ(2))
    IF(HKTEST > HMAX) THEN
        DIRECT = .FALSE.
        ! Set target Hk for inverse
        HTARG = HK1 + 0.03*(X2-X1)/T1  ! laminar
        HTARG = HK1 - 0.15*(X2-X1)/T1  ! turbulent
    ENDIF
ELSE
    ! Inverse mode (Hk prescribed at HTARG)
    VS2(4,1) = 0.
    VS2(4,2) = HK2_T2
    VS2(4,3) = HK2_D2
    VS2(4,4) = HK2_U2
    VSREZ(4) = HTARG - HK2  ! Residual: drive Hk to target
    CALL GAUSS(4,4,VS2,VSREZ,1)
    
    ! Update includes Ue change:
    UEI = UEI + RLX*VSREZ(4)  ! ← THIS IS THE KEY TO STALL PREDICTION
ENDIF
```

**FlexFoil Equivalent:**
```rust
fn march_direct_inverse(&mut self, side: usize) {
    const HLMAX: f64 = 3.8;  // Laminar separation threshold
    const HTMAX: f64 = 2.5;  // Turbulent separation threshold
    
    // Initialize at stagnation
    let mut theta = self.thwaites_initial_theta(side);
    let mut dstar = 2.2 * theta;
    let mut ctau = 0.0;
    let mut ue = self.ue_inviscid[side][1];
    
    for ibl in 1..self.nbl[side] {
        let is_laminar = ibl < self.itran[side];
        let is_wake = ibl > self.iblte[side];
        let hk_max = if is_laminar { HLMAX } else { HTMAX };
        
        // Try direct mode first
        let mut direct_mode = true;
        let ue_prescribed = self.ue_inviscid[side][ibl];
        
        for _newton_iter in 0..25 {
            // Compute BL equations
            let (residual, jacobian) = self.bl_equations_local(
                theta, dstar, ue, ctau, is_laminar, is_wake
            );
            
            if direct_mode {
                // Direct: prescribe Ue
                let (dtheta, ddstar, dctau) = solve_direct(&jacobian, &residual);
                
                // Check if separation would occur
                let hk_new = self.compute_hk(dstar + ddstar, theta + dtheta);
                
                if hk_new > hk_max && !is_wake {
                    // Switch to inverse mode
                    direct_mode = false;
                    let htarg = self.extrapolate_hk_target(ibl, side, is_laminar);
                    continue;  // Re-solve in inverse mode
                }
                
                // Apply update
                theta += dtheta;
                dstar += ddstar;
                ctau += dctau;
                ue = ue_prescribed;
                
            } else {
                // INVERSE: prescribe Hk, solve for Ue
                let htarg = self.compute_hk_target(ibl, side, is_laminar, is_wake);
                let (dtheta, ddstar, dctau, due) = solve_inverse(
                    &jacobian, &residual, htarg
                );
                
                // Apply update INCLUDING Ue change
                theta += dtheta;
                dstar += ddstar;
                ctau += dctau;
                ue += due;  // ← KEY: Ue is reduced in separated regions
            }
            
            if converged { break; }
        }
        
        // Store results
        self.theta[side][ibl] = theta;
        self.dstar[side][ibl] = dstar;
        self.ue[side][ibl] = ue;
        self.ctau[side][ibl] = ctau;
    }
}
```

---

### 3.4 `MRCHDU` (xbl.f:875)

**Purpose:** March BL in mixed mode (quasi-normal to Ue-Hk characteristic) to avoid Goldstein singularity

**Key Parameter:**
```fortran
SENSWT = 1000.0  ! Controls how far Hk can deviate from baseline
```

**Algorithm:**
```
For each station:
  1. Compute Ue-Hk characteristic slope (dUe/dHk)
  2. Prescribe a combination: Hk_target + SENSWT*(Ue_target/Ue_ref - 1)
  3. This keeps solution on a line quasi-normal to characteristic
  4. Avoids separation singularity while maintaining physicality
```

---

### 3.5 `UPDATE` (xbl.f:1253)

**Purpose:** Apply Newton deltas with under-relaxation and clamping

**Algorithm:**
```
1. Compute new Ue distribution from mass defect influence:
   Ue_new = Ue_inviscid + Σ DIJ(i,j) × MASS(j)
   
2. Compute new CL from Ue:
   CL = ∮ Cp dx  where  Cp = 1 - (Qvis/Qinf)²
   
3. Set relaxation factor:
   - Limit alpha change: |Δα| < 0.5°
   - Limit CL change: |ΔCL| < 0.5
   - Limit BL variable changes: |Δθ/θ|, |Δδ*/δ*|, |ΔUe| < 1.5
   
4. Apply relaxed updates:
   θ_new = θ_old + RLX × Δθ
   δ*_new = δ*_old + RLX × Δδ*
   Ue_new = Ue_old + RLX × ΔUe
   
5. Clamp Hk to prevent unphysical values:
   Hk_min = 1.02 (airfoil) or 1.00005 (wake)
```

**FlexFoil Equivalent:**
```rust
fn apply_update(&mut self, deltas: &BlDeltas, rlx: f64) {
    for side in 0..2 {
        for ibl in 1..self.nbl[side] {
            let d = &deltas[side][ibl];
            
            // Apply relaxed updates
            self.ctau[side][ibl] += rlx * d.ctau;
            self.theta[side][ibl] += rlx * d.theta;
            self.dstar[side][ibl] += rlx * d.dstar;
            self.ue[side][ibl] += rlx * d.ue;
            
            // Clamp Ctau
            if ibl >= self.itran[side] {
                self.ctau[side][ibl] = self.ctau[side][ibl].clamp(1e-7, 0.25);
            }
            
            // Clamp Hk (prevent Hk < 1)
            let hk_min = if ibl <= self.iblte[side] { 1.02 } else { 1.00005 };
            let hk = self.dstar[side][ibl] / self.theta[side][ibl];
            if hk < hk_min {
                self.dstar[side][ibl] = hk_min * self.theta[side][ibl];
            }
            
            // Update mass defect
            self.mass[side][ibl] = self.dstar[side][ibl] * self.ue[side][ibl];
        }
    }
}
```

---

## 4. Boundary Layer System Setup

### 4.1 `BLSYS` (xblsys.f:583)

**Purpose:** Assembles the 3-equation Newton system for one BL interval

**Equations:**
1. **Shear lag** (turbulent) or **Amplification** (laminar)
2. **Momentum integral**
3. **Shape parameter (kinetic energy)**

**Algorithm:**
```
1. Compute secondary variables: BLVAR(ityp) where ityp = 1/2/3 = lam/turb/wake
2. Compute midpoint Cf: BLMID(ityp)
3. Set up finite differences: BLDIF(ityp) or TRDIF (transition)
4. Convert to incompressible Ue and Mach
```

---

### 4.2 `BLDIF` (xblsys.f:1552)

**Purpose:** Sets up finite-difference Newton coefficients for one interval

**Equation 1 - Laminar Amplification:**
```
N₂ - N₁ = AX × (x₂ - x₁)

where AX = envelope amplification rate
```

**Equation 1 - Turbulent Shear Lag:**
```
SCC × (Cq - Ct×ALD) × Δx - δe × 2 × ln(S₂/S₁) + δe × 2 × (Uq×Δx - ln(U₂/U₁)) = 0

where:
  SCC = 5.6×1.333/(1+Us)     (shear lag constant)
  ALD = 1.0 (airfoil) or 0.9 (wake)
  Uq = equilibrium dUe/dx
```

**Equation 2 - Momentum:**
```
ln(θ₂/θ₁) + (H + 2 - M + Hw) × ln(U₂/U₁) - ln(x₂/x₁) × Cf×x/2θ = 0
```

**Equation 3 - Shape Parameter:**
```
ln(H*₂/H*₁) + (2H**/H* + 1 - H - Hw) × ln(U₂/U₁) + ln(x₂/x₁) × (Cf/2 - CD/H*) × x/θ = 0
```

**FlexFoil Equivalent:**
```rust
fn bl_equations(&self, st1: &BlStation, st2: &BlStation, turb: bool, wake: bool) 
    -> (Vec3, Mat3x5, Mat3x5) 
{
    let dx = st2.x - st1.x;
    let ulog = (st2.ue / st1.ue).ln();
    let tlog = (st2.theta / st1.theta).ln();
    let xlog = (st2.x / st1.x).ln();
    
    let mut res = Vec3::zeros();
    let mut jac1 = Mat3x5::zeros();  // d/d(station 1)
    let mut jac2 = Mat3x5::zeros();  // d/d(station 2)
    
    if turb {
        // Equation 1: Shear lag
        let scc = 5.6 * 1.333 / (1.0 + 0.5 * (st1.us + st2.us));
        let ald = if wake { 0.9 } else { 1.0 };
        let cqa = 0.5 * (st1.cq + st2.cq);
        let sa = 0.5 * (st1.ctau + st2.ctau);
        let dea = 0.5 * (st1.de + st2.de);
        
        res[0] = scc * (cqa - sa * ald) * dx 
               - dea * 2.0 * (st2.ctau / st1.ctau).ln()
               + dea * 2.0 * (self.uq(st1, st2) * dx - ulog);
        // ... linearization
    } else {
        // Equation 1: Amplification
        let ax = self.amplification_rate(st1, st2);
        res[0] = st2.ampl - st1.ampl - ax * dx;
        // ... linearization
    }
    
    // Equation 2: Momentum
    let ha = 0.5 * (st1.h + st2.h);
    let ma = 0.5 * (st1.mach + st2.mach);
    let cfx = self.cf_weighted(st1, st2);
    res[1] = tlog + (ha + 2.0 - ma) * ulog - xlog * 0.5 * cfx;
    
    // Equation 3: Shape parameter
    let hsa = 0.5 * (st1.hs + st2.hs);
    let hca = 0.5 * (st1.hc + st2.hc);
    let dix = self.cd_weighted(st1, st2);
    res[2] = (st2.hs / st1.hs).ln() 
           + (2.0 * hca / hsa + 1.0 - ha) * ulog 
           + xlog * (0.5 * cfx - dix);
    
    (res, jac1, jac2)
}
```

---

### 4.3 `TRDIF` (xblsys.f:1195)

**Purpose:** Handles transition interval where both laminar and turbulent equations apply

**Algorithm:**
```
1. Interpolate variables to transition point XT using weighting factors
2. Solve laminar equations from X1 to XT
3. Initialize turbulent Ctau at transition:
   Ct = CTRCON × exp(-CTRCEX/(Hk-1)) × Cq_eq
4. Solve turbulent equations from XT to X2
5. Sum the two contributions
```

---

### 4.4 `TESYS` (xblsys.f:664)

**Purpose:** Sets up dummy system between airfoil TE and first wake point

**Equations:** Simple continuity of variables:
```
Ct_wake - Ct_TE = 0
θ_wake - θ_TE = 0
(δ* + DW)_wake - δ*_TE = 0
```

---

## 5. Boundary Layer Equations

### 5.1 `BLPRV` (xblsys.f:701)

**Purpose:** Sets primary "2" variables from input parameters

```fortran
X2 = XSI
T2 = THI
D2 = DSI - DSWAKI  ! Subtract wake displacement
DW2 = DSWAKI

! Compressibility correction (Prandtl-Glauert-like)
U2 = UEI*(1-TKLAM) / (1 - TKLAM*(UEI/QINF)²)
```

**FlexFoil:**
```rust
fn set_station_primary(&mut self, xi: f64, theta: f64, dstar: f64, ue: f64, dswake: f64) {
    self.x = xi;
    self.theta = theta;
    self.dstar = dstar - dswake;
    self.dw = dswake;
    
    // Compressibility correction
    let tk = self.tklam;
    let ue_ratio = ue / self.qinf;
    self.ue_comp = ue * (1.0 - tk) / (1.0 - tk * ue_ratio * ue_ratio);
}
```

---

### 5.2 `BLKIN` (xblsys.f:725)

**Purpose:** Calculates turbulence-independent secondary variables

**Calculations:**
```fortran
! Edge Mach number
M2 = U2²×HSTINV / (γ₁×(1 - 0.5×U2²×HSTINV))

! Edge density (isentropic)
R2 = RST × (1 + 0.5×γ₁×M²)^(-1/γ₁)

! Shape parameter
H2 = D2/T2

! Molecular viscosity
V2 = sqrt(HERAT³) × (1+HVRAT)/(HERAT+HVRAT) / REYBL

! Kinematic shape parameter
HK2 = HKIN(H2, M2)

! Momentum thickness Reynolds number
RT2 = R2×U2×T2/V2
```

**FlexFoil:**
```rust
fn compute_secondary(&self, theta: f64, dstar: f64, ue: f64) -> BlSecondary {
    let h = dstar / theta;
    
    // Compressible quantities
    let msq = self.edge_mach_squared(ue);
    let rho = self.edge_density(msq);
    let mu = self.edge_viscosity(ue);
    
    // Kinematic shape parameter
    let hk = self.hkin(h, msq);
    
    // Reynolds number
    let ret = rho * ue * theta / mu;
    
    BlSecondary { h, hk, msq, rho, ret }
}
```

---

### 5.3 `BLVAR` (xblsys.f:784)

**Purpose:** Calculates ALL secondary "2" variables including turbulence-dependent ones

**Called with:**
- `ITYP = 1`: Laminar
- `ITYP = 2`: Turbulent
- `ITYP = 3`: Wake

**Calculates:**
```
Hk clamping:   Hk ≥ 1.05 (airfoil), Hk ≥ 1.00005 (wake)
HC (H**)   :   Density shape parameter from HCT()
HS (H*)    :   Energy shape parameter from HSL()/HST()
US         :   Normalized slip velocity
CQ         :   Equilibrium shear coefficient
CF         :   Skin friction from CFL()/CFT()
DI (2CD/H*):   Dissipation from DIL()/DIT()
DE (δ)     :   BL thickness from Green's correlation
```

---

### 5.4 `BLMID` (xblsys.f:1124)

**Purpose:** Calculates midpoint skin friction CFM for momentum equation accuracy

```fortran
HKA = 0.5*(HK1 + HK2)
RTA = 0.5*(RT1 + RT2)
MA  = 0.5*(M1  + M2)

IF(ITYP.EQ.1) CALL CFL(HKA, RTA, MA, CFM, ...)
IF(ITYP.EQ.2) CALL CFT(HKA, RTA, MA, CFM, ...)
IF(ITYP.EQ.3) CFM = 0.0  ! No friction in wake
```

---

## 6. Closure Relations

### 6.1 `HKIN` (xblsys.f:2276)

**Purpose:** Kinematic shape parameter (Whitfield)

```fortran
Hk = (H - 0.29×M²) / (1 + 0.113×M²)
```

**FlexFoil:**
```rust
fn hkin(h: f64, msq: f64) -> f64 {
    (h - 0.29 * msq) / (1.0 + 0.113 * msq)
}
```

---

### 6.2 `HSL` (xblsys.f:2327) - Laminar H*

```fortran
IF(Hk < 4.35) THEN
    H* = 0.0111×(Hk-4.35)²/(Hk+1) - 0.0278×(Hk-4.35)³/(Hk+1) + 1.528
ELSE
    H* = 0.015×(Hk-4.35)²/Hk + 1.528
ENDIF
```

---

### 6.3 `HST` (xblsys.f:2388) - Turbulent H*

```fortran
! With Re_θ dependence
IF(Rθ > 400) THEN
    Ho = 3 + 400/Rθ
ELSE
    Ho = 4
ENDIF

IF(Hk < Ho) THEN
    ! Attached branch (Swafford profiles)
    Hr = (Ho - Hk)/(Ho - 1)
    H* = (2 - H*_min - 4/Rθ)×Hr²×1.5/(Hk+0.5) + H*_min + 4/Rθ
ELSE
    ! Separated branch
    H* = (Hk - Ho)²×[0.007×ln(Rθ)/(Hk-Ho+4/ln(Rθ))² + 0.015/Hk] + H*_min + 4/Rθ
ENDIF

! Compressibility correction (Whitfield)
H* = (H* + 0.028×M²) / (1 + 0.014×M²)
```

---

### 6.4 `CFL` (xblsys.f:2354) - Laminar Cf

```fortran
! Falkner-Skan correlation
IF(Hk < 5.5) THEN
    Cf = [0.0727×(5.5-Hk)³/(Hk+1) - 0.07] / Rθ
ELSE
    Cf = [0.015×(1 - 1/(Hk-4.5))² - 0.07] / Rθ
ENDIF
```

---

### 6.5 `CFT` (xblsys.f:2483) - Turbulent Cf

```fortran
! Coles correlation
Fc = sqrt(1 + 0.5×(γ-1)×M²)
Grt = max(ln(Rθ/Fc), 3)
Gex = -1.74 - 0.31×Hk
Cfo = 0.3×exp(-1.33×Hk) × (Grt/2.3)^Gex
Cf = [Cfo + 1.1e-4×(tanh(4 - Hk/0.875) - 1)] / Fc
```

---

### 6.6 `DIL` (xblsys.f:2290) - Laminar Dissipation

```fortran
! Falkner-Skan dissipation
IF(Hk < 4) THEN
    2CD/H* = [0.00205×(4-Hk)^5.5 + 0.207] / Rθ
ELSE
    2CD/H* = [-0.0016×(Hk-4)²/(1+0.02×(Hk-4)²) + 0.207] / Rθ
ENDIF
```

---

### 6.7 `HCT` (xblsys.f:2514) - Density Shape Parameter

```fortran
H** = M² × (0.064/(Hk-0.8) + 0.251)
```

---

## 7. Transition Prediction

### 7.1 `TRCHEK` / `TRCHEK2` (xblsys.f:22, 231)

**Purpose:** Check for transition using e^N envelope method

**Algorithm (2nd order):**
```
1. Calculate amplification rate AX at both stations
2. Solve implicit equation for N₂:
   N₂ - N₁ = 0.5×(N'(XT,NT) + N'(X1,N1)) × (X₂-X₁)
3. If N₂ ≥ N_crit → transition occurs
4. Locate transition point XT by interpolation
```

---

### 7.2 `DAMPL` / `DAMPL2` (xblsys.f:1981, 2099)

**Purpose:** Envelope amplification rate dN/dx

**Algorithm:**
```
1. Calculate critical Rθ from H:
   log₁₀(Rθ_crit) = 2.492/(Hk-1)^0.43 + 0.7×(tanh(14/(Hk-1) - 9.24) + 1)
   
2. If Rθ < Rθ_crit → no amplification (AX = 0)

3. Otherwise:
   dN/dRθ = 0.028×(Hk-1) - 0.0345×exp(-(3.87/(Hk-1) - 2.52)²)
   
   θ×dRθ/dx = -0.05 + 2.7/(Hk-1) - 5.5/(Hk-1)² + 3/(Hk-1)³
   
   AX = (dN/dRθ) × (θ×dRθ/dx) / θ
```

**DAMPL2 modification for separated profiles (Hk > 4):**
```
For Hk > 4, blend to O-S maximum ai(H,Rθ) function
```

---

## 8. Newton System Solver

### 8.1 `BLSOLV` (xsolve.f:283)

**Purpose:** Custom block-elimination solver for the coupled system

**System Structure:**
```
[A  |  |  .  |  |  .  |][d1]   [R1]
[B  A  |  .  |  |  .  |][d2]   [R2]
[|  B  A  .  |  |  .  |][d3] = [R3]
[.  .  .  .  |  |  .  |][. ]   [. ]
[|  Z  |  |  B  A  .  |][dn]   [Rn]
[|  |  |  |  |  |  B  A][d ]   [R ]

where:
  A, B = 3×2 BL equation Jacobians (diagonal, off-diagonal)
  Z = 3×2 TE coupling (links wake to upper surface)
  | = 3×1 mass defect influence vectors
  d = [dCtau, dθ, dm] unknowns
  R = residuals
```

**Algorithm:**
```
Forward sweep (IV = 1 to NSYS):
  1. Invert VA(IV) block (3×3 with mass column)
  2. Eliminate VB(IV+1) block
  3. Eliminate VZ block at TE station
  4. Eliminate lower VM columns (with acceleration threshold)

Backward sweep (IV = NSYS to 2):
  1. Back-substitute upper VM columns
```

**FlexFoil:**
```rust
fn solve_bl_system(&self, system: &BlNewtonSystem) -> Vec<BlDelta> {
    // Use custom block-elimination or sparse LU
    // The key is handling the mass-influence coupling VM
    
    let mut deltas = vec![BlDelta::default(); system.nsys];
    
    // Forward elimination with block structure
    for iv in 0..system.nsys {
        // Invert 3x3 diagonal block with pivoting
        // Eliminate off-diagonal blocks
        // Handle mass-influence vectors (sparse elimination)
    }
    
    // Backward substitution
    for iv in (0..system.nsys).rev() {
        // Backsolve mass-influence columns
    }
    
    deltas
}
```

---

### 8.2 `GAUSS` (xsolve.f:22)

**Purpose:** General Gaussian elimination with partial pivoting

Used for local 4×4 systems in `MRCHUE` and `MRCHDU`.

---

## 9. Panel Method Interface

### 9.1 Key Panel Method Functions

| Function | Purpose |
|----------|---------|
| `QDCALC` | Compute dQ/dσ influence matrix DIJ |
| `QISET` | Set inviscid velocities QINV for current α |
| `UICALC` | Convert QINV to BL edge velocity UINV |
| `QVFUE` | Set QVIS from BL edge velocities UEDG |
| `GAMQV` | Update vorticity GAM from QVIS |
| `STFIND` | Locate stagnation point |
| `STMOVE` | Move stagnation point based on new CL |

### 9.2 Mass Defect Coupling

The viscous-inviscid coupling is through the mass defect:
```
σ = d(Ue × δ*)/ds = mass source/sink

Qvis(i) = Qinv(i) + Σⱼ DIJ(i,j) × σ(j)
```

---

## 10. Utility Functions

### 10.1 `XIFSET` (xbl.f:1196)
Sets forced transition x/c position to arc length coordinate.

### 10.2 `IBLSYS` (xbl.f:519)
Sets BL station → Newton system line pointers.

### 10.3 `DSLIM` (xbl.f:1564)
Limits δ* to maintain minimum Hk.

### 10.4 `BLPINI` (xbl.f:1578)
Initializes BL model parameters.

---

## 11. FlexFoil Implementation Mapping

### Module Structure

```
flexfoil/
├── src/
│   ├── lib.rs
│   ├── state.rs           # XfoilState equivalent
│   ├── parameters.rs      # BlParameters, constants
│   │
│   ├── panel/
│   │   ├── mod.rs
│   │   ├── geometry.rs    # Airfoil + wake geometry
│   │   ├── influence.rs   # AIJ, DIJ matrices
│   │   └── velocity.rs    # QINV, QVIS computation
│   │
│   ├── boundary_layer/
│   │   ├── mod.rs
│   │   ├── station.rs     # BlStation with derivatives
│   │   ├── closures.rs    # Hkin, Hs, Cf, Cd correlations
│   │   ├── equations.rs   # BLDIF, momentum, shape param
│   │   ├── transition.rs  # TRCHEK, DAMPL
│   │   └── inverse.rs     # Inverse mode (Hk prescribed)
│   │
│   ├── coupling/
│   │   ├── mod.rs
│   │   ├── newton.rs      # BLSOLV equivalent
│   │   ├── update.rs      # UPDATE with relaxation
│   │   └── stagnation.rs  # STFIND, STMOVE
│   │
│   └── solver.rs          # VISCAL main loop
```

### Implementation Priority

| Priority | Function | Critical for Stall? |
|----------|----------|---------------------|
| 1 | `HKIN`, `HSL`, `HST` | Baseline |
| 2 | `CFL`, `CFT`, `DIL` | Baseline |
| 3 | `DAMPL`, `TRCHEK2` | Transition |
| 4 | **`MRCHUE` inverse mode** | **YES - KEY** |
| 5 | `BLDIF`, `BLSYS` | Newton system |
| 6 | `BLSOLV` | System solve |
| 7 | `UPDATE` | Relaxation |
| 8 | `VISCAL` | Integration |

### Critical Implementation Notes

1. **Inverse Mode is Essential for Stall**
   - When `Hk > HLMAX (3.8 lam) or HTMAX (2.5 turb)`:
   - Prescribe `Hk_target`, solve for `Ue`
   - The reduced `Ue` feeds back to reduce lift

2. **All Derivatives Must Be Computed**
   - Every closure relation needs `∂/∂θ`, `∂/∂δ*`, `∂/∂Ue`
   - These form the Jacobian for Newton iteration

3. **Wake Handling**
   - Wake starts at TE with θ = θ_top + θ_bottom
   - Uses wake correlations (ITYP=3)
   - Dissipation length ratio = 0.9

4. **Compressibility**
   - All correlations include Mach corrections
   - Edge velocity is compressibility-corrected

---

## Appendix: Quick Reference

### Separation Thresholds
```
Laminar:   Hk > 3.8 → separation (inverse mode)
Turbulent: Hk > 2.5 → separation (inverse mode)
```

### Transition Criterion
```
N ≥ N_crit (typically 9) → transition to turbulent
```

### Key Variable Relationships
```
H = δ*/θ                    (shape parameter)
Hk = (H - 0.29M²)/(1+0.113M²)  (kinematic shape parameter)
H* = f(Hk, Rθ, M²)          (energy shape parameter)
Rθ = ρ·Ue·θ/μ               (momentum thickness Reynolds number)
mass = Ue · δ*              (mass defect for coupling)
```

---

*Document generated for FlexFoil development. Based on XFOIL 6.99 by Mark Drela.*
