# Task 09: Implement BlStation State Structure

## Objective
Create the boundary layer state structure matching XBL.INC.

## Prerequisites
- Tasks 01-08 completed

## Context
- XBL.INC: `/Users/harry/flexfoil-boundary-layer/Xfoil/src/XBL.INC`
- Contains 73+ state variables and derivatives per BL station

## Deliverables

### 1. Create src/state.rs
```rust
//! Boundary layer state at a single station
//!
//! Matches XFOIL's XBL.INC common block structure

/// Complete boundary layer state at one station
#[derive(Debug, Clone)]
pub struct BlStation {
    // === Primary Variables (Newton unknowns) ===
    /// Arc length position
    pub x: f64,
    /// Edge velocity Ue
    pub u: f64,
    /// Momentum thickness θ
    pub theta: f64,
    /// Displacement thickness δ*
    pub delta_star: f64,
    /// Shear stress coefficient √Cτ (turbulent only)
    pub ctau: f64,
    /// Amplification factor N (laminar only)
    pub ampl: f64,
    
    // === Secondary Variables (computed by blvar) ===
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Kinematic shape factor Hk
    pub hk: f64,
    /// Energy shape factor Hs
    pub hs: f64,
    /// Density shape factor Hc
    pub hc: f64,
    /// Reynolds number Rθ = Ue*θ/ν
    pub r_theta: f64,
    /// Skin friction coefficient Cf
    pub cf: f64,
    /// Dissipation coefficient CD
    pub cd: f64,
    /// Mass defect Ue*δ*
    pub mass_defect: f64,
    
    // === Mode Flags ===
    pub is_laminar: bool,
    pub is_wake: bool,
    pub is_turbulent: bool,
    
    // === Partial Derivatives ===
    /// All partial derivatives for Jacobian construction
    pub derivs: BlDerivatives,
}

/// Partial derivatives of secondary variables w.r.t. primary variables
/// 
/// These are needed to construct the Newton system Jacobian.
/// Notation: `foo_bar` means ∂foo/∂bar
#[derive(Debug, Clone, Default)]
pub struct BlDerivatives {
    // H derivatives
    pub h_theta: f64,
    pub h_delta_star: f64,
    
    // Hk derivatives  
    pub hk_h: f64,
    pub hk_msq: f64,
    
    // Hs derivatives
    pub hs_hk: f64,
    pub hs_rt: f64,
    pub hs_msq: f64,
    
    // Cf derivatives
    pub cf_hk: f64,
    pub cf_rt: f64,
    pub cf_msq: f64,
    
    // CD derivatives
    pub cd_hk: f64,
    pub cd_rt: f64,
    
    // Edge velocity derivatives
    pub u_x: f64,      // dUe/dx
    
    // ... many more for full Jacobian
}

impl BlStation {
    /// Create a new station with default values
    pub fn new() -> Self {
        Self {
            x: 0.0,
            u: 1.0,
            theta: 0.001,
            delta_star: 0.002,
            ctau: 0.0,
            ampl: 0.0,
            h: 2.0,
            hk: 2.0,
            hs: 1.5,
            hc: 0.0,
            r_theta: 1000.0,
            cf: 0.003,
            cd: 0.001,
            mass_defect: 0.002,
            is_laminar: true,
            is_wake: false,
            is_turbulent: false,
            derivs: BlDerivatives::default(),
        }
    }
    
    /// Initialize for stagnation point
    pub fn stagnation(ue: f64, re: f64) -> Self {
        let mut station = Self::new();
        station.u = ue;
        // Hiemenz stagnation point solution
        station.theta = 0.29234 * (ue.abs() / re).sqrt();
        station.delta_star = station.theta * 2.216;
        station.h = 2.216;
        station.is_laminar = true;
        station
    }
}

impl Default for BlStation {
    fn default() -> Self {
        Self::new()
    }
}
```

## Next Task
After completion, proceed to TASK_10_EQUATIONS.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
