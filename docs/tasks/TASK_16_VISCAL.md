# Task 16: Implement VISCAL Main Loop

## Objective
Implement the top-level viscous solver that coordinates all components.

## Prerequisites
- All rustfoil-bl and rustfoil-coupling tasks completed (01-15)

## Context
- VISCAL: xoper.f line 2886 - main viscous solution loop

## FORTRAN Reference

### VISCAL (xoper.f:2886-3161)
Main iteration loop:
1. Get inviscid solution for current alpha
2. Initialize BL from stagnation point
3. Build Newton system (SETBL)
4. Solve Newton system (BLSOLV)
5. Update BL variables (UPDATE)
6. Update edge velocities (UESET)
7. Check convergence, repeat if needed
8. Compute forces (CDCALC)

## Deliverables

### Create src/viscous/mod.rs in rustfoil-solver
```rust
//! Viscous solver module

pub mod viscal;
pub mod forces;
pub mod config;

pub use viscal::solve_viscous;
pub use forces::compute_forces;
pub use config::ViscousSolverConfig;
```

### Create src/viscous/config.rs
```rust
//! Viscous solver configuration

/// Configuration for viscous solver
#[derive(Debug, Clone)]
pub struct ViscousSolverConfig {
    /// Reynolds number
    pub reynolds: f64,
    /// Mach number
    pub mach: f64,
    /// Critical N factor for transition
    pub ncrit: f64,
    /// Maximum global iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Relaxation factor
    pub relaxation: f64,
}

impl Default for ViscousSolverConfig {
    fn default() -> Self {
        Self {
            reynolds: 1e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 50,
            tolerance: 1e-4,
            relaxation: 1.0,
        }
    }
}

impl ViscousSolverConfig {
    /// Create config for given Reynolds number
    pub fn with_reynolds(re: f64) -> Self {
        Self {
            reynolds: re,
            ..Default::default()
        }
    }
}
```

### Create src/viscous/viscal.rs
```rust
//! Main viscous solution loop
//!
//! XFOIL Reference: xoper.f VISCAL (line 2886)

use rustfoil_core::Body;
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::{
    dij::build_dij_matrix,
    newton::CoupledNewtonSystem,
    solve::solve_bl_system,
    march::march_fixed_ue,
    update::{update_stations, set_edge_velocities},
};
use super::config::ViscousSolverConfig;
use super::forces::{compute_forces, AeroForces};
use crate::inviscid::InviscidSolver;
use crate::SolverError;

/// Result of viscous solution
#[derive(Debug, Clone)]
pub struct ViscousResult {
    /// Angle of attack (degrees)
    pub alpha: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Moment coefficient
    pub cm: f64,
    /// Upper surface transition x/c
    pub x_tr_upper: f64,
    /// Lower surface transition x/c
    pub x_tr_lower: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Final residual
    pub residual: f64,
    /// Whether solution converged
    pub converged: bool,
}

/// Solve viscous flow for given angle of attack
///
/// # Arguments
/// * `body` - Airfoil body
/// * `alpha` - Angle of attack (degrees)
/// * `config` - Solver configuration
///
/// # Reference
/// XFOIL xoper.f VISCAL (line 2886)
pub fn solve_viscous(
    body: &Body,
    alpha: f64,
    config: &ViscousSolverConfig,
) -> Result<ViscousResult, SolverError> {
    let msq = config.mach * config.mach;
    
    // Step 1: Get inviscid solution
    let inviscid = InviscidSolver::new(body)?;
    let inv_result = inviscid.solve_alpha(alpha)?;
    
    // Extract edge velocities and arc lengths
    let n = body.panels().len();
    let x: Vec<f64> = body.panels().iter().map(|p| p.arc_length).collect();
    let ue_inviscid: Vec<f64> = inv_result.surface_velocities.clone();
    
    // Step 2: Build DIJ matrix
    let panel_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint.x).collect();
    let panel_y: Vec<f64> = body.panels().iter().map(|p| p.midpoint.y).collect();
    let dij = build_dij_matrix(&panel_x, &panel_y);
    
    // Step 3: Initialize BL from stagnation
    let march_config = rustfoil_coupling::march::MarchConfig {
        ncrit: config.ncrit,
        ..Default::default()
    };
    let mut march_result = march_fixed_ue(
        &x, &ue_inviscid, config.reynolds, msq, &march_config
    );
    let mut stations = march_result.stations;
    
    // Step 4: Global Newton iteration
    let mut converged = false;
    let mut residual = 1.0;
    let mut iter = 0;
    
    for i in 0..config.max_iterations {
        iter = i + 1;
        
        // Build Newton system
        let mut newton = CoupledNewtonSystem::new(&stations, &dij);
        newton.build(&stations, msq);
        
        // Solve
        let deltas = solve_bl_system(&newton.bl);
        
        // Update with limiting
        let update_config = rustfoil_coupling::update::UpdateConfig {
            relaxation: config.relaxation,
            ..Default::default()
        };
        let max_change = update_stations(&mut stations, &deltas, &update_config);
        
        // Update edge velocities
        set_edge_velocities(&mut stations, &ue_inviscid, &dij);
        
        // Check convergence
        residual = newton.bl.max_residual();
        if residual < config.tolerance && max_change < config.tolerance {
            converged = true;
            break;
        }
    }
    
    // Step 5: Compute forces
    let forces = compute_forces(&stations, &inv_result, config);
    
    Ok(ViscousResult {
        alpha,
        cl: forces.cl,
        cd: forces.cd,
        cm: forces.cm,
        x_tr_upper: march_result.x_transition.unwrap_or(1.0),
        x_tr_lower: march_result.x_transition.unwrap_or(1.0), // TODO: separate upper/lower
        iterations: iter,
        residual,
        converged,
    })
}

/// Solve viscous polar for multiple angles (parallel)
pub fn solve_viscous_polar_parallel(
    body: &Body,
    alphas: &[f64],
    config: &ViscousSolverConfig,
) -> Vec<Result<ViscousResult, SolverError>> {
    use rayon::prelude::*;
    
    alphas.par_iter()
        .map(|&alpha| solve_viscous(body, alpha, config))
        .collect()
}
```

### Create src/viscous/forces.rs
```rust
//! Aerodynamic force computation
//!
//! XFOIL Reference: xoper.f CDCALC

use rustfoil_bl::state::BlStation;
use super::config::ViscousSolverConfig;

/// Aerodynamic forces
#[derive(Debug, Clone)]
pub struct AeroForces {
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cd_pressure: f64,
    pub cd_friction: f64,
}

/// Compute aerodynamic forces from BL solution
///
/// # Reference
/// XFOIL xoper.f CDCALC
pub fn compute_forces(
    stations: &[BlStation],
    inviscid_result: &crate::inviscid::InviscidResult,
    config: &ViscousSolverConfig,
) -> AeroForces {
    // Friction drag from skin friction integration
    let cd_friction: f64 = stations.iter()
        .zip(stations.iter().skip(1))
        .map(|(s1, s2)| {
            let ds = s2.x - s1.x;
            let cf_avg = 0.5 * (s1.cf + s2.cf);
            cf_avg * ds
        })
        .sum();
    
    // Pressure drag from momentum deficit at wake trailing edge
    // (Squire-Young formula)
    let wake_station = stations.last().unwrap();
    let cd_pressure = 2.0 * wake_station.theta 
        * wake_station.u.powf((5.0 + wake_station.h) / 2.0);
    
    AeroForces {
        cl: inviscid_result.cl,  // Use inviscid CL (viscous correction small)
        cd: cd_friction + cd_pressure,
        cm: inviscid_result.cm,
        cd_pressure,
        cd_friction,
    }
}
```

## Verification
```bash
cargo build -p rustfoil-solver
cargo test -p rustfoil-solver viscous
```

## Next Task
After completion, proceed to TASK_17_CLI.md
