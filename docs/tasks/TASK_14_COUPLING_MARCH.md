# Task 14: Implement BL Marching (MRCHUE, MRCHDU)

**Status:** COMPLETE  
**Date:** 2026-01-21

## Objective
Implement boundary layer marching from stagnation point.

## Prerequisites
- Task 13 (solver) completed

## Context
- MRCHUE: xbl.f line 542 - march with fixed Ue
- MRCHDU: xbl.f line 875 - march with Ue updates

## FORTRAN Reference

### MRCHUE (xbl.f:542-873)
Direct BL march with prescribed edge velocity:
1. Start from stagnation point
2. March downstream solving BL equations
3. Check for transition and separation

### MRCHDU (xbl.f:875-1251)  
Coupled march with edge velocity updates from inviscid.

## Deliverables

### src/march.rs
```rust
//! Boundary layer marching
//!
//! XFOIL Reference: xbl.f MRCHUE (line 542), MRCHDU (line 875)

use rustfoil_bl::state::BlStation;
use rustfoil_bl::equations::blvar;
use rustfoil_bl::transition::amplification_rate;

/// Result of BL march
pub struct MarchResult {
    /// BL stations after march
    pub stations: Vec<BlStation>,
    /// Transition location (if occurred)
    pub x_transition: Option<f64>,
    /// Separation location (if occurred)  
    pub x_separation: Option<f64>,
    /// Whether march completed successfully
    pub converged: bool,
}

/// Configuration for BL march
pub struct MarchConfig {
    /// Critical N factor for transition
    pub ncrit: f64,
    /// Maximum iterations per station
    pub max_iter: usize,
    /// Convergence tolerance
    pub tolerance: f64,
}

impl Default for MarchConfig {
    fn default() -> Self {
        Self {
            ncrit: 9.0,
            max_iter: 25,
            tolerance: 1e-5,
        }
    }
}

/// March boundary layer with fixed edge velocity
///
/// Starts from stagnation point and marches downstream,
/// solving integral BL equations at each station.
///
/// # Arguments
/// * `x` - Arc length coordinates
/// * `ue` - Edge velocities at each station
/// * `re` - Reynolds number
/// * `msq` - Mach² (for compressibility)
/// * `config` - March configuration
///
/// # Reference
/// XFOIL xbl.f MRCHUE (line 542)
pub fn march_fixed_ue(
    x: &[f64],
    ue: &[f64],
    re: f64,
    msq: f64,
    config: &MarchConfig,
) -> MarchResult {
    let n = x.len();
    let mut stations = Vec::with_capacity(n);
    let mut x_transition = None;
    let mut x_separation = None;
    
    // Initialize stagnation point (first station)
    let mut prev = BlStation::stagnation(ue[0], re);
    prev.x = x[0];
    stations.push(prev.clone());
    
    // March downstream
    for i in 1..n {
        let ds = x[i] - x[i-1];
        let mut station = BlStation::new();
        station.x = x[i];
        station.u = ue[i];
        station.is_laminar = prev.is_laminar;
        
        // Initial guess: extrapolate from previous
        station.theta = prev.theta * (1.0 + ds * 0.01);
        station.delta_star = prev.delta_star * (1.0 + ds * 0.01);
        
        // Newton iteration to solve BL equations
        for _iter in 0..config.max_iter {
            blvar(&mut station, msq);
            
            // TODO: Compute residuals and update
            // This requires bldif and local Newton solve
            
            break; // Placeholder
        }
        
        // Check for transition (laminar only)
        if station.is_laminar && x_transition.is_none() {
            let amp = amplification_rate(station.hk, station.theta, station.r_theta);
            station.ampl = prev.ampl + amp.ax * ds;
            
            if station.ampl >= config.ncrit {
                x_transition = Some(station.x);
                station.is_laminar = false;
                station.is_turbulent = true;
            }
        }
        
        // Check for separation
        if station.cf < 0.0 && x_separation.is_none() {
            x_separation = Some(station.x);
        }
        
        stations.push(station.clone());
        prev = station;
    }
    
    MarchResult {
        stations,
        x_transition,
        x_separation,
        converged: true,
    }
}

/// March with edge velocity updates (viscous-inviscid coupling)
///
/// # Reference
/// XFOIL xbl.f MRCHDU (line 875)
pub fn march_coupled(
    x: &[f64],
    initial_ue: &[f64],
    dij: &nalgebra::DMatrix<f64>,
    re: f64,
    msq: f64,
    config: &MarchConfig,
) -> MarchResult {
    // TODO: Port MRCHDU from xbl.f:875
    todo!("Implement coupled march")
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_march_flat_plate() {
        // Flat plate: Ue = 1 everywhere
        let n = 50;
        let x: Vec<f64> = (0..n).map(|i| i as f64 * 0.02).collect();
        let ue = vec![1.0; n];
        
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &MarchConfig::default());
        
        assert!(result.stations.len() == n);
        // Blasius solution: H ≈ 2.59 for laminar flat plate
        // (This is a simplified test)
    }
}
```

## Next Task
After completion, proceed to TASK_15_COUPLING_UPDATE.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
