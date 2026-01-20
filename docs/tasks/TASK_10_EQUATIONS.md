# Task 10: Implement BLVAR and BLDIF

## Objective
Implement the core BL equation computations.

## Prerequisites
- Tasks 01-09 completed (all closures + state)

## Context
- BLVAR: xblsys.f line 784 - compute all secondary variables from primary
- BLDIF: xblsys.f line 1552 - compute residuals and Jacobian blocks

## FORTRAN Reference

### BLVAR (xblsys.f:784-1550)
~750 lines computing:
- Shape factors (H, Hk, Hs, Hc)
- Skin friction Cf
- Dissipation CD
- All partial derivatives

### BLDIF (xblsys.f:1552-1980)
~400 lines computing:
- Momentum integral residual
- Shape factor residual  
- Amplification/shear residual
- Jacobian blocks for Newton system

## Deliverables

### 1. Create src/equations.rs
```rust
//! Boundary layer integral equations
//!
//! XFOIL Reference: xblsys.f BLVAR (line 784), BLDIF (line 1552)

use crate::state::BlStation;
use crate::closures::*;

/// Compute all secondary variables from primary variables
///
/// Updates station.h, station.hk, station.hs, station.cf, etc.
/// and all partial derivatives in station.derivs
///
/// # Reference
/// XFOIL xblsys.f BLVAR (line 784)
pub fn blvar(station: &mut BlStation, msq: f64) {
    // Shape factor
    station.h = station.delta_star / station.theta;
    station.derivs.h_theta = -station.h / station.theta;
    station.derivs.h_delta_star = 1.0 / station.theta;
    
    // Kinematic shape factor
    let hkin_result = hkin(station.h, msq);
    station.hk = hkin_result.hk;
    station.derivs.hk_h = hkin_result.hk_h;
    station.derivs.hk_msq = hkin_result.hk_msq;
    
    // Reynolds number
    // Rθ = Ue * θ / ν, where ν is kinematic viscosity
    // In non-dimensional form: Rθ = Re_ref * Ue * θ
    // station.r_theta = ... (needs Reynolds number context)
    
    // Energy shape factor
    let hs_result = if station.is_laminar {
        hs_laminar(station.hk, station.r_theta, msq)
    } else {
        hs_turbulent(station.hk, station.r_theta, msq)
    };
    station.hs = hs_result.hs;
    station.derivs.hs_hk = hs_result.hs_hk;
    station.derivs.hs_rt = hs_result.hs_rt;
    station.derivs.hs_msq = hs_result.hs_msq;
    
    // Skin friction
    let cf_result = if station.is_laminar {
        cf_laminar(station.hk, station.r_theta, msq)
    } else {
        cf_turbulent(station.hk, station.r_theta, msq)
    };
    station.cf = cf_result.cf;
    station.derivs.cf_hk = cf_result.cf_hk;
    station.derivs.cf_rt = cf_result.cf_rt;
    station.derivs.cf_msq = cf_result.cf_msq;
    
    // Dissipation
    let di_result = if station.is_laminar {
        dissipation_laminar(station.hk, station.r_theta)
    } else {
        // Turbulent needs more inputs
        todo!("Turbulent dissipation")
    };
    // station.cd = ...
    
    // Density shape factor (for compressible)
    let hct_result = density_shape(station.hk, msq);
    station.hc = hct_result.hc;
}

/// Residuals of integral BL equations
#[derive(Debug, Clone)]
pub struct BlResiduals {
    /// Momentum integral residual
    pub res_mom: f64,
    /// Shape parameter residual  
    pub res_shape: f64,
    /// Third equation residual (amp or shear)
    pub res_third: f64,
}

/// Jacobian blocks for Newton system
#[derive(Debug, Clone)]
pub struct BlJacobian {
    /// Diagonal block (3x2 for upper/lower, 3x3 for wake)
    pub va: [[f64; 3]; 3],
    /// Off-diagonal block
    pub vb: [[f64; 3]; 3],
}

/// Compute BL equation residuals and Jacobian between two stations
///
/// # Arguments
/// * `s1` - Upstream station
/// * `s2` - Downstream station (current)
/// * `ds` - Arc length step
///
/// # Reference
/// XFOIL xblsys.f BLDIF (line 1552)
pub fn bldif(s1: &BlStation, s2: &BlStation, ds: f64) -> (BlResiduals, BlJacobian) {
    // TODO: Port full BLDIF from xblsys.f:1552-1980
    // This is the core of the BL solver
    todo!("Implement BLDIF - see xblsys.f line 1552")
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_blvar_laminar() {
        let mut station = BlStation::new();
        station.theta = 0.001;
        station.delta_star = 0.0026;
        station.is_laminar = true;
        
        blvar(&mut station, 0.0);
        
        // H should be δ*/θ
        assert!((station.h - 2.6).abs() < 0.01);
    }
}
```

## Next Task
After completion, proceed to TASK_11_COUPLING_DIJ.md
