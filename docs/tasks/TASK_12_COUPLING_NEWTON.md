# Task 12: Implement BLSYS Newton System Builder

## Objective
Implement the Newton system construction for viscous-inviscid coupling.

## Prerequisites
- Task 11 (QDCALC) completed

## Context
- BLSYS: xblsys.f line 583 - builds the coupled Newton system
- SETBL: xbl.f line 21 - coordinates the full system assembly

## FORTRAN Reference

### BLSYS (xblsys.f:583-782)
Builds block-tridiagonal system for BL equations coupled with inviscid.

### System Structure
At each station i, we solve for updates to:
- δ* (displacement thickness)
- θ (momentum thickness)  
- N or Cτ (amplification or shear stress)

The coupled system has form:
```
VA_i * Δx_i + VB_i * Δx_{i-1} + VC_i * ΔUe_i = -R_i
```

Where VC represents coupling to inviscid through edge velocity.

## Deliverables

### src/newton.rs
```rust
//! Newton system construction for viscous-inviscid coupling
//!
//! XFOIL Reference: xblsys.f BLSYS (line 583), xbl.f SETBL (line 21)

use nalgebra::{DMatrix, DVector};
use rustfoil_bl::state::BlStation;
use rustfoil_bl::equations::{bldif, BlResiduals, BlJacobian};

/// Block-tridiagonal system for BL Newton iteration
pub struct BlNewtonSystem {
    /// Number of stations
    pub n: usize,
    /// Diagonal blocks VA[i] (3x3)
    pub va: Vec<[[f64; 3]; 3]>,
    /// Lower diagonal blocks VB[i] (3x3)  
    pub vb: Vec<[[f64; 3]; 3]>,
    /// Right-hand side residuals
    pub rhs: Vec<[f64; 3]>,
    /// Coupling to inviscid (dUe influence)
    pub vs: Vec<[f64; 3]>,
}

impl BlNewtonSystem {
    /// Create new system for n stations
    pub fn new(n: usize) -> Self {
        Self {
            n,
            va: vec![[[0.0; 3]; 3]; n],
            vb: vec![[[0.0; 3]; 3]; n],
            rhs: vec![[0.0; 3]; n],
            vs: vec![[0.0; 3]; n],
        }
    }
    
    /// Build the Newton system from BL stations
    ///
    /// # Reference
    /// XFOIL xblsys.f BLSYS (line 583)
    pub fn build(&mut self, stations: &[BlStation], msq: f64) {
        for i in 1..self.n {
            let ds = stations[i].x - stations[i-1].x;
            let (residuals, jacobian) = bldif(&stations[i-1], &stations[i], ds);
            
            // Store residuals
            self.rhs[i] = [residuals.res_mom, residuals.res_shape, residuals.res_third];
            
            // Store Jacobian blocks
            self.va[i] = jacobian.va;
            self.vb[i] = jacobian.vb;
        }
    }
    
    /// Get maximum residual (convergence check)
    pub fn max_residual(&self) -> f64 {
        self.rhs.iter()
            .flat_map(|r| r.iter())
            .map(|&v| v.abs())
            .fold(0.0, f64::max)
    }
}

/// Full coupled Newton system including inviscid interaction
pub struct CoupledNewtonSystem {
    /// BL system (block tridiagonal)
    pub bl: BlNewtonSystem,
    /// DIJ matrix for mass defect coupling
    pub dij: DMatrix<f64>,
    /// Global system matrix (assembled)
    pub matrix: DMatrix<f64>,
    /// Global RHS vector
    pub rhs: DVector<f64>,
}

impl CoupledNewtonSystem {
    /// Assemble the full coupled system
    ///
    /// Combines BL equations with inviscid coupling through DIJ
    pub fn assemble(&mut self) {
        // TODO: Port SETBL from xbl.f:21
        // This assembles the global Newton system
        todo!("Implement full system assembly")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_system_creation() {
        let system = BlNewtonSystem::new(50);
        assert_eq!(system.n, 50);
        assert_eq!(system.va.len(), 50);
    }
}
```

## Next Task
After completion, proceed to TASK_13_COUPLING_SOLVE.md
