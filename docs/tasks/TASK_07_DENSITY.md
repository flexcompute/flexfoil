# Task 07: Implement Density Shape Factor (HCT)

## Objective
Implement the density thickness shape factor closure.

## Prerequisites
- Tasks 01-06 completed

## Context
- HCT: xblsys.f line 2514

## FORTRAN Reference (xblsys.f:2514-2525)
```fortran
      SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
C
C---- density shape parameter     (from Whitfield)
      HC     = MSQ * (0.064/(HK-0.8) + 0.251)
      HC_HK  = MSQ * (-.064/(HK-0.8)**2     )
      HC_MSQ =        0.064/(HK-0.8) + 0.251
C
      RETURN
      END
```

## Deliverables

### 1. Create src/closures/density.rs
```rust
//! Density shape factor closure (Hc)
//!
//! XFOIL Reference: xblsys.f HCT (line 2514)

#[derive(Debug, Clone, Copy)]
pub struct HctResult {
    pub hc: f64,
    pub hc_hk: f64,
    pub hc_msq: f64,
}

/// Density thickness shape parameter (Whitfield correlation)
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `msq` - Mach²
pub fn density_shape(hk: f64, msq: f64) -> HctResult {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / (hk - 0.8).powi(2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    
    HctResult { hc, hc_hk, hc_msq }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_hct_incompressible() {
        // At M=0, Hc should be 0
        let result = density_shape(2.5, 0.0);
        assert_relative_eq!(result.hc, 0.0, epsilon = 1e-15);
    }
    
    #[test]
    fn test_hct_compressible() {
        let result = density_shape(2.5, 0.25);
        // Hc = 0.25 * (0.064/(2.5-0.8) + 0.251) = 0.25 * (0.0376 + 0.251)
        let expected = 0.25 * (0.064 / 1.7 + 0.251);
        assert_relative_eq!(result.hc, expected, epsilon = 1e-10);
    }
    
    #[test]
    fn test_hct_derivatives() {
        let hk = 2.5;
        let msq = 0.25;
        let eps = 1e-7;
        
        // Check hc_hk
        let r1 = density_shape(hk - eps, msq);
        let r2 = density_shape(hk + eps, msq);
        let numerical_hk = (r2.hc - r1.hc) / (2.0 * eps);
        let result = density_shape(hk, msq);
        assert_relative_eq!(result.hc_hk, numerical_hk, epsilon = 1e-5);
        
        // Check hc_msq
        let r1 = density_shape(hk, msq - eps);
        let r2 = density_shape(hk, msq + eps);
        let numerical_msq = (r2.hc - r1.hc) / (2.0 * eps);
        assert_relative_eq!(result.hc_msq, numerical_msq, epsilon = 1e-5);
    }
}
```

## Next Task
After completion, proceed to TASK_08_TRANSITION.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
