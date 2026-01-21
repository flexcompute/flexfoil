# Task 05: Implement CFL and CFT Closures

## Objective
Implement laminar (CFL) and turbulent (CFT) skin friction closures.

## Prerequisites
- Tasks 01-04 completed

## Context
- CFL: xblsys.f line 2354 (laminar skin friction)
- CFT: xblsys.f line 2483 (turbulent skin friction)

## FORTRAN Reference

### CFL (xblsys.f:2354-2386)
```fortran
      SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
C
C---- Laminar skin friction function  ( HS = Hs(Hk,Rt,Ms) )
C     El Hady correlation
C
      IF(HK.LT.5.5) THEN
       TMP = (5.5-HK)**3 / (HK+1.0)
       CF    = ( 0.0727*TMP                           - 0.07       )/RT
       CF_HK = ( 0.0727*TMP*(-3.0/(5.5-HK) - 1.0/(HK+1.0))         )/RT
      ELSE
       TMP = 1.0 - 1.0/(HK-4.5)
       CF    = ( 0.015*TMP      - 0.07  ) / RT
       CF_HK = ( 0.015/(HK-4.5)**2      ) / RT
      ENDIF
      CF_RT  = -CF/RT
      CF_MSQ = 0.0
C
      RETURN
      END
```

### CFT (xblsys.f:2483-2512)
Turbulent uses power law with compressibility correction.

## Deliverables

### 1. Create src/closures/cf.rs
```rust
//! Skin friction closures (Cf)
//!
//! XFOIL Reference: xblsys.f CFL (line 2354), CFT (line 2483)

use crate::constants::CFFAC;

/// Result of skin friction computation
#[derive(Debug, Clone, Copy)]
pub struct CfResult {
    pub cf: f64,
    pub cf_hk: f64,
    pub cf_rt: f64,
    pub cf_msq: f64,
}

/// Laminar skin friction coefficient
///
/// Uses El Hady correlation
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Reynolds number based on θ (Rθ)
/// * `msq` - Mach² (not used for laminar)
pub fn cf_laminar(hk: f64, rt: f64, _msq: f64) -> CfResult {
    let (cf, cf_hk) = if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        let cf = (0.0727 * tmp - 0.07) / rt;
        let cf_hk = (0.0727 * tmp * (-3.0 / (5.5 - hk) - 1.0 / (hk + 1.0))) / rt;
        (cf, cf_hk)
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        let cf = (0.015 * tmp - 0.07) / rt;
        let cf_hk = (0.015 / (hk - 4.5).powi(2)) / rt;
        (cf, cf_hk)
    };
    
    CfResult {
        cf,
        cf_hk,
        cf_rt: -cf / rt,
        cf_msq: 0.0,
    }
}

/// Turbulent skin friction coefficient  
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Reynolds number based on θ (Rθ)
/// * `msq` - Mach²
pub fn cf_turbulent(hk: f64, rt: f64, msq: f64) -> CfResult {
    // TODO: Port full CFT from xblsys.f:2483-2512
    todo!("Implement CFT - see xblsys.f line 2483")
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_cfl_typical() {
        // Typical laminar values: Hk≈2.6, Rθ≈500
        let result = cf_laminar(2.6, 500.0, 0.0);
        // Cf should be small positive
        assert!(result.cf > 0.0 && result.cf < 0.01);
    }
    
    #[test]  
    fn test_cfl_derivative_rt() {
        let hk = 2.6;
        let rt = 500.0;
        let eps = 1.0;
        let r1 = cf_laminar(hk, rt - eps, 0.0);
        let r2 = cf_laminar(hk, rt + eps, 0.0);
        let numerical = (r2.cf - r1.cf) / (2.0 * eps);
        let result = cf_laminar(hk, rt, 0.0);
        assert_relative_eq!(result.cf_rt, numerical, epsilon = 1e-6);
    }
}
```

### 2. Update mod.rs
Add cf module exports.

## Next Task
After completion, proceed to TASK_06_DISSIPATION.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
