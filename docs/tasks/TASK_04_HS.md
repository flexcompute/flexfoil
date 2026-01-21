# Task 04: Implement HSL and HST Closures

## Objective
Implement laminar (HSL) and turbulent (HST) energy shape factor closures.

## Prerequisites
- Tasks 01-03 completed

## Context
- HSL: xblsys.f line 2327 (laminar energy shape factor)
- HST: xblsys.f line 2388 (turbulent energy shape factor)

## FORTRAN Reference

### HSL (xblsys.f:2327-2352)
```fortran
      SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
C---- Laminar HS correlation
      IF(HK.LT.4.35) THEN
       TMP = HK - 4.35
       HS    = 0.0111*TMP*TMP/(HK+1.0)
     &       - 0.0278*TMP*TMP*TMP/(HK+1.0)  + 1.528
     &       - 0.0002*(TMP*HK)**2
       HS_HK = 0.0111*(2.0*TMP - TMP*TMP/(HK+1.0))/(HK+1.0)
     &       - 0.0278*(3.0*TMP*TMP - TMP*TMP*TMP/(HK+1.0))/(HK+1.0)
     &       - 0.0002*2.0*TMP*HK*(HK + TMP)
      ELSE
       HS    = 0.015*    (HK-4.35)**2/HK + 1.528
       HS_HK = 0.015*2.0*(HK-4.35)   /HK
     &       - 0.015*    (HK-4.35)**2/HK**2
      ENDIF
C
      HS_RT  = 0.0
      HS_MSQ = 0.0
C
      RETURN
      END
```

### HST (xblsys.f:2388-2481)
Much longer - turbulent correlation with Rθ dependence.

## Deliverables

### 1. Create src/closures/hs.rs
```rust
//! Energy shape factor closures (Hs)
//!
//! XFOIL Reference: xblsys.f HSL (line 2327), HST (line 2388)

/// Result of energy shape factor computation
#[derive(Debug, Clone, Copy)]
pub struct HsResult {
    pub hs: f64,
    pub hs_hk: f64,
    pub hs_rt: f64,
    pub hs_msq: f64,
}

/// Laminar energy shape factor
/// 
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Reynolds number based on θ (Rθ) - not used for laminar
/// * `msq` - Mach² - not used for laminar
pub fn hs_laminar(hk: f64, _rt: f64, _msq: f64) -> HsResult {
    let (hs, hs_hk) = if hk < 4.35 {
        let tmp = hk - 4.35;
        let tmp2 = tmp * tmp;
        let tmp3 = tmp2 * tmp;
        let hk1 = hk + 1.0;
        
        let hs = 0.0111 * tmp2 / hk1
               - 0.0278 * tmp3 / hk1
               + 1.528
               - 0.0002 * (tmp * hk).powi(2);
        
        let hs_hk = 0.0111 * (2.0 * tmp - tmp2 / hk1) / hk1
                  - 0.0278 * (3.0 * tmp2 - tmp3 / hk1) / hk1
                  - 0.0002 * 2.0 * tmp * hk * (hk + tmp);
        
        (hs, hs_hk)
    } else {
        let diff = hk - 4.35;
        let hs = 0.015 * diff.powi(2) / hk + 1.528;
        let hs_hk = 0.015 * 2.0 * diff / hk - 0.015 * diff.powi(2) / hk.powi(2);
        (hs, hs_hk)
    };
    
    HsResult {
        hs,
        hs_hk,
        hs_rt: 0.0,
        hs_msq: 0.0,
    }
}

/// Turbulent energy shape factor
///
/// # Arguments
/// * `hk` - Kinematic shape factor  
/// * `rt` - Reynolds number based on θ (Rθ)
/// * `msq` - Mach²
pub fn hs_turbulent(hk: f64, rt: f64, msq: f64) -> HsResult {
    // TODO: Port full HST from xblsys.f:2388-2481
    // This is a complex correlation with multiple branches
    todo!("Implement HST - see xblsys.f line 2388")
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_hsl_low_hk() {
        let result = hs_laminar(2.5, 1000.0, 0.0);
        // Verify against known value
        assert!(result.hs > 1.0 && result.hs < 2.0);
    }
    
    #[test]
    fn test_hsl_high_hk() {
        let result = hs_laminar(5.0, 1000.0, 0.0);
        assert!(result.hs > 1.5);
    }
    
    #[test]
    fn test_hsl_derivative() {
        let hk = 3.0;
        let eps = 1e-7;
        let r1 = hs_laminar(hk - eps, 0.0, 0.0);
        let r2 = hs_laminar(hk + eps, 0.0, 0.0);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);
        let result = hs_laminar(hk, 0.0, 0.0);
        assert_relative_eq!(result.hs_hk, numerical, epsilon = 1e-5);
    }
}
```

### 2. Add FORTRAN tests for HSL and HST
Extend test_closures.f with HSL and HST test cases.

### 3. Update mod.rs
```rust
pub mod hkin;
pub mod hs;

pub use hkin::{hkin, HkinResult};
pub use hs::{hs_laminar, hs_turbulent, HsResult};
```

## Verification
```bash
cargo test -p rustfoil-bl hs
```

## Next Task
After completion, proceed to TASK_05_CF.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
