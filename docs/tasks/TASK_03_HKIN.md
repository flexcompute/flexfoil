# Task 03: Implement HKIN Closure

## Objective
Implement the kinematic shape factor transformation (HKIN) with unit tests against FORTRAN.

## Prerequisites
- Task 01 (rustfoil-testkit) completed
- Task 02 (rustfoil-bl with constants) completed

## Context
- FORTRAN source: `/Users/harry/flexfoil-boundary-layer/Xfoil/src/xblsys.f` line 2276
- HKIN converts shape factor H to kinematic shape factor Hk

## FORTRAN Reference (xblsys.f:2276-2325)
```fortran
      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )
C
C---- Calculate kinematic shape parameter (usually called H**)
C
      MSQ_H = 0.0
C
C---- incompressible limit
      HK     =    (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =     1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK)/(1.0 + 0.113*MSQ)
C
      RETURN
      END
```

## Deliverables

### 1. Create src/closures/hkin.rs
```rust
//! Kinematic shape factor transformation
//!
//! XFOIL Reference: xblsys.f HKIN (line 2276)

/// Result of HKIN computation including partial derivatives
#[derive(Debug, Clone, Copy)]
pub struct HkinResult {
    /// Kinematic shape factor Hk
    pub hk: f64,
    /// Partial derivative ∂Hk/∂H
    pub hk_h: f64,
    /// Partial derivative ∂Hk/∂M²
    pub hk_msq: f64,
}

/// Calculate kinematic shape parameter Hk from shape factor H and Mach² 
///
/// The kinematic shape factor accounts for compressibility effects on the
/// boundary layer shape. This is XFOIL's compressibility correction.
///
/// # Arguments
/// * `h` - Shape factor H = δ*/θ
/// * `msq` - Mach number squared (M²)
///
/// # Returns
/// Kinematic shape factor and its partial derivatives
///
/// # Reference
/// XFOIL xblsys.f line 2276
pub fn hkin(h: f64, msq: f64) -> HkinResult {
    let denom = 1.0 + 0.113 * msq;
    let hk = (h - 0.29 * msq) / denom;
    let hk_h = 1.0 / denom;
    let hk_msq = (-0.29 - 0.113 * hk) / denom;
    
    HkinResult { hk, hk_h, hk_msq }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_hkin_incompressible() {
        // At M=0, Hk should equal H
        let result = hkin(2.5, 0.0);
        assert_relative_eq!(result.hk, 2.5, epsilon = 1e-10);
        assert_relative_eq!(result.hk_h, 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.hk_msq, -0.29, epsilon = 1e-10);
    }
    
    #[test]
    fn test_hkin_compressible() {
        // At M=0.5 (M²=0.25)
        let result = hkin(2.5, 0.25);
        // Expected: (2.5 - 0.29*0.25) / (1 + 0.113*0.25)
        //         = (2.5 - 0.0725) / 1.02825
        //         = 2.4275 / 1.02825 ≈ 2.361
        assert_relative_eq!(result.hk, 2.4275 / 1.02825, epsilon = 1e-10);
    }
    
    #[test]
    fn test_hkin_derivative_h() {
        // Numerical derivative check
        let h = 2.5;
        let msq = 0.25;
        let eps = 1e-7;
        
        let r1 = hkin(h - eps, msq);
        let r2 = hkin(h + eps, msq);
        let numerical_deriv = (r2.hk - r1.hk) / (2.0 * eps);
        
        let result = hkin(h, msq);
        assert_relative_eq!(result.hk_h, numerical_deriv, epsilon = 1e-6);
    }
    
    #[test]
    fn test_hkin_derivative_msq() {
        // Numerical derivative check
        let h = 2.5;
        let msq = 0.25;
        let eps = 1e-7;
        
        let r1 = hkin(h, msq - eps);
        let r2 = hkin(h, msq + eps);
        let numerical_deriv = (r2.hk - r1.hk) / (2.0 * eps);
        
        let result = hkin(h, msq);
        assert_relative_eq!(result.hk_msq, numerical_deriv, epsilon = 1e-6);
    }
}
```

### 2. Add FORTRAN test to test_closures.f
Extend `/Users/harry/flexfoil/crates/rustfoil-testkit/fortran/test_closures.f` with HKIN tests (if not already done in Task 01).

### 3. Create FORTRAN comparison test
Add to `rustfoil-bl/tests/fortran_comparison.rs`:
```rust
use rustfoil_bl::closures::hkin::hkin;
use rustfoil_testkit::load_reference;
use serde::Deserialize;

#[derive(Deserialize)]
struct HkinTest {
    h: f64,
    msq: f64,
    hk: f64,
    hk_h: f64,
    hk_msq: f64,
}

#[derive(Deserialize)]
struct HkinReference {
    hkin_tests: Vec<HkinTest>,
}

#[test]
fn test_hkin_matches_fortran() {
    let reference: HkinReference = load_reference("closures_reference.json")
        .expect("Failed to load reference data");
    
    for (i, t) in reference.hkin_tests.iter().enumerate() {
        let result = hkin(t.h, t.msq);
        
        assert!(
            (result.hk - t.hk).abs() < 1e-10,
            "Test {}: hk mismatch: {} vs {} for h={}, msq={}",
            i, result.hk, t.hk, t.h, t.msq
        );
        assert!(
            (result.hk_h - t.hk_h).abs() < 1e-10,
            "Test {}: hk_h mismatch", i
        );
        assert!(
            (result.hk_msq - t.hk_msq).abs() < 1e-10,
            "Test {}: hk_msq mismatch", i
        );
    }
}
```

### 4. Update src/closures/mod.rs
```rust
pub mod hkin;
pub use hkin::{hkin, HkinResult};
```

## Verification
```bash
# Generate reference data
cd /Users/harry/flexfoil/crates/rustfoil-testkit/fortran
make run

# Run Rust tests
cd /Users/harry/flexfoil
cargo test -p rustfoil-bl
cargo test -p rustfoil-bl test_hkin_matches_fortran
```

## Next Task
After completion, proceed to TASK_04_HS.md
