# Task 08: Implement Transition Prediction

## Objective
Implement amplification rate (DAMPL) and transition checking (TRCHEK2).

## Prerequisites
- Tasks 01-07 completed

## Context
- DAMPL: xblsys.f line 1981 (amplification rate dN/ds)
- TRCHEK2: xbl.f transition logic

## FORTRAN Reference

### DAMPL (xblsys.f:1981-2097)
Computes dN/dRθ for e^n method:
```fortran
      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
C
C---- Amplification rate  dN/dRt  (= dN/ds * Theta)
C
C    IMPLICIT REAL (A-H,M,O-Z)
C
      DGRT = 0.0
C
C---- Drela-Giles correlation for Tollmien-Schlichting growth rate
C
      HMI = 1.0/(HK-1.0)
      RLN = LOG(RT)
      RTNM = RT**(-0.8)
C
      DGR = 0.01*SQRT((2.4*HK-3.7 + 2.5*TANH(1.5*HK-4.65))
     &         /HK**2 )
      ... (continues for 100+ lines)
```

## Deliverables

### 1. Create src/transition.rs
```rust
//! Transition prediction using e^n method
//!
//! XFOIL Reference: xblsys.f DAMPL, DAMPL2

use crate::constants::*;

/// Amplification rate result
#[derive(Debug, Clone, Copy)]
pub struct AmplificationResult {
    /// dN/dRθ
    pub ax: f64,
    pub ax_hk: f64,
    pub ax_th: f64,
    pub ax_rt: f64,
}

/// Compute amplification rate dN/dRθ for Tollmien-Schlichting waves
///
/// Uses Drela-Giles correlation
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `th` - Momentum thickness θ
/// * `rt` - Reynolds number Rθ
pub fn amplification_rate(hk: f64, th: f64, rt: f64) -> AmplificationResult {
    // TODO: Port full DAMPL from xblsys.f:1981-2097
    // This is ~100 lines of complex correlation
    todo!("Implement DAMPL - see xblsys.f line 1981")
}

/// Check for laminar-turbulent transition
///
/// # Arguments
/// * `ampl` - Current amplification factor N
/// * `ncrit` - Critical N factor (typically 9 for free flight)
///
/// # Returns
/// * `None` if still laminar
/// * `Some(x_tr)` if transition occurred
pub fn check_transition(ampl: f64, ncrit: f64) -> Option<f64> {
    if ampl >= ncrit {
        // Transition occurred
        // In practice, interpolate to find exact x_tr
        Some(0.0) // Placeholder
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_transition_check() {
        assert!(check_transition(8.0, 9.0).is_none());
        assert!(check_transition(10.0, 9.0).is_some());
    }
}
```

## Next Task
After completion, proceed to TASK_09_STATE.md
