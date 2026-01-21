//! Skin friction coefficient closures
//!
//! This module implements laminar and turbulent skin friction coefficient (Cf)
//! correlations from XFOIL's xblsys.f.

use crate::constants::CFFAC;

/// Result from skin friction closure calculation
///
/// Contains the skin friction coefficient and its partial derivatives
/// with respect to the input parameters.
#[derive(Debug, Clone, Copy)]
pub struct CfResult {
    /// Skin friction coefficient Cf
    pub cf: f64,
    /// ∂Cf/∂Hk (derivative w.r.t. kinematic shape factor)
    pub cf_hk: f64,
    /// ∂Cf/∂Rθ (derivative w.r.t. momentum thickness Reynolds number)
    pub cf_rt: f64,
    /// ∂Cf/∂M² (derivative w.r.t. Mach number squared)
    pub cf_msq: f64,
}

/// Laminar skin friction coefficient (Falkner-Skan correlation)
///
/// Computes Cf and its derivatives for laminar boundary layers.
/// Based on XFOIL's CFL subroutine (xblsys.f:2354-2371).
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk
/// * `rt` - Momentum thickness Reynolds number Rθ
/// * `_msq` - Mach number squared (unused, included for API consistency)
///
/// # Returns
/// `CfResult` containing Cf and partial derivatives
///
/// # Notes
/// - Uses two-branch correlation with threshold at Hk = 5.5
/// - For Hk < 5.5: attached laminar flow
/// - For Hk ≥ 5.5: separated laminar flow
pub fn cf_laminar(hk: f64, rt: f64, _msq: f64) -> CfResult {
    let (cf, cf_hk) = if hk < 5.5 {
        // Attached laminar branch
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        let cf = (0.0727 * tmp - 0.07) / rt;
        let cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
        (cf, cf_hk)
    } else {
        // Separated laminar branch
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        let cf = (0.015 * tmp.powi(2) - 0.07) / rt;
        let cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5).powi(2)) / rt;
        (cf, cf_hk)
    };

    let cf_rt = -cf / rt;
    let cf_msq = 0.0;

    CfResult {
        cf,
        cf_hk,
        cf_rt,
        cf_msq,
    }
}

/// Turbulent skin friction coefficient (Coles correlation)
///
/// Computes Cf and its derivatives for turbulent boundary layers.
/// Based on XFOIL's CFT subroutine (xblsys.f:2483-2510).
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk
/// * `rt` - Momentum thickness Reynolds number Rθ
/// * `msq` - Mach number squared M²
///
/// # Returns
/// `CfResult` containing Cf and partial derivatives
///
/// # Notes
/// - Includes compressibility correction via FC factor
/// - Uses GAM = 1.4 (ratio of specific heats, hardcoded as in XFOIL)
/// - GRT is clamped to minimum of 3.0
/// - ARG is clamped to minimum of -20.0
pub fn cf_turbulent(hk: f64, rt: f64, msq: f64) -> CfResult {
    const GAM: f64 = 1.4;
    const LN10: f64 = 2.302_585_092_994_046; // ln(10) = 2.3026...

    let gm1 = GAM - 1.0;

    // Compressibility correction factor
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();

    // Log Reynolds number (clamped)
    let grt_raw = (rt / fc).ln();
    let grt = grt_raw.max(3.0);

    // Exponent for Cf correlation
    let gex = -1.74 - 0.31 * hk;

    // Argument for exponential (clamped to avoid underflow)
    let arg = (-1.33 * hk).max(-20.0);

    // Tanh term for low-Re correction
    let thk = (4.0 - hk / 0.875).tanh();

    // Base Cf (before FC correction)
    let cfo = CFFAC * 0.3 * arg.exp() * (grt / LN10).powf(gex);

    // Final Cf with low-Re correction
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;

    // Derivative w.r.t. Hk
    // d(CFO)/d(HK) = CFO * (-1.33 + -0.31 * ln(grt/2.3026))
    // d(THK)/d(HK) = -(1 - THK^2) / 0.875
    let cf_hk = (-1.33 * cfo - 0.31 * (grt / LN10).ln() * cfo
        - 1.1e-4 * (1.0 - thk.powi(2)) / 0.875)
        / fc;

    // Derivative w.r.t. Rθ
    // Only contributes when grt > 3.0 (not clamped)
    let cf_rt = if grt_raw > 3.0 {
        gex * cfo / (fc * grt) / rt
    } else {
        0.0
    };

    // Derivative w.r.t. M²
    // Has two terms: one from CFO via grt, one from FC in denominator
    let cf_msq = if grt_raw > 3.0 {
        gex * cfo / (fc * grt) * (-0.25 * gm1 / fc.powi(2)) - 0.25 * gm1 * cf / fc.powi(2)
    } else {
        // When clamped, only the FC denominator term contributes
        -0.25 * gm1 * cf / fc.powi(2)
    };

    CfResult {
        cf,
        cf_hk,
        cf_rt,
        cf_msq,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;
    const DERIV_TOL: f64 = 1e-5;

    /// Helper to compute numerical derivative
    fn numerical_deriv<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
        (f(x + h) - f(x - h)) / (2.0 * h)
    }

    // ========== Laminar (CFL) Tests ==========

    #[test]
    fn test_cfl_attached_typical() {
        // Typical attached laminar: hk=2.6, rt=500
        let result = cf_laminar(2.6, 500.0, 0.0);

        // Cf should be positive for attached flow
        assert!(result.cf > 0.0, "Cf should be positive for attached flow");

        // Verify derivative w.r.t. Hk numerically
        let cf_hk_num =
            numerical_deriv(|h| cf_laminar(h, 500.0, 0.0).cf, 2.6, EPS);
        assert!(
            (result.cf_hk - cf_hk_num).abs() < DERIV_TOL,
            "cf_hk mismatch: analytical={}, numerical={}",
            result.cf_hk,
            cf_hk_num
        );

        // Verify derivative w.r.t. Rθ numerically
        let cf_rt_num =
            numerical_deriv(|r| cf_laminar(2.6, r, 0.0).cf, 500.0, EPS);
        assert!(
            (result.cf_rt - cf_rt_num).abs() < DERIV_TOL,
            "cf_rt mismatch: analytical={}, numerical={}",
            result.cf_rt,
            cf_rt_num
        );

        // cf_msq should be exactly 0 for laminar
        assert_eq!(result.cf_msq, 0.0);
    }

    #[test]
    fn test_cfl_branch_transition() {
        // At the transition point hk=5.5
        let result = cf_laminar(5.5, 1000.0, 0.0);

        // Both branches should give consistent values at hk=5.5
        // Test just above and below the threshold
        let below = cf_laminar(5.5 - EPS, 1000.0, 0.0);
        let above = cf_laminar(5.5 + EPS, 1000.0, 0.0);

        // Cf should be continuous across the transition
        assert!(
            (below.cf - above.cf).abs() < 1e-4,
            "Cf discontinuity at hk=5.5: below={}, above={}",
            below.cf,
            above.cf
        );

        // Verify numerical derivative at transition
        let cf_hk_num =
            numerical_deriv(|h| cf_laminar(h, 1000.0, 0.0).cf, 5.5, EPS);
        // Note: derivative may be discontinuous at branch point
        assert!(result.cf_hk.is_finite());
        println!(
            "At hk=5.5: cf_hk analytical={}, numerical={}",
            result.cf_hk, cf_hk_num
        );
    }

    #[test]
    fn test_cfl_separated() {
        // Separated laminar: hk=7.0, rt=800
        let result = cf_laminar(7.0, 800.0, 0.0);

        // In the separated branch, Cf can be negative (reversed flow)
        // Verify derivative w.r.t. Hk numerically
        let cf_hk_num =
            numerical_deriv(|h| cf_laminar(h, 800.0, 0.0).cf, 7.0, EPS);
        assert!(
            (result.cf_hk - cf_hk_num).abs() < DERIV_TOL,
            "cf_hk mismatch: analytical={}, numerical={}",
            result.cf_hk,
            cf_hk_num
        );

        // Verify derivative w.r.t. Rθ numerically
        let cf_rt_num =
            numerical_deriv(|r| cf_laminar(7.0, r, 0.0).cf, 800.0, EPS);
        assert!(
            (result.cf_rt - cf_rt_num).abs() < DERIV_TOL,
            "cf_rt mismatch: analytical={}, numerical={}",
            result.cf_rt,
            cf_rt_num
        );
    }

    // ========== Turbulent (CFT) Tests ==========

    #[test]
    fn test_cft_attached_typical() {
        // Typical attached turbulent: hk=1.4, rt=2000, msq=0.0
        let result = cf_turbulent(1.4, 2000.0, 0.0);

        // Cf should be positive
        assert!(result.cf > 0.0, "Cf should be positive for attached flow");

        // Verify derivative w.r.t. Hk numerically
        let cf_hk_num =
            numerical_deriv(|h| cf_turbulent(h, 2000.0, 0.0).cf, 1.4, EPS);
        assert!(
            (result.cf_hk - cf_hk_num).abs() < DERIV_TOL,
            "cf_hk mismatch: analytical={}, numerical={}",
            result.cf_hk,
            cf_hk_num
        );

        // Verify derivative w.r.t. Rθ numerically
        let cf_rt_num =
            numerical_deriv(|r| cf_turbulent(1.4, r, 0.0).cf, 2000.0, EPS);
        assert!(
            (result.cf_rt - cf_rt_num).abs() < DERIV_TOL,
            "cf_rt mismatch: analytical={}, numerical={}",
            result.cf_rt,
            cf_rt_num
        );

        // cf_msq should be 0 when msq=0 (no compressibility correction)
        let cf_msq_num =
            numerical_deriv(|m| cf_turbulent(1.4, 2000.0, m).cf, 0.0, EPS);
        assert!(
            (result.cf_msq - cf_msq_num).abs() < DERIV_TOL,
            "cf_msq mismatch: analytical={}, numerical={}",
            result.cf_msq,
            cf_msq_num
        );
    }

    #[test]
    fn test_cft_with_compressibility() {
        // With compressibility: hk=2.5, rt=5000, msq=0.1
        let result = cf_turbulent(2.5, 5000.0, 0.1);

        // Cf should be positive
        assert!(result.cf > 0.0);

        // Verify all derivatives numerically
        let cf_hk_num =
            numerical_deriv(|h| cf_turbulent(h, 5000.0, 0.1).cf, 2.5, EPS);
        assert!(
            (result.cf_hk - cf_hk_num).abs() < DERIV_TOL,
            "cf_hk mismatch: analytical={}, numerical={}",
            result.cf_hk,
            cf_hk_num
        );

        let cf_rt_num =
            numerical_deriv(|r| cf_turbulent(2.5, r, 0.1).cf, 5000.0, EPS);
        assert!(
            (result.cf_rt - cf_rt_num).abs() < DERIV_TOL,
            "cf_rt mismatch: analytical={}, numerical={}",
            result.cf_rt,
            cf_rt_num
        );

        let cf_msq_num =
            numerical_deriv(|m| cf_turbulent(2.5, 5000.0, m).cf, 0.1, EPS);
        assert!(
            (result.cf_msq - cf_msq_num).abs() < DERIV_TOL,
            "cf_msq mismatch: analytical={}, numerical={}",
            result.cf_msq,
            cf_msq_num
        );
    }

    #[test]
    fn test_cft_high_re_high_mach() {
        // High Re, high Mach: hk=3.0, rt=10000, msq=0.3
        let result = cf_turbulent(3.0, 10000.0, 0.3);

        // Cf should be positive
        assert!(result.cf > 0.0);

        // Verify all derivatives numerically
        let cf_hk_num =
            numerical_deriv(|h| cf_turbulent(h, 10000.0, 0.3).cf, 3.0, EPS);
        assert!(
            (result.cf_hk - cf_hk_num).abs() < DERIV_TOL,
            "cf_hk mismatch: analytical={}, numerical={}",
            result.cf_hk,
            cf_hk_num
        );

        let cf_rt_num =
            numerical_deriv(|r| cf_turbulent(3.0, r, 0.3).cf, 10000.0, EPS);
        assert!(
            (result.cf_rt - cf_rt_num).abs() < DERIV_TOL,
            "cf_rt mismatch: analytical={}, numerical={}",
            result.cf_rt,
            cf_rt_num
        );

        let cf_msq_num =
            numerical_deriv(|m| cf_turbulent(3.0, 10000.0, m).cf, 0.3, EPS);
        assert!(
            (result.cf_msq - cf_msq_num).abs() < DERIV_TOL,
            "cf_msq mismatch: analytical={}, numerical={}",
            result.cf_msq,
            cf_msq_num
        );
    }

    #[test]
    fn test_cft_grt_clamping() {
        // Test case where GRT might be clamped (low Rθ)
        // ln(rt/fc) < 3 means rt/fc < e^3 ≈ 20.09
        // With fc ≈ 1.0 (msq=0), rt < 20 triggers clamping
        let result = cf_turbulent(1.5, 15.0, 0.0);

        // Cf should still be finite
        assert!(result.cf.is_finite());
        assert!(result.cf_hk.is_finite());

        // When clamped, cf_rt should be 0
        assert_eq!(result.cf_rt, 0.0, "cf_rt should be 0 when GRT is clamped");
    }
}
