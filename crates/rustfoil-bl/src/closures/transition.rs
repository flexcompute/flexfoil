//! Transition prediction using e^n method
//!
//! This module implements the Drela-Giles correlation for Tollmien-Schlichting
//! wave amplification, used in the e^n transition prediction method.
//!
//! XFOIL Reference: xblsys.f DAMPL (line 1981)
//!
//! # Background
//!
//! The e^n method predicts laminar-turbulent transition by tracking the
//! amplification of Tollmien-Schlichting instability waves. Transition is
//! assumed to occur when the amplification factor N reaches a critical value
//! (typically Ncrit = 9 for free flight conditions).
//!
//! The amplification rate dN/dx depends on the local boundary layer state
//! (shape factor, momentum thickness, Reynolds number).

/// Smoothing half-width for critical Reynolds number transition
///
/// This value controls how smoothly the amplification "turns on" as
/// Rθ exceeds Rcrit. Larger values = smoother but less accurate transition.
const DGR: f64 = 0.08;

/// Result of amplification rate calculation
///
/// Contains the spatial amplification rate dN/dx and its partial derivatives
/// with respect to the boundary layer parameters. These derivatives are
/// essential for the Newton-Raphson solver's Jacobian matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmplificationResult {
    /// Spatial amplification rate dN/dx (1/length)
    pub ax: f64,
    /// Partial derivative ∂(dN/dx)/∂Hk
    pub ax_hk: f64,
    /// Partial derivative ∂(dN/dx)/∂θ
    pub ax_th: f64,
    /// Partial derivative ∂(dN/dx)/∂Rθ
    pub ax_rt: f64,
}

impl Default for AmplificationResult {
    fn default() -> Self {
        Self {
            ax: 0.0,
            ax_hk: 0.0,
            ax_th: 0.0,
            ax_rt: 0.0,
        }
    }
}

/// Compute amplification rate dN/dx for Tollmien-Schlichting waves
///
/// This is a direct port of XFOIL's DAMPL subroutine (xblsys.f:1981-2095),
/// implementing the Drela-Giles correlation for envelope e^n method.
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk (must be > 1.0)
/// * `th` - Momentum thickness θ (must be > 0)
/// * `rt` - Momentum thickness Reynolds number Rθ (must be > 0)
///
/// # Returns
/// `AmplificationResult` containing dN/dx and all partial derivatives
///
/// # Physics
///
/// The function returns zero amplification when Rθ is below the critical
/// Reynolds number (subcritical flow is stable). Above critical, amplification
/// increases according to the Drela-Giles correlation. A smooth cubic ramp
/// is used near the critical point to avoid numerical discontinuities.
///
/// # Reference
/// Drela, M., Giles, M., "Viscous/Inviscid Analysis of Transonic and
/// Low Reynolds Number Airfoils", AIAA Journal, Oct. 1987.
pub fn amplification_rate(hk: f64, th: f64, rt: f64) -> AmplificationResult {
    // Handle edge cases
    if hk <= 1.0 || th <= 0.0 || rt <= 0.0 {
        return AmplificationResult::default();
    }

    // Reciprocal of (Hk - 1) and its derivative
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // Log10(Critical Rθ) - H correlation for Falkner-Skan profiles
    // AA = 2.492 * (1/(Hk-1))^0.43
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    // BB = tanh(14/(Hk-1) - 9.24)
    let bb_arg = 14.0 * hmi - 9.24;
    let bb = bb_arg.tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    // Critical log10(Rθ)
    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    // Current log10(Rθ)
    let gr = rt.log10();
    let gr_rt = 1.0 / (std::f64::consts::LN_10 * rt);

    // Check if below critical (subcritical - no amplification)
    if gr < grcrit - DGR {
        return AmplificationResult::default();
    }

    // Compute smooth cubic ramp (RFAC) near critical
    // Ramp goes from 0 to 1 as log10(Rθ/Rcrit) goes from -DGR to +DGR
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);

    let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
        // Fully supercritical
        (1.0, 0.0, 0.0)
    } else {
        // In the ramp region: RFAC = 3*rnorm² - 2*rnorm³
        let rfac = 3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3);
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm.powi(2);
        (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
    };

    // Amplification envelope slope correlation for Falkner-Skan
    // ARG = 3.87/(Hk-1) - 2.52
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;

    // EX = exp(-ARG²)
    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);

    // DADR = 0.028*(Hk-1) - 0.0345*exp(-ARG²)
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;

    // New m(H) correlation (March 1991)
    // AF = -0.05 + 2.7*HMI - 5.5*HMI² + 3.0*HMI³
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi.powi(2) + 3.0 * hmi.powi(3);
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi.powi(2);
    let af_hk = af_hmi * hmi_hk;

    // Final amplification rate: AX = (AF * DADR / TH) * RFAC
    let ax_base = af * dadr / th;
    let ax = ax_base * rfac;

    // Derivatives
    // ∂AX/∂Hk has three terms: from AF, from DADR, and from RFAC
    let ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + ax_base * rfac_hk;

    // ∂AX/∂θ = -AX/θ (simple chain rule since AX ∝ 1/θ)
    let ax_th = -ax / th;

    // ∂AX/∂Rθ comes only from RFAC (the base rate doesn't depend on Rθ directly)
    let ax_rt = ax_base * rfac_rt;

    AmplificationResult {
        ax,
        ax_hk,
        ax_th,
        ax_rt,
    }
}

/// Check for laminar-turbulent transition
///
/// Compares the current amplification factor N against the critical value.
/// When N reaches Ncrit, transition to turbulence is triggered.
///
/// # Arguments
/// * `ampl` - Current amplification factor N (integrated dN/dx)
/// * `ncrit` - Critical N factor (typically 9 for free flight, 3-5 for noisy tunnels)
///
/// # Returns
/// * `None` if flow is still laminar (N < Ncrit)
/// * `Some(fraction)` if transition occurred, where fraction indicates how far
///   past the threshold (for interpolation purposes)
///
/// # Notes
/// In practice, the exact transition location is found by interpolating
/// between stations where N crosses Ncrit. This function provides the
/// basic threshold check.
pub fn check_transition(ampl: f64, ncrit: f64) -> Option<f64> {
    if ampl >= ncrit {
        // Transition occurred - return overshoot fraction for interpolation
        Some((ampl - ncrit) / ncrit.max(1e-10))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-7;
    const DERIV_TOL: f64 = 1e-5;

    /// Helper to compute numerical derivative
    fn numerical_deriv<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
        (f(x + h) - f(x - h)) / (2.0 * h)
    }

    // ========== Subcritical Tests ==========

    #[test]
    fn test_subcritical_zero() {
        // At low Reynolds number, flow is stable - no amplification
        // For Hk=2.5, critical Rθ ≈ 200-300
        let result = amplification_rate(2.5, 0.001, 50.0);

        assert_eq!(result.ax, 0.0, "AX should be zero below critical Rθ");
        assert_eq!(result.ax_hk, 0.0, "AX_HK should be zero below critical Rθ");
        assert_eq!(result.ax_th, 0.0, "AX_TH should be zero below critical Rθ");
        assert_eq!(result.ax_rt, 0.0, "AX_RT should be zero below critical Rθ");
    }

    #[test]
    fn test_subcritical_edge_cases() {
        // Invalid inputs should return zero
        assert_eq!(amplification_rate(0.5, 0.001, 1000.0).ax, 0.0, "Hk <= 1 should give zero");
        assert_eq!(amplification_rate(2.5, 0.0, 1000.0).ax, 0.0, "th = 0 should give zero");
        assert_eq!(amplification_rate(2.5, 0.001, 0.0).ax, 0.0, "rt = 0 should give zero");
        assert_eq!(amplification_rate(2.5, -0.001, 1000.0).ax, 0.0, "th < 0 should give zero");
    }

    // ========== Supercritical Tests ==========

    #[test]
    fn test_supercritical_positive() {
        // At high Reynolds number, amplification should be positive
        // For Hk=2.5, use Rθ well above critical
        let result = amplification_rate(2.5, 0.001, 5000.0);

        assert!(result.ax > 0.0, "AX should be positive above critical Rθ, got {}", result.ax);
        assert!(result.ax.is_finite(), "AX should be finite");
    }

    #[test]
    fn test_supercritical_typical_values() {
        // Test with typical boundary layer parameters
        // Hk=2.6 (attached laminar), θ=0.001, Rθ=1000
        let result = amplification_rate(2.6, 0.001, 1000.0);

        // Amplification rate should be reasonable (order of 1-100 per unit length)
        assert!(result.ax > 0.0, "Should have positive amplification");
        assert!(result.ax < 1e6, "Amplification rate should be bounded");

        // All values should be finite
        assert!(result.ax_hk.is_finite());
        assert!(result.ax_th.is_finite());
        assert!(result.ax_rt.is_finite());
    }

    // ========== Smooth Ramp Tests ==========

    #[test]
    fn test_smooth_ramp_near_critical() {
        // Test that amplification ramps smoothly near critical Rθ
        let hk = 2.5;
        let th = 0.001;

        // Find approximate critical Rθ by testing
        let mut rt_crit = 100.0;
        while amplification_rate(hk, th, rt_crit).ax == 0.0 && rt_crit < 10000.0 {
            rt_crit *= 1.1;
        }

        // Now test smoothness around this region
        let rt_below = rt_crit * 0.9;
        let rt_at = rt_crit;
        let rt_above = rt_crit * 1.1;

        let ax_below = amplification_rate(hk, th, rt_below).ax;
        let ax_at = amplification_rate(hk, th, rt_at).ax;
        let ax_above = amplification_rate(hk, th, rt_above).ax;

        // Should be monotonically increasing
        assert!(ax_at >= ax_below, "AX should increase with Rθ near critical");
        assert!(ax_above >= ax_at, "AX should increase with Rθ above critical");

        // Well above critical, check for smooth changes (ratio test only when both are non-trivial)
        let rt_well_above = rt_crit * 2.0;
        let rt_further = rt_crit * 2.2;
        let ax_well = amplification_rate(hk, th, rt_well_above).ax;
        let ax_further = amplification_rate(hk, th, rt_further).ax;

        if ax_well > 1e-10 && ax_further > 1e-10 {
            let ratio = ax_further / ax_well;
            assert!(ratio < 2.0 && ratio > 0.5, 
                "Amplification should change smoothly well above critical, ratio={}", ratio);
        }
    }

    #[test]
    fn test_no_discontinuity() {
        // Sweep through Rθ and check for continuity in the fully amplifying region
        // The ramp region (near critical) has inherently large relative changes as
        // the function transitions from zero - we only test well above critical
        let hk = 2.5;
        let th = 0.001;

        // First find where we're fully supercritical (ax_rt ≈ 0)
        let mut start_rt = 100.0;
        loop {
            let result = amplification_rate(hk, th, start_rt);
            // When ax_rt is near zero, we're past the ramp region
            if result.ax > 0.0 && result.ax_rt.abs() < 1e-12 {
                break;
            }
            start_rt *= 1.2;
            if start_rt > 100000.0 {
                panic!("Could not find supercritical region");
            }
        }

        // Now sweep and check for smooth changes in the fully supercritical region
        let mut prev_ax = amplification_rate(hk, th, start_rt).ax;
        for i in 1..50 {
            let rt = start_rt * (1.05_f64).powi(i);
            let ax = amplification_rate(hk, th, rt).ax;

            // Check for reasonable change
            if prev_ax > 0.0 && ax > 0.0 {
                let change = (ax - prev_ax).abs() / prev_ax;
                // In the supercritical region, changes should be gradual
                // (amplification rate is essentially constant w.r.t. Rθ)
                assert!(change < 0.1, "Discontinuity detected at Rθ={}: change={}", rt, change);
            }
            prev_ax = ax;
        }
    }

    // ========== Derivative Tests ==========

    #[test]
    fn test_ax_hk_derivative() {
        // Verify ∂AX/∂Hk numerically
        let hk = 2.5;
        let th = 0.001;
        let rt = 2000.0;

        let result = amplification_rate(hk, th, rt);

        // Only test if we're in the amplifying region
        if result.ax > 0.0 {
            let ax_hk_num = numerical_deriv(
                |h| amplification_rate(h, th, rt).ax,
                hk,
                EPS,
            );

            assert!(
                (result.ax_hk - ax_hk_num).abs() < DERIV_TOL * result.ax.abs().max(1.0),
                "ax_hk mismatch: analytical={}, numerical={}",
                result.ax_hk,
                ax_hk_num
            );
        }
    }

    #[test]
    fn test_ax_th_derivative() {
        // Verify ∂AX/∂θ numerically
        let hk = 2.5;
        let th = 0.001;
        let rt = 2000.0;

        let result = amplification_rate(hk, th, rt);

        // Only test if we're in the amplifying region
        if result.ax > 0.0 {
            let ax_th_num = numerical_deriv(
                |t| amplification_rate(hk, t, rt).ax,
                th,
                EPS * th, // Scale epsilon by th since it's small
            );

            assert!(
                (result.ax_th - ax_th_num).abs() < DERIV_TOL * result.ax.abs().max(1.0) / th,
                "ax_th mismatch: analytical={}, numerical={}",
                result.ax_th,
                ax_th_num
            );
        }
    }

    #[test]
    fn test_ax_rt_derivative() {
        // Verify ∂AX/∂Rθ numerically
        let hk = 2.5;
        let th = 0.001;
        let rt = 2000.0;

        let result = amplification_rate(hk, th, rt);

        // Only test if we're in the amplifying region AND in the ramp
        // (outside the ramp, ax_rt = 0)
        let ax_rt_num = numerical_deriv(
            |r| amplification_rate(hk, th, r).ax,
            rt,
            EPS * rt,
        );

        // Allow larger tolerance since this derivative can be small or zero
        let tol = DERIV_TOL * result.ax.abs().max(1.0) / rt.max(1.0);
        assert!(
            (result.ax_rt - ax_rt_num).abs() < tol.max(1e-10),
            "ax_rt mismatch: analytical={}, numerical={}",
            result.ax_rt,
            ax_rt_num
        );
    }

    #[test]
    fn test_derivatives_supercritical() {
        // Test derivatives well above critical where RFAC = 1
        let hk = 2.5;
        let th = 0.001;
        let rt = 10000.0; // Well above critical

        let result = amplification_rate(hk, th, rt);

        // ax_rt should be zero when fully supercritical (RFAC = 1, RFAC_RT = 0)
        assert!(
            result.ax_rt.abs() < 1e-10,
            "ax_rt should be ~0 when fully supercritical, got {}",
            result.ax_rt
        );

        // ax_th should equal -ax/th
        let expected_ax_th = -result.ax / th;
        assert!(
            (result.ax_th - expected_ax_th).abs() < 1e-10,
            "ax_th = -ax/th: {} vs {}",
            result.ax_th,
            expected_ax_th
        );
    }

    // ========== Transition Check Tests ==========

    #[test]
    fn test_transition_check_below() {
        // Below Ncrit - no transition
        assert!(check_transition(5.0, 9.0).is_none());
        assert!(check_transition(8.9, 9.0).is_none());
        assert!(check_transition(0.0, 9.0).is_none());
    }

    #[test]
    fn test_transition_check_at_or_above() {
        // At or above Ncrit - transition occurs
        assert!(check_transition(9.0, 9.0).is_some());
        assert!(check_transition(10.0, 9.0).is_some());
        assert!(check_transition(15.0, 9.0).is_some());
    }

    #[test]
    fn test_transition_check_fraction() {
        // Check that the returned fraction makes sense
        let result = check_transition(10.0, 9.0);
        assert!(result.is_some());
        let fraction = result.unwrap();
        // (10 - 9) / 9 ≈ 0.111
        assert!((fraction - 1.0/9.0).abs() < 1e-10);
    }

    // ========== Edge Case Tests ==========

    #[test]
    fn test_high_hk_near_separation() {
        // Near separation, Hk can be quite high (5-10)
        // The correlation should still work
        let result = amplification_rate(5.0, 0.001, 5000.0);

        assert!(result.ax.is_finite(), "AX should be finite for high Hk");
        assert!(result.ax_hk.is_finite(), "AX_HK should be finite for high Hk");
        assert!(result.ax_th.is_finite(), "AX_TH should be finite for high Hk");
        assert!(result.ax_rt.is_finite(), "AX_RT should be finite for high Hk");
    }

    #[test]
    fn test_hk_just_above_one() {
        // Hk very close to 1 (Blasius-like profile)
        let result = amplification_rate(1.01, 0.001, 10000.0);

        assert!(result.ax.is_finite(), "AX should be finite for Hk near 1");
        // With Hk close to 1, critical Rθ is very high, so we might still be subcritical
    }

    #[test]
    fn test_varying_hk() {
        // Test that amplification varies sensibly with Hk
        let th = 0.001;
        let rt = 5000.0;

        let ax_low_hk = amplification_rate(2.0, th, rt).ax;
        let ax_mid_hk = amplification_rate(2.5, th, rt).ax;
        let ax_high_hk = amplification_rate(3.5, th, rt).ax;

        // All should be finite
        assert!(ax_low_hk.is_finite());
        assert!(ax_mid_hk.is_finite());
        assert!(ax_high_hk.is_finite());

        // Generally, higher Hk = more unstable = higher amplification (if supercritical)
        // But the relationship is complex, so just check they're reasonable
        println!("AX at Hk=2.0: {}", ax_low_hk);
        println!("AX at Hk=2.5: {}", ax_mid_hk);
        println!("AX at Hk=3.5: {}", ax_high_hk);
    }
}
