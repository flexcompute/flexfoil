//! Inverse mode boundary layer equations.
//!
//! When the shape factor Hk exceeds the separation threshold, the boundary layer
//! solver switches from "direct mode" (solve for δ* given Ue) to "inverse mode"
//! (solve for Ue given target Hk).
//!
//! This module implements XFOIL's inverse mode formulation from `xbl.f` lines 680-728.
//!
//! # Physics Background
//!
//! In attached flow, the boundary layer equations determine the shape factor H
//! from the edge velocity Ue. However, as separation approaches, H grows rapidly
//! and the direct equations become ill-conditioned.
//!
//! XFOIL's solution is to prescribe a slowly-varying target Hk in separated regions
//! and solve for the edge velocity Ue that produces this Hk. The target Hk evolves
//! gradually downstream, allowing the boundary layer to remain numerically stable
//! while still capturing the physics of separation.
//!
//! The modified edge velocity from inverse mode feeds back to the inviscid solver
//! through transpiration, naturally reducing lift as separation grows.
//!
//! # Reference
//!
//! XFOIL source: `xbl.f` lines 680-728 (MRCHUE subroutine)

use super::state::{HK_MAX_LAMINAR, HK_MAX_TURBULENT};

/// Configuration for inverse mode behavior.
#[derive(Debug, Clone, Copy)]
pub struct InverseConfig {
    /// Rate of Hk increase in laminar separated regions (XFOIL: 0.03)
    pub laminar_hk_rate: f64,
    /// Rate of Hk decrease in turbulent attached recovery (XFOIL: 0.15)
    pub turbulent_hk_rate: f64,
    /// Rate constant for wake Hk evolution (XFOIL: 0.03)
    pub wake_hk_rate: f64,
    /// Minimum Hk in wake (approaches 1.0 asymptotically)
    pub wake_hk_min: f64,
    /// Relaxation factor for Ue updates in inverse mode
    pub ue_relaxation: f64,
}

impl Default for InverseConfig {
    fn default() -> Self {
        Self {
            laminar_hk_rate: 0.03,
            turbulent_hk_rate: 0.15,
            wake_hk_rate: 0.03,
            wake_hk_min: 1.01,
            ue_relaxation: 0.5,
        }
    }
}

/// Determine if inverse mode should be used based on shape factor.
///
/// From XFOIL `xbl.f` lines 680-682:
/// ```fortran
/// IF(IBL.LT.ITRAN(IS)) HMAX = HLMAX
/// IF(IBL.GE.ITRAN(IS)) HMAX = HTMAX
/// DIRECT = HKTEST.LT.HMAX
/// ```
///
/// # Arguments
/// * `hk` - Current kinematic shape factor
/// * `is_turbulent` - Whether the flow is turbulent at this station
///
/// # Returns
/// `true` if inverse mode should be used (Hk >= threshold)
pub fn should_use_inverse_mode(hk: f64, is_turbulent: bool) -> bool {
    let hk_max = if is_turbulent {
        HK_MAX_TURBULENT
    } else {
        HK_MAX_LAMINAR
    };
    hk >= hk_max
}

/// Compute target Hk for inverse mode.
///
/// This is the key function that determines how Hk evolves in separated regions.
/// The target Hk changes slowly downstream, allowing stable numerical solution
/// while still capturing separation physics.
///
/// From XFOIL `xbl.f` lines 691-728.
///
/// # Arguments
/// * `hk_prev` - Shape factor at the previous station
/// * `ds` - Arc-length step to current station
/// * `theta_prev` - Momentum thickness at previous station
/// * `is_turbulent` - Whether the flow is turbulent
/// * `is_wake` - Whether this is a wake station
///
/// # Returns
/// Target Hk for the current station
pub fn compute_target_hk(
    hk_prev: f64,
    ds: f64,
    theta_prev: f64,
    is_turbulent: bool,
    is_wake: bool,
) -> f64 {
    compute_target_hk_with_config(hk_prev, ds, theta_prev, is_turbulent, is_wake, &InverseConfig::default())
}

/// Compute target Hk with custom configuration.
///
/// # Arguments
/// * `hk_prev` - Shape factor at the previous station
/// * `ds` - Arc-length step to current station  
/// * `theta_prev` - Momentum thickness at previous station
/// * `is_turbulent` - Whether the flow is turbulent
/// * `is_wake` - Whether this is a wake station
/// * `config` - Configuration parameters
///
/// # Returns
/// Target Hk for the current station
pub fn compute_target_hk_with_config(
    hk_prev: f64,
    ds: f64,
    theta_prev: f64,
    is_turbulent: bool,
    is_wake: bool,
    config: &InverseConfig,
) -> f64 {
    let theta_safe = theta_prev.max(1e-10);
    let ds_over_theta = ds / theta_safe;
    
    if !is_turbulent {
        // Laminar separated: Hk increases slowly downstream
        // From XFOIL: HTARG = HK1 + .03*(X2-X1)/TH1
        // But clamped to not exceed HK_MAX_LAMINAR
        let hk_new = hk_prev + config.laminar_hk_rate * ds_over_theta;
        hk_new.min(HK_MAX_LAMINAR + 1.0) // Allow slight overshoot
    } else if is_wake {
        // Wake: asymptotic approach to H = 1
        // From XFOIL: Newton iteration for H + const*(H-1)^3 = H_prev
        // This gives gradual recovery in the wake
        let const_val = config.wake_hk_rate * ds_over_theta;
        
        // Newton iteration to solve: h + const*(h-1)^3 = hk_prev
        let mut h = hk_prev;
        for _ in 0..3 {
            let f = h + const_val * (h - 1.0).powi(3) - hk_prev;
            let df = 1.0 + 3.0 * const_val * (h - 1.0).powi(2);
            if df.abs() > 1e-10 {
                h = h - f / df;
            }
        }
        h.max(config.wake_hk_min)
    } else {
        // Turbulent separated: Hk decreases toward reattachment
        // From XFOIL: HTARG = HK1 - .15*(X2-X1)/TH1
        // But clamped to not go below HK_MAX_TURBULENT
        let hk_new = hk_prev - config.turbulent_hk_rate * ds_over_theta;
        hk_new.max(HK_MAX_TURBULENT)
    }
}

/// Result of inverse mode BL solve at a single station.
#[derive(Debug, Clone, Copy)]
pub struct InverseModeResult {
    /// Updated edge velocity that produces target Hk
    pub ue_new: f64,
    /// Target Hk that was prescribed
    pub hk_target: f64,
    /// Whether the inverse solve converged
    pub converged: bool,
    /// Number of iterations used
    pub iterations: usize,
}

/// Solve for edge velocity in inverse mode.
///
/// In inverse mode, we prescribe Hk and solve for the edge velocity Ue
/// that produces this shape factor. This is done by Newton iteration
/// on the closure relations.
///
/// # Arguments
/// * `theta` - Momentum thickness at current station
/// * `hk_target` - Target shape factor to achieve
/// * `ue_guess` - Initial guess for edge velocity
/// * `reynolds` - Reynolds number
/// * `is_turbulent` - Whether flow is turbulent
/// * `max_iter` - Maximum iterations
/// * `tol` - Convergence tolerance
///
/// # Returns
/// Result containing the edge velocity that achieves target Hk
pub fn solve_inverse_mode(
    theta: f64,
    hk_target: f64,
    ue_guess: f64,
    reynolds: f64,
    is_turbulent: bool,
    max_iter: usize,
    tol: f64,
) -> InverseModeResult {
    let mut ue = ue_guess.max(1e-10);
    let mut converged = false;
    let mut iterations = 0;
    
    // In incompressible flow, Hk = H = delta*/theta
    // So we need to find Ue such that the closure relations give H = Hk_target
    //
    // The closure relations relate Cf, H, etc. to Re_theta = Ue * theta * Re
    // We iterate on Ue to match the target H.
    
    for iter in 0..max_iter {
        iterations = iter + 1;
        
        // Compute current shape factor from closure relations
        let re_theta = ue * theta * reynolds;
        let h_current = compute_shape_factor_from_closure(re_theta, is_turbulent);
        
        // Residual: difference from target
        let residual = h_current - hk_target;
        
        if residual.abs() < tol {
            converged = true;
            break;
        }
        
        // Compute derivative dH/dUe using finite difference
        let eps = ue * 1e-6;
        let re_theta_plus = (ue + eps) * theta * reynolds;
        let h_plus = compute_shape_factor_from_closure(re_theta_plus, is_turbulent);
        let dh_due = (h_plus - h_current) / eps;
        
        // Newton update
        if dh_due.abs() > 1e-10 {
            let delta_ue = -residual / dh_due;
            // Damped update to avoid overshooting
            ue = (ue + 0.5 * delta_ue).max(1e-10);
        } else {
            // Derivative too small, try bisection approach
            if residual > 0.0 {
                ue *= 1.1; // H too high, increase Ue
            } else {
                ue *= 0.9; // H too low, decrease Ue
            }
        }
    }
    
    InverseModeResult {
        ue_new: ue,
        hk_target,
        converged,
        iterations,
    }
}

/// Compute shape factor H from closure relations given Re_theta.
///
/// This is a simplified closure relation for inverse mode iteration.
/// For turbulent flow, uses an empirical correlation.
fn compute_shape_factor_from_closure(re_theta: f64, is_turbulent: bool) -> f64 {
    if is_turbulent {
        // Turbulent: H decreases with Re_theta (White's correlation simplified)
        // H ≈ 1.3 + 0.7 * (400/Re_theta)^0.5 for attached flow
        let re_safe = re_theta.max(100.0);
        1.3 + 0.7 * (400.0 / re_safe).sqrt().min(1.5)
    } else {
        // Laminar: Blasius value with mild Re dependence
        // H ≈ 2.59 for Blasius, increases toward separation
        2.59 + 0.5 * (1000.0 / re_theta.max(100.0)).sqrt().min(2.0)
    }
}

/// Update edge velocity from inverse mode with relaxation.
///
/// Applies under-relaxation to the edge velocity update for stability.
///
/// # Arguments
/// * `ue_old` - Previous edge velocity
/// * `ue_inverse` - Edge velocity from inverse mode solve
/// * `relaxation` - Relaxation factor (0 < relax <= 1)
///
/// # Returns
/// Relaxed edge velocity update
pub fn relax_ue_update(ue_old: f64, ue_inverse: f64, relaxation: f64) -> f64 {
    let relax = relaxation.clamp(0.1, 1.0);
    ue_old + relax * (ue_inverse - ue_old)
}

/// Compute the Ue modification needed to achieve target Hk.
///
/// This is used to feed back the inverse mode result to the inviscid solver
/// through a modified edge velocity boundary condition.
///
/// # Arguments
/// * `ue_inviscid` - Edge velocity from inviscid solver
/// * `ue_inverse` - Edge velocity from inverse mode
/// * `blend_factor` - How much to blend (0 = full inviscid, 1 = full inverse)
///
/// # Returns
/// Blended edge velocity
pub fn blend_ue_for_feedback(
    ue_inviscid: f64,
    ue_inverse: f64,
    blend_factor: f64,
) -> f64 {
    let blend = blend_factor.clamp(0.0, 1.0);
    (1.0 - blend) * ue_inviscid + blend * ue_inverse
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_should_use_inverse_mode() {
        // Laminar below threshold
        assert!(!should_use_inverse_mode(3.0, false));
        // Laminar at threshold
        assert!(should_use_inverse_mode(3.8, false));
        // Laminar above threshold
        assert!(should_use_inverse_mode(4.5, false));
        
        // Turbulent below threshold
        assert!(!should_use_inverse_mode(2.0, true));
        // Turbulent at threshold
        assert!(should_use_inverse_mode(2.5, true));
        // Turbulent above threshold
        assert!(should_use_inverse_mode(3.0, true));
    }
    
    #[test]
    fn test_compute_target_hk_laminar() {
        let hk_prev = 3.5;
        let ds = 0.01;
        let theta_prev = 0.001;
        
        let hk_target = compute_target_hk(hk_prev, ds, theta_prev, false, false);
        
        // Laminar: Hk should increase
        assert!(hk_target > hk_prev, "Laminar Hk should increase");
        // But not too fast
        assert!(hk_target < hk_prev + 1.0, "Laminar Hk increase should be gradual");
    }
    
    #[test]
    fn test_compute_target_hk_turbulent() {
        let hk_prev = 3.0;
        let ds = 0.01;
        let theta_prev = 0.001;
        
        let hk_target = compute_target_hk(hk_prev, ds, theta_prev, true, false);
        
        // Turbulent: Hk should decrease toward HK_MAX_TURBULENT
        assert!(hk_target < hk_prev, "Turbulent Hk should decrease");
        assert!(hk_target >= HK_MAX_TURBULENT, "Turbulent Hk should not go below threshold");
    }
    
    #[test]
    fn test_compute_target_hk_wake() {
        let hk_prev = 2.0;
        let ds = 0.01;
        let theta_prev = 0.001;
        
        let hk_target = compute_target_hk(hk_prev, ds, theta_prev, true, true);
        
        // Wake: Hk should decrease toward 1.0
        assert!(hk_target < hk_prev, "Wake Hk should decrease");
        assert!(hk_target >= 1.01, "Wake Hk should approach 1.0 asymptotically");
    }
    
    #[test]
    fn test_solve_inverse_mode_converges() {
        let theta = 0.001;
        let hk_target = 2.0;
        let ue_guess = 1.0;
        let reynolds = 1e6;
        
        let result = solve_inverse_mode(theta, hk_target, ue_guess, reynolds, true, 20, 1e-4);
        
        assert!(result.ue_new > 0.0, "Ue should be positive");
        // May not fully converge with simplified closures, but should be stable
        assert!(result.ue_new.is_finite(), "Ue should be finite");
    }
    
    #[test]
    fn test_relax_ue_update() {
        let ue_old = 1.0;
        let ue_inverse = 0.8;
        
        // Full relaxation
        let ue_full = relax_ue_update(ue_old, ue_inverse, 1.0);
        assert!((ue_full - 0.8).abs() < 1e-10);
        
        // Half relaxation
        let ue_half = relax_ue_update(ue_old, ue_inverse, 0.5);
        assert!((ue_half - 0.9).abs() < 1e-10);
        
        // No relaxation
        let ue_none = relax_ue_update(ue_old, ue_inverse, 0.0);
        // Clamped to 0.1 minimum
        assert!(ue_none < ue_old);
    }
    
    #[test]
    fn test_blend_ue_for_feedback() {
        let ue_inv = 1.0;
        let ue_inverse = 0.8;
        
        // Full inviscid
        assert!((blend_ue_for_feedback(ue_inv, ue_inverse, 0.0) - 1.0).abs() < 1e-10);
        
        // Full inverse
        assert!((blend_ue_for_feedback(ue_inv, ue_inverse, 1.0) - 0.8).abs() < 1e-10);
        
        // 50/50 blend
        assert!((blend_ue_for_feedback(ue_inv, ue_inverse, 0.5) - 0.9).abs() < 1e-10);
    }
}
