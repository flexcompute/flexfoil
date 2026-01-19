//! eN transition prediction method.
//!
//! The eN method predicts laminar-to-turbulent transition by tracking
//! the amplification of Tollmien-Schlichting waves. Transition occurs
//! when the amplification factor N reaches a critical value (typically 9).
//!
//! References:
//! - van Ingen, J.L. (1956) "A suggested semi-empirical method for the
//!   calculation of the boundary layer transition region"
//! - Drela, M. & Giles, M.B. (1987) "Viscous-Inviscid Analysis of
//!   Transonic and Low Reynolds Number Airfoils", AIAA Journal

/// Information about transition on a surface.
#[derive(Debug, Clone, Default)]
pub struct TransitionInfo {
    /// X-coordinate of transition
    pub x_tr: f64,
    /// Arc-length at transition
    pub s_tr: f64,
    /// Amplification factor N at transition
    pub n_factor: f64,
    /// Re_theta at transition
    pub re_theta_tr: f64,
    /// Whether forced transition was used
    pub forced: bool,
}

/// Compute the amplification rate dN/ds using XFOIL's envelope method.
///
/// This implements the DAMPL subroutine from XFOIL (Drela & Giles, 1987).
///
/// # Arguments
/// * `h` - Shape factor (kinematic shape parameter HK)
/// * `theta` - Momentum thickness
/// * `re_theta` - Momentum thickness Reynolds number = Ue*θ*Re
///
/// # Returns
/// dN/ds (amplification rate per unit arc length)
pub fn amplification_rate(h: f64, theta: f64, re_theta: f64) -> f64 {
    let hk = h.max(1.05);
    let theta = theta.max(1e-10);
    let rt = re_theta.max(1.0);
    
    // HMI = 1/(HK-1)
    let hmi = 1.0 / (hk - 1.0).max(0.1);
    
    // Critical log10(Re_theta) correlation
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let gr_crit = aa + 0.7 * (bb + 1.0);
    
    // Current log10(Re_theta)
    let gr = rt.log10();
    
    // Ramp parameter DGR = 0.08
    let dgr = 0.08;
    
    if gr < gr_crit - dgr {
        // Below critical: no amplification
        return 0.0;
    }
    
    // Smooth ramp turn-on
    let rnorm = ((gr - (gr_crit - dgr)) / (2.0 * dgr)).clamp(0.0, 1.0);
    let rfac = if rnorm >= 1.0 {
        1.0
    } else {
        3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3)
    };
    
    // DADR = dN/dRe_theta = 0.028*(HK-1) - 0.0345*exp(-ARG^2)
    let arg = 3.87 * hmi - 2.52;
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * (-arg.powi(2)).exp();
    
    // AF = m(H) correlation
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi.powi(2) + 3.0 * hmi.powi(3);
    
    // AX = AF * DADR / theta (amplification rate per unit x)
    // This is the key: divide by theta, not multiply!
    let ax = (af * dadr / theta) * rfac;
    
    ax.max(0.0)
}

/// Critical Re_theta for instability onset (XFOIL correlation).
///
/// Returns the critical momentum-thickness Reynolds number below which
/// all disturbances are damped.
pub fn critical_re_theta(h: f64) -> f64 {
    let hk = h.max(1.05);
    let hmi = 1.0 / (hk - 1.0).max(0.1);
    
    // XFOIL correlation: log10(Re_theta_crit) = AA + 0.7*(BB+1)
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let log10_re_crit = aa + 0.7 * (bb + 1.0);
    
    10.0_f64.powf(log10_re_crit)
}

/// Compute updated amplification factor.
///
/// # Arguments
/// * `n_prev` - Previous amplification factor
/// * `h` - Current shape factor
/// * `theta` - Momentum thickness
/// * `re_theta` - Current Re_theta = Ue*θ*Re
/// * `ds` - Step size in arc-length
///
/// # Returns
/// Updated amplification factor N
pub fn compute_amplification(n_prev: f64, h: f64, theta: f64, re_theta: f64, ds: f64) -> f64 {
    let dn_ds = amplification_rate(h, theta, re_theta);
    (n_prev + dn_ds * ds).max(0.0)
}

/// Check if transition should occur.
///
/// # Arguments
/// * `n` - Current amplification factor
/// * `n_crit` - Critical N-factor (typically 9)
///
/// # Returns
/// True if transition should occur
pub fn should_transition(n: f64, n_crit: f64) -> bool {
    n >= n_crit
}

/// Michel's criterion for transition (alternative to eN).
///
/// Transition occurs when Re_theta > Re_theta_tr where
/// Re_theta_tr = 1.174 * (1 + 22400/Re_x) * Re_x^0.46
///
/// This is simpler but less accurate than eN.
pub fn michel_transition(re_x: f64, re_theta: f64) -> bool {
    if re_x < 1e3 {
        return false;
    }
    
    let re_theta_tr = 1.174 * (1.0 + 22400.0 / re_x) * re_x.powf(0.46);
    re_theta > re_theta_tr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_critical_re_theta() {
        // For Blasius profile (H ≈ 2.59), Re_theta_crit ≈ 200-300
        let re_crit = critical_re_theta(2.59);
        assert!(re_crit > 100.0 && re_crit < 500.0, "Re_crit = {}", re_crit);
    }

    #[test]
    fn test_amplification_rate() {
        let theta = 0.0005; // Typical momentum thickness
        
        // Below critical, should be zero
        let rate = amplification_rate(2.59, theta, 100.0);
        assert!(rate.abs() < 0.01, "rate below crit = {}", rate);
        
        // Above critical, should be positive and significant
        let rate = amplification_rate(2.59, theta, 500.0);
        assert!(rate > 0.1, "rate above crit = {}", rate); // Should be O(1) or larger
    }

    #[test]
    fn test_compute_amplification() {
        let theta = 0.0005;
        let n = compute_amplification(0.0, 2.59, theta, 500.0, 0.01);
        assert!(n >= 0.0);
        assert!(n < 1.0, "Single step N should be small: {}", n); // One step shouldn't give N=9
    }

    #[test]
    fn test_transition_check() {
        assert!(!should_transition(5.0, 9.0));
        assert!(should_transition(9.0, 9.0));
        assert!(should_transition(10.0, 9.0));
    }

    #[test]
    fn test_michel() {
        // At Re_x = 5e5, Michel formula gives Re_theta_tr ≈ 513
        // Test with values clearly above and below the threshold
        assert!(michel_transition(5e5, 600.0));  // Above threshold
        assert!(!michel_transition(5e5, 200.0)); // Below threshold
    }
    
    #[test]
    fn test_xfoil_transition_location() {
        // For NACA 0012 at Re=1e6, XFOIL predicts x_tr ≈ 0.69
        // This test verifies the amplification formula integrates to
        // reasonable values with a realistic Ue distribution that includes
        // adverse pressure gradient (which accelerates amplification)
        let reynolds: f64 = 1e6;
        let n_crit: f64 = 9.0;
        
        let mut n: f64 = 0.0;
        let mut x_tr: f64 = 1.0;
        let ds: f64 = 0.005; // Finer steps
        let mut theta_sq_int: f64 = 0.0;
        
        // March from LE to TE with NACA 0012-like Ue distribution
        // Key: adverse pressure gradient (due_dx < 0) increases H and amplification
        for i in 1..200 {
            let x = i as f64 * ds;
            
            // More realistic NACA 0012 Ue distribution (symmetric airfoil at α=0)
            // Ue peaks around x/c=0.1-0.2, then decelerates towards TE
            let ue: f64 = if x < 0.05 {
                1.0 + 4.0 * x  // Strong acceleration near LE
            } else if x < 0.2 {
                1.2  // Velocity peak region
            } else {
                1.2 - 0.2 * (x - 0.2) / 0.8  // Gradual deceleration
            };
            
            // Thwaites integral for laminar BL
            theta_sq_int += ds * 0.45 / reynolds * ue.powi(5);
            let theta = (theta_sq_int / ue.powi(6)).sqrt();
            let re_theta = ue * theta * reynolds;
            
            // H varies with pressure gradient
            // Adverse gradient (x > 0.2) increases H significantly
            let h: f64 = if x < 0.2 {
                2.59  // Blasius
            } else {
                2.59 + 1.5 * (x - 0.2) / 0.8  // H rises to ~4 at TE (adverse gradient effect)
            };
            
            n = compute_amplification(n, h, theta, re_theta, ds);
            
            if n >= n_crit && x_tr > 0.99 {
                x_tr = x;
                break;
            }
        }
        
        // With realistic adverse pressure gradient modeling, transition should
        // occur somewhere in 0.3-0.9 range
        assert!(x_tr > 0.2 && x_tr < 0.95, 
            "x_tr = {:.2} (expected 0.3-0.9, N={:.1})", x_tr, n);
    }
}
