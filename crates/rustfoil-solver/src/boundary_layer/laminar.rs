//! Thwaites' method for laminar boundary layers.
//!
//! Thwaites' integral method solves the momentum integral equation
//! for laminar boundary layers using empirical correlations.
//!
//! The key equation is:
//! θ²(s) = (0.45/Ue^6) * ∫₀ˢ Ue^5 ds
//!
//! Reference: Thwaites, B. (1949) "Approximate Calculation of the Laminar 
//! Boundary Layer", Aeronautical Quarterly, Vol 1.

/// Solve Thwaites' integral for momentum thickness squared.
///
/// # Arguments
/// * `s` - Arc-length from stagnation point
/// * `ue` - Edge velocity at current station
/// * `reynolds` - Chord Reynolds number
///
/// # Returns
/// θ² × 6 (for numerical stability, actual θ² = result/6)
pub fn thwaites_solve(s: f64, ue: f64, reynolds: f64) -> f64 {
    // Thwaites' equation: θ² = (0.45 * ν / Ue^6) * ∫ Ue^5 ds
    // For incompressible flow with ν = 1/Re (normalized):
    // θ² ≈ 0.45 * s / (Re * Ue)  for roughly constant Ue
    
    // This is a simplified form; full integration would use:
    // θ² = (0.45/Re) * (1/Ue^6) * ∫₀ˢ Ue^5 ds'
    
    // For a flat plate (Ue = const), this gives:
    // θ² = 0.45 * s / (Re * Ue)
    // θ = 0.671 * sqrt(s / (Re * Ue))  -- close to Blasius θ/x = 0.664/sqrt(Rex)
    
    let nu = 1.0 / reynolds;
    
    // Prevent division by zero
    let ue_safe = ue.max(1e-10);
    
    // Thwaites constant = 0.45
    let theta_sq = 0.45 * nu * s / ue_safe;
    
    // Return θ² × 6 for the closure relations
    theta_sq * 6.0
}

/// Thwaites' pressure gradient parameter λ.
///
/// λ = (θ²/ν) * (dUe/ds) = Re * θ² * (dUe/ds)
///
/// Typical values:
/// - λ = 0: flat plate (Blasius)
/// - λ > 0: favorable pressure gradient (accelerating)
/// - λ < 0: adverse pressure gradient (decelerating)
/// - λ ≈ -0.09: separation
pub fn thwaites_lambda(theta: f64, due_ds: f64, reynolds: f64, ue: f64) -> f64 {
    let nu = 1.0 / reynolds;
    theta.powi(2) * due_ds / nu / ue.max(1e-10)
}

/// Compute separation criterion.
///
/// Returns true if the laminar BL has separated (λ < -0.09).
pub fn is_separated(lambda: f64) -> bool {
    lambda < -0.09
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flat_plate() {
        // Blasius flat plate: θ/x = 0.664/sqrt(Rex)
        // At x = 0.5, Re = 1e6: Rex = 5e5
        // θ = 0.5 * 0.664 / sqrt(5e5) = 0.00047
        
        let s = 0.5;
        let ue = 1.0;
        let re = 1e6;
        
        let theta_sq_6 = thwaites_solve(s, ue, re);
        let theta = (theta_sq_6 / 6.0).sqrt();
        
        // Thwaites gives θ ≈ 0.671 * sqrt(x / Rex) vs Blasius 0.664
        let theta_blasius = 0.664 * (s / re).sqrt();
        
        // Should be within ~5% of Blasius
        assert!((theta - theta_blasius).abs() / theta_blasius < 0.1);
    }

    #[test]
    fn test_separation() {
        assert!(is_separated(-0.1));
        assert!(!is_separated(0.0));
        assert!(!is_separated(-0.05));
    }
}
