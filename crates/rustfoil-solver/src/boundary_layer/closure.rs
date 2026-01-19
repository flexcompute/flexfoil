//! Closure relations for boundary layer equations.
//!
//! These correlations relate the integral boundary layer parameters
//! (θ, δ*, H) to each other and to skin friction.
//!
//! References:
//! - Thwaites (1949) for laminar correlations
//! - Ludwieg & Tillmann (1949) for turbulent skin friction
//! - Head (1958) for entrainment shape factor


// ============================================================================
// Laminar Correlations (Thwaites)
// ============================================================================

/// Shape factor H from Thwaites' pressure gradient parameter λ.
///
/// H = f(λ) where λ = (θ²/ν)(dUe/ds)
///
/// Correlation from Thwaites (1949):
/// - λ = 0 (flat plate): H = 2.59
/// - λ = 0.09 (favorable): H ≈ 2.0
/// - λ = -0.09 (separation): H ≈ 3.5
pub fn thwaites_h(lambda: f64) -> f64 {
    // Thwaites' correlation
    if lambda >= 0.0 {
        // Favorable pressure gradient: H decreases as λ increases
        2.61 - 3.75 * lambda + 5.24 * lambda.powi(2)
    } else {
        // Adverse pressure gradient: H increases as λ becomes more negative
        // Linear fit: H ≈ 2.61 - 10.3λ for λ ∈ [-0.09, 0]
        // At λ = -0.09 (near separation): H ≈ 3.54
        let l = (-lambda).min(0.09);
        2.61 + 10.3 * l
    }
}

/// Shear correlation L(λ) from Thwaites.
///
/// Related to skin friction: Cf = 2L / Re_theta
fn thwaites_l(lambda: f64) -> f64 {
    if lambda >= 0.0 {
        0.22 + 1.57 * lambda - 1.8 * lambda.powi(2)
    } else {
        let l = (-lambda).min(0.1);
        0.22 + 1.402 * l + 0.018 * l / (0.107 + l)
    }
}

/// Skin friction from Thwaites method.
///
/// Cf = 2 * L(λ) / Re_theta
pub fn thwaites_cf(lambda: f64, theta: f64, reynolds: f64, ue: f64) -> f64 {
    let re_theta = (ue * theta * reynolds).max(1.0);
    let l = thwaites_l(lambda);
    2.0 * l / re_theta
}

/// General laminar skin friction (Blasius).
///
/// Cf = 0.664 / sqrt(Re_x)  -- flat plate
pub fn blasius_cf(re_x: f64) -> f64 {
    if re_x < 1.0 {
        return 0.0;
    }
    0.664 / re_x.sqrt()
}

// ============================================================================
// Turbulent Correlations
// ============================================================================

/// Ludwieg-Tillmann skin friction formula.
///
/// Cf = 0.246 * 10^(-0.678*H) * Re_theta^(-0.268)
///
/// Valid for equilibrium turbulent boundary layers.
pub fn ludwieg_tillmann_cf(re_theta: f64, h: f64) -> f64 {
    if re_theta < 100.0 {
        return 0.005; // Fallback for very low Re
    }
    
    let h_clamped = h.clamp(1.0, 4.0);
    0.246 * 10.0_f64.powf(-0.678 * h_clamped) * re_theta.powf(-0.268)
}

/// White-Christoph skin friction formula (alternative).
///
/// Cf = 0.3 * exp(-1.33*H) / (log10(Re_theta))^(1.74 + 0.31*H)
pub fn white_christoph_cf(re_theta: f64, h: f64) -> f64 {
    if re_theta < 100.0 {
        return 0.005;
    }
    
    let h_clamped = h.clamp(1.0, 4.0);
    let log_re = re_theta.log10();
    
    0.3 * (-1.33 * h_clamped).exp() / log_re.powf(1.74 + 0.31 * h_clamped)
}

/// Turbulent flat plate skin friction (Schlichting).
///
/// Cf = 0.0583 / Re_x^(1/5)  (for 5e5 < Re_x < 1e7)
pub fn schlichting_cf(re_x: f64) -> f64 {
    if re_x < 1e4 {
        return 0.005;
    }
    0.0583 / re_x.powf(0.2)
}

// ============================================================================
// Head's Entrainment Method
// ============================================================================

/// Head's shape factor H1 = (δ - δ*) / θ
///
/// H1 is related to H by a correlation.
/// For equilibrium flow: H1 ≈ 3.3 + 0.8234 * (H - 1.1)^(-1.287)
pub fn head_h1(h: f64) -> f64 {
    let h_safe = h.max(1.1);
    
    if h_safe < 1.6 {
        3.3 + 0.8234 * (h_safe - 1.1 + 0.01).powf(-1.287)
    } else {
        3.3 + 1.5501 * (h_safe - 0.6778).powf(-3.064)
    }
}

/// Inverse: get H from H1.
pub fn head_h_from_h1(h1: f64) -> f64 {
    // Numerical inversion using Newton's method
    let mut h = 1.4; // Initial guess
    
    for _ in 0..10 {
        let h1_calc = head_h1(h);
        let error = h1_calc - h1;
        
        if error.abs() < 1e-6 {
            break;
        }
        
        // Derivative approximation
        let dh1_dh = (head_h1(h + 0.001) - h1_calc) / 0.001;
        
        if dh1_dh.abs() > 1e-10 {
            h -= error / dh1_dh;
            h = h.clamp(1.0, 10.0);
        }
    }
    
    h
}

/// Head's entrainment coefficient Ce.
///
/// Ce = F(H1) where F is empirical.
/// dE/ds = Ce * Ue  where E = H1 * θ
pub fn head_entrainment(h1: f64) -> f64 {
    // Head's correlation
    let h1_safe = h1.max(3.0);
    0.0306 * (h1_safe - 3.0).powf(-0.6169)
}

// ============================================================================
// Shape Factor Correlations (XFOIL-style)
// ============================================================================

/// General shape factor correlations used in XFOIL.
pub fn shape_factor_correlations(hk: f64, re_theta: f64, is_turbulent: bool) -> ShapeFactors {
    if is_turbulent {
        turbulent_shape_factors(hk, re_theta)
    } else {
        laminar_shape_factors(hk)
    }
}

/// Shape factors for laminar flow.
fn laminar_shape_factors(hk: f64) -> ShapeFactors {
    let hk = hk.clamp(1.0, 10.0);
    
    // Energy thickness shape factor Hs = θ*/θ
    let hs = if hk < 4.35 {
        0.0111 * (hk - 1.0).powi(2) / (hk - 1.0 + 0.0278) + 1.528 - 0.0002 * (hk.powi(2) - 1.0)
    } else {
        0.015 * (hk - 4.35).powi(2) + 1.528
    };
    
    // Density thickness shape factor
    let hs_star = hs;
    
    ShapeFactors {
        h: hk,
        hs,
        hs_star,
    }
}

/// Shape factors for turbulent flow.
fn turbulent_shape_factors(hk: f64, re_theta: f64) -> ShapeFactors {
    let hk = hk.clamp(1.0, 10.0);
    let re_theta = re_theta.max(100.0);
    
    // Hs correlation from XFOIL
    let hs = if hk < 4.0 {
        let hk0 = if hk < 1.1 { 1.1 } else { hk };
        1.505 + 4.0 / (hk0 - 1.0) + 0.165 - 1.6 / hk0
    } else {
        1.505 + 4.0 / (hk - 1.0) + (hk - 4.0).powi(2) / (6.0 * hk)
    };
    
    // Re_theta correction
    let log_re = re_theta.ln();
    let hs_corr = hs + 0.028 * (log_re - 7.4).max(0.0);
    
    ShapeFactors {
        h: hk,
        hs: hs_corr,
        hs_star: hs_corr,
    }
}

/// Collection of shape factors.
#[derive(Debug, Clone, Copy)]
pub struct ShapeFactors {
    /// Displacement thickness shape factor H = δ*/θ
    pub h: f64,
    /// Energy thickness shape factor Hs = θ*/θ
    pub hs: f64,
    /// Density thickness shape factor
    pub hs_star: f64,
}

/// General skin friction calculation.
pub fn skin_friction(re_theta: f64, h: f64, is_turbulent: bool) -> f64 {
    if is_turbulent {
        ludwieg_tillmann_cf(re_theta, h)
    } else {
        // Laminar: approximate using Blasius relationship
        if re_theta < 1.0 {
            return 0.0;
        }
        // Cf * sqrt(Re_theta) ≈ 0.664 for Blasius
        0.664 / re_theta.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thwaites_h() {
        // Flat plate: λ = 0, H ≈ 2.59
        let h = thwaites_h(0.0);
        assert!((h - 2.59).abs() < 0.1);
        
        // Favorable: lower H
        let h_fav = thwaites_h(0.05);
        assert!(h_fav < h);
        
        // Adverse: higher H
        let h_adv = thwaites_h(-0.05);
        assert!(h_adv > h);
    }

    #[test]
    fn test_ludwieg_tillmann() {
        // At Re_theta = 1000, H = 1.4: Cf ≈ 0.003
        let cf = ludwieg_tillmann_cf(1000.0, 1.4);
        assert!(cf > 0.002 && cf < 0.005);
    }

    #[test]
    fn test_head_h1_roundtrip() {
        // H -> H1 -> H should be identity
        let h_orig = 1.5;
        let h1 = head_h1(h_orig);
        let h_back = head_h_from_h1(h1);
        assert!((h_orig - h_back).abs() < 0.01);
    }

    #[test]
    fn test_shape_factors() {
        let sf_lam = laminar_shape_factors(2.5);
        assert!(sf_lam.h > 0.0);
        assert!(sf_lam.hs > 0.0);
        
        let sf_turb = turbulent_shape_factors(1.4, 1000.0);
        assert!(sf_turb.h > 0.0);
        assert!(sf_turb.hs > 0.0);
    }
}
