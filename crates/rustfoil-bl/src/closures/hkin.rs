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
/// boundary layer shape. This is XFOIL's compressibility correction from
/// Whitfield's correlation.
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
        // At M=0: hk_msq = (-0.29 - 0.113*2.5) / 1.0 = -0.29 - 0.2825 = -0.5725
        assert_relative_eq!(result.hk_msq, -0.29 - 0.113 * 2.5, epsilon = 1e-10);
    }

    #[test]
    fn test_hkin_compressible() {
        // At M=0.5 (M²=0.25)
        let result = hkin(2.5, 0.25);
        // Expected: (2.5 - 0.29*0.25) / (1 + 0.113*0.25)
        //         = (2.5 - 0.0725) / 1.02825
        //         = 2.4275 / 1.02825
        let expected_hk = 2.4275 / 1.02825;
        assert_relative_eq!(result.hk, expected_hk, epsilon = 1e-10);
    }

    #[test]
    fn test_hkin_derivative_h() {
        // Numerical derivative check for ∂Hk/∂H
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
        // Numerical derivative check for ∂Hk/∂M²
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
