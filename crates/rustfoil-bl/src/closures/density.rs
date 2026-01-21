//! Density shape factor closure (Hc)
//!
//! XFOIL Reference: xblsys.f HCT (line 2514)

/// Result of HCT computation including partial derivatives
#[derive(Debug, Clone, Copy)]
pub struct HctResult {
    /// Density thickness shape parameter Hc
    pub hc: f64,
    /// Partial derivative ∂Hc/∂Hk
    pub hc_hk: f64,
    /// Partial derivative ∂Hc/∂M²
    pub hc_msq: f64,
}

/// Density thickness shape parameter (Whitfield correlation)
///
/// Computes the density shape factor Hc which accounts for compressibility
/// effects on the boundary layer density profile.
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk
/// * `msq` - Mach number squared (M²)
///
/// # Returns
/// Density shape parameter and its partial derivatives
///
/// # Reference
/// XFOIL xblsys.f line 2514
pub fn density_shape(hk: f64, msq: f64) -> HctResult {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / (hk - 0.8).powi(2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;

    HctResult { hc, hc_hk, hc_msq }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hct_incompressible() {
        // At M=0, Hc should be 0 (no density thickness effect)
        let result = density_shape(2.5, 0.0);
        assert_relative_eq!(result.hc, 0.0, epsilon = 1e-15);
        // hc_hk should also be 0 when msq=0
        assert_relative_eq!(result.hc_hk, 0.0, epsilon = 1e-15);
        // hc_msq should be non-zero (it's the coefficient)
        let expected_hc_msq = 0.064 / (2.5 - 0.8) + 0.251;
        assert_relative_eq!(result.hc_msq, expected_hc_msq, epsilon = 1e-10);
    }

    #[test]
    fn test_hct_compressible() {
        let result = density_shape(2.5, 0.25);
        // Hc = 0.25 * (0.064/(2.5-0.8) + 0.251) = 0.25 * (0.0376... + 0.251)
        let expected = 0.25 * (0.064 / 1.7 + 0.251);
        assert_relative_eq!(result.hc, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_hct_derivatives() {
        let hk = 2.5;
        let msq = 0.25;
        let eps = 1e-7;

        let result = density_shape(hk, msq);

        // Check hc_hk via central difference
        let r1 = density_shape(hk - eps, msq);
        let r2 = density_shape(hk + eps, msq);
        let numerical_hk = (r2.hc - r1.hc) / (2.0 * eps);
        assert_relative_eq!(result.hc_hk, numerical_hk, epsilon = 1e-5);

        // Check hc_msq via central difference
        let r1 = density_shape(hk, msq - eps);
        let r2 = density_shape(hk, msq + eps);
        let numerical_msq = (r2.hc - r1.hc) / (2.0 * eps);
        assert_relative_eq!(result.hc_msq, numerical_msq, epsilon = 1e-5);
    }
}
