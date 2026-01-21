//! Energy shape factor closures (Hs)
//!
//! The energy shape factor Hs = δ*/θ* relates the energy thickness to the
//! momentum thickness. It appears in the kinetic energy integral equation
//! and is essential for computing the dissipation coefficient.
//!
//! XFOIL Reference: xblsys.f HSL (line 2327), HST (line 2388)

/// Result of energy shape factor computation
#[derive(Debug, Clone, Copy)]
pub struct HsResult {
    /// Energy shape factor Hs
    pub hs: f64,
    /// Partial derivative ∂Hs/∂Hk
    pub hs_hk: f64,
    /// Partial derivative ∂Hs/∂Rθ
    pub hs_rt: f64,
    /// Partial derivative ∂Hs/∂M²
    pub hs_msq: f64,
}

/// Laminar energy shape factor correlation
///
/// Computes Hs for laminar boundary layers using a polynomial fit.
/// The correlation has two branches meeting at Hk = 4.35:
/// - For Hk < 4.35: cubic polynomial (typical attached flow)
/// - For Hk >= 4.35: quadratic form (near separation)
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk
/// * `_rt` - Reynolds number based on θ (Rθ) - not used for laminar
/// * `_msq` - Mach² - not used for laminar
///
/// # Returns
/// Energy shape factor and derivatives. Note: `hs_rt` and `hs_msq` are
/// always zero for laminar flow (no Rθ or Mach dependence).
///
/// # Reference
/// XFOIL xblsys.f lines 2327-2351
pub fn hs_laminar(hk: f64, _rt: f64, _msq: f64) -> HsResult {
    let (hs, hs_hk) = if hk < 4.35 {
        // Low Hk branch (typical attached laminar flow)
        let tmp = hk - 4.35;
        let tmp2 = tmp * tmp;
        let tmp3 = tmp2 * tmp;
        let hk1 = hk + 1.0;

        // HS = 0.0111*TMP²/(HK+1) - 0.0278*TMP³/(HK+1) + 1.528 - 0.0002*(TMP*HK)²
        let hs = 0.0111 * tmp2 / hk1 - 0.0278 * tmp3 / hk1 + 1.528
            - 0.0002 * (tmp * hk).powi(2);

        // HS_HK = 0.0111*(2*TMP - TMP²/(HK+1))/(HK+1)
        //       - 0.0278*(3*TMP² - TMP³/(HK+1))/(HK+1)
        //       - 0.0002*2*TMP*HK*(TMP + HK)
        let hs_hk = 0.0111 * (2.0 * tmp - tmp2 / hk1) / hk1
            - 0.0278 * (3.0 * tmp2 - tmp3 / hk1) / hk1
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);

        (hs, hs_hk)
    } else {
        // High Hk branch (near laminar separation)
        // HS2 = 0.015 (coefficient from XFOIL)
        let hs2 = 0.015;
        let diff = hk - 4.35;

        // HS = HS2*(HK-4.35)²/HK + 1.528
        let hs = hs2 * diff.powi(2) / hk + 1.528;

        // HS_HK = HS2*2*(HK-4.35)/HK - HS2*(HK-4.35)²/HK²
        let hs_hk = hs2 * 2.0 * diff / hk - hs2 * diff.powi(2) / hk.powi(2);

        (hs, hs_hk)
    };

    HsResult {
        hs,
        hs_hk,
        hs_rt: 0.0,  // No Rθ dependence for laminar
        hs_msq: 0.0, // No Mach dependence for laminar
    }
}

/// Turbulent energy shape factor correlation
///
/// Computes Hs for turbulent boundary layers using XFOIL's 1991 correlation.
/// The correlation accounts for:
/// - Reynolds number (Rθ) dependence
/// - Attached vs separated flow regimes
/// - Compressibility effects (Whitfield correction)
///
/// # Arguments
/// * `hk` - Kinematic shape factor Hk
/// * `rt` - Reynolds number based on θ (Rθ)
/// * `msq` - Mach² (M²)
///
/// # Returns
/// Energy shape factor and all partial derivatives.
///
/// # Algorithm
/// 1. Compute HO threshold (varies with Rθ for Rθ > 400)
/// 2. Clip Rθ at minimum 200 (numerical stability)
/// 3. Branch on attached (Hk < HO) vs separated (Hk >= HO)
/// 4. Apply Whitfield compressibility correction
///
/// # Reference
/// XFOIL xblsys.f lines 2388-2479
pub fn hs_turbulent(hk: f64, rt: f64, msq: f64) -> HsResult {
    // Constants from XFOIL DATA statement
    const HSMIN: f64 = 1.500;
    const DHSINF: f64 = 0.015;

    // Step 1: Compute HO (threshold between attached/separated)
    // Limited Rθ dependence for Rθ < 400
    let (ho, ho_rt) = if rt > 400.0 {
        (3.0 + 400.0 / rt, -400.0 / (rt * rt))
    } else {
        (4.0, 0.0)
    };

    // Step 2: Clip RTZ at minimum 200 (numerical stability for low Re)
    let (rtz, rtz_rt) = if rt > 200.0 {
        (rt, 1.0)
    } else {
        (200.0, 0.0)
    };

    // Step 3: Compute Hs based on attached/separated regime
    let (mut hs, mut hs_hk, mut hs_rt) = if hk < ho {
        // ATTACHED branch (new correlation from 29 Nov 1991)
        // Based on arctan(y+) + Schlichting profiles
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;

        // HS = (2-HSMIN-4/RTZ)*HR²*1.5/(HK+0.5) + HSMIN + 4/RTZ
        let coef = 2.0 - HSMIN - 4.0 / rtz;
        let hr2 = hr * hr;
        let hk_term = 1.5 / (hk + 0.5);

        let hs = coef * hr2 * hk_term + HSMIN + 4.0 / rtz;

        // HS_HK = -(2-HSMIN-4/RTZ)*HR²*1.5/(HK+0.5)²
        //       + (2-HSMIN-4/RTZ)*HR*2*1.5/(HK+0.5)*HR_HK
        let hs_hk = -coef * hr2 * 1.5 / (hk + 0.5).powi(2)
            + coef * hr * 2.0 * hk_term * hr_hk;

        // HS_RT = (2-HSMIN-4/RTZ)*HR*2*1.5/(HK+0.5)*HR_RT
        //       + (HR²*1.5/(HK+0.5) - 1)*4/RTZ²*RTZ_RT
        let hs_rt = coef * hr * 2.0 * hk_term * hr_rt
            + (hr2 * hk_term - 1.0) * 4.0 / (rtz * rtz) * rtz_rt;

        (hs, hs_hk, hs_rt)
    } else {
        // SEPARATED branch
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;

        // HTMP = 0.007*GRT/RTMP² + DHSINF/HK
        let htmp = 0.007 * grt / (rtmp * rtmp) + DHSINF / hk;

        // HTMP_HK = -0.014*GRT/RTMP³ - DHSINF/HK²
        let htmp_hk = -0.014 * grt / rtmp.powi(3) - DHSINF / (hk * hk);

        // HTMP_RT = -0.014*GRT/RTMP³ * (-HO_RT - 4/GRT²/RTZ * RTZ_RT)
        //         + 0.007/RTMP² / RTZ * RTZ_RT
        let htmp_rt = -0.014 * grt / rtmp.powi(3) * (-ho_rt - 4.0 / (grt * grt) / rtz * rtz_rt)
            + 0.007 / (rtmp * rtmp) / rtz * rtz_rt;

        // HS = HDIF²*HTMP + HSMIN + 4/RTZ
        let hdif2 = hdif * hdif;
        let hs = hdif2 * htmp + HSMIN + 4.0 / rtz;

        // HS_HK = HDIF*2*HTMP + HDIF²*HTMP_HK
        let hs_hk = hdif * 2.0 * htmp + hdif2 * htmp_hk;

        // HS_RT = HDIF²*HTMP_RT - 4/RTZ²*RTZ_RT + HDIF*2*HTMP*(-HO_RT)
        let hs_rt = hdif2 * htmp_rt - 4.0 / (rtz * rtz) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);

        (hs, hs_hk, hs_rt)
    };

    // Step 4: Apply Whitfield's compressibility correction
    // FM = 1 + 0.014*MSQ
    // HS = (HS + 0.028*MSQ) / FM
    let fm = 1.0 + 0.014 * msq;
    hs = (hs + 0.028 * msq) / fm;
    hs_hk /= fm;
    hs_rt /= fm;
    let hs_msq = 0.028 / fm - 0.014 * hs / fm;

    HsResult {
        hs,
        hs_hk,
        hs_rt,
        hs_msq,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // =========================================================================
    // HSL (Laminar) Tests
    // =========================================================================

    #[test]
    fn test_hsl_low_hk() {
        // Test primary branch (Hk < 4.35) - typical attached laminar flow
        // At Hk = 2.5, we're well within the low branch
        let result = hs_laminar(2.5, 1000.0, 0.0);

        // Hs should be physically reasonable (between 1.5 and 2.0 for laminar)
        assert!(result.hs > 1.4 && result.hs < 2.0,
            "Laminar Hs at Hk=2.5 should be between 1.4 and 2.0, got {}", result.hs);

        // At Hk = 2.5, tmp = 2.5 - 4.35 = -1.85 (negative)
        // The cubic terms contribute to raising Hs above 1.528
    }

    #[test]
    fn test_hsl_high_hk() {
        // Test secondary branch (Hk >= 4.35) - near separation
        let result = hs_laminar(5.0, 1000.0, 0.0);

        // At Hk = 5.0: diff = 0.65, hs = 0.015 * 0.65² / 5.0 + 1.528
        let expected = 0.015 * 0.65_f64.powi(2) / 5.0 + 1.528;
        assert_relative_eq!(result.hs, expected, epsilon = 1e-10);

        // Derivative should be positive (Hs increases with Hk in this branch)
        assert!(result.hs_hk > 0.0, "dHs/dHk should be positive for Hk > 4.35");
    }

    #[test]
    fn test_hsl_at_branch_point() {
        // Test continuity at Hk = 4.35 where branches meet
        let eps = 1e-6;
        let result_below = hs_laminar(4.35 - eps, 0.0, 0.0);
        let result_above = hs_laminar(4.35 + eps, 0.0, 0.0);

        // Both branches should give Hs ≈ 1.528 at the branch point
        assert_relative_eq!(result_below.hs, 1.528, epsilon = 1e-4);
        assert_relative_eq!(result_above.hs, 1.528, epsilon = 1e-4);

        // Continuity check
        assert_relative_eq!(result_below.hs, result_above.hs, epsilon = 1e-3);
    }

    #[test]
    fn test_hsl_derivative_hk() {
        // Numerical derivative check using central differences
        let hk = 3.0;
        let eps = 1e-7;

        let r1 = hs_laminar(hk - eps, 0.0, 0.0);
        let r2 = hs_laminar(hk + eps, 0.0, 0.0);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_laminar(hk, 0.0, 0.0);
        assert_relative_eq!(result.hs_hk, numerical, epsilon = 1e-5);
    }

    #[test]
    fn test_hsl_derivative_hk_high_branch() {
        // Numerical derivative check in high Hk branch
        let hk = 5.0;
        let eps = 1e-7;

        let r1 = hs_laminar(hk - eps, 0.0, 0.0);
        let r2 = hs_laminar(hk + eps, 0.0, 0.0);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_laminar(hk, 0.0, 0.0);
        assert_relative_eq!(result.hs_hk, numerical, epsilon = 1e-5);
    }

    #[test]
    fn test_hsl_rt_msq_zero() {
        // Laminar Hs has no Rθ or Mach dependence
        let result = hs_laminar(2.5, 1000.0, 0.5);

        assert_eq!(result.hs_rt, 0.0, "Laminar hs_rt should be exactly 0");
        assert_eq!(result.hs_msq, 0.0, "Laminar hs_msq should be exactly 0");

        // Verify Hs doesn't change with rt or msq
        let result2 = hs_laminar(2.5, 5000.0, 0.9);
        assert_relative_eq!(result.hs, result2.hs, epsilon = 1e-15);
    }

    // =========================================================================
    // HST (Turbulent) Tests
    // =========================================================================

    #[test]
    fn test_hst_attached() {
        // Test attached turbulent BL (Hk < HO)
        // At RT = 1000, HO = 3.0 + 400/1000 = 3.4
        // Use Hk = 2.0, which is well in attached regime
        let result = hs_turbulent(2.0, 1000.0, 0.0);

        // Hs should be physically reasonable for turbulent attached flow
        assert!(result.hs > 1.5 && result.hs < 2.5,
            "Turbulent attached Hs should be between 1.5 and 2.5, got {}", result.hs);
    }

    #[test]
    fn test_hst_separated() {
        // Test separated turbulent BL (Hk >= HO)
        // At RT = 1000, HO = 3.4, so Hk = 4.0 is separated
        let result = hs_turbulent(4.0, 1000.0, 0.0);

        // Separated flow has higher Hs
        assert!(result.hs > 1.5,
            "Turbulent separated Hs should be > 1.5, got {}", result.hs);
    }

    #[test]
    fn test_hst_low_rt() {
        // Test Rθ < 200 clipping (RTZ floored at 200)
        let result_low = hs_turbulent(2.0, 100.0, 0.0);
        let _result_at_floor = hs_turbulent(2.0, 200.0, 0.0);

        // Both should give same result since RTZ is floored at 200
        // But HO differs: RT=100 gives HO=4.0, RT=200 gives HO=4.0 (both < 400)
        // So Hs values should be similar but not identical due to HO calculation
        assert!(result_low.hs.is_finite(), "Hs should be finite even at low Rθ");
        assert!(result_low.hs > 1.0, "Hs should be > 1.0 at low Rθ");

        // hs_rt should be 0 when RTZ_RT = 0 (below floor)
        // Actually this depends on the attached/separated branch
    }

    #[test]
    fn test_hst_ho_transition() {
        // Test HO computation for different RT values
        // RT > 400: HO = 3.0 + 400/RT
        // RT <= 400: HO = 4.0

        // At RT = 400, HO should transition
        let result_below = hs_turbulent(3.5, 399.0, 0.0);  // HO = 4.0
        let result_above = hs_turbulent(3.5, 401.0, 0.0);  // HO ≈ 3.997

        // Both should give finite, reasonable values
        assert!(result_below.hs.is_finite());
        assert!(result_above.hs.is_finite());
    }

    #[test]
    fn test_hst_compressibility() {
        // Test Whitfield's compressibility correction
        let result_m0 = hs_turbulent(2.5, 1000.0, 0.0);
        let result_m05 = hs_turbulent(2.5, 1000.0, 0.25);  // M ≈ 0.5
        let result_m08 = hs_turbulent(2.5, 1000.0, 0.64);  // M ≈ 0.8

        // Compressibility should increase Hs slightly
        // (HS + 0.028*MSQ) / (1 + 0.014*MSQ) > HS when MSQ > 0
        assert!(result_m05.hs > result_m0.hs,
            "Compressibility should increase Hs");
        assert!(result_m08.hs > result_m05.hs,
            "Higher Mach should give higher Hs");
    }

    #[test]
    fn test_hst_derivative_hk() {
        // Numerical derivative check for dHs/dHk in attached regime
        let hk = 2.5;
        let rt = 1000.0;
        let msq = 0.0;
        let eps = 1e-7;

        let r1 = hs_turbulent(hk - eps, rt, msq);
        let r2 = hs_turbulent(hk + eps, rt, msq);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_turbulent(hk, rt, msq);
        assert_relative_eq!(result.hs_hk, numerical, epsilon = 1e-4);
    }

    #[test]
    fn test_hst_derivative_hk_separated() {
        // Numerical derivative check for dHs/dHk in separated regime
        let hk = 4.5;  // Above HO ≈ 3.4 at RT=1000
        let rt = 1000.0;
        let msq = 0.0;
        let eps = 1e-7;

        let r1 = hs_turbulent(hk - eps, rt, msq);
        let r2 = hs_turbulent(hk + eps, rt, msq);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_turbulent(hk, rt, msq);
        assert_relative_eq!(result.hs_hk, numerical, epsilon = 1e-4);
    }

    #[test]
    fn test_hst_derivative_rt() {
        // Numerical derivative check for dHs/dRθ
        let hk = 2.5;
        let rt = 1000.0;
        let msq = 0.0;
        let eps = 1e-4;  // Larger eps for RT since it's a larger quantity

        let r1 = hs_turbulent(hk, rt - eps, msq);
        let r2 = hs_turbulent(hk, rt + eps, msq);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_turbulent(hk, rt, msq);
        assert_relative_eq!(result.hs_rt, numerical, epsilon = 1e-5);
    }

    #[test]
    fn test_hst_derivative_msq() {
        // Numerical derivative check for dHs/dM²
        let hk = 2.5;
        let rt = 1000.0;
        let msq = 0.25;
        let eps = 1e-7;

        let r1 = hs_turbulent(hk, rt, msq - eps);
        let r2 = hs_turbulent(hk, rt, msq + eps);
        let numerical = (r2.hs - r1.hs) / (2.0 * eps);

        let result = hs_turbulent(hk, rt, msq);
        assert_relative_eq!(result.hs_msq, numerical, epsilon = 1e-5);
    }

    // =========================================================================
    // Physical Validation Tests
    // =========================================================================

    #[test]
    fn test_hs_physical_bounds() {
        // Hs should be > 1.0 for physical boundary layers
        // Energy thickness should exceed momentum thickness

        // Test various laminar conditions
        for hk in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0] {
            let result = hs_laminar(hk, 0.0, 0.0);
            assert!(result.hs > 1.0,
                "Laminar Hs should be > 1.0, got {} at Hk={}", result.hs, hk);
        }

        // Test various turbulent conditions
        for hk in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0] {
            for rt in [500.0, 1000.0, 5000.0] {
                let result = hs_turbulent(hk, rt, 0.0);
                assert!(result.hs > 1.0,
                    "Turbulent Hs should be > 1.0, got {} at Hk={}, Rθ={}", 
                    result.hs, hk, rt);
            }
        }
    }

    #[test]
    fn test_hsl_monotonic_high_branch() {
        // In high Hk branch, Hs should increase monotonically
        let mut prev_hs = hs_laminar(4.35, 0.0, 0.0).hs;

        for hk in [4.5, 5.0, 6.0, 8.0, 10.0] {
            let result = hs_laminar(hk, 0.0, 0.0);
            assert!(result.hs >= prev_hs,
                "Laminar Hs should increase with Hk in high branch");
            prev_hs = result.hs;
        }
    }

    #[test]
    fn test_hst_hs_near_hsmin() {
        // At very high Rθ and low Hk, Hs approaches HSMIN + corrections
        // At RT=100000, HO = 3.0 + 400/100000 ≈ 3.004
        // At Hk=1.5, HR = (3.004-1.5)/(3.004-1) ≈ 0.75
        // The HR² term and 4/RTZ contribute to push Hs above HSMIN=1.5
        let result = hs_turbulent(1.5, 100000.0, 0.0);

        // Hs should be reasonably close to HSMIN but with positive corrections
        // The formula gives: (2-1.5-4/RTZ)*HR²*1.5/(HK+0.5) + 1.5 + 4/RTZ
        assert!(result.hs > 1.5 && result.hs < 2.0,
            "At high Rθ, low Hk, Hs should be between HSMIN=1.5 and 2.0, got {}", result.hs);
    }
}
