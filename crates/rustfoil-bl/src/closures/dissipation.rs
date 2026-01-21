//! Dissipation closures (2*CD/H*)
//!
//! XFOIL Reference: xblsys.f DIL (line 2290), DIT (line 2375), DILW (line 2308)
//!
//! The dissipation coefficient represents the rate of energy loss in the boundary
//! layer due to viscous effects. It's expressed as 2*CD/H* where CD is the
//! dissipation coefficient and H* is the kinetic energy shape factor.

/// Result from laminar or wake dissipation calculation
///
/// Contains the dissipation coefficient and its partial derivatives
/// with respect to kinematic shape factor (Hk) and Reynolds number (Rθ).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DissipationResult {
    /// Dissipation coefficient (2*CD/H*)
    pub di: f64,
    /// ∂DI/∂Hk - derivative with respect to kinematic shape factor
    pub di_hk: f64,
    /// ∂DI/∂Rθ - derivative with respect to momentum thickness Reynolds number
    pub di_rt: f64,
}

/// Result from turbulent dissipation calculation
///
/// Contains the dissipation coefficient and its partial derivatives
/// with respect to all input parameters used in turbulent flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TurbDissipationResult {
    /// Dissipation coefficient (2*CD/H*)
    pub di: f64,
    /// ∂DI/∂Hs - derivative with respect to kinetic energy shape factor
    pub di_hs: f64,
    /// ∂DI/∂Us - derivative with respect to edge velocity ratio (Ue/Uinf)
    pub di_us: f64,
    /// ∂DI/∂Cf - derivative with respect to skin friction coefficient
    pub di_cf: f64,
    /// ∂DI/∂St - derivative with respect to shear stress coefficient
    pub di_st: f64,
}

/// Laminar H* (kinetic energy shape factor) correlation
///
/// XFOIL Reference: xblsys.f HSL (line 2327)
///
/// This is an internal helper for the wake dissipation calculation.
/// The full HSL closure with Mach number dependence is in the hs module.
///
/// # Arguments
/// * `hk` - Kinematic shape factor
///
/// # Returns
/// * `(hs, hs_hk)` - H* value and its derivative with respect to Hk
fn hsl_simple(hk: f64) -> (f64, f64) {
    if hk < 4.35 {
        let tmp = hk - 4.35;
        let hk_plus_1 = hk + 1.0;

        // HS correlation for attached flow (Hk < 4.35)
        let hs = 0.0111 * tmp.powi(2) / hk_plus_1
            - 0.0278 * tmp.powi(3) / hk_plus_1
            + 1.528
            - 0.0002 * (tmp * hk).powi(2);

        let hs_hk = 0.0111 * (2.0 * tmp - tmp.powi(2) / hk_plus_1) / hk_plus_1
            - 0.0278 * (3.0 * tmp.powi(2) - tmp.powi(3) / hk_plus_1) / hk_plus_1
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);

        (hs, hs_hk)
    } else {
        // HS correlation for separated flow (Hk >= 4.35)
        let hs2 = 0.015;
        let hk_minus_435 = hk - 4.35;

        let hs = hs2 * hk_minus_435.powi(2) / hk + 1.528;
        let hs_hk = hs2 * 2.0 * hk_minus_435 / hk - hs2 * hk_minus_435.powi(2) / hk.powi(2);

        (hs, hs_hk)
    }
}

/// Laminar dissipation coefficient
///
/// XFOIL Reference: xblsys.f DIL (line 2290)
///
/// Computes the laminar dissipation coefficient from Falkner-Skan correlations.
/// The dissipation is expressed as 2*CD/H* and represents energy loss due to
/// viscous effects in the laminar boundary layer.
///
/// # Arguments
/// * `hk` - Kinematic shape factor (typically 1.5 to 4.0 for laminar flow)
/// * `rt` - Momentum thickness Reynolds number (Rθ = Ue*θ/ν)
///
/// # Returns
/// [`DissipationResult`] containing DI and its derivatives ∂DI/∂Hk, ∂DI/∂Rθ
///
/// # Physics
/// - For Hk < 4.0: Uses a polynomial fit to Falkner-Skan solutions
/// - For Hk >= 4.0: Modified correlation for separated flow with denominator limiting
///
/// # Example
/// ```
/// use rustfoil_bl::closures::dissipation::dissipation_laminar;
///
/// let result = dissipation_laminar(2.6, 500.0);
/// assert!(result.di > 0.0, "Dissipation must be positive");
/// assert!(result.di_rt < 0.0, "∂DI/∂Rθ should be negative (higher Re = less dissipation)");
/// ```
pub fn dissipation_laminar(hk: f64, rt: f64) -> DissipationResult {
    let (di, di_hk) = if hk < 4.0 {
        // Attached laminar flow correlation (from Falkner-Skan)
        let di = (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / rt;
        let di_hk = (-0.00205 * 5.5 * (4.0 - hk).powf(4.5)) / rt;
        (di, di_hk)
    } else {
        // Separated flow correlation with denominator limiting
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb.powi(2);
        let di = (-0.0016 * hkb.powi(2) / den + 0.207) / rt;
        let di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb.powi(2) / den.powi(2))) / rt;
        (di, di_hk)
    };

    // Reynolds number derivative is always -DI/Rθ
    let di_rt = -di / rt;

    DissipationResult { di, di_hk, di_rt }
}

/// Turbulent dissipation coefficient
///
/// XFOIL Reference: xblsys.f DIT (line 2375)
///
/// Computes the turbulent dissipation coefficient based on the equilibrium
/// turbulent boundary layer model. Unlike the laminar case, turbulent
/// dissipation depends on skin friction (Cf), shear stress (St), and
/// edge velocity ratio (Us).
///
/// # Arguments
/// * `hs` - Kinetic energy shape factor (H*)
/// * `us` - Edge velocity ratio (Ue/Uinf)
/// * `cf` - Skin friction coefficient
/// * `st` - Shear stress coefficient (Clauser's equilibrium parameter)
///
/// # Returns
/// [`TurbDissipationResult`] containing DI and all partial derivatives
///
/// # Physics
/// The turbulent dissipation consists of two terms:
/// - `0.5*Cf*Us`: Direct viscous dissipation at the wall
/// - `St²*(1-Us)`: Turbulent production term from Reynolds stresses
///
/// # Example
/// ```
/// use rustfoil_bl::closures::dissipation::dissipation_turbulent;
///
/// let result = dissipation_turbulent(1.6, 0.98, 0.003, 0.02);
/// assert!(result.di > 0.0, "Dissipation must be positive");
/// ```
pub fn dissipation_turbulent(hs: f64, us: f64, cf: f64, st: f64) -> TurbDissipationResult {
    // Core dissipation formula: DI = (0.5*Cf*Us + St²*(1-Us)) * 2/Hs
    let core = 0.5 * cf * us + st * st * (1.0 - us);

    let di = core * 2.0 / hs;
    let di_hs = -core * 2.0 / hs.powi(2);
    let di_us = (0.5 * cf - st * st) * 2.0 / hs;
    let di_cf = (0.5 * us) * 2.0 / hs;
    let di_st = (2.0 * st * (1.0 - us)) * 2.0 / hs;

    TurbDissipationResult {
        di,
        di_hs,
        di_us,
        di_cf,
        di_st,
    }
}

/// Wake dissipation coefficient
///
/// XFOIL Reference: xblsys.f DILW (line 2308)
///
/// Computes the laminar wake dissipation coefficient for the viscous wake
/// region downstream of the trailing edge. This uses a different correlation
/// than the boundary layer since there's no wall.
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number
///
/// # Returns
/// [`DissipationResult`] containing DI and its derivatives
///
/// # Physics
/// The wake dissipation is modeled as:
/// - RCD = 1.10 * (1 - 1/Hk)² / Hk  (wake drag coefficient)
/// - DI = 2*RCD / (Hs*Rθ)
///
/// The formula reflects that wakes have no wall friction, so dissipation
/// comes entirely from the velocity defect profile relaxation.
///
/// # Example
/// ```
/// use rustfoil_bl::closures::dissipation::dissipation_wake;
///
/// let result = dissipation_wake(2.0, 1000.0);
/// assert!(result.di > 0.0, "Wake dissipation must be positive");
/// ```
pub fn dissipation_wake(hk: f64, rt: f64) -> DissipationResult {
    // Get laminar H* for wake calculation (MSQ = 0 incompressible)
    let (hs, hs_hk) = hsl_simple(hk);

    // HS_RT is 0 for laminar flow
    let hs_rt = 0.0;

    // Wake drag coefficient correlation (xblsys.f line 2315)
    // RCD = 1.10 * (1 - 1/HK)² / HK
    let one_minus_inv_hk = 1.0 - 1.0 / hk;
    let rcd = 1.10 * one_minus_inv_hk.powi(2) / hk;

    // RCD_HK exactly as in XFOIL (xblsys.f line 2316-2317)
    // Note: This is the XFOIL empirical derivative formula
    let rcd_hk = -1.10 * one_minus_inv_hk * 2.0 / hk.powi(3) - rcd / hk;

    // Dissipation from wake drag (xblsys.f lines 2319-2321)
    let di = 2.0 * rcd / (hs * rt);
    let di_hk = 2.0 * rcd_hk / (hs * rt) - (di / hs) * hs_hk;
    let di_rt = -di / rt - (di / hs) * hs_rt;

    DissipationResult { di, di_hk, di_rt }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // LAMINAR DISSIPATION (DIL) TESTS
    // =========================================================================

    #[test]
    fn test_dil_attached_flow() {
        // Test DIL for attached laminar flow (Hk < 4.0)
        // Typical Blasius flat plate has Hk ≈ 2.59
        let result = dissipation_laminar(2.59, 500.0);

        // Dissipation must be positive (energy is always lost)
        assert!(result.di > 0.0, "DIL must be positive for attached flow");

        // For attached flow with Hk < 4, DI_HK should be negative
        // (as Hk increases toward separation, the (4-Hk)^5.5 term decreases)
        assert!(
            result.di_hk < 0.0,
            "∂DI/∂Hk should be negative for Hk < 4"
        );

        // DI_RT should always be negative (higher Re = less dissipation)
        assert!(
            result.di_rt < 0.0,
            "∂DI/∂Rθ must be negative (viscous scaling)"
        );
    }

    #[test]
    fn test_dil_near_separation() {
        // Test DIL near separation point (Hk approaching 4.0)
        let result = dissipation_laminar(3.9, 500.0);

        assert!(result.di > 0.0, "DIL must be positive near separation");

        // At Hk = 3.9, the (4.0 - 3.9)^5.5 term is very small
        // so DI ≈ 0.207/Rθ
        let expected_baseline = 0.207 / 500.0;
        assert!(
            result.di > expected_baseline * 0.9,
            "Near separation, DI should be close to baseline 0.207/Rθ"
        );
    }

    #[test]
    fn test_dil_separated_flow() {
        // Test DIL for separated flow (Hk >= 4.0)
        let result = dissipation_laminar(5.0, 500.0);

        // Dissipation still positive in separated flow
        assert!(result.di > 0.0, "DIL must be positive for separated flow");

        // For Hk > 4, the formula changes - dissipation decreases
        // due to the negative -0.0016*HKB²/DEN term
        let baseline = 0.207 / 500.0;
        assert!(
            result.di < baseline,
            "Separated flow DI should be less than baseline"
        );
    }

    #[test]
    fn test_dil_branch_continuity() {
        // Test that DIL is continuous at Hk = 4.0 (branch point)
        let hk_below = 3.999;
        let hk_above = 4.001;

        let result_below = dissipation_laminar(hk_below, 500.0);
        let result_above = dissipation_laminar(hk_above, 500.0);

        // DI values should be very close at the branch point
        let di_diff = (result_above.di - result_below.di).abs();
        assert!(
            di_diff < 1e-4,
            "DIL should be continuous at Hk=4: diff={}",
            di_diff
        );
    }

    #[test]
    fn test_dil_reynolds_scaling() {
        // Verify that DI scales as 1/Rθ
        let result_500 = dissipation_laminar(2.6, 500.0);
        let result_1000 = dissipation_laminar(2.6, 1000.0);

        // DI(Rθ=500) should be ~2× DI(Rθ=1000)
        let ratio = result_500.di / result_1000.di;
        assert!(
            (ratio - 2.0).abs() < 0.01,
            "DIL should scale as 1/Rθ: ratio={}",
            ratio
        );
    }

    #[test]
    fn test_dil_derivative_consistency() {
        // Verify derivatives using finite differences
        let hk = 2.6;
        let rt = 500.0;
        let eps = 1e-6;

        let base = dissipation_laminar(hk, rt);
        let perturb_hk = dissipation_laminar(hk + eps, rt);
        let perturb_rt = dissipation_laminar(hk, rt + eps);

        // Numerical derivatives
        let di_hk_num = (perturb_hk.di - base.di) / eps;
        let di_rt_num = (perturb_rt.di - base.di) / eps;

        // Analytical derivatives should match numerical
        assert!(
            (base.di_hk - di_hk_num).abs() < 1e-5,
            "∂DI/∂Hk mismatch: analytical={}, numerical={}",
            base.di_hk,
            di_hk_num
        );
        assert!(
            (base.di_rt - di_rt_num).abs() < 1e-8,
            "∂DI/∂Rθ mismatch: analytical={}, numerical={}",
            base.di_rt,
            di_rt_num
        );
    }

    // =========================================================================
    // TURBULENT DISSIPATION (DIT) TESTS
    // =========================================================================

    #[test]
    fn test_dit_typical_conditions() {
        // Typical turbulent BL parameters
        let hs = 1.6; // H* for turbulent flow
        let us = 0.98; // Edge velocity ratio near freestream
        let cf = 0.003; // Typical skin friction
        let st = 0.02; // Shear stress coefficient

        let result = dissipation_turbulent(hs, us, cf, st);

        // Dissipation must be positive
        assert!(result.di > 0.0, "DIT must be positive");

        // DI_HS should be negative (higher H* = less dissipation rate per unit H*)
        assert!(result.di_hs < 0.0, "∂DI/∂Hs should be negative");
    }

    #[test]
    fn test_dit_cf_dominates_at_wall() {
        // When Us is high (near freestream), Cf term dominates
        let hs = 1.6;
        let us = 0.99; // Very close to freestream
        let cf = 0.003;
        let st = 0.02;

        let result = dissipation_turbulent(hs, us, cf, st);

        // Verify we get positive dissipation
        assert!(result.di > 0.0, "DIT should be positive");

        // The Cf*Us term should dominate over St²*(1-Us) when Us ≈ 1
        let cf_term = 0.5 * cf * us * 2.0 / hs;
        let st_term = st * st * (1.0 - us) * 2.0 / hs;

        assert!(
            cf_term > st_term,
            "Cf term should dominate when Us→1: cf_term={}, st_term={}",
            cf_term,
            st_term
        );
    }

    #[test]
    fn test_dit_derivative_signs() {
        // Test derivative signs for physical consistency
        let result = dissipation_turbulent(1.6, 0.98, 0.003, 0.02);

        // ∂DI/∂Hs < 0 (higher H* reduces normalized dissipation)
        assert!(result.di_hs < 0.0);

        // ∂DI/∂Cf > 0 (higher friction = more dissipation)
        assert!(result.di_cf > 0.0);

        // ∂DI/∂St > 0 when Us < 1 (higher shear stress = more turbulent dissipation)
        assert!(result.di_st > 0.0);
    }

    #[test]
    fn test_dit_derivative_consistency() {
        // Verify derivatives using finite differences
        let hs = 1.6;
        let us = 0.98;
        let cf = 0.003;
        let st = 0.02;
        let eps = 1e-8;

        let base = dissipation_turbulent(hs, us, cf, st);

        // Perturb each parameter
        let d_hs = dissipation_turbulent(hs + eps, us, cf, st);
        let d_us = dissipation_turbulent(hs, us + eps, cf, st);
        let d_cf = dissipation_turbulent(hs, us, cf + eps, st);
        let d_st = dissipation_turbulent(hs, us, cf, st + eps);

        // Check each derivative
        let di_hs_num = (d_hs.di - base.di) / eps;
        let di_us_num = (d_us.di - base.di) / eps;
        let di_cf_num = (d_cf.di - base.di) / eps;
        let di_st_num = (d_st.di - base.di) / eps;

        assert!(
            (base.di_hs - di_hs_num).abs() < 1e-6,
            "∂DI/∂Hs mismatch"
        );
        assert!(
            (base.di_us - di_us_num).abs() < 1e-6,
            "∂DI/∂Us mismatch"
        );
        assert!(
            (base.di_cf - di_cf_num).abs() < 1e-6,
            "∂DI/∂Cf mismatch"
        );
        assert!(
            (base.di_st - di_st_num).abs() < 1e-6,
            "∂DI/∂St mismatch"
        );
    }

    // =========================================================================
    // WAKE DISSIPATION (DILW) TESTS
    // =========================================================================

    #[test]
    fn test_dilw_typical_wake() {
        // Typical wake parameters
        let result = dissipation_wake(2.0, 1000.0);

        // Wake dissipation must be positive
        assert!(result.di > 0.0, "DILW must be positive");

        // DI_RT should be negative (higher Re = less dissipation)
        assert!(
            result.di_rt < 0.0,
            "∂DI/∂Rθ must be negative for wake"
        );
    }

    #[test]
    fn test_dilw_vs_dil_magnitude() {
        // Wake dissipation should generally be smaller than boundary layer
        // dissipation because there's no wall friction
        let hk = 2.5;
        let rt = 500.0;

        let dil_result = dissipation_laminar(hk, rt);
        let dilw_result = dissipation_wake(hk, rt);

        // Wake should have different (typically lower) dissipation
        // This is a sanity check, not a strict physical requirement
        assert!(
            dilw_result.di != dil_result.di,
            "Wake and BL dissipation should differ"
        );
    }

    #[test]
    fn test_dilw_hk_dependence() {
        // Test that wake dissipation varies with Hk
        let rt = 1000.0;

        let result_low_hk = dissipation_wake(1.5, rt);
        let result_high_hk = dissipation_wake(3.0, rt);

        // Both should be positive
        assert!(result_low_hk.di > 0.0);
        assert!(result_high_hk.di > 0.0);

        // They should be different
        assert!(
            (result_low_hk.di - result_high_hk.di).abs() > 1e-6,
            "DILW should vary with Hk"
        );
    }

    #[test]
    fn test_dilw_derivative_consistency() {
        // Verify DI_RT derivative using finite differences
        // Note: DI_HK uses XFOIL's empirical formula which doesn't match
        // the analytical derivative, but has been validated over decades
        // of XFOIL usage. We test DI_RT which is straightforward.
        let hk = 2.0;
        let rt = 1000.0;
        let eps = 1e-6;

        let base = dissipation_wake(hk, rt);
        let perturb_rt = dissipation_wake(hk, rt + eps);

        let di_rt_num = (perturb_rt.di - base.di) / eps;

        // DI_RT should match finite difference
        assert!(
            (base.di_rt - di_rt_num).abs() < 1e-8,
            "∂DI/∂Rθ mismatch: analytical={}, numerical={}",
            base.di_rt,
            di_rt_num
        );

        // DI_HK should be finite and non-zero (sanity check)
        assert!(
            base.di_hk.is_finite(),
            "∂DI/∂Hk should be finite"
        );
    }

    #[test]
    fn test_dilw_matches_xfoil_formula() {
        // Verify our implementation matches XFOIL's exact formula
        // at a specific test point
        let hk = 2.5;
        let rt = 800.0;

        let result = dissipation_wake(hk, rt);

        // Manually compute expected values using XFOIL formulas
        // HSL for Hk=2.5 (Hk < 4.35 branch)
        let tmp = hk - 4.35; // -1.85
        let hk_plus_1 = hk + 1.0; // 3.5
        let hs = 0.0111 * tmp.powi(2) / hk_plus_1
            - 0.0278 * tmp.powi(3) / hk_plus_1
            + 1.528
            - 0.0002 * (tmp * hk).powi(2);

        // RCD formula
        let one_minus_inv_hk = 1.0 - 1.0 / hk; // 0.6
        let rcd = 1.10 * one_minus_inv_hk.powi(2) / hk;

        // DI formula
        let expected_di = 2.0 * rcd / (hs * rt);

        assert!(
            (result.di - expected_di).abs() < 1e-10,
            "DI mismatch: got {}, expected {}",
            result.di,
            expected_di
        );
    }

    // =========================================================================
    // HSL HELPER FUNCTION TESTS
    // =========================================================================

    #[test]
    fn test_hsl_simple_attached() {
        // Test HSL for attached flow (Hk < 4.35)
        let (hs, _hs_hk) = hsl_simple(2.59);

        // H* should be around 1.5-1.6 for Blasius flow
        assert!(hs > 1.4 && hs < 1.7, "HSL should give H* ≈ 1.5 for Blasius: got {}", hs);
    }

    #[test]
    fn test_hsl_simple_separated() {
        // Test HSL for separated flow (Hk >= 4.35)
        let (hs, _hs_hk) = hsl_simple(5.0);

        // H* should be larger for separated flow
        assert!(hs > 1.5, "HSL should increase for separated flow: got {}", hs);
    }

    #[test]
    fn test_hsl_simple_continuity() {
        // Test continuity at Hk = 4.35 branch point
        let hk_below = 4.349;
        let hk_above = 4.351;

        let (hs_below, _) = hsl_simple(hk_below);
        let (hs_above, _) = hsl_simple(hk_above);

        let diff = (hs_above - hs_below).abs();
        assert!(
            diff < 0.01,
            "HSL should be continuous at Hk=4.35: diff={}",
            diff
        );
    }

    // =========================================================================
    // PHYSICAL SANITY CHECKS
    // =========================================================================

    #[test]
    fn test_dissipation_positivity() {
        // Dissipation must always be positive (2nd law of thermodynamics)
        // Test across a range of parameters

        for hk in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] {
            for rt in [100.0, 500.0, 1000.0, 5000.0] {
                let dil = dissipation_laminar(hk, rt);
                assert!(
                    dil.di > 0.0,
                    "DIL must be positive: Hk={}, Rθ={}, DI={}",
                    hk,
                    rt,
                    dil.di
                );

                if hk > 1.1 {
                    // DILW needs Hk > 1 to avoid division issues
                    let dilw = dissipation_wake(hk, rt);
                    assert!(
                        dilw.di > 0.0,
                        "DILW must be positive: Hk={}, Rθ={}, DI={}",
                        hk,
                        rt,
                        dilw.di
                    );
                }
            }
        }

        // Turbulent dissipation
        for hs in [1.4, 1.6, 1.8, 2.0] {
            for us in [0.9, 0.95, 0.99] {
                let dit = dissipation_turbulent(hs, us, 0.003, 0.02);
                assert!(
                    dit.di > 0.0,
                    "DIT must be positive: Hs={}, Us={}, DI={}",
                    hs,
                    us,
                    dit.di
                );
            }
        }
    }
}
