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

use super::hkin::hkin;

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

/// Result of AXSET RMS averaging computation
///
/// Contains the averaged amplification rate and intermediate values for debugging.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxsetResult {
    /// Individual amplification rate at station 1 (from DAMPL)
    pub ax1: f64,
    /// Individual amplification rate at station 2 (from DAMPL)
    pub ax2: f64,
    /// RMS average of AX1 and AX2
    pub axa_rms: f64,
    /// Near-Ncrit correction term
    pub dax: f64,
    /// Final averaged amplification rate (AXA + DAX)
    pub ax: f64,
}

/// Extended AXSET result with all derivatives for TRDIF
///
/// Contains the averaged amplification rate and its partial derivatives
/// with respect to all station variables. Used by TRCHEK2 and TRDIF.
///
/// # Reference
/// XFOIL xblsys.f AXSET (lines 35-144)
#[derive(Debug, Clone, Copy, Default)]
pub struct AxsetFullResult {
    /// Final averaged amplification rate
    pub ax: f64,
    
    // Derivatives w.r.t. station 1 variables
    /// ∂AX/∂Hk₁
    pub ax_hk1: f64,
    /// ∂AX/∂θ₁
    pub ax_t1: f64,
    /// ∂AX/∂Rθ₁
    pub ax_rt1: f64,
    /// ∂AX/∂N₁ (amplification)
    pub ax_a1: f64,
    
    // Derivatives w.r.t. station 2 variables
    /// ∂AX/∂Hk₂
    pub ax_hk2: f64,
    /// ∂AX/∂θ₂
    pub ax_t2: f64,
    /// ∂AX/∂Rθ₂
    pub ax_rt2: f64,
    /// ∂AX/∂N₂ (amplification)
    pub ax_a2: f64,
}

impl Default for AxsetResult {
    fn default() -> Self {
        Self {
            ax1: 0.0,
            ax2: 0.0,
            axa_rms: 0.0,
            dax: 0.0,
            ax: 0.0,
        }
    }
}

/// Compute averaged amplification rate over interval using RMS averaging
///
/// This is a direct port of XFOIL's AXSET subroutine (xblsys.f:35-144),
/// which computes an RMS-averaged amplification rate between two stations.
///
/// XFOIL uses this instead of simple midpoint averaging because RMS averaging
/// performs better on coarse grids and provides smoother transition prediction.
///
/// # Arguments
/// * `hk1`, `th1`, `rt1` - Boundary layer state at station 1
/// * `ampl1` - Current N-factor at station 1
/// * `hk2`, `th2`, `rt2` - Boundary layer state at station 2
/// * `ampl2` - Current N-factor at station 2
/// * `ncrit` - Critical N-factor for transition
///
/// # Returns
/// `AxsetResult` containing the averaged amplification rate and intermediate values
///
/// # Physics
///
/// The function computes:
/// 1. Individual amplification rates AX1 and AX2 at each station using DAMPL
/// 2. RMS average: AXA = sqrt(0.5*(AX1² + AX2²))
/// 3. A small correction DAX to ensure positive dN/dx near Ncrit:
///    DAX = exp(-20*(Ncrit - 0.5*(N1+N2))) * 0.002/(θ1+θ2)
/// 4. Final rate: AX = AXA + DAX
///
/// The DAX term prevents the solver from getting stuck when N is very close
/// to Ncrit and both AX1 and AX2 are near zero.
///
/// # Reference
/// XFOIL xblsys.f lines 35-144
pub fn axset(
    hk1: f64,
    th1: f64,
    rt1: f64,
    ampl1: f64,
    hk2: f64,
    th2: f64,
    rt2: f64,
    ampl2: f64,
    ncrit: f64,
) -> AxsetResult {
    // Compute individual amplification rates at each station
    let result1 = amplification_rate(hk1, th1, rt1);
    let result2 = amplification_rate(hk2, th2, rt2);

    let ax1 = result1.ax;
    let ax2 = result2.ax;

    // RMS average (XFOIL style - better on coarse grids)
    // AXSQ = 0.5*(AX1**2 + AX2**2)
    // AXA = SQRT(AXSQ)
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let axa_rms = if axsq > 0.0 { axsq.sqrt() } else { 0.0 };

    // Small additional term to ensure dN/dx > 0 near N = Ncrit
    // This prevents the solver from stalling when we're very close to transition
    //
    // ARG = MIN(20.0*(ACRIT - 0.5*(A1+A2)), 20.0)
    // EXN = EXP(-ARG)
    // DAX = EXN * 0.002/(T1+T2)
    let avg_ampl = 0.5 * (ampl1 + ampl2);
    let arg = (20.0 * (ncrit - avg_ampl)).min(20.0);
    let exn = if arg <= 0.0 {
        1.0 // When N >= Ncrit, exn = 1
    } else {
        (-arg).exp()
    };
    
    let th_sum = th1 + th2;
    let dax = if th_sum > 0.0 {
        exn * 0.002 / th_sum
    } else {
        0.0
    };

    // Final averaged amplification rate
    let ax = axa_rms + dax;

    AxsetResult {
        ax1,
        ax2,
        axa_rms,
        dax,
        ax,
    }
}

/// Compute averaged amplification rate with full derivatives
///
/// Extended version of `axset()` that returns all partial derivatives
/// needed for TRCHEK2 and TRDIF Jacobian calculations.
///
/// # Reference
/// XFOIL xblsys.f AXSET (lines 35-144)
pub fn axset_full(
    hk1: f64,
    th1: f64,
    rt1: f64,
    ampl1: f64,
    hk2: f64,
    th2: f64,
    rt2: f64,
    ampl2: f64,
    ncrit: f64,
) -> AxsetFullResult {
    // Compute individual amplification rates and their derivatives at each station
    let result1 = amplification_rate(hk1, th1, rt1);
    let result2 = amplification_rate(hk2, th2, rt2);

    let ax1 = result1.ax;
    let ax2 = result2.ax;
    let ax1_hk1 = result1.ax_hk;
    let ax1_t1 = result1.ax_th;
    let ax1_rt1 = result1.ax_rt;
    let ax2_hk2 = result2.ax_hk;
    let ax2_t2 = result2.ax_th;
    let ax2_rt2 = result2.ax_rt;

    // RMS average: AXA = sqrt(0.5*(AX1² + AX2²))
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let (axa, axa_ax1, axa_ax2) = if axsq > 0.0 {
        let axa = axsq.sqrt();
        (axa, 0.5 * ax1 / axa, 0.5 * ax2 / axa)
    } else {
        (0.0, 0.0, 0.0)
    };

    // Small additional term to ensure dN/dx > 0 near N = Ncrit
    // ARG = MIN(20.0*(ACRIT - 0.5*(A1+A2)), 20.0)
    let avg_ampl = 0.5 * (ampl1 + ampl2);
    let arg = (20.0 * (ncrit - avg_ampl)).min(20.0);
    
    let (exn, exn_a1, exn_a2) = if arg <= 0.0 {
        (1.0, 0.0, 0.0)
    } else {
        let exn = (-arg).exp();
        // dEXN/dA1 = d/dA1 exp(-20*(Ncrit - 0.5*(A1+A2))) = 20*0.5*EXN
        (exn, 20.0 * 0.5 * exn, 20.0 * 0.5 * exn)
    };

    let th_sum = th1 + th2;
    let (dax, dax_a1, dax_a2, dax_t1, dax_t2) = if th_sum > 0.0 {
        let dax = exn * 0.002 / th_sum;
        (
            dax,
            exn_a1 * 0.002 / th_sum,
            exn_a2 * 0.002 / th_sum,
            -dax / th_sum,
            -dax / th_sum,
        )
    } else {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    };

    // Final AX = AXA + DAX
    let ax = axa + dax;

    // Chain rule for final derivatives
    // AX_HK1 = AXA_AX1 * AX1_HK1
    // AX_T1  = AXA_AX1 * AX1_T1  + DAX_T1
    // AX_RT1 = AXA_AX1 * AX1_RT1
    // AX_A1  =                    DAX_A1
    let ax_hk1 = axa_ax1 * ax1_hk1;
    let ax_t1 = axa_ax1 * ax1_t1 + dax_t1;
    let ax_rt1 = axa_ax1 * ax1_rt1;
    let ax_a1 = dax_a1;

    let ax_hk2 = axa_ax2 * ax2_hk2;
    let ax_t2 = axa_ax2 * ax2_t2 + dax_t2;
    let ax_rt2 = axa_ax2 * ax2_rt2;
    let ax_a2 = dax_a2;

    AxsetFullResult {
        ax,
        ax_hk1,
        ax_t1,
        ax_rt1,
        ax_a1,
        ax_hk2,
        ax_t2,
        ax_rt2,
        ax_a2,
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

/// Result of TRCHEK2 transition check
///
/// Contains the computed N-factor at station 2, transition status,
/// and transition location if transition occurred.
#[derive(Debug, Clone, Copy)]
pub struct Trchek2Result {
    /// Final N-factor at station 2 (may be > Ncrit if transition occurred)
    pub ampl2: f64,
    /// Whether transition was detected in this interval
    pub transition: bool,
    /// Transition location as arc length (if transition occurred)
    pub xt: Option<f64>,
    /// Final amplification rate used (for debugging)
    pub ax: f64,
    /// Number of iterations used
    pub iterations: usize,
    /// Whether the implicit iteration converged
    pub converged: bool,
}

/// Full TRCHEK2 result with all derivatives for TRDIF
///
/// Contains the transition point XT, weighting factors WF1/WF2,
/// and all derivatives needed for the TRDIF Jacobian transformation.
///
/// # Variable naming convention
/// - `xt_*1` = ∂XT/∂(variable at station 1)
/// - `xt_*2` = ∂XT/∂(variable at station 2)
///
/// # Reference
/// XFOIL xblsys.f TRCHEK2 (lines 231-580)
#[derive(Debug, Clone, Copy, Default)]
pub struct Trchek2FullResult {
    /// Critical N-factor used for this transition solve
    pub ncrit: f64,
    /// Transition x-location
    pub xt: f64,
    /// Final N-factor at station 2
    pub ampl2: f64,
    /// Final averaged amplification rate (AX)
    pub ax: f64,
    /// Whether transition occurred
    pub transition: bool,
    /// Whether the selected transition is forced rather than free e^N
    pub forced: bool,
    /// Whether iteration converged
    pub converged: bool,
    /// Number of TRCHEK2 iterations used
    pub iterations: usize,
    
    /// Weighting factor: 1 - fraction from X1 to XT
    pub wf1: f64,
    /// Weighting factor: fraction from X1 to XT
    pub wf2: f64,
    
    // XT derivatives w.r.t. station 1 variables
    /// ∂XT/∂N₁ (amplification)
    pub xt_a1: f64,
    /// ∂XT/∂X₁
    pub xt_x1: f64,
    /// ∂XT/∂θ₁
    pub xt_t1: f64,
    /// ∂XT/∂δ*₁
    pub xt_d1: f64,
    /// ∂XT/∂Ue₁
    pub xt_u1: f64,
    
    // XT derivatives w.r.t. station 2 variables
    /// ∂XT/∂X₂
    pub xt_x2: f64,
    /// ∂XT/∂θ₂
    pub xt_t2: f64,
    /// ∂XT/∂δ*₂
    pub xt_d2: f64,
    /// ∂XT/∂Ue₂
    pub xt_u2: f64,
    
    // XT derivatives w.r.t. global parameters
    /// ∂XT/∂M² (Mach number squared)
    pub xt_ms: f64,
    /// ∂XT/∂Re (Reynolds number)
    pub xt_re: f64,
    
    // Interpolation derivatives for TT (theta at XT)
    /// ∂TT/∂N₁
    pub tt_a1: f64,
    /// ∂TT/∂X₁
    pub tt_x1: f64,
    /// ∂TT/∂X₂
    pub tt_x2: f64,
    /// ∂TT/∂θ₁
    pub tt_t1: f64,
    /// ∂TT/∂θ₂
    pub tt_t2: f64,
    /// ∂TT/∂δ*₁
    pub tt_d1: f64,
    /// ∂TT/∂δ*₂
    pub tt_d2: f64,
    /// ∂TT/∂Ue₁
    pub tt_u1: f64,
    /// ∂TT/∂Ue₂
    pub tt_u2: f64,
    /// ∂TT/∂M²
    pub tt_ms: f64,
    /// ∂TT/∂Re
    pub tt_re: f64,
    
    // Interpolation derivatives for DT (delta_star at XT)
    /// ∂DT/∂N₁
    pub dt_a1: f64,
    /// ∂DT/∂X₁
    pub dt_x1: f64,
    /// ∂DT/∂X₂
    pub dt_x2: f64,
    /// ∂DT/∂θ₁
    pub dt_t1: f64,
    /// ∂DT/∂θ₂
    pub dt_t2: f64,
    /// ∂DT/∂δ*₁
    pub dt_d1: f64,
    /// ∂DT/∂δ*₂
    pub dt_d2: f64,
    /// ∂DT/∂Ue₁
    pub dt_u1: f64,
    /// ∂DT/∂Ue₂
    pub dt_u2: f64,
    /// ∂DT/∂M²
    pub dt_ms: f64,
    /// ∂DT/∂Re
    pub dt_re: f64,
    
    // Interpolation derivatives for UT (Ue at XT)
    /// ∂UT/∂N₁
    pub ut_a1: f64,
    /// ∂UT/∂X₁
    pub ut_x1: f64,
    /// ∂UT/∂X₂
    pub ut_x2: f64,
    /// ∂UT/∂θ₁
    pub ut_t1: f64,
    /// ∂UT/∂θ₂
    pub ut_t2: f64,
    /// ∂UT/∂δ*₁
    pub ut_d1: f64,
    /// ∂UT/∂δ*₂
    pub ut_d2: f64,
    /// ∂UT/∂Ue₁
    pub ut_u1: f64,
    /// ∂UT/∂Ue₂
    pub ut_u2: f64,
    /// ∂UT/∂M²
    pub ut_ms: f64,
    /// ∂UT/∂Re
    pub ut_re: f64,
}

/// XFOIL's TRCHEK2 - 2nd order implicit N-factor integration
///
/// This is a direct port of XFOIL's TRCHEK2 subroutine (xblsys.f:235-466),
/// which uses a 2nd-order implicit scheme to compute N-factor evolution
/// and detect transition.
///
/// The implicit equation solved is:
/// ```text
///     N2 - N1     N'(XT,NT) + N'(X1,N1)
///     -------  =  ---------------------
///     X2 - X1               2
/// ```
///
/// Where XT,NT are either (X2,N2) if no transition, or interpolated values
/// at the transition point if N2 >= Ncrit.
///
/// # Arguments
/// * `x1`, `x2` - Arc length positions at station 1 and 2
/// * `hk1`, `hk2` - Kinematic shape factors
/// * `t1`, `t2` - Momentum thickness θ
/// * `rt1`, `rt2` - Momentum thickness Reynolds number Rθ
/// * `d1`, `d2` - Displacement thickness δ*
/// * `u1`, `u2` - Edge velocity
/// * `ampl1` - N-factor at station 1
/// * `ncrit` - Critical N-factor for transition
/// * `x_forced` - Forced transition location (None for free transition)
/// * `re` - Reynolds number (for computing Rθ at interpolated point)
///
/// # Returns
/// `Trchek2Result` containing the computed N-factor, transition status, and location
///
/// # Reference
/// XFOIL xblsys.f TRCHEK2 subroutine (lines 235-466)
pub fn trchek2(
    x1: f64,
    x2: f64,
    hk1: f64,
    hk2: f64,
    t1: f64,
    t2: f64,
    rt1: f64,
    rt2: f64,
    d1: f64,
    d2: f64,
    u1: f64,
    u2: f64,
    ampl1: f64,
    ncrit: f64,
    x_forced: Option<f64>,
    re: f64,
) -> Trchek2Result {
    const DAEPS: f64 = 5.0e-5; // Convergence tolerance
    const MAX_ITER: usize = 30;
    
    let dx = x2 - x1;
    if dx <= 0.0 {
        // Invalid interval
        return Trchek2Result {
            ampl2: ampl1,
            transition: false,
            xt: None,
            ax: 0.0,
            iterations: 0,
            converged: true,
        };
    }
    
    // Calculate initial average amplification rate
    let ax_result = axset(hk1, t1, rt1, ampl1, hk2, t2, rt2, ampl1, ncrit);
    
    // Initial guess for N2
    let mut ampl2 = ampl1 + ax_result.ax * dx;
    
    let mut ax = ax_result.ax;
    let mut xt = x2;
    let mut converged = false;
    let mut iterations = 0;
    
    // Implicit iteration for N2
    for itam in 0..MAX_ITER {
        iterations = itam + 1;
        
        // Define weighting factors based on whether transition occurred
        let (sfa, amplt) = if ampl2 <= ncrit {
            // No transition yet - "T" is the same as "2"
            (1.0, ampl2)
        } else {
            // Transition in interval - "T" is at Ncrit
            let sfa = (ncrit - ampl1) / (ampl2 - ampl1).max(1e-20);
            (sfa, ncrit)
        };
        
        // Check for forced transition
        let sfx = match x_forced {
            Some(xf) if xf < x2 => (xf - x1) / dx,
            _ => 1.0,
        };
        
        // Use the smaller of free or forced transition factor
        let wf2 = sfa.min(sfx);
        let wf1 = 1.0 - wf2;
        
        // Interpolate BL variables to transition point XT
        xt = x1 * wf1 + x2 * wf2;
        let tt = t1 * wf1 + t2 * wf2;
        let dt = d1 * wf1 + d2 * wf2;
        let ut = u1 * wf1 + u2 * wf2;
        
        // Compute Hk at the interpolated point using proper H = δ*/θ ratio
        // This is more accurate than linear interpolation of Hk because
        // the ratio of interpolated values ≠ interpolation of ratios
        let ht = if tt > 1e-20 { dt / tt } else { hk1 * wf1 + hk2 * wf2 };
        let hkt = hkin(ht, 0.0).hk;  // Assume incompressible (msq = 0)
        
        // Rθ at transition point (using interpolated values)
        let rtt = ut * tt * re;
        
        // Recompute amplification rate over X1..XT interval
        let ax_result = axset(hk1, t1, rt1, ampl1, hkt, tt, rtt, amplt, ncrit);
        ax = ax_result.ax;
        
        // Exit early if no amplification
        if ax <= 0.0 {
            converged = true;
            break;
        }
        
        // Residual for implicit N2 definition
        // Note: In XFOIL, they use (X2-X1) in the residual, not (XT-X1)
        // because XT depends on AMPL2
        let res = ampl2 - ampl1 - ax * dx;
        
        // Simple derivative approximation: RES_A2 ≈ 1 (neglecting AX dependence on A2)
        // This is simplified from XFOIL which computes full sensitivities
        let res_a2 = 1.0;
        
        // Newton update
        let da2 = -res / res_a2;
        
        // Check convergence
        if da2.abs() < DAEPS {
            converged = true;
            break;
        }
        
        // Relaxation to prevent oscillation
        let mut rlx = 1.0;
        if da2.abs() > 1.0 {
            rlx = 1.0 / da2.abs();
        }
        
        // Prevent crossing Ncrit boundary in a single step
        let ampl2_new = ampl2 + rlx * da2;
        if (ampl2 > ncrit && ampl2_new < ncrit) || (ampl2 < ncrit && ampl2_new > ncrit) {
            ampl2 = ncrit;
        } else {
            ampl2 = ampl2_new;
        }
    }
    
    // Test for free or forced transition
    let mut tr_free = ampl2 >= ncrit;
    let mut tr_forced = x_forced.is_some_and(|xf| xf > x1 && xf <= x2);
    let transition = tr_free || tr_forced;
    
    // Resolve if both forced and free transition
    let final_xt = if transition {
        if tr_forced && x_forced.is_some() {
            let xf = x_forced.unwrap();
            if !tr_free || xf < xt {
                Some(xf)
            } else {
                Some(xt)
            }
        } else {
            Some(xt)
        }
    } else {
        None
    };
    
    Trchek2Result {
        ampl2,
        transition,
        xt: final_xt,
        ax,
        iterations,
        converged,
    }
}

/// Simplified TRCHEK for use in fixed-Ue march
///
/// This is a simpler version that works with BlStation data directly,
/// avoiding the need to pass all individual variables.
///
/// # Arguments
/// * `x1`, `x2` - Arc length positions
/// * `hk1`, `hk2` - Kinematic shape factors
/// * `t1`, `t2` - Momentum thickness θ
/// * `rt1`, `rt2` - Momentum thickness Reynolds number Rθ
/// * `d1`, `d2` - Displacement thickness δ* (for proper Hk interpolation)
/// * `ampl1` - N-factor at station 1
/// * `ncrit` - Critical N-factor for transition
pub fn trchek2_stations(
    x1: f64,
    x2: f64,
    hk1: f64,
    t1: f64,
    rt1: f64,
    d1: f64,
    u1: f64,
    ampl1: f64,
    hk2: f64,
    t2: f64,
    rt2: f64,
    d2: f64,
    u2: f64,
    ncrit: f64,
    msq: f64,
    re: f64,
) -> Trchek2Result {
    const DAEPS: f64 = 5.0e-5;
    const MAX_ITER: usize = 30;
    
    let dx = x2 - x1;
    if dx <= 0.0 {
        return Trchek2Result {
            ampl2: ampl1,
            transition: false,
            xt: None,
            ax: 0.0,
            iterations: 0,
            converged: true,
        };
    }
    
    // Initial amplification rate (XFOIL xblsys.f:270-278)
    let ax_init = axset(hk1, t1, rt1, ampl1, hk2, t2, rt2, ampl1, ncrit);
    let mut ampl2 = ampl1 + ax_init.ax * dx;
    
    let mut ax = ax_init.ax;
    let mut xt = x2;
    let mut converged = false;
    let mut iterations = 0;
    
    for itam in 0..MAX_ITER {
        iterations = itam + 1;
        
        // Weighting factors (XFOIL xblsys.f:316-321)
        let (wf2, amplt, amplt_a2) = if ampl2 <= ncrit {
            (1.0, ampl2, 1.0)
        } else {
            let da = (ampl2 - ampl1).max(1e-20);
            let sfa = (ncrit - ampl1) / da;
            (sfa, ncrit, 0.0)
        };
        let wf1 = 1.0 - wf2;
        
        // dSFA/dA2 for interpolation sensitivities (XFOIL xblsys.f:322-330)
        let sfa_a2 = if ampl2 > ncrit {
            let da = (ampl2 - ampl1).max(1e-20);
            -(ncrit - ampl1) / (da * da)
        } else {
            0.0
        };
        
        // Interpolate to transition point
        xt = x1 * wf1 + x2 * wf2;
        let tt = t1 * wf1 + t2 * wf2;
        let dt = d1 * wf1 + d2 * wf2;
        let ut = u1 * wf1 + u2 * wf2;
        let rtt = re * ut * tt;
        
        // Sensitivities through interpolation (XFOIL xblsys.f:370-390)
        let tt_a2 = (t2 - t1) * sfa_a2;
        let dt_a2 = (d2 - d1) * sfa_a2;
        let ut_a2 = (u2 - u1) * sfa_a2;
        let rtt_tt = re * ut;
        let rtt_ut = re * tt;
        
        // Compute Hk at interpolated point
        let ht = if tt > 1e-20 { dt / tt } else { hk1 * wf1 + hk2 * wf2 };
        let hk_result = hkin(ht, msq);
        let hkt = hk_result.hk;
        let hkt_h = hk_result.hk_h;
        let hkt_tt = if tt > 1e-20 { hkt_h * (-dt / (tt * tt)) } else { 0.0 };
        let hkt_dt = if tt > 1e-20 { hkt_h * (1.0 / tt) } else { 0.0 };
        let hkt_ut = 0.0; // Ue doesn't affect H directly
        
        // Recompute AX with full derivatives at transition point
        let ax_result = axset_full(hk1, t1, rt1, ampl1, hkt, tt, rtt, amplt, ncrit);
        ax = ax_result.ax;
        
        if ax <= 0.0 {
            converged = true;
            break;
        }
        
        // Compute AX_A2 through the interpolation chain (XFOIL xblsys.f:396-400)
        let ax_hkt = ax_result.ax_hk2;
        let ax_tt = ax_result.ax_t2;
        let ax_rtt = ax_result.ax_rt2;
        let ax_at = ax_result.ax_a2;
        
        let ax_a2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2
            + (ax_hkt * hkt_dt) * dt_a2
            + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2
            + ax_at * amplt_a2;
        
        // Residual and proper Newton update (XFOIL xblsys.f:402-406)
        let res = ampl2 - ampl1 - ax * dx;
        let res_a2 = 1.0 - ax_a2 * dx;
        let da2 = if res_a2.abs() > 1e-20 { -res / res_a2 } else { -res };
        
        if da2.abs() < DAEPS {
            converged = true;
            break;
        }
        
        // Relaxed update
        let rlx = if da2.abs() > 1.0 { 1.0 / da2.abs() } else { 1.0 };
        let ampl2_new = ampl2 + rlx * da2;
        
        // Prevent crossing Ncrit boundary
        if (ampl2 > ncrit && ampl2_new < ncrit) || (ampl2 < ncrit && ampl2_new > ncrit) {
            ampl2 = ncrit;
        } else {
            ampl2 = ampl2_new;
        }
    }
    
    let transition = ampl2 >= ncrit;
    
    Trchek2Result {
        ampl2,
        transition,
        xt: if transition { Some(xt) } else { None },
        ax,
        iterations,
        converged,
    }
}

/// Full TRCHEK2 with all derivatives for TRDIF
///
/// This is the complete version of TRCHEK2 that returns all the derivatives
/// needed for the TRDIF Jacobian transformation. It computes:
/// - XT (transition location) and its derivatives
/// - WF1, WF2 (weighting factors) and their derivatives  
/// - TT, DT, UT (interpolated variables) and their derivatives
///
/// # Arguments
/// * `x1`, `x2` - Arc length positions
/// * `t1`, `t2` - Momentum thickness θ
/// * `d1`, `d2` - Displacement thickness δ*
/// * `u1`, `u2` - Edge velocity Ue
/// * `hk1`, `hk2` - Kinematic shape factors (for Hk derivatives)
/// * `rt1`, `rt2` - Momentum thickness Reynolds number Rθ
/// * `ampl1` - N-factor at station 1
/// * `ncrit` - Critical N-factor for transition
/// * `msq` - Mach number squared (for compressibility effects)
/// * `re` - Reynolds number
///
/// # Reference
/// XFOIL xblsys.f TRCHEK2 (lines 231-580)
pub fn trchek2_full(
    x1: f64,
    x2: f64,
    t1: f64,
    t2: f64,
    d1: f64,
    d2: f64,
    u1: f64,
    u2: f64,
    hk1: f64,
    hk2: f64,
    rt1: f64,
    rt2: f64,
    ampl1: f64,
    ncrit: f64,
    x_forced: Option<f64>,
    msq: f64,
    re: f64,
) -> Trchek2FullResult {
    const DAEPS: f64 = 5.0e-5;
    const MAX_ITER: usize = 30;

    let dx = x2 - x1;
    if dx <= 1e-20 {
        return Trchek2FullResult::default();
    }

    // Hk derivatives from hkin
    let hkin1 = hkin(d1 / t1.max(1e-20), msq);
    let hkin2 = hkin(d2 / t2.max(1e-20), msq);
    
    // H = δ*/θ, so dH/dθ = -H/θ, dH/dδ* = 1/θ
    let h1 = d1 / t1.max(1e-20);
    let h2 = d2 / t2.max(1e-20);
    let h1_t1 = -h1 / t1.max(1e-20);
    let h1_d1 = 1.0 / t1.max(1e-20);
    let h2_t2 = -h2 / t2.max(1e-20);
    let h2_d2 = 1.0 / t2.max(1e-20);
    
    let hk1_t1 = hkin1.hk_h * h1_t1;
    let hk1_d1 = hkin1.hk_h * h1_d1;
    let hk1_ms = hkin1.hk_msq;
    
    let hk2_t2 = hkin2.hk_h * h2_t2;
    let hk2_d2 = hkin2.hk_h * h2_d2;
    let hk2_ms = hkin2.hk_msq;
    
    let reybl = re;
    // Rθ = RE * U * θ, matching XFOIL's RT2 = RE * U2 * T2.
    let rt1_t1 = reybl * u1;
    let rt1_u1 = reybl * t1;
    let rt1_re = u1 * t1;
    let rt1_ms = 0.0;
    
    let rt2_t2 = reybl * u2;
    let rt2_u2 = reybl * t2;
    let rt2_re = u2 * t2;
    let rt2_ms = 0.0;

    // Initial amplification rate with derivatives
    let ax_init = axset_full(hk1, t1, rt1, ampl1, hk2, t2, rt2, ampl1, ncrit);
    let mut ampl2 = ampl1 + ax_init.ax * dx;

    let mut ax = ax_init.ax;
    let mut xt = x2;
    let mut converged = false;
    let mut iterations = 0usize;
    
    // Weighting factor derivatives (will be computed after convergence)
    let mut wf2 = 1.0;
    let mut wf1 = 0.0;
    let mut wf2_a1 = 0.0;
    let mut wf2_a2 = 0.0;
    let mut wf2_x1 = 0.0;
    let mut wf2_x2 = 0.0;
    let mut amplt_a2 = 1.0;
    
    // Interpolated variables at transition
    let mut tt = t2;
    let mut dt = d2;
    let mut ut = u2;
    
    // XT_A2 for implicit function theorem
    let mut xt_a2 = 0.0;
    let mut tt_a2 = 0.0;
    let mut dt_a2 = 0.0;
    let mut ut_a2 = 0.0;

    // === Implicit iteration for AMPL2 ===
    for itam in 0..MAX_ITER {
        iterations = itam + 1;
        // Define weighting factors WF1, WF2 for interpolation to transition point
        let (sfa, sfa_a1, sfa_a2, amplt_local, amplt_a2_local) = if ampl2 <= ncrit {
            // No transition yet, "T" is the same as "2"
            (1.0, 0.0, 0.0, ampl2, 1.0)
        } else {
            // Transition in interval, "T" is at Ncrit
            let sfa = (ncrit - ampl1) / (ampl2 - ampl1).max(1e-20);
            let denom = (ampl2 - ampl1).max(1e-20);
            (sfa, (sfa - 1.0) / denom, -sfa / denom, ncrit, 0.0)
        };

        amplt_a2 = amplt_a2_local;

        let (sfx, sfx_x1, sfx_x2) = match x_forced {
            Some(xf) if xf < x2 => {
                let sfx = (xf - x1) / dx;
                (sfx, (sfx - 1.0) / dx, -sfx / dx)
            }
            _ => (1.0, 0.0, 0.0),
        };

        if sfa < sfx {
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
        } else {
            wf2 = sfx;
            wf2_a1 = 0.0;
            wf2_a2 = 0.0;
            wf2_x1 = sfx_x1;
            wf2_x2 = sfx_x2;
        }
        
        wf1 = 1.0 - wf2;
        
        // Interpolate BL variables to XT
        xt = x1 * wf1 + x2 * wf2;
        tt = t1 * wf1 + t2 * wf2;
        dt = d1 * wf1 + d2 * wf2;
        ut = u1 * wf1 + u2 * wf2;
        
        // Derivatives w.r.t. A2 for iteration
        xt_a2 = x1 * (-wf2_a2) + x2 * wf2_a2;
        tt_a2 = t1 * (-wf2_a2) + t2 * wf2_a2;
        dt_a2 = d1 * (-wf2_a2) + d2 * wf2_a2;
        ut_a2 = u1 * (-wf2_a2) + u2 * wf2_a2;
        
        // Compute Hk at transition point
        let ht = if tt > 1e-20 { dt / tt } else { hk1 * wf1 + hk2 * wf2 };
        let hkt_result = hkin(ht, msq);
        let hkt = hkt_result.hk;
        
        // HKT derivatives for AX_A2 calculation
        let hkt_tt = hkt_result.hk_h * (-ht / tt.max(1e-20));
        let hkt_dt = hkt_result.hk_h * (1.0 / tt.max(1e-20));
        let hkt_ut = 0.0;
        let hkt_ms = hkt_result.hk_msq;
        
        // Rθ at transition point
        let rtt = reybl * ut * tt;
        let rtt_tt = reybl * ut;
        let rtt_ut = reybl * tt;
        let rtt_re = ut * tt;
        let rtt_ms = 0.0;
        
        // Compute averaged amplification rate with derivatives
        let amplt = amplt_local;
        let ax_result = axset_full(hk1, t1, rt1, ampl1, hkt, tt, rtt, amplt, ncrit);
        ax = ax_result.ax;
        
        if ax <= 0.0 {
            converged = true;
            break;
        }
        
        // AX_A2 = chain rule through HKT, TT, RTT, AMPLT
        let ax_a2 = (ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt) * tt_a2
                  + (ax_result.ax_hk2 * hkt_dt) * dt_a2
                  + (ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut) * ut_a2
                  + ax_result.ax_a2 * amplt_a2;
        
        // Residual: RES = AMPL2 - AMPL1 - AX*(X2-X1)
        let res = ampl2 - ampl1 - ax * dx;
        let res_a2 = 1.0 - ax_a2 * dx;
        
        let da2 = -res / res_a2.max(1e-20);
        
        // Relaxation
        let dxt = xt_a2 * da2;
        let mut rlx = 1.0;
        if (rlx * dxt.abs() / dx) > 0.05 {
            rlx = 0.05 * dx / dxt.abs().max(1e-20);
        }
        if rlx * da2.abs() > 1.0 {
            rlx = 1.0 / da2.abs();
        }
        
        // Check convergence
        if da2.abs() < DAEPS {
            converged = true;
            break;
        }
        
        // Prevent crossing Ncrit boundary
        let ampl2_new = ampl2 + rlx * da2;
        if (ampl2 > ncrit && ampl2_new < ncrit) || (ampl2 < ncrit && ampl2_new > ncrit) {
            ampl2 = ncrit;
        } else {
            ampl2 = ampl2_new;
        }
    }

    // Test for free or forced transition
    let tr_free = ampl2 >= ncrit;
    let mut tr_forced = x_forced.is_some_and(|xf| xf > x1 && xf <= x2);
    let transition = tr_free || tr_forced;
    
    if !transition {
        // No transition - return minimal result
        return Trchek2FullResult {
            ncrit,
            xt,
            ampl2,
            ax,
            transition: false,
            forced: false,
            converged,
            iterations,
            wf1,
            wf2,
            ..Default::default()
        };
    }

    // Match XFOIL's forced-transition exit: XT is prescribed, so XT sensitivities
    // are zero and the interpolated T/D/U terms are rebuilt from the fixed XT.
    if tr_free && tr_forced {
        if let Some(xf) = x_forced {
            tr_forced = xf < xt;
        }
    }

    if tr_forced {
        let xf = x_forced.expect("forced transition requires x_forced");
        let xt = xf;
        let wf2 = (xt - x1) / dx;
        let wf1 = 1.0 - wf2;

        let wf2_x1 = (wf2 - 1.0) / dx;
        let wf2_x2 = -wf2 / dx;
        let wf1_x1 = -wf2_x1;
        let wf1_x2 = -wf2_x2;

        let tt_x1 = t1 * wf1_x1 + t2 * wf2_x1;
        let tt_x2 = t1 * wf1_x2 + t2 * wf2_x2;
        let dt_x1 = d1 * wf1_x1 + d2 * wf2_x1;
        let dt_x2 = d1 * wf1_x2 + d2 * wf2_x2;
        let ut_x1 = u1 * wf1_x1 + u2 * wf2_x1;
        let ut_x2 = u1 * wf1_x2 + u2 * wf2_x2;

        return Trchek2FullResult {
            ncrit,
            xt,
            ampl2,
            ax,
            transition: true,
            forced: true,
            converged,
            iterations,
            wf1,
            wf2,
            xt_a1: 0.0,
            xt_x1: 0.0,
            xt_t1: 0.0,
            xt_d1: 0.0,
            xt_u1: 0.0,
            xt_x2: 0.0,
            xt_t2: 0.0,
            xt_d2: 0.0,
            xt_u2: 0.0,
            xt_ms: 0.0,
            xt_re: 0.0,
            tt_a1: 0.0,
            tt_x1,
            tt_x2,
            tt_t1: wf1,
            tt_t2: wf2,
            tt_d1: 0.0,
            tt_d2: 0.0,
            tt_u1: 0.0,
            tt_u2: 0.0,
            tt_ms: 0.0,
            tt_re: 0.0,
            dt_a1: 0.0,
            dt_x1,
            dt_x2,
            dt_t1: 0.0,
            dt_t2: 0.0,
            dt_d1: wf1,
            dt_d2: wf2,
            dt_u1: 0.0,
            dt_u2: 0.0,
            dt_ms: 0.0,
            dt_re: 0.0,
            ut_a1: 0.0,
            ut_x1,
            ut_x2,
            ut_t1: 0.0,
            ut_t2: 0.0,
            ut_d1: 0.0,
            ut_d2: 0.0,
            ut_u1: wf1,
            ut_u2: wf2,
            ut_ms: 0.0,
            ut_re: 0.0,
            ..Default::default()
        };
    }

    // === Compute full derivatives for free transition ===
    // At this point we have converged values and need to compute all sensitivities
    
    // Recompute interpolated variables with final wf1, wf2
    tt = t1 * wf1 + t2 * wf2;
    dt = d1 * wf1 + d2 * wf2;
    ut = u1 * wf1 + u2 * wf2;
    
    // Compute Hk and Rθ at transition with derivatives
    let ht = if tt > 1e-20 { dt / tt } else { hk1 * wf1 + hk2 * wf2 };
    let hkt_result = hkin(ht, msq);
    let hkt = hkt_result.hk;
    
    let hkt_tt = hkt_result.hk_h * (-ht / tt.max(1e-20));
    let hkt_dt = hkt_result.hk_h * (1.0 / tt.max(1e-20));
    let hkt_ut = 0.0;
    let hkt_ms = hkt_result.hk_msq;
    
    let rtt = reybl * ut * tt;
    let rtt_tt = reybl * ut;
    let rtt_ut = reybl * tt;
    let rtt_ms = 0.0;
    let rtt_re = ut * tt;
    
    let amplt = if ampl2 <= ncrit { ampl2 } else { ncrit };
    let ax_result = axset_full(hk1, t1, rt1, ampl1, hkt, tt, rtt, amplt, ncrit);
    
    // === Basic interpolation derivatives ===
    // XT = X1*WF1 + X2*WF2
    let xt_x1_base = wf1;
    let xt_x2_base = wf2;
    let tt_t1_base = wf1;
    let tt_t2_base = wf2;
    let dt_d1_base = wf1;
    let dt_d2_base = wf2;
    let ut_u1_base = wf1;
    let ut_u2_base = wf2;
    
    // Derivatives through WF (via A1)
    let wf1_a1 = -wf2_a1;
    let xt_a1_via_wf = x1 * wf1_a1 + x2 * wf2_a1;
    let tt_a1_via_wf = t1 * wf1_a1 + t2 * wf2_a1;
    let dt_a1_via_wf = d1 * wf1_a1 + d2 * wf2_a1;
    let ut_a1_via_wf = u1 * wf1_a1 + u2 * wf2_a1;
    
    // Derivatives through WF (via X1, X2 for forced transition - all zero for free)
    let wf1_x1 = -wf2_x1;
    let wf1_x2 = -wf2_x2;
    let xt_x1_via_wf = x1 * wf1_x1 + x2 * wf2_x1;
    let xt_x2_via_wf = x1 * wf1_x2 + x2 * wf2_x2;
    let tt_x1 = t1 * wf1_x1 + t2 * wf2_x1;
    let tt_x2 = t1 * wf1_x2 + t2 * wf2_x2;
    let dt_x1 = d1 * wf1_x1 + d2 * wf2_x1;
    let dt_x2 = d1 * wf1_x2 + d2 * wf2_x2;
    let ut_x1 = u1 * wf1_x1 + u2 * wf2_x1;
    let ut_x2 = u1 * wf1_x2 + u2 * wf2_x2;
    
    // === Full AX sensitivities w.r.t. T1, D1, U1, A1, T2, D2, U2, X1, X2, MS, RE ===
    // AX = AX(HK1, T1, RT1, A1, HKT, TT, RTT, AT)
    // Need to chain through HKT(TT,DT), RTT(TT,UT), and interpolation
    
    let ax_tt = ax_result.ax_hk2 * hkt_tt + ax_result.ax_t2 + ax_result.ax_rt2 * rtt_tt;
    let ax_dt = ax_result.ax_hk2 * hkt_dt;
    let ax_ut = ax_result.ax_hk2 * hkt_ut + ax_result.ax_rt2 * rtt_ut;
    let ax_at = ax_result.ax_a2;
    
    // AX sensitivities w.r.t. original variables
    // Station 1 direct derivatives
    let ax_t1 = ax_result.ax_hk1 * hk1_t1 + ax_result.ax_t1 + ax_result.ax_rt1 * rt1_t1
              + ax_tt * tt_t1_base;
    let ax_d1 = ax_result.ax_hk1 * hk1_d1 + ax_dt * dt_d1_base;
    let ax_u1 = ax_result.ax_rt1 * rt1_u1 + ax_ut * ut_u1_base;
    let ax_a1 = ax_result.ax_a1
              + ax_tt * tt_a1_via_wf + ax_dt * dt_a1_via_wf + ax_ut * ut_a1_via_wf;
    let ax_x1 = ax_tt * tt_x1 + ax_dt * dt_x1 + ax_ut * ut_x1;
    
    // Station 2 derivatives  
    let ax_t2 = ax_tt * tt_t2_base;
    let ax_d2 = ax_dt * dt_d2_base;
    let ax_u2 = ax_ut * ut_u2_base;
    let ax_a2 = ax_tt * tt_a2 + ax_dt * dt_a2 + ax_ut * ut_a2 + ax_at * amplt_a2;
    let ax_x2 = ax_tt * tt_x2 + ax_dt * dt_x2 + ax_ut * ut_x2;
    
    // Global parameters
    let ax_ms = ax_result.ax_hk1 * hk1_ms + ax_result.ax_hk2 * hkt_ms
              + ax_result.ax_rt1 * rt1_ms + ax_result.ax_rt2 * rtt_ms;
    let ax_re = ax_result.ax_rt1 * rt1_re + ax_result.ax_rt2 * rtt_re;
    
    // === Residual sensitivities Z_* ===
    // RES = AMPL2 - AMPL1 - AX*(X2-X1)
    let z_ax = -dx;
    
    let z_a1 = z_ax * ax_a1 - 1.0;
    let z_t1 = z_ax * ax_t1;
    let z_d1 = z_ax * ax_d1;
    let z_u1 = z_ax * ax_u1;
    let z_x1 = z_ax * ax_x1 + ax;
    
    let z_a2 = z_ax * ax_a2 + 1.0;
    let z_t2 = z_ax * ax_t2;
    let z_d2 = z_ax * ax_d2;
    let z_u2 = z_ax * ax_u2;
    let z_x2 = z_ax * ax_x2 - ax;
    
    let z_ms = z_ax * ax_ms;
    let z_re = z_ax * ax_re;
    
    // === Final XT sensitivities using implicit function theorem ===
    // XT_* = XT_*_basic - (XT_A2/Z_A2)*Z_*
    // This ensures that when we change a variable, RES remains = 0
    
    let xt_a2_over_z_a2 = if z_a2.abs() > 1e-20 { xt_a2 / z_a2 } else { 0.0 };
    
    let xt_a1 = xt_a1_via_wf - xt_a2_over_z_a2 * z_a1;
    let xt_t1 = -xt_a2_over_z_a2 * z_t1;
    let xt_d1 = -xt_a2_over_z_a2 * z_d1;
    let xt_u1 = -xt_a2_over_z_a2 * z_u1;
    let xt_x1 = xt_x1_base + xt_x1_via_wf - xt_a2_over_z_a2 * z_x1;
    
    let xt_t2 = -xt_a2_over_z_a2 * z_t2;
    let xt_d2 = -xt_a2_over_z_a2 * z_d2;
    let xt_u2 = -xt_a2_over_z_a2 * z_u2;
    let xt_x2 = xt_x2_base + xt_x2_via_wf - xt_a2_over_z_a2 * z_x2;
    
    let xt_ms = -xt_a2_over_z_a2 * z_ms;
    let xt_re = -xt_a2_over_z_a2 * z_re;
    
    // === Compute WF derivatives from XT ===
    // WF2 = (XT - X1)/(X2 - X1) = XT/dx - X1/dx
    // WF2_XT = 1/dx
    let wf2_xt = 1.0 / dx;
    
    // Full WF2 derivatives
    let wf2_a1_full = wf2_xt * xt_a1;
    let wf2_x1_full = wf2_xt * xt_x1 + (wf2 - 1.0) / dx;
    let wf2_x2_full = wf2_xt * xt_x2 - wf2 / dx;
    let wf2_t1 = wf2_xt * xt_t1;
    let wf2_t2 = wf2_xt * xt_t2;
    let wf2_d1 = wf2_xt * xt_d1;
    let wf2_d2 = wf2_xt * xt_d2;
    let wf2_u1 = wf2_xt * xt_u1;
    let wf2_u2 = wf2_xt * xt_u2;
    let wf2_ms = wf2_xt * xt_ms;
    let wf2_re = wf2_xt * xt_re;
    
    let wf1_a1_full = -wf2_a1_full;
    let wf1_x1_full = -wf2_x1_full;
    let wf1_x2_full = -wf2_x2_full;
    let wf1_t1 = -wf2_t1;
    let wf1_t2 = -wf2_t2;
    let wf1_d1 = -wf2_d1;
    let wf1_d2 = -wf2_d2;
    let wf1_u1 = -wf2_u1;
    let wf1_u2 = -wf2_u2;
    let wf1_ms = -wf2_ms;
    let wf1_re = -wf2_re;
    
    // === Full interpolation derivatives for TT, DT, UT ===
    // TT = T1*WF1 + T2*WF2
    let tt_a1 = t1 * wf1_a1_full + t2 * wf2_a1_full;
    let tt_x1 = t1 * wf1_x1_full + t2 * wf2_x1_full;
    let tt_x2 = t1 * wf1_x2_full + t2 * wf2_x2_full;
    let tt_t1 = t1 * wf1_t1 + t2 * wf2_t1 + wf1;
    let tt_t2 = t1 * wf1_t2 + t2 * wf2_t2 + wf2;
    let tt_d1 = t1 * wf1_d1 + t2 * wf2_d1;
    let tt_d2 = t1 * wf1_d2 + t2 * wf2_d2;
    let tt_u1 = t1 * wf1_u1 + t2 * wf2_u1;
    let tt_u2 = t1 * wf1_u2 + t2 * wf2_u2;
    let tt_ms = t1 * wf1_ms + t2 * wf2_ms;
    let tt_re = t1 * wf1_re + t2 * wf2_re;
    
    // DT = D1*WF1 + D2*WF2
    let dt_a1 = d1 * wf1_a1_full + d2 * wf2_a1_full;
    let dt_x1 = d1 * wf1_x1_full + d2 * wf2_x1_full;
    let dt_x2 = d1 * wf1_x2_full + d2 * wf2_x2_full;
    let dt_t1 = d1 * wf1_t1 + d2 * wf2_t1;
    let dt_t2 = d1 * wf1_t2 + d2 * wf2_t2;
    let dt_d1 = d1 * wf1_d1 + d2 * wf2_d1 + wf1;
    let dt_d2 = d1 * wf1_d2 + d2 * wf2_d2 + wf2;
    let dt_u1 = d1 * wf1_u1 + d2 * wf2_u1;
    let dt_u2 = d1 * wf1_u2 + d2 * wf2_u2;
    let dt_ms = d1 * wf1_ms + d2 * wf2_ms;
    let dt_re = d1 * wf1_re + d2 * wf2_re;
    
    // UT = U1*WF1 + U2*WF2
    let ut_a1 = u1 * wf1_a1_full + u2 * wf2_a1_full;
    let ut_x1 = u1 * wf1_x1_full + u2 * wf2_x1_full;
    let ut_x2 = u1 * wf1_x2_full + u2 * wf2_x2_full;
    let ut_t1 = u1 * wf1_t1 + u2 * wf2_t1;
    let ut_t2 = u1 * wf1_t2 + u2 * wf2_t2;
    let ut_d1 = u1 * wf1_d1 + u2 * wf2_d1;
    let ut_d2 = u1 * wf1_d2 + u2 * wf2_d2;
    let ut_u1 = u1 * wf1_u1 + u2 * wf2_u1 + wf1;
    let ut_u2 = u1 * wf1_u2 + u2 * wf2_u2 + wf2;
    let ut_ms = u1 * wf1_ms + u2 * wf2_ms;
    let ut_re = u1 * wf1_re + u2 * wf2_re;

    Trchek2FullResult {
        ncrit,
        xt,
        ampl2,
        ax,
        transition: true,
        forced: false,
        converged,
        iterations,
        wf1,
        wf2,
        
        xt_a1,
        xt_x1,
        xt_t1,
        xt_d1,
        xt_u1,
        xt_x2,
        xt_t2,
        xt_d2,
        xt_u2,
        xt_ms,
        xt_re,
        
        tt_a1,
        tt_x1,
        tt_x2,
        tt_t1,
        tt_t2,
        tt_d1,
        tt_d2,
        tt_u1,
        tt_u2,
        tt_ms,
        tt_re,
        
        dt_a1,
        dt_x1,
        dt_x2,
        dt_t1,
        dt_t2,
        dt_d1,
        dt_d2,
        dt_u1,
        dt_u2,
        dt_ms,
        dt_re,
        
        ut_a1,
        ut_x1,
        ut_x2,
        ut_t1,
        ut_t2,
        ut_d1,
        ut_d2,
        ut_u1,
        ut_u2,
        ut_ms,
        ut_re,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-7;
    const DERIV_TOL: f64 = 1e-5;

    // ========== AXSET Tests ==========

    #[test]
    fn test_axset_zero_amplification() {
        // When both stations are subcritical, all amplification should be ~0
        // (except for the small DAX term near Ncrit)
        let result = axset(2.5, 0.001, 50.0, 0.0, 2.5, 0.001, 60.0, 0.0, 9.0);
        
        assert_eq!(result.ax1, 0.0, "AX1 should be 0 below critical");
        assert_eq!(result.ax2, 0.0, "AX2 should be 0 below critical");
        assert_eq!(result.axa_rms, 0.0, "RMS should be 0 when both are 0");
        assert!(result.dax > 0.0, "DAX should be positive (near-Ncrit term)");
        assert!(result.ax > 0.0, "Final AX should be positive due to DAX");
    }

    #[test]
    fn test_axset_supercritical() {
        // When both stations are supercritical, RMS should be positive
        let result = axset(
            2.6, 0.001, 1000.0, 3.0,  // station 1: supercritical
            2.7, 0.001, 1200.0, 4.0,  // station 2: supercritical
            9.0,                       // Ncrit
        );
        
        assert!(result.ax1 > 0.0, "AX1 should be positive");
        assert!(result.ax2 > 0.0, "AX2 should be positive");
        assert!(result.axa_rms > 0.0, "RMS average should be positive");
        
        // RMS should be between min and max of individual values
        let min_ax = result.ax1.min(result.ax2);
        let max_ax = result.ax1.max(result.ax2);
        assert!(result.axa_rms >= min_ax * 0.9, "RMS >= 0.9*min");
        assert!(result.axa_rms <= max_ax * 1.1, "RMS <= 1.1*max");
    }

    #[test]
    fn test_axset_rms_averaging() {
        // Verify RMS formula: sqrt(0.5*(a² + b²))
        let result = axset(
            2.6, 0.001, 2000.0, 5.0,  // station 1
            2.6, 0.001, 2000.0, 5.0,  // station 2 (same)
            9.0,
        );
        
        // When both stations are identical, RMS = individual value
        assert!(
            (result.axa_rms - result.ax1).abs() < 1e-10,
            "When stations are identical, RMS should equal individual: {} vs {}",
            result.axa_rms,
            result.ax1
        );
    }

    #[test]
    fn test_axset_dax_near_ncrit() {
        // DAX term should be larger when close to Ncrit
        let result_close = axset(
            2.6, 0.001, 1000.0, 8.5,  // N close to Ncrit=9
            2.6, 0.001, 1000.0, 8.5,
            9.0,
        );
        
        let result_far = axset(
            2.6, 0.001, 1000.0, 2.0,  // N far from Ncrit=9
            2.6, 0.001, 1000.0, 2.0,
            9.0,
        );
        
        assert!(
            result_close.dax > result_far.dax,
            "DAX should be larger near Ncrit: {} vs {}",
            result_close.dax,
            result_far.dax
        );
    }

    #[test]
    fn test_axset_dax_at_ncrit() {
        // When N >= Ncrit, DAX should be at maximum (exn = 1)
        let result = axset(
            2.6, 0.001, 1000.0, 10.0,  // N > Ncrit
            2.6, 0.001, 1000.0, 10.0,
            9.0,
        );
        
        // DAX = 0.002 / (th1 + th2) when exn = 1
        let expected_dax = 0.002 / (0.001 + 0.001);
        assert!(
            (result.dax - expected_dax).abs() < 1e-10,
            "DAX at N>=Ncrit: {} vs {}",
            result.dax,
            expected_dax
        );
    }

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
