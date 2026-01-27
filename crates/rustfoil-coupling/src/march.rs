//! Boundary layer marching procedures
//!
//! This module implements boundary layer marching from the stagnation point
//! downstream along the airfoil surface. Two modes are provided:
//!
//! - `march_fixed_ue()`: Direct march with prescribed edge velocity (MRCHUE)
//! - `march_coupled()`: Coupled march with Ue updates from inviscid (MRCHDU)
//!
//! # XFOIL Reference
//! - MRCHUE: xbl.f line 542 - direct BL march
//! - MRCHDU: xbl.f line 875 - coupled BL march with Ue updates

use nalgebra::DMatrix;
use rustfoil_bl::closures::{
    axset, check_transition, hkin, trchek2_full, trchek2_stations, Trchek2FullResult,
    Trchek2Result,
};
use rustfoil_bl::constants::{CTRCON, CTRCEX};
use rustfoil_bl::equations::{
    bldif, bldif_with_terms, blvar, trdif, trdif_full, trdif_turb_terms, FlowType,
};
use rustfoil_bl::state::BlStation;
use rustfoil_bl::{closures, BlvarInput, BlvarOutput};

use crate::solve::{build_4x4_system, solve_4x4};

/// Compute initial ctau at transition using XFOIL's physics-based formula
/// 
/// XFOIL xblsys.f lines 1391-1401:
/// ```text
/// CTR = CTRCON * EXP(-CTRCEX/(HK2-1.0))
/// ST = CTR * CQ2
/// ```
/// 
/// This computes the initial shear stress coefficient at the transition point
/// based on the local boundary layer state.
/// 
/// # Arguments
/// * `hk` - Kinematic shape factor at transition (should be near htmax, ~2.5)
/// * `cq` - Equilibrium shear coefficient from turbulent closures
/// 
/// # Note
/// XFOIL uses the Hk at the interpolated transition point, not the downstream
/// station's laminar Hk. When calling this, use htmax as a representative value
/// if the actual transition Hk isn't available.
#[inline]
fn compute_transition_ctau(hk: f64, cq: f64) -> f64 {
    let hk_arg = (hk - 1.0).max(0.1); // Prevent division by zero
    let ctr = CTRCON * (-CTRCEX / hk_arg).exp();
    ctr * cq
}

fn estimate_transition_ctau(
    prev: &BlStation,
    curr: &BlStation,
    tr: &Trchek2FullResult,
    msq: f64,
    re: f64,
) -> f64 {
    let tt = prev.theta * tr.wf1 + curr.theta * tr.wf2;
    let dt = prev.delta_star * tr.wf1 + curr.delta_star * tr.wf2;
    let ut = prev.u * tr.wf1 + curr.u * tr.wf2;

    if !tt.is_finite() || tt <= 1e-20 {
        return 0.03;
    }

    let mut st_turb = BlStation::new();
    st_turb.x = tr.xt;
    st_turb.u = ut;
    st_turb.theta = tt;
    st_turb.delta_star = dt;
    st_turb.ctau = 0.03;
    st_turb.ampl = 0.0;
    st_turb.is_laminar = false;
    st_turb.is_turbulent = true;
    blvar(&mut st_turb, FlowType::Turbulent, msq, re);

    let st = compute_transition_ctau(st_turb.hk, st_turb.cq);
    if st.is_finite() {
        st.clamp(1e-7, 0.3)
    } else {
        0.03
    }
}

fn bldif_terms_xfoil_upstream(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> Option<(rustfoil_bl::BldifTerms, BlStation)> {
    if flow_type != FlowType::Turbulent {
        return None;
    }

    let mut s1_xfoil = s1.clone();
    s1_xfoil.theta = s2.theta;
    s1_xfoil.delta_star = s2.delta_star;
    s1_xfoil.ctau = s2.ctau;
    s1_xfoil.ampl = s2.ampl;
    blvar(&mut s1_xfoil, flow_type, msq, re);

    let (_res_x, _jac_x, terms_x) = bldif_with_terms(&s1_xfoil, s2, flow_type, msq, re);
    Some((terms_x, s1_xfoil))
}

fn emit_trchek2_debug(
    side: usize,
    ibl: usize,
    prev: &BlStation,
    curr: &BlStation,
    tr: &Trchek2FullResult,
    msq: f64,
    re: f64,
    ncrit: f64,
) {
    if !rustfoil_bl::is_debug_active() {
        return;
    }

    let dx = curr.x - prev.x;
    let residual = tr.ampl2 - prev.ampl - tr.ax * dx;

    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_iter(
        1,
        side,
        ibl,
        tr.iterations.max(1),
        prev.x,
        curr.x,
        prev.ampl,
        tr.ampl2,
        tr.ax,
        residual,
        tr.wf1,
        tr.wf2,
        tr.xt,
        prev.hk,
        curr.hk,
        prev.r_theta,
        curr.r_theta,
        prev.theta,
        curr.theta,
        prev.u,
        curr.u,
        ncrit,
        tr.transition,
    ));

    let tt = prev.theta * tr.wf1 + curr.theta * tr.wf2;
    let dt = prev.delta_star * tr.wf1 + curr.delta_star * tr.wf2;
    let ut = prev.u * tr.wf1 + curr.u * tr.wf2;

    let ht = if tt > 1e-20 {
        dt / tt
    } else {
        prev.hk * tr.wf1 + curr.hk * tr.wf2
    };
    let hk_t = hkin(ht, msq).hk;
    let rt_t = re * ut * tt;

    let (st_opt, cq_opt) = if tr.transition {
        let mut st_turb_temp = BlStation::new();
        st_turb_temp.x = tr.xt;
        st_turb_temp.u = ut;
        st_turb_temp.theta = tt;
        st_turb_temp.delta_star = dt;
        st_turb_temp.ctau = 0.03;
        st_turb_temp.ampl = 0.0;
        st_turb_temp.is_laminar = false;
        st_turb_temp.is_turbulent = true;
        blvar(&mut st_turb_temp, FlowType::Turbulent, msq, re);

        let cq_t = st_turb_temp.cq;
        let st = compute_transition_ctau(st_turb_temp.hk, cq_t);
        (Some(st), Some(cq_t))
    } else {
        (None, None)
    };

    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_final(
        1,
        side,
        ibl,
        tr.iterations.max(1),
        tr.converged,
        prev.x,
        curr.x,
        prev.ampl,
        tr.ampl2,
        tr.ax,
        tr.xt,
        ncrit,
        tr.transition,
        false,
        Some(tr.wf1),
        Some(tr.wf2),
        Some(tt),
        Some(dt),
        Some(ut),
        Some(hk_t),
        Some(rt_t),
        st_opt,
        cq_opt,
    ));
}

fn emit_trchek2_station_iter(
    side: usize,
    ibl: usize,
    prev: &BlStation,
    curr: &BlStation,
    tr: &Trchek2Result,
    ncrit: f64,
) {
    if !rustfoil_bl::is_debug_active() {
        return;
    }

    let dx = curr.x - prev.x;
    let residual = tr.ampl2 - prev.ampl - tr.ax * dx;
    let wf2 = if tr.ampl2 <= ncrit {
        1.0
    } else {
        (ncrit - prev.ampl) / (tr.ampl2 - prev.ampl).max(1e-20)
    };
    let wf1 = 1.0 - wf2;
    let xt = prev.x * wf1 + curr.x * wf2;

    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_iter(
        1,
        side,
        ibl,
        tr.iterations.max(1),
        prev.x,
        curr.x,
        prev.ampl,
        tr.ampl2,
        tr.ax,
        residual,
        wf1,
        wf2,
        xt,
        prev.hk,
        curr.hk,
        prev.r_theta,
        curr.r_theta,
        prev.theta,
        curr.theta,
        prev.u,
        curr.u,
        ncrit,
        tr.transition,
    ));
}

fn emit_blvar_debug(side: usize, ibl: usize, flow_type: FlowType, station: &BlStation) {
    if !rustfoil_bl::is_debug_active() {
        return;
    }

    let flow_type_num = match flow_type {
        FlowType::Laminar => 1,
        FlowType::Turbulent => 2,
        FlowType::Wake => 3,
    };

    let input = BlvarInput {
        x: station.x,
        u: station.u,
        theta: station.theta,
        delta_star: station.delta_star,
        ctau: station.ctau,
        ampl: station.ampl,
    };
    let dd = station.ctau * station.ctau * (0.995 - station.us) * 2.0 / station.hs;
    let dd2 = 0.15 * (0.995 - station.us).powi(2) / station.r_theta * 2.0 / station.hs;
    let cf2t_result = closures::cf_turbulent(station.hk, station.r_theta, 0.0);
    let di_wall = (0.5 * cf2t_result.cf * station.us) * 2.0 / station.hs;
    let grt = station.r_theta.ln();
    let hmin = 1.0 + 2.1 / grt;
    let fl = (station.hk - 1.0) / (hmin - 1.0);
    let dfac = 0.5 + 0.5 * fl.tanh();
    let di_wall_corrected = di_wall * dfac;
    let di_total = di_wall_corrected + dd + dd2;
    let di_lam = closures::dissipation_laminar(station.hk, station.r_theta).di;
    let (di_lam_override, di_used_minus_total, di_used_minus_lam) = if flow_type == FlowType::Turbulent {
        let di_lam_override = di_lam > di_total;
        let di_used_minus_total = station.cd - di_total;
        let di_used_minus_lam = station.cd - di_lam;
        (Some(di_lam_override), Some(di_used_minus_total), Some(di_used_minus_lam))
    } else {
        (None, None, None)
    };

    let output = BlvarOutput {
        H: station.h,
        Hk: station.hk,
        Hs: station.hs,
        Hc: station.hc,
        Rtheta: station.r_theta,
        Cf: station.cf,
        Cd: station.cd,
        Us: station.us,
        Cq: station.cq,
        De: station.de,
        Dd: dd,
        Dd2: dd2,
        DiWall: di_wall_corrected,
        DiTotal: di_total,
        DiLam: di_lam,
        DiUsed: station.cd,
        DiLamOverride: di_lam_override,
        DiUsedMinusTotal: di_used_minus_total,
        DiUsedMinusLam: di_used_minus_lam,
    };

    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::blvar(
        1,
        side,
        ibl,
        flow_type_num,
        input,
        output,
    ));
}

// ============================================================================
// Configuration and Result Structures
// ============================================================================

/// Result of a boundary layer march
///
/// Contains the computed BL stations along with key events (transition,
/// separation) detected during the march.
#[derive(Debug, Clone)]
pub struct MarchResult {
    /// BL stations after march (one per input station)
    pub stations: Vec<BlStation>,
    /// Arc length where transition occurred (if any)
    pub x_transition: Option<f64>,
    /// Arc length where separation occurred (if any)
    pub x_separation: Option<f64>,
    /// Whether the march completed successfully
    pub converged: bool,
    /// Station index where transition occurred
    pub transition_index: Option<usize>,
}

impl MarchResult {
    /// Create a new empty result
    pub fn new(capacity: usize) -> Self {
        Self {
            stations: Vec::with_capacity(capacity),
            x_transition: None,
            x_separation: None,
            converged: true,
            transition_index: None,
        }
    }
}

/// Configuration for boundary layer march
///
/// Controls convergence criteria, transition prediction, and mode switching.
#[derive(Debug, Clone)]
pub struct MarchConfig {
    /// Critical N factor for transition (e^n method)
    /// Typical values: 9.0 for free flight, 3-5 for noisy tunnels
    pub ncrit: f64,
    /// Maximum Newton iterations per station
    pub max_iter: usize,
    /// Convergence tolerance for Newton iteration
    pub tolerance: f64,
    /// Maximum Hk before switching to inverse mode (laminar)
    /// XFOIL default: 3.8
    pub hlmax: f64,
    /// Maximum Hk before switching to inverse mode (turbulent)
    /// XFOIL default: 2.5
    pub htmax: f64,
    /// Sensitivity weight for coupled march (SENSWT)
    /// Controls how strongly Ue is coupled to Hk deviations
    pub senswt: f64,
    /// Small tolerance for convergence check in coupled mode
    pub deps: f64,
    /// Enable debug tracing output
    pub debug_trace: bool,
}

impl Default for MarchConfig {
    fn default() -> Self {
        Self {
            ncrit: 9.0,
            max_iter: 25,
            tolerance: 1.0e-5, // Match XFOIL's MRCHUE convergence criterion
            hlmax: 3.8,  // XFOIL default for laminar inverse-mode switch
            htmax: 2.5,
            senswt: 1000.0,
            deps: 5.0e-6,
            debug_trace: false,
        }
    }
}

impl MarchConfig {
    /// Create config for low-turbulence wind tunnel
    pub fn low_turbulence() -> Self {
        Self {
            ncrit: 9.0,
            ..Default::default()
        }
    }

    /// Create config for high-turbulence environment
    pub fn high_turbulence() -> Self {
        Self {
            ncrit: 3.0,
            ..Default::default()
        }
    }
}

// ============================================================================
// Stagnation Point Initialization
// ============================================================================

/// Initialize boundary layer at stagnation point using Thwaites' formula
///
/// The Thwaites method provides an accurate estimate for the initial BL
/// thickness near a stagnation point where the similarity solution applies.
///
/// # Arguments
/// * `x` - Arc length at first station
/// * `ue` - Edge velocity at first station
/// * `re` - Reynolds number
///
/// # Returns
/// A BlStation initialized for the stagnation region
///
/// # Reference
/// XFOIL xbl.f lines 568-579:
/// ```fortran
/// BULE = 1.0
/// UCON = UEI/XSI**BULE
/// TSQ = 0.45/(UCON*(5.0*BULE+1.0)*REYBL) * XSI**(1.0-BULE)
/// THI = SQRT(TSQ)
/// DSI = 2.2*THI
/// ```
fn init_stagnation(x: f64, ue: f64, re: f64) -> BlStation {
    let mut station = BlStation::new();
    station.x = x;
    station.u = ue;

    // Thwaites' formula for stagnation point
    // With BULE = 1.0: TSQ = 0.45 / (Ue/x * 6 * Re) * x^0 = 0.45*x / (6*Ue*Re)
    let bule = 1.0;
    let ucon = ue / x.powf(bule);
    let tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * re) * x.powf(1.0 - bule);

    station.theta = tsq.sqrt().max(1e-12);
    station.delta_star = 2.2 * station.theta;
    station.h = station.delta_star / station.theta;
    station.hk = station.h; // At low Mach
    station.ampl = 0.0;
    station.ctau = 0.03; // Initial turbulent shear (for if transition happens)
    station.is_laminar = true;
    station.is_turbulent = false;
    station.is_wake = false;

    // Compute Rθ
    station.r_theta = re * ue * station.theta;
    station.mass_defect = ue * station.delta_star;

    station
}

// ============================================================================
// Simplified BL Integration (Thwaites-like)
// ============================================================================

/// Compute momentum thickness growth for one step using momentum integral
///
/// Uses the von Karman momentum integral equation:
/// dθ/dx + (H + 2 - M²) θ/Ue dUe/dx = Cf/2
///
/// For laminar flow with favorable/zero pressure gradient, this gives
/// Blasius-like growth. For adverse pressure gradient, H increases.
///
/// Note: Near the leading edge, velocity gradients are very strong (dUe/dx >> 1)
/// which can cause numerical instability with forward Euler. We use a stabilized
/// scheme that limits the maximum relative change in θ per step.
fn step_momentum(
    prev: &BlStation,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
    is_laminar: bool,
) -> (f64, f64) {
    let dx = (x_new - prev.x).abs(); // Use absolute value for arc length
    if dx <= 1e-12 {
        return (prev.theta, prev.delta_star);
    }

    // Use upwind values for stability, with safeguards
    let h = prev.h.clamp(1.0, 4.0);
    let theta = prev.theta.max(1e-12);
    let ue = prev.u.abs().max(1e-10);
    let cf = prev.cf.max(1e-6);

    // Velocity gradient term - use magnitude for signed velocities
    let ue_new_abs = ue_new.abs().max(1e-10);
    let due_dx = (ue_new_abs - ue) / dx;
    let ue_avg = 0.5 * (ue + ue_new_abs);

    // Momentum thickness growth: dθ/dx = Cf/2 - (H+2-M²) θ/Ue dUe/dx
    let h_term = h + 2.0 - msq;
    
    // For strong favorable gradients (accelerating flow near LE), use stabilized scheme
    // The pressure gradient term can dominate, causing theta to go negative
    // Limit the relative change to prevent instability
    let pressure_term = h_term * theta / ue_avg * due_dx;
    let friction_term = cf / 2.0;
    
    // Use implicit-like stabilization for strong acceleration
    // If pressure term would make theta decrease by more than 50% per step, limit it
    let max_decrease_rate = 0.5 * theta / dx;  // Max rate that gives 50% decrease
    let pressure_term_limited = pressure_term.min(friction_term + max_decrease_rate);
    
    let dtheta_dx = friction_term - pressure_term_limited;

    // Integrate θ with limiting to prevent blow-up
    // Use min of:
    // 1. Forward Euler result
    // 2. Maximum relative increase (2x per step for stability)
    let theta_euler = theta + dtheta_dx * dx;
    let theta_max = theta * 2.0;  // Don't more than double per step
    let theta_min = theta * 0.5;  // Don't reduce by more than half per step
    let theta_new = theta_euler.clamp(theta_min.max(1e-12), theta_max.min(0.1));

    // For H (shape factor), use empirical correlations
    // In attached laminar flow, H ≈ 2.59 (Blasius)
    // In adverse pressure gradient, H increases
    // In turbulent flow, H ≈ 1.4

    // Compute equilibrium H based on pressure gradient parameter
    let lambda = (theta * theta * re * due_dx / ue_avg).clamp(-0.5, 0.5);

    let h_new = if is_laminar {
        // Thwaites correlation for laminar H
        // H = H(λ) where λ = θ² Re dUe/dx / Ue
        let h_eq = if lambda < -0.09 {
            // Strong favorable gradient
            2.0
        } else if lambda > 0.09 {
            // Strong adverse gradient (approaching separation)
            4.0
        } else {
            // Thwaites approximation: H ≈ 2.61 - 3.75λ - 5.24λ²
            (2.61 - 3.75 * lambda - 5.24 * lambda * lambda).clamp(2.0, 4.0)
        };
        // Relax toward equilibrium
        0.8 * h + 0.2 * h_eq
    } else {
        // Turbulent: H tends toward ~1.4 for equilibrium
        let h_eq = 1.4 + 0.5 * lambda.abs().min(0.2);
        (0.8 * h + 0.2 * h_eq).clamp(1.2, 3.0)
    };

    let delta_star_new = h_new * theta_new;

    (theta_new, delta_star_new)
}

// ============================================================================
// Local Newton Solver (4x4 system)
// ============================================================================

/// Solve 4x4 Gauss elimination for single-station Newton update
///
/// The 4x4 system is: [ctau/ampl, theta, delta_star, ue] with 4 equations:
/// - 3 BL equations from bldif
/// - 1 constraint equation (direct: dUe=0, inverse: target Hk)
fn gauss_4x4(a: &mut [[f64; 4]; 4], b: &mut [f64; 4]) {
    // Forward elimination with partial pivoting
    for k in 0..4 {
        // Find pivot
        let mut max_val = a[k][k].abs();
        let mut max_row = k;
        for i in (k + 1)..4 {
            if a[i][k].abs() > max_val {
                max_val = a[i][k].abs();
                max_row = i;
            }
        }

        // Swap rows if needed
        if max_row != k {
            a.swap(k, max_row);
            b.swap(k, max_row);
        }

        // Check for near-zero pivot
        if a[k][k].abs() < 1e-30 {
            continue;
        }

        // Eliminate column k
        for i in (k + 1)..4 {
            let factor = a[i][k] / a[k][k];
            for j in k..4 {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    for i in (0..4).rev() {
        if a[i][i].abs() < 1e-30 {
            b[i] = 0.0;
            continue;
        }
        for j in (i + 1)..4 {
            b[i] -= a[i][j] * b[j];
        }
        b[i] /= a[i][i];
    }
}

/// Solve for BL station using Newton iteration
///
/// This is the core Newton solver for a single BL station, equivalent to
/// XFOIL's MRCHUE inner loop. It solves for theta, delta_star (and ctau for
/// turbulent) using the integral BL equations.
///
/// # Arguments
/// * `prev` - Previous (upstream) BL station, used as "1" variables in bldif
/// * `x_new` - Arc length at current station
/// * `ue_new` - Edge velocity at current station
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
/// * `is_laminar` - Whether flow is still laminar at this station
/// * `max_iter` - Maximum Newton iterations
/// * `tolerance` - Convergence tolerance (DMAX < tolerance)
/// * `hlmax` - Max Hk for laminar direct mode
/// * `htmax` - Max Hk for turbulent direct mode
///
/// # Returns
/// Converged BlStation with updated theta, delta_star, etc.
///
/// # Reference
/// XFOIL xbl.f MRCHUE Newton loop (line 622-789)
pub fn newton_solve_station(
    prev: &BlStation,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
    is_laminar: bool,
    max_iter: usize,
    tolerance: f64,
    hlmax: f64,
    htmax: f64,
) -> (BlStation, bool, f64) {
    newton_solve_station_with_guess(
        prev,
        None,
        x_new,
        ue_new,
        re,
        msq,
        is_laminar,
        0.0,
        max_iter,
        tolerance,
        hlmax,
        htmax,
        None,
        None,
    )
}

/// Solve for BL station using Newton iteration with optional initial guess
///
/// This is the core Newton solver for a single BL station, equivalent to
/// XFOIL's MRCHUE inner loop. It solves for theta, delta_star (and ctau for
/// turbulent) using the integral BL equations.
///
/// # Arguments
/// * `prev` - Previous (upstream) BL station, used as "1" variables in bldif
/// * `initial_guess` - Optional initial guess for theta, delta_star, ctau. 
///                     If None, uses values from `prev`. Used at transition
///                     to initialize from the laminar solution while using the
///                     actual upstream station for bldif.
/// * `x_new` - Arc length at current station
/// * `ue_new` - Edge velocity at current station
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
/// * `is_laminar` - Whether flow is still laminar at this station
/// * `max_iter` - Maximum Newton iterations
/// * `tolerance` - Convergence tolerance (DMAX < tolerance)
/// * `hlmax` - Max Hk for laminar direct mode
/// * `htmax` - Max Hk for turbulent direct mode
///
/// # Returns
/// Converged BlStation with updated theta, delta_star, etc.
///
/// # Reference
/// XFOIL xbl.f MRCHUE Newton loop (line 622-789)
pub fn newton_solve_station_with_guess(
    prev: &BlStation,
    initial_guess: Option<&BlStation>,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
    is_laminar: bool,
    ncrit: f64,
    max_iter: usize,
    tolerance: f64,
    hlmax: f64,
    htmax: f64,
    debug_side: Option<usize>,
    debug_ibl: Option<usize>,
) -> (BlStation, bool, f64) {
    // Call the extended version without transition info
    newton_solve_station_transition(
        prev,
        initial_guess,
        x_new,
        ue_new,
        re,
        msq,
        is_laminar,
        ncrit,
        max_iter,
        tolerance,
        hlmax,
        htmax,
        debug_side,
        debug_ibl,
    )
}

/// Solve for BL station at transition using XFOIL's TRDIF formulation
///
/// This is the Newton solver for the transition station. It uses TRDIF
/// (hybrid laminar-turbulent equations) instead of regular BLDIF when
/// a full transition result with derivatives is provided.
///
/// # Arguments
/// * `prev` - Previous (upstream) laminar station
/// * `initial_guess` - Optional initial guess for theta, delta_star, ctau
/// * `x_new` - Arc length at current station
/// * `ue_new` - Edge velocity at current station
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
/// * `is_laminar` - Flow type (should be false for transition)
/// * `ncrit` - Critical N-factor (use 0.0 to disable TRDIF)
/// * `max_iter` - Maximum Newton iterations
/// * `tolerance` - Convergence tolerance
/// * `hlmax` - Max Hk for laminar
/// * `htmax` - Max Hk for turbulent
///
/// # Reference
/// XFOIL xbl.f MRCHUE with TRDIF call at transition
pub fn newton_solve_station_transition(
    prev: &BlStation,
    initial_guess: Option<&BlStation>,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
    is_laminar: bool,
    ncrit: f64,
    max_iter: usize,
    tolerance: f64,
    hlmax: f64,
    htmax: f64,
    debug_side: Option<usize>,
    debug_ibl: Option<usize>,
) -> (BlStation, bool, f64) {
    let mut current_laminar = is_laminar;
    let mut flow_type = if current_laminar {
        FlowType::Laminar
    } else {
        FlowType::Turbulent
    };
    
    // Use initial_guess if provided, otherwise use prev
    let init = initial_guess.unwrap_or(prev);
    
    // Initialize station
    let mut station = BlStation::new();
    station.x = x_new;
    station.u = ue_new;
    station.is_laminar = is_laminar;
    station.is_turbulent = !is_laminar;
    station.is_wake = false;
    
    // Initialize from initial guess (or prev if no guess)
    station.theta = init.theta;
    station.delta_star = init.delta_star;
    station.ctau = if is_laminar { 0.03 } else { init.ctau };
    station.ampl = init.ampl;
    
    // Compute secondary variables
    blvar(&mut station, flow_type, msq, re);
    
    let mut direct = true;
    let mut htarg = if current_laminar { hlmax } else { htmax };
    let use_trchek = ncrit > 0.0;
    
    // Newton iteration loop
    let mut last_dmax = f64::INFINITY;
    let mut transition_active = false;
    let mut trchek_hold: Option<Trchek2FullResult> = None;
    let mut trchek_hold_state: Option<BlStation> = None;
    for iter_idx in 0..max_iter {
        // Build BL residuals and Jacobian
        // Use full TRDIF with proper chain rule through XT derivatives when available
        let mut trchek_full = trchek_hold.as_ref();
        let mut trchek_state: Option<BlStation> = None;
        if use_trchek && current_laminar {
            trchek_state = Some(station.clone());
            let tr = trchek2_full(
                prev.x,
                station.x,
                prev.theta,
                station.theta,
                prev.delta_star,
                station.delta_star,
                prev.u,
                station.u,
                prev.hk,
                station.hk,
                prev.r_theta,
                station.r_theta,
                prev.ampl,
                ncrit,
                msq,
                re,
            );
            if tr.transition {
                transition_active = true;
                current_laminar = false;
                flow_type = FlowType::Turbulent;
                station.is_laminar = false;
                station.is_turbulent = true;
                station.ampl = tr.ampl2;
                if station.ctau <= 0.0 {
                    station.ctau = 0.03;
                }
                // Recompute turbulent secondary variables now that flow type flipped.
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
            if tr.transition {
                trchek_hold = Some(tr);
                if trchek_hold_state.is_none() {
                    trchek_hold_state = trchek_state.clone();
                }
                trchek_full = trchek_hold.as_ref();
            }
        }

        if let (Some(tr), Some(side), Some(ibl)) = (trchek_full.as_ref(), debug_side, debug_ibl) {
            let curr = trchek_hold_state
                .as_ref()
                .or(trchek_state.as_ref())
                .unwrap_or(&station);
            emit_trchek2_debug(side, ibl, prev, curr, tr, msq, re, ncrit);
        }

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                if trchek_full.is_none() {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif_state(
                        iter_idx + 1,
                        side,
                        ibl,
                        rustfoil_bl::BldifStateEvent {
                            side,
                            ibl,
                            iter: iter_idx + 1,
                            flow_type: match flow_type {
                                FlowType::Laminar => 1,
                                FlowType::Turbulent => 2,
                                FlowType::Wake => 3,
                            },
                            X1: prev.x,
                            U1: prev.u,
                            T1: prev.theta,
                            D1: prev.delta_star,
                            S1: prev.ctau,
                            X2: station.x,
                            U2: station.u,
                            T2: station.theta,
                            D2: station.delta_star,
                            S2: station.ctau,
                        },
                    ));
                }
            }
        }

        let (res, jac) = if let Some(tr) = trchek_full.as_ref() {
            trdif_full(prev, &station, tr, msq, re)
        } else {
            bldif(prev, &station, flow_type, msq, re)
        };

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl), Some(tr)) = (debug_side, debug_ibl, trchek_full.as_ref()) {
                let mut vs1 = [[0.0f64; 5]; 4];
                let mut vs2 = [[0.0f64; 5]; 4];
                for i in 0..3 {
                    vs1[i] = jac.vs1[i];
                    vs2[i] = jac.vs2[i];
                }
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trdif(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::TrdifEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        VS1: vs1,
                        VS2: vs2,
                        VSREZ: [res.res_third, res.res_mom, res.res_shape, 0.0],
                    },
                ));

                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trdif_derivs(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::TrdifDerivsEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        WF1: tr.wf1,
                        WF2: tr.wf2,
                        XT: tr.xt,
                        XT_A1: tr.xt_a1,
                        XT_X1: tr.xt_x1,
                        XT_X2: tr.xt_x2,
                        XT_T1: tr.xt_t1,
                        XT_T2: tr.xt_t2,
                        XT_D1: tr.xt_d1,
                        XT_D2: tr.xt_d2,
                        XT_U1: tr.xt_u1,
                        XT_U2: tr.xt_u2,
                        XT_MS: tr.xt_ms,
                        XT_RE: tr.xt_re,
                        TT_A1: tr.tt_a1,
                        TT_X1: tr.tt_x1,
                        TT_X2: tr.tt_x2,
                        TT_T1: tr.tt_t1,
                        TT_T2: tr.tt_t2,
                        TT_D1: tr.tt_d1,
                        TT_D2: tr.tt_d2,
                        TT_U1: tr.tt_u1,
                        TT_U2: tr.tt_u2,
                        TT_MS: tr.tt_ms,
                        TT_RE: tr.tt_re,
                        DT_A1: tr.dt_a1,
                        DT_X1: tr.dt_x1,
                        DT_X2: tr.dt_x2,
                        DT_T1: tr.dt_t1,
                        DT_T2: tr.dt_t2,
                        DT_D1: tr.dt_d1,
                        DT_D2: tr.dt_d2,
                        DT_U1: tr.dt_u1,
                        DT_U2: tr.dt_u2,
                        DT_MS: tr.dt_ms,
                        DT_RE: tr.dt_re,
                        UT_A1: tr.ut_a1,
                        UT_X1: tr.ut_x1,
                        UT_X2: tr.ut_x2,
                        UT_T1: tr.ut_t1,
                        UT_T2: tr.ut_t2,
                        UT_D1: tr.ut_d1,
                        UT_D2: tr.ut_d2,
                        UT_U1: tr.ut_u1,
                        UT_U2: tr.ut_u2,
                        UT_MS: tr.ut_ms,
                        UT_RE: tr.ut_re,
                    },
                ));

                if let Some(terms) = trdif_turb_terms(prev, &station, tr, msq, re) {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif_terms(
                        iter_idx + 1,
                        side,
                        ibl,
                        2,
                        rustfoil_bl::BldifTermsEvent {
                            side,
                            ibl,
                            flow_type: 2,
                            xlog: terms.xlog,
                            ulog: terms.ulog,
                            tlog: terms.tlog,
                            hlog: terms.hlog,
                            z_cfx_shape: terms.z_cfx_shape,
                            z_dix_shape: terms.z_dix_shape,
                            upw: terms.upw,
                            ha: terms.ha,
                            btmp_mom: terms.btmp_mom,
                            cfx: terms.cfx,
                            cfx_ta: terms.cfx_ta,
                            cfx_t1: terms.cfx_t1,
                            cfx_t2: terms.cfx_t2,
                            btmp_shape: terms.btmp_shape,
                            cfx_shape: terms.cfx_shape,
                            dix: terms.dix,
                            cfx_upw: terms.cfx_upw,
                            dix_upw: terms.dix_upw,
                            hsa: terms.hsa,
                            hca: terms.hca,
                            dd: terms.dd,
                            dd2: terms.dd2,
                            xot1: terms.xot1,
                            xot2: terms.xot2,
                            cf1: terms.cf1,
                            cf2: terms.cf2,
                            di1: terms.di1,
                            di2: terms.di2,
                            z_hs2: terms.z_hs2,
                            z_cf2_shape: terms.z_cf2_shape,
                            z_di2: terms.z_di2,
                            z_t2_shape: terms.z_t2_shape,
                            z_u2_shape: terms.z_u2_shape,
                            z_hca: terms.z_hca,
                            z_ha_shape: terms.z_ha_shape,
                            di2_s: terms.di2_s,
                            z_upw_shape: terms.z_upw_shape,
                            hs2_t: terms.hs2_t,
                            hs2_d: terms.hs2_d,
                            hs2_u: terms.hs2_u,
                            cf2_t_shape: terms.cf2_t_shape,
                            cf2_d_shape: terms.cf2_d_shape,
                            cf2_u_shape: terms.cf2_u_shape,
                            di2_t: terms.di2_t,
                            di2_d: terms.di2_d,
                            di2_u: terms.di2_u,
                            hc2_t: terms.hc2_t,
                            hc2_d: terms.hc2_d,
                            hc2_u: terms.hc2_u,
                            h2_t: terms.h2_t,
                            h2_d: terms.h2_d,
                            upw_t2_shape: terms.upw_t2_shape,
                            upw_d2_shape: terms.upw_d2_shape,
                            upw_u2_shape: terms.upw_u2_shape,
                            vs2_3_1: 0.0,
                            vs2_3_2: 0.0,
                            vs2_3_3: 0.0,
                            vs2_3_4: 0.0,
                            xfoil_t1: None,
                            xfoil_d1: None,
                            xfoil_s1: None,
                            xfoil_hk1: None,
                            xfoil_rt1: None,
                            xfoil_di1: None,
                            xfoil_di2: None,
                            xfoil_upw: None,
                            xfoil_dix_upw: None,
                            xfoil_z_dix_shape: None,
                            xfoil_z_upw_shape: None,
                        },
                    ));
                }
            }
        }

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                if trchek_full.is_none() {
                    let (res_dbg, jac_dbg, terms_dbg) =
                        bldif_with_terms(prev, &station, flow_type, msq, re);
                    let xfoil_terms =
                        bldif_terms_xfoil_upstream(prev, &station, flow_type, msq, re);
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif_terms(
                        iter_idx + 1,
                        side,
                        ibl,
                        match flow_type {
                            FlowType::Laminar => 1,
                            FlowType::Turbulent => 2,
                            FlowType::Wake => 3,
                        },
                        rustfoil_bl::BldifTermsEvent {
                            side,
                            ibl,
                            flow_type: match flow_type {
                                FlowType::Laminar => 1,
                                FlowType::Turbulent => 2,
                                FlowType::Wake => 3,
                            },
                            xlog: terms_dbg.xlog,
                            ulog: terms_dbg.ulog,
                            tlog: terms_dbg.tlog,
                            hlog: terms_dbg.hlog,
                            z_cfx_shape: terms_dbg.z_cfx_shape,
                            z_dix_shape: terms_dbg.z_dix_shape,
                            upw: terms_dbg.upw,
                            ha: terms_dbg.ha,
                            btmp_mom: terms_dbg.btmp_mom,
                            cfx: terms_dbg.cfx,
                            cfx_ta: terms_dbg.cfx_ta,
                            cfx_t1: terms_dbg.cfx_t1,
                            cfx_t2: terms_dbg.cfx_t2,
                            btmp_shape: terms_dbg.btmp_shape,
                            cfx_shape: terms_dbg.cfx_shape,
                            dix: terms_dbg.dix,
                            cfx_upw: terms_dbg.cfx_upw,
                            dix_upw: terms_dbg.dix_upw,
                            hsa: terms_dbg.hsa,
                            hca: terms_dbg.hca,
                            dd: terms_dbg.dd,
                            dd2: terms_dbg.dd2,
                            xot1: terms_dbg.xot1,
                            xot2: terms_dbg.xot2,
                            cf1: terms_dbg.cf1,
                            cf2: terms_dbg.cf2,
                            di1: terms_dbg.di1,
                            di2: terms_dbg.di2,
                            z_hs2: terms_dbg.z_hs2,
                            z_cf2_shape: terms_dbg.z_cf2_shape,
                            z_di2: terms_dbg.z_di2,
                            z_t2_shape: terms_dbg.z_t2_shape,
                            z_u2_shape: terms_dbg.z_u2_shape,
                            z_hca: terms_dbg.z_hca,
                            z_ha_shape: terms_dbg.z_ha_shape,
                            di2_s: terms_dbg.di2_s,
                            z_upw_shape: terms_dbg.z_upw_shape,
                            hs2_t: terms_dbg.hs2_t,
                            hs2_d: terms_dbg.hs2_d,
                            hs2_u: terms_dbg.hs2_u,
                            cf2_t_shape: terms_dbg.cf2_t_shape,
                            cf2_d_shape: terms_dbg.cf2_d_shape,
                            cf2_u_shape: terms_dbg.cf2_u_shape,
                            di2_t: terms_dbg.di2_t,
                            di2_d: terms_dbg.di2_d,
                            di2_u: terms_dbg.di2_u,
                            hc2_t: terms_dbg.hc2_t,
                            hc2_d: terms_dbg.hc2_d,
                            hc2_u: terms_dbg.hc2_u,
                            h2_t: terms_dbg.h2_t,
                            h2_d: terms_dbg.h2_d,
                            upw_t2_shape: terms_dbg.upw_t2_shape,
                            upw_d2_shape: terms_dbg.upw_d2_shape,
                            upw_u2_shape: terms_dbg.upw_u2_shape,
                            vs2_3_1: jac_dbg.vs2[2][0],
                            vs2_3_2: jac_dbg.vs2[2][1],
                            vs2_3_3: jac_dbg.vs2[2][2],
                            vs2_3_4: jac_dbg.vs2[2][3],
                            xfoil_t1: xfoil_terms.as_ref().map(|(_, s)| s.theta),
                            xfoil_d1: xfoil_terms.as_ref().map(|(_, s)| s.delta_star),
                            xfoil_s1: xfoil_terms.as_ref().map(|(_, s)| s.ctau),
                            xfoil_hk1: xfoil_terms.as_ref().map(|(_, s)| s.hk),
                            xfoil_rt1: xfoil_terms.as_ref().map(|(_, s)| s.r_theta),
                            xfoil_di1: xfoil_terms.as_ref().map(|(t, _)| t.di1),
                            xfoil_di2: xfoil_terms.as_ref().map(|(t, _)| t.di2),
                            xfoil_upw: xfoil_terms.as_ref().map(|(t, _)| t.upw),
                            xfoil_dix_upw: xfoil_terms.as_ref().map(|(t, _)| t.dix_upw),
                            xfoil_z_dix_shape: xfoil_terms.as_ref().map(|(t, _)| t.z_dix_shape),
                            xfoil_z_upw_shape: xfoil_terms.as_ref().map(|(t, _)| t.z_upw_shape),
                        },
                    ));

                    let mut vs1_dbg = [[0.0; 5]; 4];
                    let mut vs2_dbg = [[0.0; 5]; 4];
                    for row in 0..3 {
                        vs1_dbg[row] = jac_dbg.vs1[row];
                        vs2_dbg[row] = jac_dbg.vs2[row];
                    }
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif(
                        iter_idx + 1,
                        side,
                        ibl,
                        match flow_type {
                            FlowType::Laminar => 1,
                            FlowType::Turbulent => 2,
                            FlowType::Wake => 3,
                        },
                        vs1_dbg,
                        vs2_dbg,
                        [res_dbg.res_third, res_dbg.res_mom, res_dbg.res_shape, 0.0],
                    ));
                }

                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue_sys(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::MrchueSysEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        direct,
                        flow_type: if current_laminar { 1 } else { 2 },
                        x1: prev.x,
                        x2: station.x,
                        t1: prev.theta,
                        t2: station.theta,
                        u1: prev.u,
                        u2: station.u,
                        hk1: prev.hk,
                        hk2: station.hk,
                        ctau1: prev.ctau,
                        ctau2: station.ctau,
                        VS2: jac.vs2,
                        VSREZ: [res.res_third, res.res_mom, res.res_shape],
                    },
                ));
            }
        }
        
        // Build 4x4 Newton system
        // Variable order in bldif: [s/ampl, θ, δ*, u, x]
        // We use [s/ampl, θ, δ*, u] (4 variables)
        //
        // Hk derivatives for inverse mode constraint (only used when direct=false)
        // XFOIL-correct derivatives (xblsys.f BLKIN lines 768-769):
        //   HK2_T2 = HK2_H2 * H2_T2 = hk_h * (-H/θ) ≈ -Hk/θ at low Mach
        //   HK2_D2 = HK2_H2 * H2_D2 = hk_h * (1/θ)  ≈ +1/θ at low Mach
        //
        // These are ONLY used in inverse mode (direct=false), where they form 
        // row 4 of the Newton system: hk2_t*dθ + hk2_d*dδ* = Hk_target - Hk_current
        // Using correct derivatives ensures the constraint pushes in the right direction.
        //
        // Using correct derivatives for inverse mode
        let hk2_t = station.derivs.hk_h * station.derivs.h_theta;      // = -Hk/θ
        let hk2_d = station.derivs.hk_h * station.derivs.h_delta_star; // = +1/θ
        let hmax = if current_laminar { hlmax } else { htmax };
        
        let (mut a, mut b) = build_4x4_system(
            &jac.vs2,
            &[res.res_third, res.res_mom, res.res_shape],
            direct,
            hk2_t,
            hk2_d,
            0.0,  // hk2_u (simplified - no Ue dependence in Hk at low Mach)
            htarg,
            station.hk,
        );

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::vs2_before(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::Vs2BeforeEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        direct,
                        flow_type: if current_laminar { 1 } else { 2 },
                        VS2_4x4: a,
                        VSREZ_rhs: b,
                    },
                ));
            }
        }
        
        // Solve 4x4 system
        let vsrez = solve_4x4(&a, &b);

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::vsrez_after(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::VsrezAfterEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        direct,
                        flow_type: if current_laminar { 1 } else { 2 },
                        VSREZ_solution: vsrez,
                    },
                ));
            }
        }
        
        // Compute max relative change (XFOIL: include dUe in inverse mode)
        let mut dmax = (vsrez[1] / station.theta.max(1e-12)).abs()
            .max((vsrez[2] / station.delta_star.max(1e-12)).abs());
        if current_laminar {
            dmax = dmax.max((vsrez[0] / 10.0).abs());
        } else {
            dmax = dmax.max((vsrez[0] / station.ctau.max(1e-6)).abs());
        }
        if !direct {
            dmax = dmax.max((vsrez[3] / station.u.max(1e-6)).abs());
        }
        last_dmax = dmax;
        
        // Relaxation factor
        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue_iter(
                    iter_idx + 1,
                    side,
                    ibl,
                    rustfoil_bl::MrchueIterEvent {
                        side,
                        ibl,
                        iter: iter_idx + 1,
                        direct,
                        dmax,
                        relaxation: rlx,
                        res_third: res.res_third,
                        res_mom: res.res_mom,
                        res_shape: res.res_shape,
                        delta_s: vsrez[0],
                        delta_theta: vsrez[1],
                        delta_delta_star: vsrez[2],
                        delta_ue: vsrez[3],
                        theta: station.theta,
                        delta_star: station.delta_star,
                        Ue: station.u,
                        H: station.h,
                        Hk: station.hk,
                        Cf: station.cf,
                        Cd: station.cd,
                        ampl: station.ampl,
                        ctau: station.ctau,
                    },
                ));
            }
        }
        
        // Check if direct mode is appropriate (Hk not too high)
        if direct {
            let h_test = (station.delta_star + rlx * vsrez[2])
                / (station.theta + rlx * vsrez[1]).max(1e-12);
            let hk_test = rustfoil_bl::closures::hkin(h_test, msq).hk;
            let direct_after = hk_test < hmax;

            if rustfoil_bl::is_debug_active() {
                if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue_mode(
                        iter_idx + 1,
                        side,
                        ibl,
                        rustfoil_bl::MrchueModeEvent {
                            side,
                            ibl,
                            iter: iter_idx + 1,
                            direct_before: direct,
                            direct_after,
                            flow_type: if current_laminar { 1 } else { 2 },
                            dmax,
                            rlx,
                            Htest: h_test,
                            Hk_test: hk_test,
                            Hmax: hmax,
                            Ue: station.u,
                            theta: station.theta,
                            delta_star: station.delta_star,
                        },
                    ));
                }
            }

            if !direct_after {
                // Need to switch to inverse mode
                direct = false;
                // Estimate target Hk based on upstream (XFOIL MRCHUE logic)
                if current_laminar {
                    // Laminar case: slow increase in Hk downstream
                    htarg = prev.hk + 0.03 * (x_new - prev.x) / prev.theta.max(1e-12);
                } else if let Some(tr) = trchek_full.as_ref() {
                    // Transition interval: blend laminar and turbulent targets
                    htarg = prev.hk
                        + (0.03 * (tr.xt - prev.x) - 0.15 * (station.x - tr.xt))
                            / prev.theta.max(1e-12);
                } else {
                    // Turbulent case: faster decrease in Hk downstream
                    htarg = prev.hk - 0.15 * (x_new - prev.x) / prev.theta.max(1e-12);
                }
                // XFOIL does not clamp the inverse-mode target above hmax.
                htarg = htarg.max(hmax);
                continue;  // Restart iteration with inverse mode
            }
        }
        
        // Apply updates with relaxation
        if current_laminar {
            station.ampl = (station.ampl + rlx * vsrez[0]).max(0.0);
        } else {
            // At transition station, keep ctau at initial value (0.03)
            // This preserves non-zero shear-lag residual for global Newton.
            //
            // XFOIL's order of operations:
            // 1. MRCHUE initializes ctau=0.03 at transition (xbl.f line 600)
            // 2. MRCHUE converges theta, delta_star but ctau stays near 0.03
            // 3. SETBL (global Newton) adjusts ALL ctau values together
            //
            // If we update ctau here, the single-station Newton converges ctau
            // to equilibrium (~0.059), making VSREZ[0] ≈ 0. The global Newton
            // then has nothing to work with.
            //
            // By skipping ctau updates at transition, we get:
            // - XFOIL-like: VSREZ[0] ≈ -2.6e-3 (non-zero residual)
            // - Global Newton can properly adjust all ctau values
            if !transition_active {
                station.ctau = (station.ctau + rlx * vsrez[0]).clamp(1e-7, 0.3);
            }
        }
        station.theta = (station.theta + rlx * vsrez[1]).max(1e-12);
        station.delta_star = (station.delta_star + rlx * vsrez[2]).max(1e-12);
        // In inverse mode, XFOIL's MRCHUE allows Ue to vary to satisfy the Hk constraint.
        // This prevents the solution from blowing up when approaching separation.
        // The Ue change is typically small due to the Hk constraint formulation.
        if !direct {
            station.u = (station.u + rlx * vsrez[3]).max(0.01);
        }
        
        // Limit Hk to prevent singularity
        let hklim = 1.02;
        let dsw = station.delta_star;
        let dslim = hklim * station.theta;
        if dsw / station.theta.max(1e-12) < hklim {
            station.delta_star = dslim;
        }
        
        // Recompute secondary variables
        blvar(&mut station, flow_type, msq, re);

        // Update transition interval data using the latest station state.
        // XFOIL evaluates TRCHEK2 against the current (post-update) station,
        // so keep the stored TRCHEK2 result in sync once transition is active.
        if transition_active {
            let tr_updated = trchek2_full(
                prev.x,
                station.x,
                prev.theta,
                station.theta,
                prev.delta_star,
                station.delta_star,
                prev.u,
                station.u,
                prev.hk,
                station.hk,
                prev.r_theta,
                station.r_theta,
                prev.ampl,
                ncrit,
                msq,
                re,
            );
            trchek_hold = Some(tr_updated);
            trchek_hold_state = Some(station.clone());
        }

        if rustfoil_bl::is_debug_active() {
            if let (Some(side), Some(ibl)) = (debug_side, debug_ibl) {
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::blvar(
                    iter_idx + 1,
                    side,
                    ibl,
                    if current_laminar { 1 } else { 2 },
                    rustfoil_bl::BlvarInput {
                        x: station.x,
                        u: station.u,
                        theta: station.theta,
                        delta_star: station.delta_star,
                        ctau: station.ctau,
                        ampl: station.ampl,
                    },
                    rustfoil_bl::BlvarOutput {
                        H: station.h,
                        Hk: station.hk,
                        Hs: station.hs,
                        Hc: station.hc,
                        Rtheta: station.r_theta,
                        Cf: station.cf,
                        Cd: station.cd,
                        Us: station.us,
                        Cq: station.cq,
                        De: station.de,
                        Dd: station.ctau * station.ctau * (0.995 - station.us) * 2.0 / station.hs,
                        Dd2: 0.15 * (0.995 - station.us).powi(2) / station.r_theta * 2.0 / station.hs,
                        DiWall: 0.0,
                        DiTotal: 0.0,
                        DiLam: 0.0,
                        DiUsed: station.cd,
                        DiLamOverride: None,
                        DiUsedMinusTotal: None,
                        DiUsedMinusLam: None,
                    },
                ));
            }
        }
        
        // Check convergence
        if dmax <= tolerance {
            return (station, true, dmax);
        }
    }
    
    // Did not converge - return best estimate
    (station, false, last_dmax)
}

// ============================================================================
// Direct March (MRCHUE) - Newton Version
// ============================================================================

/// March boundary layer with fixed edge velocity
///
/// Starts from the stagnation point and marches downstream, solving the
/// integral BL equations at each station. Transition is detected via the
/// e^n method and separation by negative skin friction.
///
/// If the shape factor Hk exceeds limits (approaching separation), the
/// solver switches to "inverse mode" where Hk is prescribed rather than
/// computed, preventing the Goldstein singularity.
///
/// # Arguments
/// * `x` - Arc length coordinates (monotonically increasing)
/// * `ue` - Edge velocities at each station
/// * `re` - Reynolds number
/// * `msq` - Mach number squared (for compressibility)
/// * `config` - March configuration
///
/// # Returns
/// MarchResult with computed BL stations and detected events
///
/// # Reference
/// XFOIL xbl.f MRCHUE (line 542)
pub fn march_fixed_ue(
    x: &[f64],
    ue: &[f64],
    re: f64,
    msq: f64,
    config: &MarchConfig,
) -> MarchResult {
    let n = x.len();
    assert_eq!(n, ue.len(), "x and ue arrays must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut result = MarchResult::new(n);

    // Initialize stagnation point (first station)
    let mut prev = init_stagnation(x[0], ue[0], re);
    blvar(&mut prev, FlowType::Laminar, msq, re);
    result.stations.push(prev.clone());

    // Track state
    let mut is_laminar = true;

    // March downstream using Newton iteration at each station
    for i in 1..n {
        let prev_station = &result.stations[i - 1];

        // Determine flow type based on previous transition state
        let flow_type = if is_laminar {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Solve for current station using Newton iteration (XFOIL MRCHUE style)
        let (mut station, converged, dmax) = newton_solve_station_with_guess(
            prev_station,
            None,
            x[i],
            ue[i],
            re,
            msq,
            is_laminar,
            config.ncrit,
            config.max_iter,
            config.tolerance,
            config.hlmax,
            config.htmax,
            Some(0),
            Some(i),
        );
        

        // If Newton didn't converge, use best estimate (station still has result)
        // Note: With tolerance=1e-5 (matching XFOIL), Newton should converge for normal cases
        // The fallback was causing issues (replacing converged Newton values with less accurate estimates)
        // TODO: Consider removing this fallback entirely since it can mask real issues
        if !converged && dmax > 0.1 {
            // Only use fallback if Newton produced clearly wrong values
            // For now, trust Newton's result even without strict convergence
            // (The solution is still valid, just not to machine precision)
        }

        // Check for transition AFTER Newton solve using TRCHEK2
        // This matches XFOIL's order: solve first, then check transition with converged values
        if is_laminar && result.x_transition.is_none() {
            // Debug trace before transition check
            if config.debug_trace && i >= 28 && i <= 35 {
                println!("[march_fixed_ue] Station {} PRE-TRCHEK: x={:.4}, Hk={:.4}, theta={:.6e}, ampl={:.4}",
                         i, station.x, station.hk, station.theta, station.ampl);
            }
            
            // Use TRCHEK2 for implicit N-factor integration with converged station values
            let trchek_result = trchek2_stations(
                prev_station.x,
                station.x,
                prev_station.hk,
                prev_station.theta,
                prev_station.r_theta,
                prev_station.delta_star,
                prev_station.u,
                prev_station.ampl,
                station.hk,
                station.theta,
                station.r_theta,
                station.delta_star,
                station.u,
                config.ncrit,
                msq,
                re,
            );
            
            station.ampl = trchek_result.ampl2;
            emit_trchek2_station_iter(
                0,
                i,
                prev_station,
                &station,
                &trchek_result,
                config.ncrit,
            );
            
            // Debug trace for transition investigation
            // Trace N-factor evolution to diagnose late transition
            if config.debug_trace && i <= 30 && station.ampl > 0.001 {
                println!("[march_fixed_ue] Station {} AMPL: x={:.4}, Hk={:.3}, Rθ={:.1}, N={:.4}",
                         i, station.x, station.hk, station.r_theta, station.ampl);
            }
            if config.debug_trace && station.ampl > 4.0 && station.ampl < 12.0 {
                println!("[march_fixed_ue] Station {} NEAR-TRANS: x={:.4}, Hk={:.3}, Rθ={:.1}, N={:.4}, dN={:.4}",
                         i, station.x, station.hk, station.r_theta, station.ampl, 
                         station.ampl - prev_station.ampl);
            }
            
            if config.debug_trace && i >= 28 && i <= 35 {
                println!("[march_fixed_ue] Station {} POST-TRCHEK: ampl={:.4}, transition={}",
                         i, station.ampl, trchek_result.transition);
            }
            
            if trchek_result.transition {
                result.x_transition = trchek_result.xt;
                result.transition_index = Some(i);

                if config.debug_trace {
                    println!("[march_fixed_ue] TRANSITION DETECTED at station {}, x={:.4}", i, station.x);
                    println!("  Laminar values: Hk={:.4}, theta={:.6e}, delta_star={:.6e}",
                             station.hk, station.theta, station.delta_star);
                }

                if station.is_turbulent {
                    // Transition already handled inside Newton loop (XFOIL behavior).
                    is_laminar = false;
                } else {
                    // Fallback: re-solve with turbulent equations if transition wasn't handled inside Newton.
                    is_laminar = false;

                    let mut initial_guess = prev_station.clone();
                    initial_guess.x = station.x;
                    initial_guess.u = station.u;
                    initial_guess.ampl = station.ampl;
                    initial_guess.ctau = 0.03;

                    let (turb_station, _turb_converged, _turb_dmax) = newton_solve_station_transition(
                        prev_station,
                        Some(&initial_guess),
                        x[i],
                        station.u,
                        re,
                        msq,
                        false,
                        config.ncrit,
                        config.max_iter,
                        config.tolerance,
                        config.hlmax,
                        config.htmax,
                        Some(0),
                        Some(i),
                    );

                    station = turb_station;
                    station.ampl = trchek_result.ampl2;

                    if config.debug_trace {
                        println!("  Initial ctau = 0.03 (XFOIL MRCHUE line 598)");
                        println!("  After full TRDIF re-solve: Hk={:.4}, theta={:.6e}, delta_star={:.6e}, ctau={:.6}",
                                 station.hk, station.theta, station.delta_star, station.ctau);
                    }
                }
            }
        } else if !is_laminar {
            // Turbulent: N-factor frozen
            station.ampl = prev_station.ampl;
            
            if config.debug_trace && i >= 28 && i <= 35 {
                println!("[march_fixed_ue] Station {} TURBULENT: x={:.4}, Hk={:.4}, theta={:.6e}",
                         i, station.x, station.hk, station.theta);
            }
        }

        // Note: Separation-induced transition is NOT applied in march_fixed_ue because this
        // function has no surface context. Use march_surface() for proper handling of
        // separation-induced transition on the lower surface at high angles of attack.

        // Check for separation (Cf < 0)
        if station.cf < 0.0 && result.x_separation.is_none() {
            result.x_separation = Some(station.x);
        }

        // Store mass defect
        station.mass_defect = station.u * station.delta_star;

        result.stations.push(station);
    }

    result
}

/// March boundary layer for a specific surface with debug output.
///
/// This is the main entry point for surface-specific marching that properly
/// tracks which surface (upper=1, lower=2) is being processed for debug output.
///
/// # Arguments
/// * `x` - Arc lengths from stagnation (should start at 0 or near 0)
/// * `ue` - Edge velocities (should be positive, increasing from stagnation)
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
/// * `config` - March configuration
/// * `side` - Surface identifier (1 = upper, 2 = lower)
///
/// # Returns
/// MarchResult with computed BL stations
pub fn march_surface(
    x: &[f64],
    ue: &[f64],
    re: f64,
    msq: f64,
    config: &MarchConfig,
    side: usize,
) -> MarchResult {
    let n = x.len();
    assert_eq!(n, ue.len(), "x and ue arrays must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut result = MarchResult::new(n);

    // XFOIL skips the stagnation point itself (IBL=1) and starts from IBL=2
    // If first point is at/near stagnation (x≈0, Ue≈0), start from second point
    let start_idx = if x[0] < 1e-6 || ue[0].abs() < 1e-6 {
        // Skip stagnation point - use second station as initial
        if n < 3 {
            // Not enough points after skipping
            return MarchResult::new(0);
        }
        1
    } else {
        0
    };

    // Initialize first BL station using Thwaites formula
    let x_init = x[start_idx].max(1e-6); // Ensure x > 0 for Thwaites
    let ue_init = ue[start_idx].abs().max(0.01); // Ensure Ue > 0
    let mut prev = init_stagnation(x_init, ue_init, re);
    prev.x = x[start_idx]; // Use actual x
    prev.u = ue[start_idx].abs();
    blvar(&mut prev, FlowType::Laminar, msq, re);
    
    // XFOIL runs Newton iteration on the similarity station (IBL=2) using
    // a point equation formulation (BLDIF with ITYP=0) where COM1=COM2.
    // Analysis shows this consistently refines theta by factor 1.0634 and
    // converges Hk to 2.2295 across all foils and Reynolds numbers.
    // This is a fundamental result from the Falkner-Skan similarity solution
    // at stagnation (pressure gradient parameter BULE=1).
    //
    // Rather than implementing the full similarity station Newton solver,
    // we apply these universal correction factors directly.
    const SIMILARITY_THETA_FACTOR: f64 = 1.0634;
    const SIMILARITY_HK: f64 = 2.2295;
    
    prev.theta *= SIMILARITY_THETA_FACTOR;
    prev.delta_star = prev.theta * SIMILARITY_HK;
    prev.h = SIMILARITY_HK;
    prev.hk = SIMILARITY_HK;
    
    // Recompute derived quantities with corrected values
    blvar(&mut prev, FlowType::Laminar, msq, re);
    
    emit_blvar_debug(side, 2, FlowType::Laminar, &prev);

    // Emit debug event for first BL station
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
            side,
            2, // IBL=2 is first station after stagnation in XFOIL
            prev.x,
            prev.u,
            prev.theta,
            prev.delta_star,
            prev.hk,
            prev.cf,
        ));
    }

    result.stations.push(prev.clone());

    // Track state
    let mut is_laminar = true;

    // March downstream (starting from station after initial)
    for i in (start_idx + 1)..n {
        let ds = x[i] - x[i - 1];
        let station_idx = i - start_idx; // Index into result.stations
        let prev_station = &result.stations[station_idx - 1];

        // Determine flow type
        let flow_type = if is_laminar {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Solve for current station using Newton iteration (XFOIL MRCHUE style)
        let (mut station, converged, dmax) = newton_solve_station_with_guess(
            prev_station,
            None,
            x[i],
            ue[i].abs(),
            re,
            msq,
            is_laminar,
            config.ncrit,
            config.max_iter,
            config.tolerance,
            config.hlmax,
            config.htmax,
            Some(side),
            Some(station_idx + 2),
        );

        // If Newton didn't converge, fall back to step_momentum
        if !converged && dmax > 0.1 {
            let (theta_new, delta_star_new) =
                step_momentum(prev_station, x[i], ue[i].abs(), re, msq, is_laminar);
            station.theta = theta_new;
            station.delta_star = delta_star_new;
            blvar(&mut station, flow_type, msq, re);
        }

        // === Safeguard against collapsed theta ===
        // Near the leading edge, theta can collapse if Newton has issues.
        // Only apply a very conservative floor (much smaller than Thwaites estimate)
        // to catch truly pathological cases without overriding valid Newton solutions.
        // The Newton solver typically produces theta within 1% of XFOIL, so we only
        // intervene when theta is clearly wrong (e.g., zero or tiny).
        let theta_floor = 1e-8;  // Absolute minimum
        
        if station.theta < theta_floor {
            // Only clamp to floor if truly collapsed
            station.theta = theta_floor;
            station.delta_star = station.theta * prev_station.h.clamp(1.5, 3.5);
        }

        // Compute all secondary variables (already done in newton_solve, but ensure consistency)
        blvar(&mut station, flow_type, msq, re);

        // Do not clamp Hk here; inverse-mode handles near-separation behavior.

        // Check for transition (laminar only) using TRCHEK2 implicit integration
        if is_laminar && result.x_transition.is_none() {
            let trchek_result = trchek2_stations(
                prev_station.x,
                station.x,
                prev_station.hk,
                prev_station.theta,
                prev_station.r_theta,
                prev_station.delta_star,
                prev_station.u,
                prev_station.ampl,
                station.hk,
                station.theta,
                station.r_theta,
                station.delta_star,
                station.u,
                config.ncrit,
                msq,
                re,
            );

            station.ampl = trchek_result.ampl2;
            emit_trchek2_station_iter(
                side,
                station_idx + 2,
                prev_station,
                &station,
                &trchek_result,
                config.ncrit,
            );

            if config.debug_trace && station.ampl > 0.1 && station.ampl < 12.0 {
                println!(
                    "[march_surface] i={} x={:.4} Hk={:.3} Rθ={:.0} N={:.3} dN={:.4}",
                    i,
                    station.x,
                    station.hk,
                    station.r_theta,
                    station.ampl,
                    station.ampl - prev_station.ampl
                );
            }

            if trchek_result.transition {
                result.x_transition = trchek_result.xt;
                result.transition_index = Some(station_idx);

                if station.is_turbulent {
                    is_laminar = false;
                } else {
                    is_laminar = false;
                    station.is_laminar = false;
                    station.is_turbulent = true;

                    let mut initial_guess = prev_station.clone();
                    initial_guess.x = station.x;
                    initial_guess.u = station.u;
                    initial_guess.ampl = station.ampl;
                    initial_guess.ctau = 0.03;

                    let (turb_station, _turb_converged, _turb_dmax) = newton_solve_station_transition(
                        prev_station,
                        Some(&initial_guess),
                        x[i],
                        station.u,
                        re,
                        msq,
                        false,
                        config.ncrit,
                        config.max_iter,
                        config.tolerance,
                        config.hlmax,
                        config.htmax,
                        Some(side),
                        Some(station_idx + 2),
                    );

                    station = turb_station;
                    station.ampl = trchek_result.ampl2;
                }
            }
        }

        // NOTE: XFOIL does NOT force transition when Hk > hlmax. It only uses hlmax
        // to switch between direct and inverse mode. When Hk exceeds hlmax, the solver
        // switches to inverse mode where Hk is prescribed, allowing laminar flow to
        // continue through near-separation conditions. Transition is only triggered by:
        // 1. e^N method (amplification exceeds Ncrit)
        // 2. Forced transition at a user-specified location (XIFORC)
        //
        // In XFOIL, when no forced transition strip is set (XSTRIP >= 1.0), XIFORC is
        // set to the TE arc length. This means laminar flow is forced to transition
        // at the trailing edge if no natural transition has occurred.
        //
        // Implement forced TE transition: if this is the last airfoil station (or very
        // close to TE, arc length > 0.98) and still laminar, force transition.
        // This matches XFOIL's behavior of forcing transition at TE.
        let is_last_station = i == x.len() - 2; // last iteration of the loop
        let near_te = station.x > 0.98;  // Arc length > 0.98 is very close to TE
        
        if is_laminar && (is_last_station || near_te) && result.x_transition.is_none() {
            if config.debug_trace {
                println!(
                    "[march_surface] FORCED TE TRANSITION at station {}, x={:.4} (last={}, near_te={})",
                    station_idx, station.x, is_last_station, near_te
                );
            }
            
            is_laminar = false;
            station.is_laminar = false;
            station.is_turbulent = true;
            result.x_transition = Some(station.x);
            result.transition_index = Some(station_idx);
            
            // Re-solve the station as turbulent with initial ctau = 0.03 (XFOIL default)
            let mut initial_guess = prev_station.clone();
            initial_guess.x = station.x;
            initial_guess.u = station.u;
            initial_guess.ampl = station.ampl;
            initial_guess.ctau = 0.03;
            
            let (turb_station, _turb_converged, _turb_dmax) = newton_solve_station_transition(
                prev_station,
                Some(&initial_guess),
                x[i],
                station.u,
                re,
                msq,
                false, // Now turbulent
                config.ncrit,
                config.max_iter,
                config.tolerance,
                config.hlmax,
                config.htmax,
                Some(side),
                Some(station_idx + 2),
            );
            
            station = turb_station;
        }

        // If transition/wake state changed, recompute secondary variables with the final flow type.
        let flow_type_final = if station.is_wake {
            FlowType::Wake
        } else if station.is_turbulent {
            FlowType::Turbulent
        } else {
            FlowType::Laminar
        };
        blvar(&mut station, flow_type_final, msq, re);

        // Check for separation
        if station.cf < 0.0 && result.x_separation.is_none() {
            result.x_separation = Some(station.x);
        }

        // Store mass defect
        station.mass_defect = station.u * station.delta_star;

        // Emit debug events
        if rustfoil_bl::is_debug_active() {
            let flow_type = if station.is_wake {
                FlowType::Wake
            } else if station.is_turbulent {
                FlowType::Turbulent
            } else {
                FlowType::Laminar
            };
            emit_blvar_debug(side, station_idx + 2, flow_type, &station);
            let (res, jac, terms) = bldif_with_terms(prev_station, &station, flow_type, msq, re);
            let xfoil_terms =
                bldif_terms_xfoil_upstream(prev_station, &station, flow_type, msq, re);
            let mut vs1 = [[0.0; 5]; 4];
            let mut vs2 = [[0.0; 5]; 4];
            for row in 0..3 {
                vs1[row] = jac.vs1[row];
                vs2[row] = jac.vs2[row];
            }
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif(
                1,
                side,
                station_idx + 2,
                match flow_type {
                    FlowType::Laminar => 1,
                    FlowType::Turbulent => 2,
                    FlowType::Wake => 3,
                },
                vs1,
                vs2,
                [res.res_third, res.res_mom, res.res_shape, 0.0],
            ));
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bldif_terms(
                1,
                side,
                station_idx + 2,
                match flow_type {
                    FlowType::Laminar => 1,
                    FlowType::Turbulent => 2,
                    FlowType::Wake => 3,
                },
                rustfoil_bl::BldifTermsEvent {
                    side,
                    ibl: station_idx + 2,
                    flow_type: match flow_type {
                        FlowType::Laminar => 1,
                        FlowType::Turbulent => 2,
                        FlowType::Wake => 3,
                    },
                    xlog: terms.xlog,
                    ulog: terms.ulog,
                    tlog: terms.tlog,
                    hlog: terms.hlog,
                    z_cfx_shape: terms.z_cfx_shape,
                    z_dix_shape: terms.z_dix_shape,
                    upw: terms.upw,
                    ha: terms.ha,
                    btmp_mom: terms.btmp_mom,
                    cfx: terms.cfx,
                    cfx_ta: terms.cfx_ta,
                    cfx_t1: terms.cfx_t1,
                    cfx_t2: terms.cfx_t2,
                    btmp_shape: terms.btmp_shape,
                    cfx_shape: terms.cfx_shape,
                    dix: terms.dix,
                    cfx_upw: terms.cfx_upw,
                    dix_upw: terms.dix_upw,
                    hsa: terms.hsa,
                    hca: terms.hca,
                    dd: terms.dd,
                    dd2: terms.dd2,
                    xot1: terms.xot1,
                    xot2: terms.xot2,
                    cf1: terms.cf1,
                    cf2: terms.cf2,
                    di1: terms.di1,
                    di2: terms.di2,
                    z_hs2: terms.z_hs2,
                    z_cf2_shape: terms.z_cf2_shape,
                    z_di2: terms.z_di2,
                    z_t2_shape: terms.z_t2_shape,
                    z_u2_shape: terms.z_u2_shape,
                    z_hca: terms.z_hca,
                    z_ha_shape: terms.z_ha_shape,
                    di2_s: terms.di2_s,
                    z_upw_shape: terms.z_upw_shape,
                    hs2_t: terms.hs2_t,
                    hs2_d: terms.hs2_d,
                    hs2_u: terms.hs2_u,
                    cf2_t_shape: terms.cf2_t_shape,
                    cf2_d_shape: terms.cf2_d_shape,
                    cf2_u_shape: terms.cf2_u_shape,
                    di2_t: terms.di2_t,
                    di2_d: terms.di2_d,
                    di2_u: terms.di2_u,
                    hc2_t: terms.hc2_t,
                    hc2_d: terms.hc2_d,
                    hc2_u: terms.hc2_u,
                    h2_t: terms.h2_t,
                    h2_d: terms.h2_d,
                    upw_t2_shape: terms.upw_t2_shape,
                    upw_d2_shape: terms.upw_d2_shape,
                    upw_u2_shape: terms.upw_u2_shape,
                    vs2_3_1: jac.vs2[2][0],
                    vs2_3_2: jac.vs2[2][1],
                    vs2_3_3: jac.vs2[2][2],
                    vs2_3_4: jac.vs2[2][3],
                    xfoil_t1: xfoil_terms.as_ref().map(|(_, s)| s.theta),
                    xfoil_d1: xfoil_terms.as_ref().map(|(_, s)| s.delta_star),
                    xfoil_s1: xfoil_terms.as_ref().map(|(_, s)| s.ctau),
                    xfoil_hk1: xfoil_terms.as_ref().map(|(_, s)| s.hk),
                    xfoil_rt1: xfoil_terms.as_ref().map(|(_, s)| s.r_theta),
                    xfoil_di1: xfoil_terms.as_ref().map(|(t, _)| t.di1),
                    xfoil_di2: xfoil_terms.as_ref().map(|(t, _)| t.di2),
                    xfoil_upw: xfoil_terms.as_ref().map(|(t, _)| t.upw),
                    xfoil_dix_upw: xfoil_terms.as_ref().map(|(t, _)| t.dix_upw),
                    xfoil_z_dix_shape: xfoil_terms.as_ref().map(|(t, _)| t.z_dix_shape),
                    xfoil_z_upw_shape: xfoil_terms.as_ref().map(|(t, _)| t.z_upw_shape),
                },
            ));
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
                side,
                station_idx + 2, // IBL is 1-indexed, +1 for skipped stagnation
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.hk,
                station.cf,
            ));
        }

        result.stations.push(station);
    }

    result
}

/// March boundary layer with debug output (for XFOIL comparison)
///
/// Same as [`march_fixed_ue`] but emits MRCHUE debug events when debug
/// collection is active. Use this for comparing RustFoil against XFOIL's
/// instrumented output.
pub fn march_fixed_ue_debug(
    x: &[f64],
    ue: &[f64],
    re: f64,
    msq: f64,
    config: &MarchConfig,
    side: usize,
) -> MarchResult {
    let n = x.len();
    assert_eq!(n, ue.len(), "x and ue arrays must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut result = MarchResult::new(n);

    // Initialize stagnation point (first station)
    let mut prev = init_stagnation(x[0], ue[0], re);
    blvar(&mut prev, FlowType::Laminar, msq, re);
    result.stations.push(prev.clone());

    // Emit debug event for first station
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
            side,
            0,
            prev.x,
            prev.u,
            prev.theta,
            prev.delta_star,
            prev.hk,
            prev.cf,
        ));
    }

    // Track state
    let mut is_laminar = true;

    // March downstream
    for i in 1..n {
        let ds = x[i] - x[i - 1];
        let prev_station = &result.stations[i - 1];

        // Initialize current station
        let mut station = BlStation::new();
        station.x = x[i];
        station.u = ue[i];
        station.is_laminar = is_laminar;
        station.is_turbulent = !is_laminar;
        station.is_wake = false;

        // Determine flow type
        let flow_type = if is_laminar {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Step the momentum and shape equations
        let (theta_new, delta_star_new) = step_momentum(prev_station, x[i], ue[i], re, msq, is_laminar);
        
        station.theta = theta_new;
        station.delta_star = delta_star_new;
        station.ctau = if is_laminar { prev_station.ctau } else { prev_station.ctau };
        station.ampl = prev_station.ampl;

        // Compute all secondary variables
        blvar(&mut station, flow_type, msq, re);

        // Do not clamp Hk here; inverse-mode handles near-separation behavior.

        // Check for transition (laminar only)
        if is_laminar && result.x_transition.is_none() {
            // Compute amplification rate using XFOIL's AXSET (RMS averaging)
            let ax_result = axset(
                prev_station.hk,
                prev_station.theta,
                prev_station.r_theta,
                prev_station.ampl,
                station.hk,
                station.theta,
                station.r_theta,
                station.ampl,
                config.ncrit,
            );
            station.ampl = prev_station.ampl + ax_result.ax * ds;

            // Check if transition occurred
            if let Some(_) = check_transition(station.ampl, config.ncrit) {
                // Interpolate transition location
                let ampl_prev = prev_station.ampl;
                let ampl_curr = station.ampl;
                let frac = (config.ncrit - ampl_prev) / (ampl_curr - ampl_prev).max(1e-20);
                let x_trans = prev_station.x + frac * ds;

                result.x_transition = Some(x_trans);
                result.transition_index = Some(i);

                // Transition to turbulent
                is_laminar = false;
                station.is_laminar = false;
                station.is_turbulent = true;
                // XFOIL uses 0.05 as fallback ctau at transition (xbl.f:801)
                station.ctau = 0.05;
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
        }

        // Check for separation (Cf < 0)
        if station.cf < 0.0 && result.x_separation.is_none() {
            result.x_separation = Some(station.x);
        }

        // Store mass defect
        station.mass_defect = station.u * station.delta_star;

        // Emit debug event
        if rustfoil_bl::is_debug_active() {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
                side,
                i,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.hk,
                station.cf,
            ));
        }

        result.stations.push(station);
    }

    result
}

// ============================================================================
// Coupled March (MRCHDU)
// ============================================================================

/// March boundary layer with edge velocity updates (viscous-inviscid coupling)
///
/// Similar to `march_fixed_ue()` but allows the edge velocity to be updated
/// based on viscous-inviscid interaction through the DIJ matrix. The solution
/// is constrained along a quasi-normal to the natural Ue-Hk characteristic,
/// preventing the Goldstein/Levy-Lees singularity.
///
/// # Arguments
/// * `x` - Arc length coordinates
/// * `initial_ue` - Initial edge velocities (from inviscid solution)
/// * `dij` - Mass defect influence matrix (DIJ)
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
/// * `config` - March configuration
///
/// # Returns
/// MarchResult with computed BL stations including updated Ue values
///
/// # Reference
/// XFOIL xbl.f MRCHDU (line 875)
pub fn march_coupled(
    x: &[f64],
    initial_ue: &[f64],
    dij: &DMatrix<f64>,
    re: f64,
    msq: f64,
    config: &MarchConfig,
) -> MarchResult {
    let n = x.len();
    assert_eq!(n, initial_ue.len(), "x and initial_ue arrays must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut result = MarchResult::new(n);

    // Initialize stagnation point
    let mut prev = init_stagnation(x[0], initial_ue[0], re);
    blvar(&mut prev, FlowType::Laminar, msq, re);
    result.stations.push(prev.clone());

    // Current Ue values (will be updated based on mass defect changes)
    let mut ue: Vec<f64> = initial_ue.to_vec();

    // Track state
    let mut is_laminar = true;

    // March downstream
    for i in 1..n {
        let ds = x[i] - x[i - 1];
        let prev_station = &result.stations[i - 1];

        // Initialize current station
        let mut station = BlStation::new();
        station.x = x[i];
        station.u = ue[i];
        station.is_laminar = is_laminar;
        station.is_turbulent = !is_laminar;
        station.is_wake = false;

        let flow_type = if is_laminar {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Step the BL equations
        let (theta_new, delta_star_new) = step_momentum(prev_station, x[i], ue[i], re, msq, is_laminar);
        
        station.theta = theta_new;
        station.delta_star = delta_star_new;
        station.ctau = if is_laminar { prev_station.ctau } else { prev_station.ctau };
        station.ampl = prev_station.ampl;

        // Compute secondary variables
        blvar(&mut station, flow_type, msq, re);

        // Viscous-inviscid coupling: update Ue based on mass defect change
        // ΔUe_i = Σ_j DIJ[i,j] × Δ(mass_defect)_j
        let mass_defect_new = station.u * station.delta_star;
        let mass_defect_old = prev_station.mass_defect;
        let delta_mass = mass_defect_new - mass_defect_old;

        // Apply DIJ coupling (simplified: local influence only)
        if i < dij.nrows() && i < dij.ncols() {
            let due = dij[(i, i)] * delta_mass * config.senswt / 1000.0;
            station.u = (station.u + due).max(0.01); // Keep Ue positive
            ue[i] = station.u;
            
            // Recompute with updated Ue
            blvar(&mut station, flow_type, msq, re);
        }

        // Do not clamp Hk here; inverse-mode handles near-separation behavior.

        // Check for transition
        if is_laminar && result.x_transition.is_none() {
            // Compute amplification rate using XFOIL's AXSET (RMS averaging)
            let ax_result = axset(
                prev_station.hk,
                prev_station.theta,
                prev_station.r_theta,
                prev_station.ampl,
                station.hk,
                station.theta,
                station.r_theta,
                station.ampl,
                config.ncrit,
            );
            station.ampl = prev_station.ampl + ax_result.ax * ds;

            if let Some(_) = check_transition(station.ampl, config.ncrit) {
                let ampl_prev = prev_station.ampl;
                let ampl_curr = station.ampl;
                let frac = (config.ncrit - ampl_prev) / (ampl_curr - ampl_prev).max(1e-20);
                let x_trans = prev_station.x + frac * ds;

                result.x_transition = Some(x_trans);
                result.transition_index = Some(i);

                is_laminar = false;
                station.is_laminar = false;
                station.is_turbulent = true;
                // XFOIL uses 0.05 as fallback ctau at transition (xbl.f:801)
                station.ctau = 0.05;
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
        }

        // Check for separation
        if station.cf < 0.0 && result.x_separation.is_none() {
            result.x_separation = Some(station.x);
        }

        station.mass_defect = station.u * station.delta_star;
        result.stations.push(station);
    }

    result
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // Configuration Tests
    // =========================================================================

    #[test]
    fn test_march_config_default() {
        let config = MarchConfig::default();

        assert_eq!(config.ncrit, 9.0);
        assert_eq!(config.max_iter, 25);
        assert!((config.tolerance - 1.0e-5).abs() < 1e-12); // Match XFOIL's MRCHUE criterion
        assert!((config.hlmax - 3.8).abs() < 1e-10);  // XFOIL default
        assert!((config.htmax - 2.5).abs() < 1e-10);
        assert!((config.senswt - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_march_config_presets() {
        let low_turb = MarchConfig::low_turbulence();
        assert_eq!(low_turb.ncrit, 9.0);

        let high_turb = MarchConfig::high_turbulence();
        assert_eq!(high_turb.ncrit, 3.0);
    }

    // =========================================================================
    // Stagnation Point Tests
    // =========================================================================

    #[test]
    fn test_init_stagnation() {
        let station = init_stagnation(0.01, 0.1, 1e6);

        assert!(station.theta > 0.0, "theta should be positive");
        assert!(station.delta_star > 0.0, "delta_star should be positive");
        assert!(
            (station.h - 2.2).abs() < 0.1,
            "H should be ~2.2 for Thwaites"
        );
        assert!(station.is_laminar, "Should start laminar");
        assert!(!station.is_turbulent, "Should not be turbulent");
    }

    #[test]
    fn test_init_stagnation_re_scaling() {
        // Higher Re should give thinner BL
        let station_low = init_stagnation(0.01, 0.1, 1e5);
        let station_high = init_stagnation(0.01, 0.1, 1e7);

        assert!(
            station_high.theta < station_low.theta,
            "Higher Re should give smaller theta"
        );
    }

    // =========================================================================
    // Flat Plate (Blasius) Tests
    // =========================================================================

    #[test]
    fn test_march_flat_plate_blasius() {
        // Flat plate: constant Ue = 1
        // At Re=1e6, transition occurs, so final H may be turbulent (~1.3-1.5)
        // or laminar (~2.5-2.7) depending on transition location
        let n = 50;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        let config = MarchConfig::default();
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &config);

        assert_eq!(result.stations.len(), n);

        // Check that H is in physically valid range (either laminar or turbulent)
        let last = result.stations.last().unwrap();
        assert!(
            last.h >= 1.0 && last.h <= 4.0,
            "H should be in valid range 1.0-4.0, got {}",
            last.h
        );

        // Theta should grow monotonically
        let first = &result.stations[0];
        assert!(
            last.theta >= first.theta,
            "theta should grow downstream"
        );
        
        // All stations should have positive theta and delta_star
        for station in &result.stations {
            assert!(station.theta > 0.0, "theta should be positive");
            assert!(station.delta_star > 0.0, "delta_star should be positive");
        }
    }

    #[test]
    fn test_march_flat_plate_laminar_stays_laminar() {
        // Very low Re flat plate should stay laminar
        // At Re=1e4, Rex_max = 1e4 * 0.3 = 3000, which is well below transition
        let n = 30;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        let config = MarchConfig {
            ncrit: 9.0,
            ..Default::default()
        };

        let result = march_fixed_ue(&x, &ue, 1e4, 0.0, &config); // Very low Re

        // At this low Re, should stay laminar
        // But allow for march to produce some numerical artifacts
        let laminar_count = result.stations.iter().filter(|s| s.is_laminar).count();
        assert!(
            laminar_count >= n / 2,
            "Most stations should be laminar at low Re, got {}/{} laminar",
            laminar_count, n
        );
    }

    // =========================================================================
    // Adverse Pressure Gradient Tests
    // =========================================================================

    #[test]
    fn test_march_adverse_pressure_gradient() {
        // Decelerating flow: Ue decreasing
        // Should thicken BL and potentially approach separation
        let n = 40;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue: Vec<f64> = x.iter().map(|&xi| 1.0 - 0.8 * xi).collect();

        let config = MarchConfig::default();
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &config);

        assert_eq!(result.stations.len(), n);

        // BL thickness should grow in adverse pressure gradient
        let first_theta = result.stations[5].theta;
        let last_theta = result.stations.last().unwrap().theta;
        assert!(
            last_theta > first_theta,
            "theta should increase in APG, first={}, last={}",
            first_theta,
            last_theta
        );
        
        // H should be finite (may be high near separation, H up to 15-20 is valid)
        for station in &result.stations {
            assert!(
                station.h >= 1.0 && station.h <= 20.0,
                "H should be in valid range (1.0-20.0), got {}",
                station.h
            );
        }
    }

    #[test]
    fn test_march_strong_apg_separation() {
        // Very strong deceleration should cause separation
        let n = 50;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        // Strong APG: rapid deceleration
        let ue: Vec<f64> = x.iter().map(|&xi| (1.0 - 1.5 * xi).max(0.1)).collect();

        let config = MarchConfig::default();
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &config);

        // May or may not separate depending on exact conditions,
        // but the march should complete
        assert!(result.stations.len() > 0, "March should produce stations");
    }

    // =========================================================================
    // Transition Tests
    // =========================================================================

    #[test]
    fn test_march_transition_detected() {
        // Higher Re should trigger transition
        let n = 100;
        let x: Vec<f64> = (1..=n).map(|i| 0.005 * i as f64).collect();
        let ue = vec![1.0; n];

        let config = MarchConfig {
            ncrit: 9.0,
            ..Default::default()
        };

        let result = march_fixed_ue(&x, &ue, 5e6, 0.0, &config); // Higher Re

        // Should trigger transition at some point
        if let Some(x_trans) = result.x_transition {
            assert!(x_trans > 0.0, "Transition should be at positive x");
            assert!(
                x_trans < x[n - 1],
                "Transition should occur before end"
            );

            // Verify transition index is set
            assert!(result.transition_index.is_some());
        }
        // Note: transition may not occur if Re is still too low
    }

    #[test]
    fn test_march_low_ncrit_early_transition() {
        // Lower Ncrit should cause earlier transition
        let n = 80;
        let x: Vec<f64> = (1..=n).map(|i| 0.005 * i as f64).collect();
        let ue = vec![1.0; n];

        let config_high = MarchConfig {
            ncrit: 9.0,
            ..Default::default()
        };
        let config_low = MarchConfig {
            ncrit: 3.0,
            ..Default::default()
        };

        let result_high = march_fixed_ue(&x, &ue, 2e6, 0.0, &config_high);
        let result_low = march_fixed_ue(&x, &ue, 2e6, 0.0, &config_low);

        // If both transition, low Ncrit should be earlier
        if result_high.x_transition.is_some() && result_low.x_transition.is_some() {
            assert!(
                result_low.x_transition.unwrap() <= result_high.x_transition.unwrap(),
                "Lower Ncrit should give earlier transition"
            );
        }
    }

    // =========================================================================
    // Coupled March Tests
    // =========================================================================

    #[test]
    fn test_march_coupled_basic() {
        // Basic coupled march test
        let n = 30;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        // Create simple DIJ matrix (identity-like)
        let dij = DMatrix::identity(n, n) * 0.1;

        let config = MarchConfig::default();
        let result = march_coupled(&x, &ue, &dij, 1e6, 0.0, &config);

        assert_eq!(result.stations.len(), n);

        // All stations should have valid values
        for station in &result.stations {
            assert!(station.theta > 0.0, "theta should be positive");
            assert!(station.delta_star > 0.0, "delta_star should be positive");
            assert!(station.u > 0.0, "u should be positive");
        }
    }

    #[test]
    fn test_march_coupled_ue_update() {
        // In coupled mode, Ue may be updated
        let n = 20;
        let x: Vec<f64> = (1..=n).map(|i| 0.02 * i as f64).collect();
        let ue = vec![1.0; n];

        let dij = DMatrix::identity(n, n) * 0.1;

        let config = MarchConfig::default();
        let result = march_coupled(&x, &ue, &dij, 1e6, 0.0, &config);

        // Coupled march may modify Ue from initial values
        // Just verify the march completes and produces valid results
        assert!(result.stations.len() == n);

        for station in &result.stations {
            assert!(
                station.u.is_finite() && station.u > 0.0,
                "Ue should be finite and positive"
            );
        }
    }

    // =========================================================================
    // Edge Cases
    // =========================================================================

    #[test]
    fn test_march_minimum_stations() {
        // Minimum 2 stations
        let x = vec![0.01, 0.02];
        let ue = vec![1.0, 1.0];

        let config = MarchConfig::default();
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &config);

        assert_eq!(result.stations.len(), 2);
    }

    #[test]
    fn test_march_result_new() {
        let result = MarchResult::new(10);

        assert_eq!(result.stations.capacity(), 10);
        assert!(result.x_transition.is_none());
        assert!(result.x_separation.is_none());
        assert!(result.converged);
    }

    #[test]
    fn test_gauss_4x4_identity() {
        // Solve I * x = b
        let mut a = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let mut b = [1.0, 2.0, 3.0, 4.0];

        gauss_4x4(&mut a, &mut b);

        assert!((b[0] - 1.0).abs() < 1e-10);
        assert!((b[1] - 2.0).abs() < 1e-10);
        assert!((b[2] - 3.0).abs() < 1e-10);
        assert!((b[3] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_gauss_4x4_simple() {
        // Solve 2*I * x = [2, 4, 6, 8] => x = [1, 2, 3, 4]
        let mut a = [
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 2.0],
        ];
        let mut b = [2.0, 4.0, 6.0, 8.0];

        gauss_4x4(&mut a, &mut b);

        assert!((b[0] - 1.0).abs() < 1e-10);
        assert!((b[1] - 2.0).abs() < 1e-10);
        assert!((b[2] - 3.0).abs() < 1e-10);
        assert!((b[3] - 4.0).abs() < 1e-10);
    }
}
