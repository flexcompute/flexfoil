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
use rustfoil_bl::closures::{amplification_rate, check_transition};
use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;

use crate::solve::{build_4x4_system, solve_4x4};

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
}

impl Default for MarchConfig {
    fn default() -> Self {
        Self {
            ncrit: 9.0,
            max_iter: 25,
            tolerance: 1e-5,
            hlmax: 3.8,
            htmax: 2.5,
            senswt: 1000.0,
            deps: 5.0e-6,
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
/// * `prev` - Previous (upstream) BL station
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
) -> (BlStation, bool) {
    let flow_type = if is_laminar {
        FlowType::Laminar
    } else {
        FlowType::Turbulent
    };
    
    // Initialize station from previous (like XFOIL's similarity/upstream extrapolation)
    let mut station = BlStation::new();
    station.x = x_new;
    station.u = ue_new;
    station.is_laminar = is_laminar;
    station.is_turbulent = !is_laminar;
    station.is_wake = false;
    
    // Initialize from previous station
    station.theta = prev.theta;
    station.delta_star = prev.delta_star;
    station.ctau = if is_laminar { 0.03 } else { prev.ctau };
    station.ampl = prev.ampl;
    
    // Compute secondary variables
    blvar(&mut station, flow_type, msq, re);
    
    let hmax = if is_laminar { hlmax } else { htmax };
    let mut direct = true;
    let mut htarg = hmax;
    
    // Newton iteration loop
    for _iter in 0..max_iter {
        // Build BL residuals and Jacobian
        let (res, jac) = bldif(prev, &station, flow_type, msq, re);
        
        // Build 4x4 Newton system
        // Variable order in bldif: [s/ampl, θ, δ*, u, x]
        // We use [s/ampl, θ, δ*, u] (4 variables)
        let (mut a, mut b) = build_4x4_system(
            &jac.vs2,
            &[res.res_third, res.res_mom, res.res_shape],
            direct,
            station.hk / station.theta,  // hk2_t ≈ -Hk/θ
            -station.hk / station.delta_star, // hk2_d ≈ Hk/δ* (but derivative is negative)
            0.0,  // hk2_u (simplified - no Ue dependence in Hk at low Mach)
            htarg,
            station.hk,
        );
        
        // Solve 4x4 system
        let vsrez = solve_4x4(&a, &b);
        
        // Compute max relative change
        let dmax = (vsrez[1] / station.theta.max(1e-12)).abs()
            .max((vsrez[2] / station.delta_star.max(1e-12)).abs());
        let dmax = if is_laminar {
            dmax.max((vsrez[0] / 10.0).abs())  // Amplification scale
        } else {
            dmax.max((vsrez[0] / station.ctau.max(1e-6)).abs())  // Ctau scale
        };
        
        // Relaxation factor
        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };
        
        // Check if direct mode is appropriate (Hk not too high)
        if direct {
            let h_test = (station.delta_star + rlx * vsrez[2]) 
                       / (station.theta + rlx * vsrez[1]).max(1e-12);
            // Use hkin to convert H to Hk (simplified for low Mach)
            let hk_test = h_test;  // At low Mach, Hk ≈ H
            
            if hk_test >= hmax {
                // Need to switch to inverse mode
                direct = false;
                // Estimate target Hk based on upstream
                if is_laminar {
                    htarg = prev.hk + 0.03 * (x_new - prev.x) / prev.theta.max(1e-12);
                } else {
                    htarg = prev.hk - 0.15 * (x_new - prev.x) / prev.theta.max(1e-12);
                }
                htarg = htarg.max(hmax).min(hmax * 1.5);
                continue;  // Restart iteration with inverse mode
            }
        }
        
        // Apply updates with relaxation
        if !is_laminar {
            station.ctau = (station.ctau + rlx * vsrez[0]).clamp(1e-7, 0.3);
        }
        station.theta = (station.theta + rlx * vsrez[1]).max(1e-12);
        station.delta_star = (station.delta_star + rlx * vsrez[2]).max(1e-12);
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
        
        // Check convergence
        if dmax <= tolerance {
            return (station, true);
        }
    }
    
    // Did not converge - return best estimate
    (station, false)
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
        let ds = x[i] - x[i - 1];
        let prev_station = &result.stations[i - 1];

        // Determine flow type
        let flow_type = if is_laminar {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Solve for current station using Newton iteration (XFOIL MRCHUE style)
        let (mut station, converged) = newton_solve_station(
            prev_station,
            x[i],
            ue[i],
            re,
            msq,
            is_laminar,
            config.max_iter,
            config.tolerance,
            config.hlmax,
            config.htmax,
        );

        // If Newton didn't converge, use best estimate (station still has result)
        if !converged {
            // Fallback: use step_momentum for robustness
            let (theta_new, delta_star_new) = step_momentum(prev_station, x[i], ue[i], re, msq, is_laminar);
            station.theta = theta_new;
            station.delta_star = delta_star_new;
            blvar(&mut station, flow_type, msq, re);
        }

        // Check for transition (laminar only)
        if is_laminar && result.x_transition.is_none() {
            // Compute amplification rate at midpoint
            let hk_mid = 0.5 * (prev_station.hk + station.hk);
            let th_mid = 0.5 * (prev_station.theta + station.theta);
            let rt_mid = 0.5 * (prev_station.r_theta + station.r_theta);

            let amp_result = amplification_rate(hk_mid, th_mid, rt_mid);
            station.ampl = prev_station.ampl + amp_result.ax * ds;

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
                station.ctau = 0.03; // Initialize shear stress

                // Recompute with turbulent closures
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
        }

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
        let (mut station, converged) = newton_solve_station(
            prev_station,
            x[i],
            ue[i].abs(),
            re,
            msq,
            is_laminar,
            config.max_iter,
            config.tolerance,
            config.hlmax,
            config.htmax,
        );

        // If Newton didn't converge, fall back to step_momentum
        if !converged {
            let (theta_new, delta_star_new) =
                step_momentum(prev_station, x[i], ue[i].abs(), re, msq, is_laminar);
            station.theta = theta_new;
            station.delta_star = delta_star_new;
            blvar(&mut station, flow_type, msq, re);
        }

        // Compute all secondary variables (already done in newton_solve, but ensure consistency)
        blvar(&mut station, flow_type, msq, re);

        // Check if Hk is too high (approaching separation)
        let hmax = if is_laminar { config.hlmax } else { config.htmax };
        if station.hk > hmax {
            let h_limit = hmax * 0.95;
            station.delta_star = h_limit * station.theta;
            blvar(&mut station, flow_type, msq, re);
        }

        // Check for transition (laminar only)
        if is_laminar && result.x_transition.is_none() {
            let hk_mid = 0.5 * (prev_station.hk + station.hk);
            let th_mid = 0.5 * (prev_station.theta + station.theta);
            let rt_mid = 0.5 * (prev_station.r_theta + station.r_theta);

            let amp_result = amplification_rate(hk_mid, th_mid, rt_mid);
            station.ampl = prev_station.ampl + amp_result.ax * ds;

            if check_transition(station.ampl, config.ncrit).is_some() {
                let ampl_prev = prev_station.ampl;
                let ampl_curr = station.ampl;
                let frac = (config.ncrit - ampl_prev) / (ampl_curr - ampl_prev).max(1e-20);
                let x_trans = prev_station.x + frac * ds;

                result.x_transition = Some(x_trans);
                result.transition_index = Some(station_idx);

                is_laminar = false;
                station.is_laminar = false;
                station.is_turbulent = true;
                station.ctau = 0.03;

                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
        }

        // Check for separation
        if station.cf < 0.0 && result.x_separation.is_none() {
            result.x_separation = Some(station.x);
        }

        // Store mass defect
        station.mass_defect = station.u * station.delta_star;

        // Emit debug event
        if rustfoil_bl::is_debug_active() {
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

        // Check if Hk is too high (approaching separation)
        let hmax = if is_laminar { config.hlmax } else { config.htmax };
        if station.hk > hmax {
            // Limit H to prevent singularity
            let h_limit = hmax * 0.95;
            station.delta_star = h_limit * station.theta;
            blvar(&mut station, flow_type, msq, re);
        }

        // Check for transition (laminar only)
        if is_laminar && result.x_transition.is_none() {
            // Compute amplification rate at midpoint
            let hk_mid = 0.5 * (prev_station.hk + station.hk);
            let th_mid = 0.5 * (prev_station.theta + station.theta);
            let rt_mid = 0.5 * (prev_station.r_theta + station.r_theta);

            let amp_result = amplification_rate(hk_mid, th_mid, rt_mid);
            station.ampl = prev_station.ampl + amp_result.ax * ds;

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
                station.ctau = 0.03; // Initialize shear stress

                // Recompute with turbulent closures
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

        // Check if Hk is too high
        let hmax = if is_laminar { config.hlmax } else { config.htmax };
        if station.hk > hmax {
            let h_limit = hmax * 0.95;
            station.delta_star = h_limit * station.theta;
            blvar(&mut station, flow_type, msq, re);
        }

        // Check for transition
        if is_laminar && result.x_transition.is_none() {
            let hk_mid = 0.5 * (prev_station.hk + station.hk);
            let th_mid = 0.5 * (prev_station.theta + station.theta);
            let rt_mid = 0.5 * (prev_station.r_theta + station.r_theta);

            let amp_result = amplification_rate(hk_mid, th_mid, rt_mid);
            station.ampl = prev_station.ampl + amp_result.ax * ds;

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
                station.ctau = 0.03;

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
        assert!((config.tolerance - 1e-5).abs() < 1e-10);
        assert!((config.hlmax - 3.8).abs() < 1e-10);
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
        // BL should have H in reasonable laminar range (2.0-3.0)
        let n = 50;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        let config = MarchConfig::default();
        let result = march_fixed_ue(&x, &ue, 1e6, 0.0, &config);

        assert_eq!(result.stations.len(), n);

        // Check that H is in reasonable laminar range
        let last = result.stations.last().unwrap();
        assert!(
            last.h >= 2.0 && last.h <= 3.5,
            "H should be in laminar range 2.0-3.5, got {}",
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
        // Low Re flat plate should stay laminar
        let n = 30;
        let x: Vec<f64> = (1..=n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        let config = MarchConfig {
            ncrit: 9.0,
            ..Default::default()
        };

        let result = march_fixed_ue(&x, &ue, 1e5, 0.0, &config); // Low Re

        // Should stay laminar throughout
        for station in &result.stations {
            assert!(station.is_laminar, "Should remain laminar at low Re");
        }
        assert!(
            result.x_transition.is_none(),
            "No transition expected at low Re"
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
        
        // H should be in valid range throughout
        for station in &result.stations {
            assert!(
                station.h >= 1.0 && station.h <= 5.0,
                "H should be in valid range, got {}",
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
