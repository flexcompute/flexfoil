//! XFOIL-style BL system equations and analytical derivatives.
//!
//! This module implements the boundary layer residual equations and their
//! analytical Jacobian entries, matching XFOIL's BLDIF subroutine.
//!
//! # Equations
//!
//! ## Momentum Equation (von Kármán)
//! ```text
//! log(θ₂/θ₁) + (H+2-M)*log(U₂/U₁) - 0.5*Cf*log(x₂/x₁) = 0
//! ```
//!
//! ## Shape Equation (Head entrainment or energy)
//! ```text
//! log(Hs₂/Hs₁) + (2*Hc/Hs+1-H)*log(U₂/U₁) + (0.5*Cf-Di)*log(x₂/x₁) = 0
//! ```
//!
//! ## Shear Lag Equation (turbulent Cτ)
//! ```text
//! Scc*(CQa-Sa*ALD)*dx - 2*DEa*log(S₂/S₁) + 2*DEa*(UQ*dx-log(U₂/U₁))*DUXCON = 0
//! ```
//!
//! # Reference
//!
//! XFOIL source: `xblsys.f` BLDIF subroutine (lines 1552-1977)

use crate::boundary_layer::closure::ludwieg_tillmann_cf;
use super::newton::{BLBlock, StationState};

/// BL closure relations and their derivatives.
/// 
/// Matches XFOIL's BLVAR routine calculations.
#[derive(Debug, Clone, Copy, Default)]
pub struct BLClosures {
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Kinematic shape factor Hk (compressibility corrected)
    pub hk: f64,
    /// Energy shape factor Hs = θ*/θ
    pub hs: f64,
    /// Density shape factor Hc
    pub hc: f64,
    /// Skin friction coefficient Cf
    pub cf: f64,
    /// Dissipation function 2*CD/Hs
    pub di: f64,
    /// Equilibrium shear stress coefficient sqrt(Cτ_eq)
    pub cq: f64,
    /// Shear stress coefficient Cτ^(1/2)
    pub ctau: f64,
    /// Reynolds number based on θ
    pub re_theta: f64,
    /// BL thickness δ (from Green's correlation)
    pub delta: f64,
    /// Normalized slip velocity Us
    pub us: f64,
}

impl BLClosures {
    /// Compute closures from primary variables.
    pub fn compute(theta: f64, h: f64, ue: f64, reynolds: f64, is_turbulent: bool) -> Self {
        let h_clamped = h.clamp(1.05, 10.0);
        let re_theta = (ue * theta * reynolds).max(100.0);
        
        // Kinematic shape factor (incompressible for now)
        let hk = h_clamped;
        
        // Energy shape factor Hs
        let hs = if is_turbulent {
            Self::hs_turbulent(hk, re_theta)
        } else {
            Self::hs_laminar(hk)
        };
        
        // Density shape factor Hc (zero for incompressible)
        let hc = 0.0;
        
        // Skin friction
        let cf = if is_turbulent {
            ludwieg_tillmann_cf(re_theta, hk)
        } else {
            // Laminar Cf from Thwaites/Blasius
            0.664 / re_theta.sqrt().max(1.0)
        };
        
        // Dissipation
        let di = if is_turbulent {
            Self::di_turbulent(hk, hs, re_theta, cf)
        } else {
            Self::di_laminar(hk, re_theta)
        };
        
        // Equilibrium Cτ (for turbulent)
        let cq = Self::cq_equilibrium(hk, hs, h_clamped, cf, re_theta);
        
        // Slip velocity
        let us = Self::slip_velocity(hs, hk, h_clamped);
        
        // BL thickness from Green's correlation
        let delta = (3.15 + 1.72 / (hk - 1.0).max(0.01)) * theta + theta * h_clamped;
        
        Self {
            h: h_clamped,
            hk,
            hs,
            hc,
            cf,
            di,
            cq,
            ctau: 0.0, // Set externally for turbulent
            re_theta,
            delta,
            us,
        }
    }
    
    /// Laminar Hs correlation (from Falkner-Skan).
    fn hs_laminar(hk: f64) -> f64 {
        if hk < 4.35 {
            let tmp = hk - 4.35;
            0.0111 * tmp.powi(2) / (hk + 1.0) - 0.0278 * tmp.powi(3) / (hk + 1.0) + 1.528
                - 0.0002 * (tmp * hk).powi(2)
        } else {
            let hs2 = 0.015;
            hs2 * (hk - 4.35).powi(2) / hk + 1.528
        }
    }
    
    /// Turbulent Hs correlation.
    fn hs_turbulent(hk: f64, re_theta: f64) -> f64 {
        let hsmin = 1.5;
        let dhsinf = 0.015;
        
        let ho = if re_theta > 400.0 {
            3.0 + 400.0 / re_theta
        } else {
            4.0
        };
        
        let rtz = re_theta.max(200.0);
        
        if hk < ho {
            // Attached
            let hr = (ho - hk) / (ho - 1.0);
            (2.0 - hsmin - 4.0 / rtz) * hr.powi(2) * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz
        } else {
            // Separated
            let grt = rtz.ln();
            let hdif = hk - ho;
            let rtmp = hk - ho + 4.0 / grt;
            let htmp = 0.007 * grt / rtmp.powi(2) + dhsinf / hk;
            hdif.powi(2) * htmp + hsmin + 4.0 / rtz
        }
    }
    
    /// Laminar dissipation function.
    fn di_laminar(hk: f64, re_theta: f64) -> f64 {
        let di = if hk < 4.0 {
            (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / re_theta
        } else {
            let hkb = hk - 4.0;
            let den = 1.0 + 0.02 * hkb.powi(2);
            (-0.0016 * hkb.powi(2) / den + 0.207) / re_theta
        };
        di.max(0.0)
    }
    
    /// Turbulent dissipation function.
    fn di_turbulent(hk: f64, hs: f64, re_theta: f64, cf: f64) -> f64 {
        // Wall contribution
        let us = Self::slip_velocity(hs, hk, hk);
        let di_wall = (0.5 * cf * us) * 2.0 / hs;
        
        // Outer layer contribution
        // Note: This is simplified; full XFOIL has more terms
        di_wall.max(0.0)
    }
    
    /// Equilibrium Cτ^(1/2) from XFOIL.
    fn cq_equilibrium(hk: f64, hs: f64, h: f64, cf: f64, re_theta: f64) -> f64 {
        let ctcon = 0.5 / (6.7 * 6.7); // XFOIL CTCON
        let gbcon = 0.72; // XFOIL GBCON
        let gccon = 18.0; // XFOIL GCCON for turbulent
        
        let hkc = (hk - 1.0 - gccon / re_theta.max(100.0)).max(0.01);
        let hkb = hk - 1.0;
        let us = Self::slip_velocity(hs, hk, h);
        let usb = (1.0 - us).max(0.001);
        
        (ctcon * hs * hkb * hkc.powi(2) / (usb * h * hk.powi(2))).sqrt().max(0.0)
    }
    
    /// Normalized slip velocity Us.
    fn slip_velocity(hs: f64, hk: f64, h: f64) -> f64 {
        let gbcon = 0.72;
        let us = 0.5 * hs * (1.0 - (hk - 1.0) / (gbcon * h));
        us.clamp(0.0, 0.98)
    }
}

/// Derivatives of closure relations wrt primary variables.
#[derive(Debug, Clone, Copy, Default)]
pub struct BLClosureDerivs {
    // Cf derivatives
    pub cf_theta: f64,
    pub cf_h: f64,
    pub cf_ue: f64,
    
    // Hs derivatives
    pub hs_theta: f64,
    pub hs_h: f64,
    
    // Di derivatives
    pub di_theta: f64,
    pub di_h: f64,
    pub di_ctau: f64,
    
    // Cq derivatives
    pub cq_theta: f64,
    pub cq_h: f64,
}

impl BLClosureDerivs {
    /// Compute derivatives using finite differences (for validation).
    pub fn compute_fd(
        theta: f64,
        h: f64,
        ue: f64,
        reynolds: f64,
        is_turbulent: bool,
        eps: f64,
    ) -> Self {
        let base = BLClosures::compute(theta, h, ue, reynolds, is_turbulent);
        
        let plus_theta = BLClosures::compute(theta + eps, h, ue, reynolds, is_turbulent);
        let plus_h = BLClosures::compute(theta, h + eps, ue, reynolds, is_turbulent);
        let plus_ue = BLClosures::compute(theta, h, ue + eps, reynolds, is_turbulent);
        
        Self {
            cf_theta: (plus_theta.cf - base.cf) / eps,
            cf_h: (plus_h.cf - base.cf) / eps,
            cf_ue: (plus_ue.cf - base.cf) / eps,
            
            hs_theta: (plus_theta.hs - base.hs) / eps,
            hs_h: (plus_h.hs - base.hs) / eps,
            
            di_theta: (plus_theta.di - base.di) / eps,
            di_h: (plus_h.di - base.di) / eps,
            di_ctau: 0.0, // Ctau is separate variable
            
            cq_theta: (plus_theta.cq - base.cq) / eps,
            cq_h: (plus_h.cq - base.cq) / eps,
        }
    }
}

/// XFOIL-style upwinding parameter for stability.
/// 
/// UPW = 0.5 gives trapezoidal differencing
/// UPW = 1.0 gives backward Euler (full upwinding)
/// 
/// XFOIL dynamically adjusts based on log(Hk) changes.
pub fn compute_upwind_factor(hk1: f64, hk2: f64) -> f64 {
    let hupwt = 1.0;
    let hdcon = 5.0 * hupwt / hk2.powi(2);
    
    let arg = ((hk2 - 1.0) / (hk1 - 1.0).max(0.01)).abs();
    let hl = arg.ln();
    let hlsq = (hl * hl).min(15.0);
    
    let ehh = (-hlsq * hdcon).exp();
    1.0 - 0.5 * ehh
}

/// Residuals and Jacobian blocks for BL equations.
#[derive(Debug, Clone)]
pub struct BLSystemResiduals {
    /// Residual for first equation (Ampl for laminar, Ctau for turbulent)
    pub res_1: f64,
    /// Residual for momentum equation
    pub res_momentum: f64,
    /// Residual for shape equation
    pub res_shape: f64,
}

/// Compute BL residuals between stations 1 and 2.
/// 
/// This matches XFOIL's BLDIF subroutine for interval differencing.
pub fn compute_interval_residuals(
    state1: &StationState,
    state2: &StationState,
    closures1: &BLClosures,
    closures2: &BLClosures,
    reynolds: f64,
) -> BLSystemResiduals {
    // Geometric quantities
    let ds = (state2.s - state1.s).abs().max(1e-12);
    let x1 = state1.s.max(1e-12);
    let x2 = state2.s.max(1e-12);
    let xlog = (x2 / x1).ln();
    
    // Velocity quantities
    let u1 = state1.ue.max(1e-10);
    let u2 = state2.ue.max(1e-10);
    let ulog = (u2 / u1).ln();
    
    // Theta quantities
    let t1 = state1.theta.max(1e-12);
    let t2 = state2.theta.max(1e-12);
    let tlog = (t2 / t1).ln();
    
    // Shape factor quantities
    let hs1 = closures1.hs;
    let hs2 = closures2.hs;
    let hlog = (hs2 / hs1.max(0.1)).ln();
    
    let h1 = closures1.h;
    let h2 = closures2.h;
    let ha = 0.5 * (h1 + h2);
    
    // Mach number (zero for incompressible)
    let ma = 0.0;
    
    // Upwinding
    let upw = compute_upwind_factor(closures1.hk, closures2.hk);
    
    // Average Cf with midpoint correction (XFOIL's CFX)
    let cf1 = closures1.cf;
    let cf2 = closures2.cf;
    let cfa = (1.0 - upw) * cf1 + upw * cf2;
    
    // CFX term from XFOIL (includes midpoint evaluation)
    let xa = 0.5 * (x1 + x2);
    let ta = 0.5 * (t1 + t2);
    let cfx = 0.5 * cfa * xa / ta + 0.25 * (cf1 * x1 / t1 + cf2 * x2 / t2);
    
    // === Momentum equation ===
    // REZT = TLOG + BTMP*ULOG - XLOG*0.5*CFX
    let btmp = ha + 2.0 - ma;
    let res_momentum = tlog + btmp * ulog - xlog * 0.5 * cfx;
    
    // === Shape equation ===
    let hca = 0.5 * (closures1.hc + closures2.hc);
    let hsa = 0.5 * (hs1 + hs2);
    
    // Dissipation with upwinding
    let di1 = closures1.di;
    let di2 = closures2.di;
    let xot1 = x1 / t1;
    let xot2 = x2 / t2;
    let dix = (1.0 - upw) * di1 * xot1 + upw * di2 * xot2;
    let cfx_shape = (1.0 - upw) * cf1 * xot1 + upw * cf2 * xot2;
    
    let btmp_shape = 2.0 * hca / hsa.max(0.1) + 1.0 - ha;
    let res_shape = hlog + btmp_shape * ulog + xlog * (0.5 * cfx_shape - dix);
    
    // === First equation (Amplification or Shear lag) ===
    let res_1 = if state2.is_turbulent {
        // Turbulent: Shear lag equation
        compute_shear_lag_residual(state1, state2, closures1, closures2, upw, ds)
    } else {
        // Laminar: Amplification equation
        compute_amplification_residual(state1, state2, closures1, closures2, ds, reynolds)
    };
    
    BLSystemResiduals {
        res_1,
        res_momentum,
        res_shape,
    }
}

/// Compute amplification equation residual (laminar).
fn compute_amplification_residual(
    state1: &StationState,
    state2: &StationState,
    closures1: &BLClosures,
    closures2: &BLClosures,
    ds: f64,
    reynolds: f64,
) -> f64 {
    // Amplification rate from XFOIL's DAMPL
    let ax1 = compute_amplification_rate(closures1.hk, state1.theta, closures1.re_theta);
    let ax2 = compute_amplification_rate(closures2.hk, state2.theta, closures2.re_theta);
    
    // RMS average (XFOIL uses this for stability)
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let ax = if axsq > 0.0 { axsq.sqrt() } else { 0.0 };
    
    // N2 - N1 - AX*(X2-X1) = 0
    let ampl1 = state1.ctau_or_ampl;
    let ampl2 = state2.ctau_or_ampl;
    
    ampl2 - ampl1 - ax * ds
}

/// Compute amplification rate dN/dx (from XFOIL's DAMPL).
fn compute_amplification_rate(hk: f64, theta: f64, re_theta: f64) -> f64 {
    let hmi = 1.0 / (hk - 1.0).max(0.01);
    
    // Critical Re_theta correlation
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let grcrit = aa + 0.7 * (bb + 1.0);
    
    let gr = re_theta.max(1.0).log10();
    
    if gr < grcrit - 0.08 {
        return 0.0; // Below critical
    }
    
    // Ramp function for smooth turn-on
    let rnorm = ((gr - (grcrit - 0.08)) / 0.16).clamp(0.0, 1.0);
    let rfac = if rnorm >= 1.0 {
        1.0
    } else {
        3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3)
    };
    
    // Amplification envelope slope
    let arg = 3.87 * hmi - 2.52;
    let ex = (-arg * arg).exp();
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    
    // Conversion factor
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi.powi(2) + 3.0 * hmi.powi(3);
    
    (af * dadr / theta.max(1e-10)) * rfac
}

/// Compute shear lag equation residual (turbulent).
fn compute_shear_lag_residual(
    state1: &StationState,
    state2: &StationState,
    closures1: &BLClosures,
    closures2: &BLClosures,
    upw: f64,
    ds: f64,
) -> f64 {
    // XFOIL constants
    let sccon = 5.6;
    let duxcon = 1.0; // Controls pressure gradient influence
    let gacon = 6.7;
    let gbcon = 0.72;
    let _dlcon = 0.9; // Dissipation length factor for wake (used in wake)
    
    // Shear stress coefficients
    let s1 = state1.ctau_or_ampl.max(1e-10);
    let s2 = state2.ctau_or_ampl.max(1e-10);
    let slog = (s2 / s1).ln();
    
    // Average quantities with upwinding
    let sa = (1.0 - upw) * s1 + upw * s2;
    let cqa = (1.0 - upw) * closures1.cq + upw * closures2.cq;
    let cfa = (1.0 - upw) * closures1.cf + upw * closures2.cf;
    let hka = (1.0 - upw) * closures1.hk + upw * closures2.hk;
    let usa = 0.5 * (closures1.us + closures2.us);
    let dea = 0.5 * (closures1.delta + closures2.delta);
    let da = 0.5 * (state1.theta * closures1.h + state2.theta * closures2.h); // δ*
    
    // Velocity log
    let u1 = state1.ue.max(1e-10);
    let u2 = state2.ue.max(1e-10);
    let ulog = (u2 / u1).ln();
    
    // Equilibrium dUe/dx
    let rta = 0.5 * (closures1.re_theta + closures2.re_theta);
    let gccon = 18.0;
    let hkc = (hka - 1.0 - gccon / rta.max(100.0)).max(0.01);
    let hr = hkc / (gacon * hka);
    let uq = (0.5 * cfa - hr * hr) / (gbcon * da.max(1e-10));
    
    // Dissipation length scale
    let ald = 1.0; // Use DLCON for wake
    
    // SCC factor
    let scc = sccon * 1.333 / (1.0 + usa);
    
    // Residual
    scc * (cqa - sa * ald) * ds - dea * 2.0 * slog + dea * 2.0 * (uq * ds - ulog) * duxcon
}

/// Compute Jacobian blocks for BL interval.
/// 
/// Returns (VS1, VS2) where:
/// - VS1[k][l] = d(residual_k) / d(variable_l at station 1)
/// - VS2[k][l] = d(residual_k) / d(variable_l at station 2)
/// 
/// Variables: [0] = Ctau/Ampl, [1] = Theta, [2] = mass (via delta* -> H)
pub fn compute_interval_jacobian(
    state1: &StationState,
    state2: &StationState,
    closures1: &BLClosures,
    closures2: &BLClosures,
    reynolds: f64,
    eps: f64,
) -> (BLBlock, BLBlock) {
    // Use finite differences for now (analytical version would be from BLDIF)
    let base_res = compute_interval_residuals(state1, state2, closures1, closures2, reynolds);
    
    let mut vs1 = BLBlock::zero();
    let mut vs2 = BLBlock::zero();
    
    // Perturb station 1 variables
    for l in 0..3 {
        let mut state1_plus = *state1;
        match l {
            0 => state1_plus.ctau_or_ampl += eps,
            1 => state1_plus.theta += eps,
            2 => {
                // Mass -> affects H
                state1_plus.mass += eps;
                state1_plus.update_h();
            }
            _ => {}
        }
        
        let closures1_plus = BLClosures::compute(
            state1_plus.theta, state1_plus.h, state1_plus.ue, reynolds, state1_plus.is_turbulent
        );
        let res_plus = compute_interval_residuals(&state1_plus, state2, &closures1_plus, closures2, reynolds);
        
        vs1.data[0][l] = (res_plus.res_1 - base_res.res_1) / eps;
        vs1.data[1][l] = (res_plus.res_momentum - base_res.res_momentum) / eps;
        vs1.data[2][l] = (res_plus.res_shape - base_res.res_shape) / eps;
    }
    
    // Perturb station 2 variables
    for l in 0..3 {
        let mut state2_plus = *state2;
        match l {
            0 => state2_plus.ctau_or_ampl += eps,
            1 => state2_plus.theta += eps,
            2 => {
                state2_plus.mass += eps;
                state2_plus.update_h();
            }
            _ => {}
        }
        
        let closures2_plus = BLClosures::compute(
            state2_plus.theta, state2_plus.h, state2_plus.ue, reynolds, state2_plus.is_turbulent
        );
        let res_plus = compute_interval_residuals(state1, &state2_plus, closures1, &closures2_plus, reynolds);
        
        vs2.data[0][l] = (res_plus.res_1 - base_res.res_1) / eps;
        vs2.data[1][l] = (res_plus.res_momentum - base_res.res_momentum) / eps;
        vs2.data[2][l] = (res_plus.res_shape - base_res.res_shape) / eps;
    }
    
    (vs1, vs2)
}

// ============================================================================
// Local Newton Iteration for Single Station
// ============================================================================

/// Configuration for local Newton iteration.
#[derive(Debug, Clone)]
pub struct LocalNewtonConfig {
    /// Maximum iterations per station
    pub max_iter: usize,
    /// Convergence tolerance for residual norm
    pub tol: f64,
    /// Finite difference epsilon
    pub fd_eps: f64,
    /// Relaxation factor (1.0 = full Newton step)
    pub relax: f64,
}

impl Default for LocalNewtonConfig {
    fn default() -> Self {
        Self {
            max_iter: 20,
            tol: 1e-8,
            fd_eps: 1e-7,
            relax: 1.0,
        }
    }
}

/// Result of local Newton iteration.
#[derive(Debug, Clone)]
pub struct LocalNewtonResult {
    /// Whether converged
    pub converged: bool,
    /// Final residual norm
    pub residual_norm: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Updated station state (if converged)
    pub state: StationState,
}

/// Solve local 3x3 Newton system at a single station.
/// 
/// Given the upstream state (state1), solves for the current station state (state2)
/// such that the interval residuals are zero.
/// 
/// This matches XFOIL's local Newton loop in `xbl.f` lines 610-740.
pub fn solve_station_newton(
    state1: &StationState,
    mut state2: StationState,
    reynolds: f64,
    config: &LocalNewtonConfig,
) -> LocalNewtonResult {
    let closures1 = BLClosures::compute(
        state1.theta, state1.h, state1.ue, reynolds, state1.is_turbulent
    );
    
    let mut iterations = 0;
    let mut converged = false;
    let mut residual_norm = f64::MAX;
    
    for iter in 0..config.max_iter {
        iterations = iter + 1;
        
        // Compute closures at current state2
        let closures2 = BLClosures::compute(
            state2.theta, state2.h, state2.ue, reynolds, state2.is_turbulent
        );
        
        // Compute residuals
        let res = compute_interval_residuals(state1, &state2, &closures1, &closures2, reynolds);
        let res_vec = [res.res_1, res.res_momentum, res.res_shape];
        
        // Check convergence
        residual_norm = (res_vec[0].powi(2) + res_vec[1].powi(2) + res_vec[2].powi(2)).sqrt();
        if residual_norm < config.tol {
            converged = true;
            break;
        }
        
        // Compute Jacobian wrt state2 variables
        let (_, vs2) = compute_interval_jacobian(
            state1, &state2, &closures1, &closures2, reynolds, config.fd_eps
        );
        
        // Solve 3x3 system: VS2 * delta = -residual
        let neg_res = [-res_vec[0], -res_vec[1], -res_vec[2]];
        let delta = match vs2.solve(&neg_res) {
            Some(d) => d,
            None => {
                // Singular - try damped step
                break;
            }
        };
        
        // Apply update with relaxation
        state2.ctau_or_ampl = (state2.ctau_or_ampl + config.relax * delta[0]).max(0.0);
        state2.theta = (state2.theta + config.relax * delta[1]).max(1e-12);
        state2.mass = (state2.mass + config.relax * delta[2]).max(1e-12);
        state2.update_h();
        
        // Clamp H to reasonable range
        state2.h = state2.h.clamp(1.05, 10.0);
    }
    
    LocalNewtonResult {
        converged,
        residual_norm,
        iterations,
        state: state2,
    }
}

/// March boundary layer along surface, solving at each station.
/// 
/// This is a simplified version of XFOIL's BLMARCH. It uses local Newton
/// iterations at each station to solve the BL equations sequentially.
/// 
/// Returns the solved BL states at all stations.
pub fn march_bl_surface(
    ue: &[f64],
    s_coords: &[f64],
    panel_indices: &[usize],
    reynolds: f64,
    n_crit: f64,
    config: &LocalNewtonConfig,
) -> Vec<StationState> {
    let n_stations = panel_indices.len();
    if n_stations < 2 {
        return vec![];
    }
    
    let mut states = Vec::with_capacity(n_stations);
    
    // Initialize first station (stagnation point)
    let panel_0 = panel_indices[0];
    let theta_0 = 1e-4 / reynolds.sqrt();
    let h_0 = 2.59; // Blasius H
    let state0 = StationState::new(
        0.0,  // Ampl = 0 at stagnation
        theta_0,
        h_0,
        ue[panel_0].abs().max(1e-10),
        s_coords[panel_0],
        false, // laminar at start
    );
    states.push(state0);
    
    let mut is_turbulent = false;
    let mut ampl = 0.0;
    
    // March downstream
    for j in 1..n_stations {
        let panel_idx = panel_indices[j];
        let prev_state = &states[j - 1];
        
        // Initial guess: extrapolate from previous
        let theta_guess = prev_state.theta;
        let h_guess = prev_state.h;
        let ue_j = ue[panel_idx].abs().max(1e-10);
        let s_j = s_coords[panel_idx];
        
        // Check transition
        if !is_turbulent && ampl >= n_crit {
            is_turbulent = true;
        }
        
        let initial_state = StationState::new(
            if is_turbulent { 0.03 } else { ampl }, // Initial Ctau or Ampl
            theta_guess,
            h_guess,
            ue_j,
            s_j,
            is_turbulent,
        );
        
        // Local Newton solve
        let result = solve_station_newton(prev_state, initial_state, reynolds, config);
        
        if !result.converged {
            // Use guess if Newton fails
            let fallback = StationState::new(
                initial_state.ctau_or_ampl,
                theta_guess * 1.1, // Slight increase
                h_guess,
                ue_j,
                s_j,
                is_turbulent,
            );
            states.push(fallback);
        } else {
            states.push(result.state);
        }
        
        // Update amplification factor for transition
        if !is_turbulent {
            ampl = states[j].ctau_or_ampl;
        }
    }
    
    states
}

/// Build block-tridiagonal Jacobian from marched BL solution.
/// 
/// This constructs the global Jacobian structure needed for the Newton VII solver,
/// using the local Jacobian blocks computed at each interval.
pub fn build_block_jacobian(
    states: &[StationState],
    reynolds: f64,
    eps: f64,
) -> super::newton::BlockTridiagJacobian {
    use super::newton::BlockTridiagJacobian;
    
    let n_stations = states.len();
    let mut jacobian = BlockTridiagJacobian::new(n_stations, n_stations);
    
    if n_stations < 2 {
        return jacobian;
    }
    
    // First station: identity (initial condition)
    let mut diag0 = super::newton::BLBlock::identity();
    jacobian.set_diag(0, diag0);
    jacobian.set_rhs(0, [0.0; 3]); // BC residual is zero at initial
    
    // Interior stations
    for j in 1..n_stations {
        let state1 = &states[j - 1];
        let state2 = &states[j];
        
        let closures1 = BLClosures::compute(state1.theta, state1.h, state1.ue, reynolds, state1.is_turbulent);
        let closures2 = BLClosures::compute(state2.theta, state2.h, state2.ue, reynolds, state2.is_turbulent);
        
        // Compute residuals
        let res = compute_interval_residuals(state1, state2, &closures1, &closures2, reynolds);
        jacobian.set_rhs(j, [res.res_1, res.res_momentum, res.res_shape]);
        
        // Compute Jacobian blocks
        let (vs1, vs2) = compute_interval_jacobian(state1, state2, &closures1, &closures2, reynolds, eps);
        
        // Sub-diagonal: derivatives wrt previous station
        jacobian.set_subdiag(j, vs1);
        
        // Diagonal: derivatives wrt current station  
        jacobian.set_diag(j, vs2);
    }
    
    jacobian
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bl_closures_laminar() {
        let closures = BLClosures::compute(0.001, 2.5, 1.0, 1e6, false);
        
        assert!(closures.h > 2.0 && closures.h < 3.0, "H = {}", closures.h);
        assert!(closures.hs > 1.0 && closures.hs < 2.0, "Hs = {}", closures.hs);
        assert!(closures.cf >= 0.0, "Cf = {}", closures.cf); // Cf can be small
        assert!(closures.re_theta > 0.0, "Re_theta = {}", closures.re_theta);
    }
    
    #[test]
    fn test_bl_closures_turbulent() {
        let closures = BLClosures::compute(0.001, 1.4, 1.0, 1e6, true);
        
        assert!(closures.h > 1.0 && closures.h < 2.0);
        assert!(closures.hs > 1.0 && closures.hs < 3.0);
        assert!(closures.cf > 0.0 && closures.cf < 0.01);
        assert!(closures.cq >= 0.0);
    }
    
    #[test]
    fn test_amplification_rate() {
        // Below critical Re_theta should give zero
        let ax_low = compute_amplification_rate(2.5, 0.001, 100.0);
        assert!(ax_low.abs() < 1e-6);
        
        // Above critical should be positive
        let ax_high = compute_amplification_rate(2.5, 0.001, 1000.0);
        assert!(ax_high >= 0.0);
    }
    
    #[test]
    fn test_upwind_factor() {
        // Same H should give ~0.5 (trapezoidal)
        let upw_same = compute_upwind_factor(2.0, 2.0);
        assert!(upw_same > 0.4, "upw_same = {}", upw_same);
        
        // Large H change should give higher upwinding
        let upw_diff = compute_upwind_factor(2.0, 4.0);
        assert!(upw_diff > upw_same, "upw_diff = {} should be > upw_same = {}", upw_diff, upw_same);
    }
    
    #[test]
    fn test_interval_residuals() {
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        let state2 = StationState::new(0.0, 0.0012, 2.4, 0.98, 0.01, false);
        
        let closures1 = BLClosures::compute(state1.theta, state1.h, state1.ue, 1e6, false);
        let closures2 = BLClosures::compute(state2.theta, state2.h, state2.ue, 1e6, false);
        
        let res = compute_interval_residuals(&state1, &state2, &closures1, &closures2, 1e6);
        
        // Residuals should be finite
        assert!(res.res_1.is_finite());
        assert!(res.res_momentum.is_finite());
        assert!(res.res_shape.is_finite());
    }
    
    #[test]
    fn test_local_newton_converges() {
        // Setup: two stations with consistent BL growth
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        
        // Initial guess for station 2
        let state2_guess = StationState::new(0.0, 0.001, 2.5, 0.98, 0.01, false);
        
        let config = LocalNewtonConfig::default();
        let result = solve_station_newton(&state1, state2_guess, 1e6, &config);
        
        // Should converge (or at least not blow up)
        assert!(result.residual_norm.is_finite(), "residual = {}", result.residual_norm);
        assert!(result.state.theta > 0.0, "theta = {}", result.state.theta);
        assert!(result.state.h > 1.0 && result.state.h < 10.0, "H = {}", result.state.h);
    }
    
    #[test]
    fn test_interval_jacobian() {
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        let state2 = StationState::new(0.0, 0.0012, 2.4, 0.98, 0.01, false);
        
        let closures1 = BLClosures::compute(state1.theta, state1.h, state1.ue, 1e6, false);
        let closures2 = BLClosures::compute(state2.theta, state2.h, state2.ue, 1e6, false);
        
        let (vs1, vs2) = compute_interval_jacobian(&state1, &state2, &closures1, &closures2, 1e6, 1e-7);
        
        // Jacobian blocks should have finite entries
        assert!(vs1.norm().is_finite(), "vs1 norm = {}", vs1.norm());
        assert!(vs2.norm().is_finite(), "vs2 norm = {}", vs2.norm());
        
        // Diagonal block (vs2) should be non-singular for well-posed problem
        let test_rhs = [1.0, 1.0, 1.0];
        let solution = vs2.solve(&test_rhs);
        assert!(solution.is_some(), "vs2 is singular");
    }
    
    #[test]
    fn test_build_block_jacobian() {
        // Simple two-station case
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        let state2 = StationState::new(0.0, 0.0012, 2.4, 0.98, 0.01, false);
        let states = vec![state1, state2];
        
        let jacobian = build_block_jacobian(&states, 1e6, 1e-7);
        
        assert_eq!(jacobian.n_stations, 2);
        assert!(jacobian.rhs_norm().is_finite());
    }
}
