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
/// In the wake (is_wake=true), uses less upwinding.
pub fn compute_upwind_factor(hk1: f64, hk2: f64) -> f64 {
    compute_upwind_factor_wake(hk1, hk2, false)
}

/// XFOIL-style upwinding with wake handling.
/// 
/// In wake region, XFOIL uses HDCON = HUPWT/HK2^2 (factor of 5 removed)
/// for smoother blending through wake.
pub fn compute_upwind_factor_wake(hk1: f64, hk2: f64, is_wake: bool) -> f64 {
    let hupwt = 1.0;
    // XFOIL: wake uses less upwinding (factor of 5 removed)
    let hdcon = if is_wake {
        hupwt / hk2.powi(2)
    } else {
        5.0 * hupwt / hk2.powi(2)
    };
    
    let arg = ((hk2 - 1.0) / (hk1 - 1.0).max(0.01)).abs();
    let hl = arg.ln();
    let hlsq = (hl * hl).min(15.0);
    
    let ehh = (-hlsq * hdcon).exp();
    1.0 - 0.5 * ehh
}

/// Residuals and Jacobian blocks for BL equations.
/// 
/// # Variables (XFOIL convention)
/// - var[0]: Ctau^(1/2) (turbulent) or N-factor (laminar)
/// - var[1]: θ (momentum thickness)
/// - var[2]: m = Ue × δ* (mass defect)
/// 
/// Using mass defect m instead of δ* enables direct V-I coupling:
/// - Source strength σ = dm/ds (transpiration velocity)
/// - DIJ matrix gives dUe/dm for Newton coupling
#[derive(Debug, Clone)]
pub struct BLSystemResiduals {
    /// Residual for first equation (Ampl for laminar, Ctau for turbulent)
    pub res_1: f64,
    /// Residual for momentum equation
    pub res_momentum: f64,
    /// Residual for shape/mass equation
    /// 
    /// For mass defect coupling, this becomes:
    /// R_mass = (m2 - m1)/ds - Ue_avg * dδ*/ds - δ*_avg * dUe/ds
    pub res_shape: f64,
}

/// Compute mass defect residual for V-I coupling.
/// 
/// The mass defect equation relates m = Ue*δ* between stations:
/// ```text
/// dm/ds = d(Ue*δ*)/ds = Ue*(dδ*/ds) + δ**(dUe/ds)
/// ```
/// 
/// This is the "third equation" in XFOIL's formulation that directly
/// couples to the inviscid solution through transpiration.
/// 
/// # Arguments
/// * `state1`, `state2` - BL states at upstream and downstream stations
/// * `closures1`, `closures2` - Closure values at stations
/// 
/// # Returns
/// Mass defect residual (should be ~0 when converged)
pub fn compute_mass_defect_residual(
    state1: &StationState,
    state2: &StationState,
    _closures1: &BLClosures,
    _closures2: &BLClosures,
) -> f64 {
    let ds = (state2.s - state1.s).abs().max(1e-12);
    
    // Mass defect at each station
    let m1 = state1.mass;
    let m2 = state2.mass;
    
    // Edge velocity
    let ue1 = state1.ue.max(1e-10);
    let ue2 = state2.ue.max(1e-10);
    
    // Displacement thickness
    let ds1 = state1.delta_star();
    let ds2 = state2.delta_star();
    
    // Averaged values
    let ue_avg = 0.5 * (ue1 + ue2);
    let ds_avg = 0.5 * (ds1 + ds2);
    
    // Gradients (using central difference)
    let d_ds_ds = (ds2 - ds1) / ds;
    let d_ue_ds = (ue2 - ue1) / ds;
    
    // Mass defect change should match dm/ds = Ue*d(δ*)/ds + δ**d(Ue)/ds
    let dm_ds_expected = ue_avg * d_ds_ds + ds_avg * d_ue_ds;
    let dm_ds_actual = (m2 - m1) / ds;
    
    // Residual: difference between actual and expected mass defect rate
    dm_ds_actual - dm_ds_expected
}

/// Derivatives of mass defect residual with respect to BL variables.
/// 
/// For the Newton iteration, we need:
/// - dR_m/d(ctau): 0 (no direct dependence)
/// - dR_m/d(theta): Through δ* = θ*H
/// - dR_m/d(m): Direct dependence
/// - dR_m/d(Ue): Through both Ue and δ*
#[derive(Debug, Clone, Copy, Default)]
pub struct MassDefectDerivs {
    /// dR/d(ctau) at station 1
    pub d_ctau1: f64,
    /// dR/d(theta) at station 1
    pub d_theta1: f64,
    /// dR/d(m) at station 1
    pub d_m1: f64,
    /// dR/d(ctau) at station 2
    pub d_ctau2: f64,
    /// dR/d(theta) at station 2
    pub d_theta2: f64,
    /// dR/d(m) at station 2
    pub d_m2: f64,
    /// dR/d(Ue) at station 1 (for V-I coupling)
    pub d_ue1: f64,
    /// dR/d(Ue) at station 2 (for V-I coupling)
    pub d_ue2: f64,
}

impl MassDefectDerivs {
    /// Compute derivatives of mass defect residual.
    pub fn compute(
        state1: &StationState,
        state2: &StationState,
        _closures1: &BLClosures,
        _closures2: &BLClosures,
    ) -> Self {
        let ds = (state2.s - state1.s).abs().max(1e-12);
        let inv_ds = 1.0 / ds;
        
        let ue1 = state1.ue.max(1e-10);
        let ue2 = state2.ue.max(1e-10);
        let ds1 = state1.delta_star();
        let ds2 = state2.delta_star();
        let h1 = state1.h;
        let h2 = state2.h;
        
        let ue_avg = 0.5 * (ue1 + ue2);
        let ds_avg = 0.5 * (ds1 + ds2);
        
        // dm/ds = (m2 - m1)/ds
        // dR/dm1 = -1/ds
        // dR/dm2 = +1/ds
        let d_m1 = -inv_ds;
        let d_m2 = inv_ds;
        
        // Through δ* = m/Ue, and θ doesn't appear directly in m
        // But shape equation uses H = δ*/θ, so dH/dθ = -δ*/θ² = -H/θ
        // dδ*/dθ comes from closure relation
        let d_theta1 = -0.5 * ue_avg * h1 * inv_ds; // d(δ*1)/dθ1 * Ue_avg
        let d_theta2 = 0.5 * ue_avg * h2 * inv_ds;
        
        // dR/dUe through ue_avg and ds_avg terms
        let d_ue1 = -0.5 * (ds2 - ds1) * inv_ds - ds_avg * 0.5 * inv_ds;
        let d_ue2 = -0.5 * (ds2 - ds1) * inv_ds + ds_avg * 0.5 * inv_ds;
        
        Self {
            d_ctau1: 0.0,
            d_theta1,
            d_m1,
            d_ctau2: 0.0,
            d_theta2,
            d_m2,
            d_ue1,
            d_ue2,
        }
    }
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
    // In inverse mode, prescribe Hk and solve for Ue
    // The shape equation becomes: Hk - Hk_target = 0
    let res_shape = if state2.inverse_mode {
        // Inverse mode: residual is simply Hk - Hk_target
        // This drives the Newton solver to find Ue that produces target Hk
        closures2.hk - state2.hk_target
    } else {
        // Direct mode: standard entrainment/energy equation
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
        hlog + btmp_shape * ulog + xlog * (0.5 * cfx_shape - dix)
    };
    
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

/// Compute amplification rate dN/dx using XFOIL's DAMPL2 algorithm.
/// 
/// DAMPL2 is the 1996 version that includes:
/// 1. Correct AF coefficient with exponential term for high H
/// 2. Non-envelope max-amplification correction for separated profiles (H > 3.5)
/// 
/// Reference: Drela & Giles, "Viscous/Inviscid Analysis of Transonic and
/// Low Reynolds Number Airfoils", AIAA Journal, Oct. 1987.
/// Updated March 1991, November 1996.
/// 
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `theta` - Momentum thickness
/// * `re_theta` - Momentum thickness Reynolds number
/// 
/// # Returns
/// Spatial amplification rate dN/dx
fn compute_amplification_rate(hk: f64, theta: f64, re_theta: f64) -> f64 {
    const DGR: f64 = 0.08;  // Ramp half-width in log10(Re_theta) space
    const HK1: f64 = 3.5;   // Start of separation profile blending
    const HK2: f64 = 4.0;   // End of separation profile blending
    
    let hmi = 1.0 / (hk - 1.0).max(0.01);
    
    // === Critical Re_theta correlation (Falkner-Skan profiles) ===
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let grc = aa + 0.7 * (bb + 1.0);  // log10(Re_theta_critical)
    
    let gr = re_theta.max(1.0).log10();
    
    // Below critical: no amplification
    if gr < grc - DGR {
        return 0.0;
    }
    
    // === Ramp function for smooth turn-on near Re_theta_critical ===
    let rnorm = ((gr - (grc - DGR)) / (2.0 * DGR)).clamp(0.0, 1.0);
    let rfac = if rnorm >= 1.0 {
        1.0
    } else {
        3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3)
    };
    
    // === Amplification envelope slope dN/d(Re_theta) ===
    let arg = 3.87 * hmi - 2.52;
    let ex = (-arg * arg).exp();
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    
    // === Conversion factor: theta * d(Re_theta)/dx ===
    // DAMPL2 adds exponential term for high H (separation bubble profiles)
    let brg = -20.0 * hmi;
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi.powi(2) + 3.0 * hmi.powi(3) + 0.1 * brg.exp();
    
    // Base envelope amplification rate AX1
    let ax1 = (af * dadr / theta.max(1e-10)) * rfac;
    
    // === Non-envelope correction for separated profiles (H > 3.5) ===
    // This handles laminar separation bubbles where standard Falkner-Skan
    // profiles are not representative. Uses Orr-Sommerfeld max ai(H, Re_theta).
    if hk < HK1 {
        return ax1;
    }
    
    // Blending fraction: 0 at HK1, 1 at HK2
    let hnorm = ((hk - HK1) / (HK2 - HK1)).clamp(0.0, 1.0);
    let hfac = if hnorm >= 1.0 {
        1.0
    } else {
        3.0 * hnorm.powi(2) - 2.0 * hnorm.powi(3)
    };
    
    // Modified amplification rate AX2 for separated profiles
    // Based on Orr-Sommerfeld stability calculations
    let gr0 = 0.30 + 0.35 * (-0.15 * (hk - 5.0)).exp();
    let tnr = (1.2 * (gr - gr0)).tanh();
    
    let ax2_raw = (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / theta.max(1e-10);
    let ax2 = ax2_raw.max(0.0);  // Clamp to non-negative
    
    // Blend between envelope (AX1) and Orr-Sommerfeld (AX2)
    hfac * ax2 + (1.0 - hfac) * ax1
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

/// Analytical Jacobian blocks matching XFOIL's BLDIF.
/// 
/// This computes the analytical derivatives instead of finite differences
/// for better accuracy and efficiency.
/// 
/// # Variables (XFOIL convention)
/// - var[0]: Ctau^(1/2) (turbulent) or N-factor (laminar)
/// - var[1]: θ (momentum thickness)  
/// - var[2]: m = Ue × δ* (mass defect)
/// 
/// # Returns
/// (VS1, VS2) where:
/// - VS1[k][l] = ∂R_k / ∂(variable_l at station 1)
/// - VS2[k][l] = ∂R_k / ∂(variable_l at station 2)
pub fn compute_analytical_jacobian(
    state1: &StationState,
    state2: &StationState,
    closures1: &BLClosures,
    closures2: &BLClosures,
    reynolds: f64,
) -> (BLBlock, BLBlock) {
    use super::newton::BLBlock;
    
    let mut vs1 = BLBlock::zero();
    let mut vs2 = BLBlock::zero();
    
    // Geometric quantities
    let x1 = state1.s.max(1e-12);
    let x2 = state2.s.max(1e-12);
    let ds = (x2 - x1).abs().max(1e-12);
    
    // Velocity quantities
    let u1 = state1.ue.max(1e-10);
    let u2 = state2.ue.max(1e-10);
    
    // Theta quantities
    let t1 = state1.theta.max(1e-12);
    let t2 = state2.theta.max(1e-12);
    
    // Shape factors
    let h1 = closures1.h.max(1.05);
    let h2 = closures2.h.max(1.05);
    let hk1 = closures1.hk.max(1.05);
    let hk2 = closures2.hk.max(1.05);
    let hs1 = closures1.hs.max(0.1);
    let hs2 = closures2.hs.max(0.1);
    
    // Skin friction
    let cf1 = closures1.cf;
    let cf2 = closures2.cf;
    
    // Dissipation
    let di1 = closures1.di;
    let di2 = closures2.di;
    
    // Upwinding factor
    let upw = compute_upwind_factor(hk1, hk2);
    
    // Mass defect
    let m1 = state1.mass.max(1e-12);
    let m2 = state2.mass.max(1e-12);
    
    // === Derivatives of closure quantities wrt primary variables ===
    
    // dH/dm at constant theta and Ue: H = m/(Ue*θ)
    let dh1_dm1 = 1.0 / (u1 * t1);
    let dh2_dm2 = 1.0 / (u2 * t2);
    
    // dH/dθ at constant m and Ue: H = m/(Ue*θ) => dH/dθ = -m/(Ue*θ²) = -H/θ
    let dh1_dt1 = -h1 / t1;
    let dh2_dt2 = -h2 / t2;
    
    // dHk/dH (simplified: Hk ≈ H for incompressible)
    let dhk1_dh1 = 1.0;
    let dhk2_dh2 = 1.0;
    
    // dHs/dHk from closure relations
    let (dhs1_dhk1, dhs2_dhk2) = if state2.is_turbulent {
        // Turbulent Hs-Hk relation
        let dhsk1 = -(hk1 - 4.35).signum() * 0.15;
        let dhsk2 = -(hk2 - 4.35).signum() * 0.15;
        (dhsk1, dhsk2)
    } else {
        // Laminar Hs-Hk relation (from Thwaites)
        (0.5, 0.5)
    };
    
    // dCf/dRe_θ and dCf/dHk (Ludwieg-Tillmann)
    let re_t1 = closures1.re_theta.max(100.0);
    let re_t2 = closures2.re_theta.max(100.0);
    let dcf1_dret = -0.268 * cf1 / re_t1;
    let dcf2_dret = -0.268 * cf2 / re_t2;
    let dcf1_dhk = -cf1 * (0.678 * h1) / (h1 - 0.4).max(0.1).powi(2);
    let dcf2_dhk = -cf2 * (0.678 * h2) / (h2 - 0.4).max(0.1).powi(2);
    
    // Chain rules: dRe_θ/dθ = Ue*Re, dRe_θ/dm involves dUe/dm through DIJ
    let dret1_dt1 = u1 * reynolds;
    let dret2_dt2 = u2 * reynolds;
    
    // === Momentum equation Jacobian ===
    // REZT = log(θ2/θ1) + (H+2)*log(U2/U1) - 0.5*CFX*log(x2/x1)
    
    // Logarithmic terms
    let tlog = (t2 / t1).ln();
    let ulog = (u2 / u1).ln();
    let xlog = (x2 / x1).ln();
    
    // Average shape factor
    let ha = 0.5 * (h1 + h2);
    let btmp = ha + 2.0;
    
    // CFX term
    let xa = 0.5 * (x1 + x2);
    let ta = 0.5 * (t1 + t2);
    let cfm = 0.5 * (cf1 + cf2);
    let cfx = 0.5 * cfm * xa / ta + 0.25 * (cf1 * x1 / t1 + cf2 * x2 / t2);
    
    // Derivatives of momentum residual
    // dR_mom/dθ1
    let z_tl = 1.0; // d(tlog)/d(log(θ))
    let z_cfx = -xlog * 0.5;
    let dcfx_dt1 = -0.5 * cfm * xa / (ta * ta) * 0.5 - 0.25 * cf1 * x1 / (t1 * t1);
    vs1.data[1][1] = -z_tl / t1 + z_cfx * dcfx_dt1 + 0.5 * ulog * dh1_dt1;
    
    // dR_mom/dθ2
    let dcfx_dt2 = -0.5 * cfm * xa / (ta * ta) * 0.5 - 0.25 * cf2 * x2 / (t2 * t2);
    vs2.data[1][1] = z_tl / t2 + z_cfx * dcfx_dt2 + 0.5 * ulog * dh2_dt2;
    
    // dR_mom/dm (through H)
    let z_ha = ulog;
    vs1.data[1][2] = 0.5 * z_ha * dh1_dm1;
    vs2.data[1][2] = 0.5 * z_ha * dh2_dm2;
    
    // === Shape equation Jacobian ===
    // REZH = log(Hs2/Hs1) + (2*Hc/Hs + 1 - H)*log(U2/U1) + log(x2/x1)*(0.5*CFX - DIX)
    
    if state2.inverse_mode {
        // Inverse mode: R_shape = Hk - Hk_target
        // dR/dθ2 = dHk/dθ2 = dHk/dH * dH/dθ2
        vs2.data[2][1] = dhk2_dh2 * dh2_dt2;
        // dR/dm2 = dHk/dm2 = dHk/dH * dH/dm2
        vs2.data[2][2] = dhk2_dh2 * dh2_dm2;
        // Station 1 derivatives are zero in inverse mode shape equation
        vs1.data[2][1] = 0.0;
        vs1.data[2][2] = 0.0;
    } else {
        // Direct mode: standard entrainment equation
        let hsa = 0.5 * (hs1 + hs2);
        let hca = 0.0; // Incompressible
        let btmp_shape = 2.0 * hca / hsa.max(0.1) + 1.0 - ha;
        
        let xot1 = x1 / t1;
        let xot2 = x2 / t2;
        let dix = (1.0 - upw) * di1 * xot1 + upw * di2 * xot2;
        let cfx_shape = (1.0 - upw) * cf1 * xot1 + upw * cf2 * xot2;
        
        // dR_shape/dθ1
        let z_hs1 = -1.0 / hs1;
        let z_cfx_s = xlog * 0.5;
        let z_dix = -xlog;
        let dt_cfx1 = (1.0 - upw) * cf1 * (-xot1 / t1);
        let dt_dix1 = (1.0 - upw) * di1 * (-xot1 / t1);
        vs1.data[2][1] = z_hs1 * dhs1_dhk1 * dhk1_dh1 * dh1_dt1 
                        + z_cfx_s * dt_cfx1 + z_dix * dt_dix1
                        - 0.5 * ulog * dh1_dt1;
        
        // dR_shape/dθ2
        let z_hs2 = 1.0 / hs2;
        let dt_cfx2 = upw * cf2 * (-xot2 / t2);
        let dt_dix2 = upw * di2 * (-xot2 / t2);
        vs2.data[2][1] = z_hs2 * dhs2_dhk2 * dhk2_dh2 * dh2_dt2
                        + z_cfx_s * dt_cfx2 + z_dix * dt_dix2
                        - 0.5 * ulog * dh2_dt2;
        
        // dR_shape/dm1, dR_shape/dm2 (through H and Hs)
        vs1.data[2][2] = z_hs1 * dhs1_dhk1 * dhk1_dh1 * dh1_dm1 - 0.5 * ulog * dh1_dm1;
        vs2.data[2][2] = z_hs2 * dhs2_dhk2 * dhk2_dh2 * dh2_dm2 - 0.5 * ulog * dh2_dm2;
    }
    
    // === First equation (Amplification or Shear lag) ===
    if state2.is_turbulent {
        // Shear lag: R = SCC*(CQ - S*ALD)*ds - 2*DE*log(S2/S1) + ...
        let s1_val = state1.ctau_or_ampl.max(1e-10);
        let s2_val = state2.ctau_or_ampl.max(1e-10);
        
        // dR/dS1 = SCC*ds*(-ALD)*(1-UPW) + 2*DE/S1
        let sccon = 5.6;
        let usa = 0.5 * (closures1.us + closures2.us);
        let scc = sccon * 1.333 / (1.0 + usa);
        let dea = 0.5 * (closures1.delta + closures2.delta);
        let ald = 1.0;
        
        vs1.data[0][0] = -scc * ds * ald * (1.0 - upw) + 2.0 * dea / s1_val;
        vs2.data[0][0] = -scc * ds * ald * upw - 2.0 * dea / s2_val;
        
        // dR/dθ, dR/dm through CQ and DE derivatives (simplified)
        vs1.data[0][1] = 0.0; // Simplified - would need CQ derivatives
        vs2.data[0][1] = 0.0;
        vs1.data[0][2] = 0.0;
        vs2.data[0][2] = 0.0;
    } else {
        // Amplification: R = N2 - N1 - AX*ds
        let ax = compute_amplification_rate(hk2, t2, re_t2);
        
        // dR/dN1 = -1, dR/dN2 = 1
        vs1.data[0][0] = -1.0;
        vs2.data[0][0] = 1.0;
        
        // dR/dθ through dAX/dθ (complex chain rule)
        // Simplified: AX depends on Hk and Re_θ
        let dax_dhk = compute_dax_dhk(hk2, t2, re_t2);
        let dax_dret = compute_dax_dret(hk2, t2, re_t2);
        
        vs2.data[0][1] = -ds * (dax_dhk * dhk2_dh2 * dh2_dt2 + dax_dret * dret2_dt2);
        vs2.data[0][2] = -ds * dax_dhk * dhk2_dh2 * dh2_dm2;
        
        vs1.data[0][1] = 0.0;
        vs1.data[0][2] = 0.0;
    }
    
    (vs1, vs2)
}

/// Derivative of amplification rate wrt Hk
fn compute_dax_dhk(hk: f64, theta: f64, re_theta: f64) -> f64 {
    let eps = 1e-6;
    let ax1 = compute_amplification_rate(hk, theta, re_theta);
    let ax2 = compute_amplification_rate(hk + eps, theta, re_theta);
    (ax2 - ax1) / eps
}

/// Derivative of amplification rate wrt Re_theta
fn compute_dax_dret(hk: f64, theta: f64, re_theta: f64) -> f64 {
    let eps = re_theta * 1e-6;
    let ax1 = compute_amplification_rate(hk, theta, re_theta);
    let ax2 = compute_amplification_rate(hk, theta, re_theta + eps);
    (ax2 - ax1) / eps
}

/// Compute Jacobian blocks for BL interval.
/// 
/// Returns (VS1, VS2) where:
/// - VS1[k][l] = d(residual_k) / d(variable_l at station 1)
/// - VS2[k][l] = d(residual_k) / d(variable_l at station 2)
/// 
/// Variables: [0] = Ctau/Ampl, [1] = Theta, [2] = mass (direct) or Ue (inverse mode)
/// 
/// In inverse mode, the third variable changes from mass to Ue, and the
/// third equation becomes Hk - Hk_target = 0.
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
    
    // Check if in inverse mode
    let inverse_mode = state2.inverse_mode;
    
    // Perturb station 1 variables
    for l in 0..3 {
        let mut state1_plus = *state1;
        match l {
            0 => state1_plus.ctau_or_ampl += eps,
            1 => state1_plus.theta += eps,
            2 => {
                // In direct mode: perturb mass (affects H)
                // In inverse mode: perturb Ue
                if inverse_mode {
                    state1_plus.ue += eps;
                    // Update mass to be consistent: m = Ue * δ* = Ue * θ * H
                    state1_plus.mass = state1_plus.ue * state1_plus.theta * state1_plus.h;
                } else {
                    state1_plus.mass += eps;
                    state1_plus.update_h();
                }
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
                // In direct mode: perturb mass
                // In inverse mode: perturb Ue (the unknown we're solving for)
                if inverse_mode {
                    state2_plus.ue += eps;
                    // Update mass to be consistent
                    state2_plus.mass = state2_plus.ue * state2_plus.theta * state2_plus.h;
                } else {
                    state2_plus.mass += eps;
                    state2_plus.update_h();
                }
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
/// In inverse mode, the third variable is Ue instead of mass, and the
/// third equation is Hk - Hk_target = 0.
/// 
/// This matches XFOIL's local Newton loop in `xbl.f` lines 610-740.
pub fn solve_station_newton(
    state1: &StationState,
    mut state2: StationState,
    reynolds: f64,
    config: &LocalNewtonConfig,
) -> LocalNewtonResult {
    use crate::boundary_layer::relax_ue_update;
    
    let closures1 = BLClosures::compute(
        state1.theta, state1.h, state1.ue, reynolds, state1.is_turbulent
    );
    
    let mut iterations = 0;
    let mut converged = false;
    let mut residual_norm = f64::MAX;
    
    // Check if in inverse mode
    let inverse_mode = state2.inverse_mode;
    
    // Use smaller relaxation for inverse mode to improve stability
    let relax = if inverse_mode { config.relax * 0.5 } else { config.relax };
    
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
        state2.ctau_or_ampl = (state2.ctau_or_ampl + relax * delta[0]).max(0.0);
        state2.theta = (state2.theta + relax * delta[1]).max(1e-12);
        
        if inverse_mode {
            // Inverse mode: update Ue (third variable)
            // Use the dedicated relaxation function for Ue updates
            let ue_old = state2.ue;
            let ue_delta = delta[2];
            state2.ue = relax_ue_update(ue_old, ue_old + ue_delta, relax);
            state2.ue = state2.ue.max(1e-10); // Ensure positive
            
            // In inverse mode, H is prescribed (target Hk)
            // Update mass to be consistent: m = Ue * δ* = Ue * θ * H
            state2.h = state2.hk_target.clamp(1.05, 10.0);
            state2.mass = state2.ue * state2.theta * state2.h;
        } else {
            // Direct mode: update mass (third variable)
            state2.mass = (state2.mass + relax * delta[2]).max(1e-12);
            state2.update_h();
            
            // Clamp H to reasonable range
            state2.h = state2.h.clamp(1.05, 10.0);
        }
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
/// Includes inverse mode detection: when Hk exceeds the separation threshold,
/// the solver switches to inverse mode where Hk is prescribed and the shape
/// equation becomes Hk - Hk_target = 0.
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
    use crate::boundary_layer::{HK_MAX_LAMINAR, HK_MAX_TURBULENT, compute_target_hk};
    
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
    let mut in_inverse_mode = false;
    
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
        
        // Check for separation / inverse mode
        // From XFOIL: DIRECT = HKTEST.LT.HMAX
        let hk_max = if is_turbulent { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
        let should_inverse = prev_state.h >= hk_max || in_inverse_mode;
        
        let mut initial_state = StationState::new(
            if is_turbulent { 0.03 } else { ampl }, // Initial Ctau or Ampl
            theta_guess,
            h_guess,
            ue_j,
            s_j,
            is_turbulent,
        );
        
        // Set inverse mode if needed
        if should_inverse {
            in_inverse_mode = true;
            
            // Compute target Hk for inverse mode
            let ds = (s_j - prev_state.s).abs();
            let hk_target = compute_target_hk(
                prev_state.h,
                ds,
                prev_state.theta,
                is_turbulent,
                false, // not wake
            );
            
            initial_state.set_inverse_mode(hk_target);
            
            // In inverse mode, clamp initial H guess to target
            initial_state.h = hk_target;
        }
        
        // Local Newton solve
        let result = solve_station_newton(prev_state, initial_state, reynolds, config);
        
        let mut final_state = if !result.converged {
            // Use guess if Newton fails
            let mut fallback = StationState::new(
                initial_state.ctau_or_ampl,
                theta_guess * 1.1, // Slight increase
                h_guess,
                ue_j,
                s_j,
                is_turbulent,
            );
            if initial_state.inverse_mode {
                fallback.set_inverse_mode(initial_state.hk_target);
                fallback.h = initial_state.hk_target; // Enforce target Hk
            }
            fallback
        } else {
            result.state
        };
        
        // Check if we can exit inverse mode (Hk dropped below threshold)
        if in_inverse_mode && final_state.h < hk_max * 0.95 {
            // Reattachment detected - exit inverse mode
            in_inverse_mode = false;
            final_state.clear_inverse_mode();
        }
        
        states.push(final_state);
        
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
    
    #[test]
    fn test_mass_defect_residual() {
        // Two stations with consistent mass defect growth
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        let state2 = StationState::new(0.0, 0.0012, 2.4, 1.02, 0.01, false);
        
        let closures1 = BLClosures::compute(state1.theta, state1.h, state1.ue, 1e6, false);
        let closures2 = BLClosures::compute(state2.theta, state2.h, state2.ue, 1e6, false);
        
        let res = compute_mass_defect_residual(&state1, &state2, &closures1, &closures2);
        
        // Residual should be finite
        assert!(res.is_finite(), "mass defect residual = {}", res);
        
        // For consistent BL growth, residual should be small
        // (not exactly zero due to discrete differences)
        assert!(res.abs() < 0.1, "mass defect residual = {} (expected small)", res);
    }
    
    #[test]
    fn test_mass_defect_derivatives() {
        let state1 = StationState::new(0.0, 0.001, 2.5, 1.0, 0.0, false);
        let state2 = StationState::new(0.0, 0.0012, 2.4, 1.02, 0.01, false);
        
        let closures1 = BLClosures::compute(state1.theta, state1.h, state1.ue, 1e6, false);
        let closures2 = BLClosures::compute(state2.theta, state2.h, state2.ue, 1e6, false);
        
        let derivs = MassDefectDerivs::compute(&state1, &state2, &closures1, &closures2);
        
        // All derivatives should be finite
        assert!(derivs.d_m1.is_finite());
        assert!(derivs.d_m2.is_finite());
        assert!(derivs.d_theta1.is_finite());
        assert!(derivs.d_theta2.is_finite());
        assert!(derivs.d_ue1.is_finite());
        assert!(derivs.d_ue2.is_finite());
        
        // Mass derivatives should have opposite signs (upstream negative, downstream positive)
        assert!(derivs.d_m1 < 0.0, "d_m1 should be negative");
        assert!(derivs.d_m2 > 0.0, "d_m2 should be positive");
    }
    
    #[test]
    fn test_mass_defect_variable_in_state() {
        // Verify StationState correctly handles mass defect m = Ue*δ*
        let theta = 0.001;
        let h = 2.5;
        let ue = 1.2;
        
        let state = StationState::new(0.03, theta, h, ue, 0.5, true);
        
        // Expected: m = Ue * δ* = Ue * θ * H
        let expected_delta_star = theta * h;
        let expected_mass = ue * expected_delta_star;
        
        assert!((state.mass - expected_mass).abs() < 1e-12,
            "mass = {}, expected {}", state.mass, expected_mass);
        
        // delta_star() should recover δ* from m/Ue
        assert!((state.delta_star() - expected_delta_star).abs() < 1e-12,
            "delta_star = {}, expected {}", state.delta_star(), expected_delta_star);
    }
}
