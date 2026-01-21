//! Boundary layer integral equations
//!
//! This module implements the core BL equation computations from XFOIL:
//! - BLVAR: Compute all secondary variables from primary variables
//! - BLDIF: Compute residuals and Jacobian blocks for Newton system
//!
//! # XFOIL Reference
//! - BLVAR: xblsys.f line 784
//! - BLDIF: xblsys.f line 1552

use crate::closures::{
    amplification_rate, cf_laminar, cf_turbulent, density_shape, dissipation_laminar,
    dissipation_wake, hkin, hs_laminar, hs_turbulent,
};
use crate::constants::{CTCON, DLCON, DUXCON, GACON, GBCON, GCCON, SCCON};
use crate::state::BlStation;

/// Flow type for boundary layer calculations
///
/// Maps to XFOIL's ITYP parameter in BLVAR/BLDIF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowType {
    /// Laminar boundary layer (ITYP = 1)
    Laminar,
    /// Turbulent boundary layer (ITYP = 2)
    Turbulent,
    /// Turbulent wake (ITYP = 3)
    Wake,
}

/// Residuals of the integral BL equations
///
/// These residuals are driven to zero by the Newton solver.
/// The three equations are:
/// 1. Amplification (laminar) or shear-lag (turbulent/wake)
/// 2. Momentum integral
/// 3. Shape parameter (kinetic energy)
#[derive(Debug, Clone, Default)]
pub struct BlResiduals {
    /// First equation residual: amplification (laminar) or shear-lag (turbulent)
    pub res_third: f64,
    /// Momentum integral equation residual
    pub res_mom: f64,
    /// Shape parameter equation residual
    pub res_shape: f64,
}

/// Jacobian blocks for Newton system
///
/// The BL equations relate two adjacent stations (upstream 1, downstream 2).
/// The Jacobian has blocks for each station plus sensitivities to global parameters.
///
/// Variable ordering: [s/ampl, θ, δ*, u, x] (5 variables per station)
/// Equation ordering: [third, momentum, shape] (3 equations per interval)
#[derive(Debug, Clone)]
pub struct BlJacobian {
    /// Jacobian w.r.t. upstream station variables [3 equations × 5 variables]
    pub vs1: [[f64; 5]; 3],
    /// Jacobian w.r.t. downstream station variables [3 equations × 5 variables]
    pub vs2: [[f64; 5]; 3],
    /// Sensitivity to Mach number (via M² derivatives)
    pub vsm: [f64; 3],
    /// Sensitivity to Reynolds number
    pub vsr: [f64; 3],
}

impl Default for BlJacobian {
    fn default() -> Self {
        Self {
            vs1: [[0.0; 5]; 3],
            vs2: [[0.0; 5]; 3],
            vsm: [0.0; 3],
            vsr: [0.0; 3],
        }
    }
}

/// Extended derivatives for BLVAR computations
///
/// These track derivatives of secondary variables w.r.t. primary variables
/// and global parameters (Mach, Reynolds). Used internally by blvar().
/// Note: Kept for future use when full Jacobian propagation is needed.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
struct BlvarDerivatives {
    // Hk derivatives (from hkin)
    hk_u: f64,
    hk_t: f64,
    hk_d: f64,
    hk_ms: f64,

    // Hc derivatives
    hc_u: f64,
    hc_t: f64,
    hc_d: f64,
    hc_ms: f64,

    // Hs derivatives
    hs_u: f64,
    hs_t: f64,
    hs_d: f64,
    hs_ms: f64,
    hs_re: f64,

    // Us derivatives
    us_u: f64,
    us_t: f64,
    us_d: f64,
    us_ms: f64,
    us_re: f64,

    // CQ derivatives
    cq_u: f64,
    cq_t: f64,
    cq_d: f64,
    cq_ms: f64,
    cq_re: f64,

    // Cf derivatives
    cf_u: f64,
    cf_t: f64,
    cf_d: f64,
    cf_ms: f64,
    cf_re: f64,

    // DI derivatives
    di_u: f64,
    di_t: f64,
    di_d: f64,
    di_s: f64,
    di_ms: f64,
    di_re: f64,

    // DE derivatives
    de_u: f64,
    de_t: f64,
    de_d: f64,
    de_ms: f64,
}

/// Compute all secondary BL variables from primary variables
///
/// This is the Rust port of XFOIL's BLVAR subroutine (xblsys.f:784-1121).
/// Given primary variables (x, u, θ, δ*, s, ampl), it computes all
/// secondary variables (H, Hk, Hs, Hc, Cf, CD, etc.) and their derivatives.
///
/// # Arguments
/// * `station` - BlStation containing primary variables (modified in place)
/// * `flow_type` - Laminar, Turbulent, or Wake
/// * `msq` - Mach number squared (M²)
/// * `re` - Reference Reynolds number
///
/// # Reference
/// XFOIL xblsys.f BLVAR (line 784)
pub fn blvar(station: &mut BlStation, flow_type: FlowType, msq: f64, re: f64) {
    // Compute shape factor H = δ*/θ
    station.h = station.delta_star / station.theta;
    let h_t = -station.h / station.theta; // ∂H/∂θ
    let h_d = 1.0 / station.theta; // ∂H/∂δ*

    // Kinematic shape factor Hk with compressibility correction
    let hkin_result = hkin(station.h, msq);
    let mut hk = hkin_result.hk;
    let hk_h = hkin_result.hk_h;
    let hk_msq = hkin_result.hk_msq;

    // Clamp Hk to minimum values (xblsys.f:798-799)
    let hk_min = if flow_type == FlowType::Wake {
        1.00005
    } else {
        1.05
    };
    if hk < hk_min {
        hk = hk_min;
    }
    station.hk = hk;

    // Chain rule: Hk derivatives w.r.t. primary variables
    // Hk depends on H and M², H depends on θ and δ*
    let hk_t = hk_h * h_t;
    let hk_d = hk_h * h_d;
    let hk_u = 0.0; // No direct U dependence through H
    let hk_ms = hk_msq; // M² sensitivity

    // Store Hk derivatives
    station.derivs.hk_h = hk_h;
    station.derivs.hk_msq = hk_msq;

    // Compute Reynolds number Rθ = Ue * θ * Re
    // In XFOIL: RT2 = RE * U2 * T2
    station.r_theta = re * station.u * station.theta;
    let rt_t = re * station.u; // ∂Rθ/∂θ
    let rt_u = re * station.theta; // ∂Rθ/∂u
    let rt_re = station.u * station.theta; // ∂Rθ/∂Re
    let rt_ms = 0.0; // No Mach dependence

    // Compute M² derivatives w.r.t. U (M² = U²/a² where a is speed of sound)
    // In incompressible limit, M² doesn't depend on U through this route
    // but we track it for compressible extensions
    let m_u = 0.0; // Simplified for now
    let m_ms = 1.0; // Direct sensitivity

    // === Density shape factor Hc (xblsys.f:802-806) ===
    let hct_result = density_shape(hk, msq);
    station.hc = hct_result.hc;
    let hc_hk = hct_result.hc_hk;
    let hc_msq = hct_result.hc_msq;

    // Chain rule for Hc
    let hc_u = hc_hk * hk_u + hc_msq * m_u;
    let hc_t = hc_hk * hk_t;
    let hc_d = hc_hk * hk_d;
    let hc_ms = hc_hk * hk_ms + hc_msq * m_ms;

    // === Energy shape factor Hs (xblsys.f:809-819) ===
    let hs_result = match flow_type {
        FlowType::Laminar => hs_laminar(hk, station.r_theta, msq),
        FlowType::Turbulent | FlowType::Wake => hs_turbulent(hk, station.r_theta, msq),
    };
    station.hs = hs_result.hs;
    let hs_hk = hs_result.hs_hk;
    let hs_rt = hs_result.hs_rt;
    let hs_msq = hs_result.hs_msq;

    // Chain rule for Hs
    let hs_u = hs_hk * hk_u + hs_rt * rt_u + hs_msq * m_u;
    let hs_t = hs_hk * hk_t + hs_rt * rt_t;
    let hs_d = hs_hk * hk_d;
    let hs_ms = hs_hk * hk_ms + hs_rt * rt_ms + hs_msq * m_ms;
    let hs_re = hs_rt * rt_re;

    // Store Hs derivatives
    station.derivs.hs_hk = hs_hk;
    station.derivs.hs_rt = hs_rt;
    station.derivs.hs_msq = hs_msq;

    // === Normalized slip velocity Us (xblsys.f:821-851) ===
    // Us = 0.5 * Hs * (1 - (Hk-1)/(GBCON*H))
    let mut us = 0.5 * station.hs * (1.0 - (hk - 1.0) / (GBCON * station.h));
    let us_hs = 0.5 * (1.0 - (hk - 1.0) / (GBCON * station.h));
    let us_hk = 0.5 * station.hs * (-1.0 / (GBCON * station.h));
    let us_h = 0.5 * station.hs * (hk - 1.0) / (GBCON * station.h * station.h);

    // Chain rule for Us
    let mut us_u = us_hs * hs_u + us_hk * hk_u;
    let mut us_t = us_hs * hs_t + us_hk * hk_t + us_h * h_t;
    let mut us_d = us_hs * hs_d + us_hk * hk_d + us_h * h_d;
    let mut us_ms = us_hs * hs_ms + us_hk * hk_ms;
    let mut us_re = us_hs * hs_re;

    // Clamp Us to prevent instability (xblsys.f:833-851)
    let us_max = if flow_type == FlowType::Wake {
        0.99995
    } else {
        0.95
    };
    if us > us_max {
        us = if flow_type == FlowType::Wake {
            0.99995
        } else {
            0.98
        };
        us_u = 0.0;
        us_t = 0.0;
        us_d = 0.0;
        us_ms = 0.0;
        us_re = 0.0;
    }

    // === Equilibrium shear coefficient CQ (xblsys.f:853-895) ===
    // CQ = sqrt(CTCON * Hs * HKB * HKC² / (USB * H * Hk²))
    let gcc = if flow_type == FlowType::Turbulent {
        GCCON
    } else {
        0.0
    };

    let mut hkc = hk - 1.0 - gcc / station.r_theta.max(1.0);
    let mut hkc_hk = 1.0;
    let mut hkc_rt = gcc / (station.r_theta.max(1.0) * station.r_theta.max(1.0));

    if flow_type == FlowType::Turbulent && hkc < 0.01 {
        hkc = 0.01;
        hkc_hk = 0.0;
        hkc_rt = 0.0;
    }
    if flow_type != FlowType::Turbulent {
        hkc = hk - 1.0;
        hkc_hk = 1.0;
        hkc_rt = 0.0;
    }

    let hkb = hk - 1.0;
    let usb = 1.0 - us;

    // Compute CQ
    let cq_arg = CTCON * station.hs * hkb * hkc * hkc / (usb * station.h * hk * hk);
    let cq = if cq_arg > 0.0 { cq_arg.sqrt() } else { 0.0 };

    // CQ derivatives (chain rule through the square root)
    let (cq_hs, cq_us, cq_hk, cq_rt, cq_h) = if cq > 1e-12 {
        let cq_arg_hs = CTCON * hkb * hkc * hkc / (usb * station.h * hk * hk);
        let cq_arg_us = CTCON * station.hs * hkb * hkc * hkc / (usb * usb * station.h * hk * hk);
        let cq_arg_hk = CTCON * station.hs * hkc * hkc / (usb * station.h * hk * hk)
            - 2.0 * CTCON * station.hs * hkb * hkc * hkc / (usb * station.h * hk * hk * hk)
            + 2.0 * CTCON * station.hs * hkb * hkc * hkc_hk / (usb * station.h * hk * hk);
        let cq_arg_rt = 2.0 * CTCON * station.hs * hkb * hkc * hkc_rt / (usb * station.h * hk * hk);
        let cq_arg_h =
            -CTCON * station.hs * hkb * hkc * hkc / (usb * station.h * station.h * hk * hk);

        (
            0.5 * cq_arg_hs / cq,
            0.5 * cq_arg_us / cq,
            0.5 * cq_arg_hk / cq,
            0.5 * cq_arg_rt / cq,
            0.5 * cq_arg_h / cq,
        )
    } else {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    };

    // Chain rule for CQ w.r.t. primary variables
    let cq_u = cq_hs * hs_u + cq_us * us_u + cq_hk * hk_u + cq_rt * rt_u;
    let cq_t = cq_hs * hs_t + cq_us * us_t + cq_hk * hk_t + cq_rt * rt_t + cq_h * h_t;
    let cq_d = cq_hs * hs_d + cq_us * us_d + cq_hk * hk_d + cq_h * h_d;
    let cq_ms_val = cq_hs * hs_ms + cq_us * us_ms + cq_hk * hk_ms;
    let cq_re_val = cq_hs * hs_re + cq_us * us_re + cq_rt * rt_re;

    // === Skin friction Cf (xblsys.f:899-927) ===
    let (cf, cf_hk, cf_rt, cf_msq) = match flow_type {
        FlowType::Wake => (0.0, 0.0, 0.0, 0.0), // Wake has no wall friction
        FlowType::Laminar => {
            let cf_result = cf_laminar(hk, station.r_theta, msq);
            (cf_result.cf, cf_result.cf_hk, cf_result.cf_rt, cf_result.cf_msq)
        }
        FlowType::Turbulent => {
            let cf_turb = cf_turbulent(hk, station.r_theta, msq);
            let cf_lam = cf_laminar(hk, station.r_theta, msq);
            // Use laminar if it gives higher Cf (low Rθ case)
            if cf_lam.cf > cf_turb.cf {
                (cf_lam.cf, cf_lam.cf_hk, cf_lam.cf_rt, cf_lam.cf_msq)
            } else {
                (cf_turb.cf, cf_turb.cf_hk, cf_turb.cf_rt, cf_turb.cf_msq)
            }
        }
    };
    station.cf = cf;

    // Chain rule for Cf
    let cf_u = cf_hk * hk_u + cf_rt * rt_u + cf_msq * m_u;
    let cf_t = cf_hk * hk_t + cf_rt * rt_t;
    let cf_d = cf_hk * hk_d;
    let cf_ms = cf_hk * hk_ms + cf_rt * rt_ms + cf_msq * m_ms;
    let cf_re = cf_rt * rt_re;

    // Store Cf derivatives
    station.derivs.cf_hk = cf_hk;
    station.derivs.cf_rt = cf_rt;
    station.derivs.cf_msq = cf_msq;

    // === Dissipation coefficient DI (xblsys.f:930-1098) ===
    let (mut di, mut di_s, mut di_u, mut di_t, mut di_d, mut di_ms, mut di_re) = match flow_type {
        FlowType::Laminar => {
            let di_result = dissipation_laminar(hk, station.r_theta);
            let di_hk = di_result.di_hk;
            let di_rt = di_result.di_rt;
            // Chain rule
            let di_u = di_hk * hk_u + di_rt * rt_u;
            let di_t = di_hk * hk_t + di_rt * rt_t;
            let di_d = di_hk * hk_d;
            let di_ms = di_hk * hk_ms;
            let di_re = di_rt * rt_re;
            (di_result.di, 0.0, di_u, di_t, di_d, di_ms, di_re)
        }
        FlowType::Turbulent => {
            // Turbulent wall contribution (xblsys.f:947-965)
            let cf2t_result = cf_turbulent(hk, station.r_theta, msq);
            let cf2t = cf2t_result.cf;
            let cf2t_hk = cf2t_result.cf_hk;
            let cf2t_rt = cf2t_result.cf_rt;
            let cf2t_msq = cf2t_result.cf_msq;

            // DI = (0.5*CF2T*US) * 2.0/HS
            let di_wall = (0.5 * cf2t * us) * 2.0 / station.hs;
            let di_hs_wall = -(0.5 * cf2t * us) * 2.0 / (station.hs * station.hs);
            let di_us_wall = (0.5 * cf2t) * 2.0 / station.hs;
            let di_cf2t = (0.5 * us) * 2.0 / station.hs;

            // Chain rule for wall contribution
            let cf2t_u = cf2t_hk * hk_u + cf2t_rt * rt_u + cf2t_msq * m_u;
            let cf2t_t = cf2t_hk * hk_t + cf2t_rt * rt_t;
            let cf2t_d = cf2t_hk * hk_d;
            let cf2t_ms = cf2t_hk * hk_ms + cf2t_rt * rt_ms + cf2t_msq * m_ms;
            let cf2t_re = cf2t_rt * rt_re;

            let di_u = di_hs_wall * hs_u + di_us_wall * us_u + di_cf2t * cf2t_u;
            let di_t = di_hs_wall * hs_t + di_us_wall * us_t + di_cf2t * cf2t_t;
            let di_d = di_hs_wall * hs_d + di_us_wall * us_d + di_cf2t * cf2t_d;
            let di_ms = di_hs_wall * hs_ms + di_us_wall * us_ms + di_cf2t * cf2t_ms;
            let di_re = di_hs_wall * hs_re + di_us_wall * us_re + di_cf2t * cf2t_re;

            // DFAC correction for very low Hk (xblsys.f:969-991)
            let grt = station.r_theta.ln();
            let hmin = 1.0 + 2.1 / grt;
            let fl = (hk - 1.0) / (hmin - 1.0);
            let tfl = fl.tanh();
            let dfac = 0.5 + 0.5 * tfl;

            let di_wall_corrected = di_wall * dfac;
            let di_u_corrected = di_u * dfac;
            let di_t_corrected = di_t * dfac;
            let di_d_corrected = di_d * dfac;
            let di_ms_corrected = di_ms * dfac;
            let di_re_corrected = di_re * dfac;

            // Outer layer contribution (xblsys.f:1008-1036)
            // DD = S2² * (0.995 - US) * 2.0/HS
            let s = station.ctau;
            let dd = s * s * (0.995 - us) * 2.0 / station.hs;
            let dd_hs = -s * s * (0.995 - us) * 2.0 / (station.hs * station.hs);
            let dd_us = -s * s * 2.0 / station.hs;
            let dd_s = s * 2.0 * (0.995 - us) * 2.0 / station.hs;

            // Laminar stress contribution (xblsys.f:1024-1035)
            let dd2 = 0.15 * (0.995 - us) * (0.995 - us) / station.r_theta * 2.0 / station.hs;
            let dd2_us = -0.15 * (0.995 - us) * 2.0 / station.r_theta * 2.0 / station.hs;
            let dd2_hs = -dd2 / station.hs;
            let dd2_rt = -dd2 / station.r_theta;

            // Total DI
            let di = di_wall_corrected + dd + dd2;
            let di_s = dd_s;
            let di_u = di_u_corrected + dd_hs * hs_u + dd_us * us_u + dd2_hs * hs_u + dd2_us * us_u + dd2_rt * rt_u;
            let di_t = di_t_corrected + dd_hs * hs_t + dd_us * us_t + dd2_hs * hs_t + dd2_us * us_t + dd2_rt * rt_t;
            let di_d = di_d_corrected + dd_hs * hs_d + dd_us * us_d + dd2_hs * hs_d + dd2_us * us_d;
            let di_ms = di_ms_corrected + dd_hs * hs_ms + dd_us * us_ms + dd2_hs * hs_ms + dd2_us * us_ms;
            let di_re = di_re_corrected + dd_hs * hs_re + dd_us * us_re + dd2_hs * hs_re + dd2_us * us_re + dd2_rt * rt_re;

            // Check if laminar DI is higher (low Rθ case)
            let di_lam = dissipation_laminar(hk, station.r_theta);
            if di_lam.di > di {
                let di_hk = di_lam.di_hk;
                let di_rt = di_lam.di_rt;
                let di_u = di_hk * hk_u + di_rt * rt_u;
                let di_t = di_hk * hk_t + di_rt * rt_t;
                let di_d = di_hk * hk_d;
                let di_ms = di_hk * hk_ms;
                let di_re = di_rt * rt_re;
                (di_lam.di, 0.0, di_u, di_t, di_d, di_ms, di_re)
            } else {
                (di, di_s, di_u, di_t, di_d, di_ms, di_re)
            }
        }
        FlowType::Wake => {
            // Wake: use DILW (xblsys.f:1070-1086)
            let di_result = dissipation_wake(hk, station.r_theta);
            let di_hk = di_result.di_hk;
            let di_rt = di_result.di_rt;
            // Chain rule
            let di_u = di_hk * hk_u + di_rt * rt_u;
            let di_t = di_hk * hk_t + di_rt * rt_t;
            let di_d = di_hk * hk_d;
            let di_ms = di_hk * hk_ms;
            let di_re = di_rt * rt_re;
            // Double dissipation for wake (two halves)
            (
                di_result.di * 2.0,
                0.0,
                di_u * 2.0,
                di_t * 2.0,
                di_d * 2.0,
                di_ms * 2.0,
                di_re * 2.0,
            )
        }
    };

    station.cd = di;
    station.derivs.cd_hk = di_u; // Store for later use
    station.derivs.cd_rt = di_t;

    // === BL thickness DE from Green's correlation (xblsys.f:1100-1118) ===
    // DE = (3.15 + 1.72/(HK-1)) * T + D
    let de = (3.15 + 1.72 / (hk - 1.0)) * station.theta + station.delta_star;
    let de_hk = -1.72 / ((hk - 1.0) * (hk - 1.0)) * station.theta;

    let de_u = de_hk * hk_u;
    let de_t = de_hk * hk_t + 3.15 + 1.72 / (hk - 1.0);
    let de_d = de_hk * hk_d + 1.0;
    let de_ms = de_hk * hk_ms;

    // Clamp DE to maximum HDMAX * T
    let hdmax = 12.0;
    let de_final = if de > hdmax * station.theta {
        hdmax * station.theta
    } else {
        de
    };

    // Store mass defect
    station.mass_defect = station.u * station.delta_star;

    // Store H derivative (for compatibility)
    station.derivs.h_theta = h_t;
    station.derivs.h_delta_star = h_d;
}

/// Compute midpoint skin friction coefficient
///
/// XFOIL Reference: xblsys.f BLMID (lines 1124-1182)
///
/// Computes CFM at the midpoint between two stations for use in the
/// momentum equation. This improves accuracy of the momentum integral.
fn blmid(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    // Average values at midpoint
    let hka = 0.5 * (s1.hk + s2.hk);
    let rta = 0.5 * (s1.r_theta + s2.r_theta);
    let ma = 0.5 * msq; // Simplified - would need proper M at each station

    // Midpoint Cf
    let (cfm, cfm_hka, cfm_rta, cfm_ma) = match flow_type {
        FlowType::Wake => (0.0, 0.0, 0.0, 0.0),
        FlowType::Laminar => {
            let cf_result = cf_laminar(hka, rta, ma);
            (cf_result.cf, cf_result.cf_hk, cf_result.cf_rt, cf_result.cf_msq)
        }
        FlowType::Turbulent => {
            let cf_turb = cf_turbulent(hka, rta, ma);
            let cf_lam = cf_laminar(hka, rta, ma);
            if cf_lam.cf > cf_turb.cf {
                (cf_lam.cf, cf_lam.cf_hk, cf_lam.cf_rt, cf_lam.cf_msq)
            } else {
                (cf_turb.cf, cf_turb.cf_hk, cf_turb.cf_rt, cf_turb.cf_msq)
            }
        }
    };

    // Derivatives w.r.t station 1 and 2 (via midpoint)
    let cfm_hk1 = 0.5 * cfm_hka;
    let cfm_hk2 = 0.5 * cfm_hka;
    let cfm_rt1 = 0.5 * cfm_rta;
    let cfm_rt2 = 0.5 * cfm_rta;

    (cfm, cfm_hk1, cfm_hk2, cfm_rt1, cfm_rt2, cfm_ma)
}

/// Compute BL equation residuals and Jacobian between two stations
///
/// This is the Rust port of XFOIL's BLDIF subroutine (xblsys.f:1552-1977).
/// It sets up the Newton system for the boundary layer equations.
///
/// # Arguments
/// * `s1` - Upstream station (station 1)
/// * `s2` - Downstream station (station 2)
/// * `flow_type` - Laminar, Turbulent, or Wake
/// * `msq` - Mach number squared
/// * `re` - Reference Reynolds number
///
/// # Returns
/// * `(BlResiduals, BlJacobian)` - Residuals and Jacobian blocks
///
/// # Reference
/// XFOIL xblsys.f BLDIF (line 1552)
pub fn bldif(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    let mut res = BlResiduals::default();
    let mut jac = BlJacobian::default();

    // === Logarithmic differences (xblsys.f:1574-1585) ===
    let xlog = (s2.x / s1.x.max(1e-20)).ln();
    let ulog = (s2.u / s1.u.max(1e-20)).ln();
    let tlog = (s2.theta / s1.theta.max(1e-20)).ln();
    let hlog = (s2.hs / s1.hs.max(1e-20)).ln();

    // === Local upwinding parameter UPW (xblsys.f:1598-1644) ===
    // Based on log(Hk-1) changes
    let hupwt = 1.0;
    let hdcon = if flow_type == FlowType::Wake {
        hupwt / (s2.hk * s2.hk)
    } else {
        5.0 * hupwt / (s2.hk * s2.hk)
    };

    let arg = ((s2.hk - 1.0) / (s1.hk - 1.0).max(1e-6)).abs();
    let hl = arg.ln();
    let hlsq = (hl * hl).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;

    // UPW derivatives (simplified)
    let upw_hl = ehh * hl * hdcon;
    let hl_hk1 = -1.0 / (s1.hk - 1.0).max(1e-6);
    let hl_hk2 = 1.0 / (s2.hk - 1.0).max(1e-6);

    // === First equation: Amplification (laminar) or Shear-lag (turbulent/wake) ===
    match flow_type {
        FlowType::Laminar => {
            // Amplification equation (xblsys.f:1654-1682)
            // REZC = AMPL2 - AMPL1 - AX*(X2-X1)
            let ax_result = amplification_rate(
                0.5 * (s1.hk + s2.hk),
                0.5 * (s1.theta + s2.theta),
                0.5 * (s1.r_theta + s2.r_theta),
            );
            let ax = ax_result.ax;
            let dx = s2.x - s1.x;

            res.res_third = -(s2.ampl - s1.ampl - ax * dx);

            // Jacobian entries for amplification equation
            // Variable order: [ampl, θ, δ*, u, x]
            jac.vs1[0][0] = -1.0; // -∂/∂AMPL1
            jac.vs1[0][4] = ax; // ∂/∂X1
            jac.vs2[0][0] = 1.0; // ∂/∂AMPL2
            jac.vs2[0][4] = -ax; // -∂/∂X2
        }
        FlowType::Turbulent | FlowType::Wake => {
            // Shear-lag equation (xblsys.f:1684-1841)
            let sa = (1.0 - upw) * s1.ctau + upw * s2.ctau;
            let cqa = (1.0 - upw) * s1.hk + upw * s2.hk; // Simplified CQ average
            let cfa = (1.0 - upw) * s1.cf + upw * s2.cf;
            let hka = (1.0 - upw) * s1.hk + upw * s2.hk;

            let usa = 0.5 * (0.5 * s1.hs + 0.5 * s2.hs); // Simplified Us average
            let rta = 0.5 * (s1.r_theta + s2.r_theta);
            let dea = 0.5
                * ((3.15 + 1.72 / (s1.hk - 1.0).max(0.01)) * s1.theta
                    + (3.15 + 1.72 / (s2.hk - 1.0).max(0.01)) * s2.theta);
            let da = 0.5 * (s1.delta_star + s2.delta_star);

            let ald = if flow_type == FlowType::Wake {
                DLCON
            } else {
                1.0
            };

            // Equilibrium dUe/dx
            let gcc = if flow_type == FlowType::Turbulent {
                GCCON
            } else {
                0.0
            };
            let hkc = (hka - 1.0 - gcc / rta.max(1.0)).max(0.01);
            let hr = hkc / (GACON * ald * hka);
            let uq = (0.5 * cfa - hr * hr) / (GBCON * da.max(1e-10));

            let scc = SCCON * 1.333 / (1.0 + usa);
            let slog = (s2.ctau / s1.ctau.max(1e-20)).ln();
            let dxi = s2.x - s1.x;

            // Residual (simplified form)
            res.res_third = -(scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog
                + dea * 2.0 * (uq * dxi - ulog) * DUXCON);

            // Jacobian entries for shear-lag equation
            jac.vs1[0][0] = (1.0 - upw) * (-scc * ald * dxi) + dea * 2.0 / s1.ctau.max(1e-20);
            jac.vs2[0][0] = upw * (-scc * ald * dxi) - dea * 2.0 / s2.ctau.max(1e-20);
            jac.vs1[0][4] = -scc * (cqa - sa * ald) - dea * 2.0 * uq * DUXCON;
            jac.vs2[0][4] = scc * (cqa - sa * ald) + dea * 2.0 * uq * DUXCON;
        }
    }

    // === Momentum equation (xblsys.f:1843-1899) ===
    let ha = 0.5 * (s1.h + s2.h);
    let ma_avg = 0.5 * msq;
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);

    // Get midpoint Cf
    let (cfm, _, _, _, _, _) = blmid(s1, s2, flow_type, msq);

    // CFX term
    let cfx = 0.5 * cfm * xa / ta + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);

    let btmp = ha + 2.0 - ma_avg;

    // Momentum residual: d(ln θ)/dx + (H+2-M²) * d(ln Ue)/dx = Cf/2 * x/θ
    res.res_mom = -(tlog + btmp * ulog - xlog * 0.5 * cfx);

    // Jacobian entries for momentum equation
    // ∂/∂T1, ∂/∂T2 (via tlog)
    jac.vs1[1][1] = 1.0 / s1.theta;
    jac.vs2[1][1] = -1.0 / s2.theta;
    // ∂/∂U1, ∂/∂U2 (via ulog)
    jac.vs1[1][3] = btmp / s1.u;
    jac.vs2[1][3] = -btmp / s2.u;
    // ∂/∂X1, ∂/∂X2 (via xlog)
    jac.vs1[1][4] = 0.5 * cfx / s1.x;
    jac.vs2[1][4] = -0.5 * cfx / s2.x;

    // === Shape parameter equation (xblsys.f:1901-1974) ===
    let hsa = 0.5 * (s1.hs + s2.hs);
    let hca = 0.5 * (s1.hc + s2.hc);

    let xot1 = s1.x / s1.theta;
    let xot2 = s2.x / s2.theta;

    // Upwinded DI and CF
    let dix = (1.0 - upw) * s1.cd * xot1 + upw * s2.cd * xot2;
    let cfx_shape = (1.0 - upw) * s1.cf * xot1 + upw * s2.cf * xot2;

    let btmp_shape = 2.0 * hca / hsa + 1.0 - ha;

    // Shape parameter residual
    res.res_shape = -(hlog + btmp_shape * ulog + xlog * (0.5 * cfx_shape - dix));

    // Jacobian entries for shape parameter equation
    // ∂/∂Hs1, ∂/∂Hs2 (via hlog)
    jac.vs1[2][1] = -hlog / s1.theta; // Simplified - Hs depends on θ
    jac.vs2[2][1] = hlog / s2.theta;
    // ∂/∂U1, ∂/∂U2 (via ulog)
    jac.vs1[2][3] = btmp_shape / s1.u;
    jac.vs2[2][3] = -btmp_shape / s2.u;
    // ∂/∂X1, ∂/∂X2
    jac.vs1[2][4] = -(0.5 * cfx_shape - dix) / s1.x;
    jac.vs2[2][4] = (0.5 * cfx_shape - dix) / s2.x;

    // DI contribution to Jacobian
    jac.vs1[2][0] = (1.0 - upw) * xot1 * xlog; // ∂/∂S1 (ctau affects DI for turbulent)
    jac.vs2[2][0] = upw * xot2 * xlog;

    (res, jac)
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // FlowType Tests
    // =========================================================================

    #[test]
    fn test_flow_type_equality() {
        assert_eq!(FlowType::Laminar, FlowType::Laminar);
        assert_ne!(FlowType::Laminar, FlowType::Turbulent);
        assert_ne!(FlowType::Turbulent, FlowType::Wake);
    }

    // =========================================================================
    // blvar() Tests
    // =========================================================================

    #[test]
    fn test_blvar_laminar_basic() {
        let mut station = BlStation::new();
        station.theta = 0.001;
        station.delta_star = 0.0026;
        station.u = 1.0;

        blvar(&mut station, FlowType::Laminar, 0.0, 1e6);

        // H should be δ*/θ = 2.6
        assert!(
            (station.h - 2.6).abs() < 0.01,
            "H should be δ*/θ, got {}",
            station.h
        );

        // Hk should equal H at M=0
        assert!(
            (station.hk - station.h).abs() < 0.1,
            "Hk should ≈ H at M=0, got Hk={}",
            station.hk
        );

        // Cf should be positive for laminar flow
        assert!(station.cf > 0.0, "Laminar Cf should be positive");

        // CD should be positive
        assert!(station.cd > 0.0, "CD should be positive");

        // Hs should be between 1.4 and 2.0 for laminar
        assert!(
            station.hs > 1.4 && station.hs < 2.0,
            "Hs should be reasonable, got {}",
            station.hs
        );
    }

    #[test]
    fn test_blvar_turbulent_basic() {
        let mut station = BlStation::new();
        station.theta = 0.002;
        station.delta_star = 0.003;
        station.u = 1.0;
        station.ctau = 0.1; // Non-zero shear stress for turbulent

        blvar(&mut station, FlowType::Turbulent, 0.0, 1e6);

        // H should be δ*/θ = 1.5
        assert!(
            (station.h - 1.5).abs() < 0.01,
            "H should be δ*/θ, got {}",
            station.h
        );

        // Cf should be positive for turbulent flow
        assert!(station.cf > 0.0, "Turbulent Cf should be positive");

        // CD should be positive and larger than laminar (outer layer contribution)
        assert!(station.cd > 0.0, "CD should be positive for turbulent");
    }

    #[test]
    fn test_blvar_wake_cf_zero() {
        let mut station = BlStation::new();
        station.theta = 0.003;
        station.delta_star = 0.005;
        station.u = 0.9;
        station.ctau = 0.05;

        blvar(&mut station, FlowType::Wake, 0.0, 1e6);

        // Wake should have Cf = 0 (no wall)
        assert_eq!(station.cf, 0.0, "Wake Cf should be exactly 0");

        // CD should still be positive (wake dissipation)
        assert!(station.cd > 0.0, "Wake CD should be positive");
    }

    #[test]
    fn test_blvar_hk_clamping() {
        let mut station = BlStation::new();
        // Set up conditions that would give Hk < 1.05
        station.theta = 0.01;
        station.delta_star = 0.01; // H = 1.0, which gives Hk ≈ 1.0
        station.u = 1.0;

        blvar(&mut station, FlowType::Laminar, 0.0, 1e6);

        // Hk should be clamped to minimum 1.05
        assert!(
            station.hk >= 1.05,
            "Hk should be clamped to >= 1.05, got {}",
            station.hk
        );
    }

    #[test]
    fn test_blvar_compressibility() {
        let mut station1 = BlStation::new();
        station1.theta = 0.001;
        station1.delta_star = 0.0026;
        station1.u = 1.0;

        let mut station2 = station1.clone();

        // Incompressible
        blvar(&mut station1, FlowType::Laminar, 0.0, 1e6);

        // Compressible (M ≈ 0.5)
        blvar(&mut station2, FlowType::Laminar, 0.25, 1e6);

        // Hk should differ due to compressibility correction
        assert!(
            (station1.hk - station2.hk).abs() > 0.01,
            "Hk should differ with compressibility"
        );

        // Hc should be 0 at M=0 and positive at M>0
        assert!(station1.hc.abs() < 1e-10, "Hc should be ~0 at M=0");
        assert!(station2.hc > 0.0, "Hc should be positive at M>0");
    }

    #[test]
    fn test_blvar_reynolds_dependence() {
        let mut station_low_re = BlStation::new();
        station_low_re.theta = 0.001;
        station_low_re.delta_star = 0.0026;
        station_low_re.u = 1.0;

        let mut station_high_re = station_low_re.clone();

        blvar(&mut station_low_re, FlowType::Laminar, 0.0, 1e5);
        blvar(&mut station_high_re, FlowType::Laminar, 0.0, 1e7);

        // Rθ should differ
        assert!(
            station_high_re.r_theta > station_low_re.r_theta,
            "Higher Re should give higher Rθ"
        );

        // Cf should be lower at higher Re (laminar scales as 1/sqrt(Re))
        // Actually laminar Cf scales as 1/Rθ, so higher Rθ = lower Cf
        assert!(
            station_high_re.cf < station_low_re.cf,
            "Higher Re should give lower laminar Cf"
        );
    }

    // =========================================================================
    // bldif() Tests
    // =========================================================================

    #[test]
    fn test_bldif_laminar_residuals() {
        let mut s1 = BlStation::new();
        s1.x = 0.1;
        s1.u = 1.0;
        s1.theta = 0.001;
        s1.delta_star = 0.0026;
        s1.ampl = 0.0;

        let mut s2 = BlStation::new();
        s2.x = 0.15;
        s2.u = 0.98;
        s2.theta = 0.0012;
        s2.delta_star = 0.0032;
        s2.ampl = 1.0;

        // Compute secondary variables
        blvar(&mut s1, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s2, FlowType::Laminar, 0.0, 1e6);

        let (res, jac) = bldif(&s1, &s2, FlowType::Laminar, 0.0, 1e6);

        // Residuals should be finite
        assert!(res.res_third.is_finite(), "res_third should be finite");
        assert!(res.res_mom.is_finite(), "res_mom should be finite");
        assert!(res.res_shape.is_finite(), "res_shape should be finite");

        // Jacobian should have non-zero entries
        let jac_sum: f64 = jac.vs1.iter().flatten().map(|x| x.abs()).sum::<f64>()
            + jac.vs2.iter().flatten().map(|x| x.abs()).sum::<f64>();
        assert!(jac_sum > 0.0, "Jacobian should have non-zero entries");
    }

    #[test]
    fn test_bldif_turbulent_residuals() {
        let mut s1 = BlStation::new();
        s1.x = 0.3;
        s1.u = 0.9;
        s1.theta = 0.003;
        s1.delta_star = 0.005;
        s1.ctau = 0.1;
        s1.is_laminar = false;
        s1.is_turbulent = true;

        let mut s2 = BlStation::new();
        s2.x = 0.35;
        s2.u = 0.85;
        s2.theta = 0.0035;
        s2.delta_star = 0.006;
        s2.ctau = 0.11;
        s2.is_laminar = false;
        s2.is_turbulent = true;

        // Compute secondary variables
        blvar(&mut s1, FlowType::Turbulent, 0.0, 1e6);
        blvar(&mut s2, FlowType::Turbulent, 0.0, 1e6);

        let (res, jac) = bldif(&s1, &s2, FlowType::Turbulent, 0.0, 1e6);

        // Residuals should be finite
        assert!(res.res_third.is_finite(), "res_third should be finite");
        assert!(res.res_mom.is_finite(), "res_mom should be finite");
        assert!(res.res_shape.is_finite(), "res_shape should be finite");

        // Jacobian should have entries for ctau (s variable)
        assert!(
            jac.vs1[0][0].abs() > 0.0 || jac.vs2[0][0].abs() > 0.0,
            "Jacobian should have ctau dependencies for turbulent"
        );
    }

    #[test]
    fn test_bldif_wake_residuals() {
        let mut s1 = BlStation::new();
        s1.x = 1.0;
        s1.u = 0.5;
        s1.theta = 0.01;
        s1.delta_star = 0.02;
        s1.ctau = 0.05;
        s1.is_wake = true;

        let mut s2 = BlStation::new();
        s2.x = 1.1;
        s2.u = 0.52;
        s2.theta = 0.011;
        s2.delta_star = 0.021;
        s2.ctau = 0.04;
        s2.is_wake = true;

        // Compute secondary variables
        blvar(&mut s1, FlowType::Wake, 0.0, 1e6);
        blvar(&mut s2, FlowType::Wake, 0.0, 1e6);

        let (res, _jac) = bldif(&s1, &s2, FlowType::Wake, 0.0, 1e6);

        // Residuals should be finite
        assert!(res.res_third.is_finite(), "res_third should be finite");
        assert!(res.res_mom.is_finite(), "res_mom should be finite");
        assert!(res.res_shape.is_finite(), "res_shape should be finite");
    }

    #[test]
    fn test_bldif_momentum_equation() {
        // Test that momentum equation behaves correctly
        // For flat plate with constant Ue, d(ln θ)/dx ≈ Cf/2

        let mut s1 = BlStation::new();
        s1.x = 0.1;
        s1.u = 1.0;
        s1.theta = 0.001;
        s1.delta_star = 0.0026;

        let mut s2 = BlStation::new();
        s2.x = 0.11;
        s2.u = 1.0; // Constant velocity
        s2.theta = 0.00105; // Slight growth
        s2.delta_star = 0.00273;

        blvar(&mut s1, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s2, FlowType::Laminar, 0.0, 1e6);

        let (res, _) = bldif(&s1, &s2, FlowType::Laminar, 0.0, 1e6);

        // With constant velocity, the ulog term vanishes
        // Residual should primarily reflect θ growth vs Cf
        assert!(
            res.res_mom.abs() < 10.0,
            "Momentum residual should be bounded, got {}",
            res.res_mom
        );
    }

    #[test]
    fn test_bldif_jacobian_structure() {
        let mut s1 = BlStation::new();
        s1.x = 0.1;
        s1.u = 1.0;
        s1.theta = 0.001;
        s1.delta_star = 0.0026;

        let mut s2 = BlStation::new();
        s2.x = 0.15;
        s2.u = 0.98;
        s2.theta = 0.0012;
        s2.delta_star = 0.0032;

        blvar(&mut s1, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s2, FlowType::Laminar, 0.0, 1e6);

        let (_, jac) = bldif(&s1, &s2, FlowType::Laminar, 0.0, 1e6);

        // Check that Jacobian has expected structure
        // For laminar, vs1[0][0] should be -1 (∂res/∂ampl1)
        assert!(
            (jac.vs1[0][0] - (-1.0)).abs() < 0.01,
            "vs1[0][0] should be -1 for ampl equation"
        );

        // vs2[0][0] should be +1 (∂res/∂ampl2)
        assert!(
            (jac.vs2[0][0] - 1.0).abs() < 0.01,
            "vs2[0][0] should be +1 for ampl equation"
        );
    }

    // =========================================================================
    // BlResiduals Tests
    // =========================================================================

    #[test]
    fn test_bl_residuals_default() {
        let res = BlResiduals::default();
        assert_eq!(res.res_third, 0.0);
        assert_eq!(res.res_mom, 0.0);
        assert_eq!(res.res_shape, 0.0);
    }

    // =========================================================================
    // BlJacobian Tests
    // =========================================================================

    #[test]
    fn test_bl_jacobian_default() {
        let jac = BlJacobian::default();

        // All entries should be zero
        for row in &jac.vs1 {
            for &val in row {
                assert_eq!(val, 0.0);
            }
        }
        for row in &jac.vs2 {
            for &val in row {
                assert_eq!(val, 0.0);
            }
        }
        for &val in &jac.vsm {
            assert_eq!(val, 0.0);
        }
        for &val in &jac.vsr {
            assert_eq!(val, 0.0);
        }
    }

    // =========================================================================
    // Numerical Derivative Verification
    // =========================================================================

    #[test]
    fn test_blvar_derivative_chain_h() {
        // Verify ∂H/∂θ and ∂H/∂δ* numerically
        let eps = 1e-7;

        let mut station = BlStation::new();
        station.theta = 0.001;
        station.delta_star = 0.0026;
        station.u = 1.0;

        // Perturb theta
        let mut s_plus = station.clone();
        s_plus.theta = station.theta + eps;
        let mut s_minus = station.clone();
        s_minus.theta = station.theta - eps;

        blvar(&mut station, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s_plus, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s_minus, FlowType::Laminar, 0.0, 1e6);

        let dh_dtheta_num = (s_plus.h - s_minus.h) / (2.0 * eps);
        let dh_dtheta_ana = station.derivs.h_theta;

        assert!(
            (dh_dtheta_num - dh_dtheta_ana).abs() < 1e-4,
            "∂H/∂θ: numerical={}, analytical={}",
            dh_dtheta_num,
            dh_dtheta_ana
        );
    }

    #[test]
    fn test_blvar_derivative_chain_delta_star() {
        let eps = 1e-7;

        let mut station = BlStation::new();
        station.theta = 0.001;
        station.delta_star = 0.0026;
        station.u = 1.0;

        // Perturb delta_star
        let mut s_plus = station.clone();
        s_plus.delta_star = station.delta_star + eps;
        let mut s_minus = station.clone();
        s_minus.delta_star = station.delta_star - eps;

        blvar(&mut station, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s_plus, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s_minus, FlowType::Laminar, 0.0, 1e6);

        let dh_dd_num = (s_plus.h - s_minus.h) / (2.0 * eps);
        let dh_dd_ana = station.derivs.h_delta_star;

        assert!(
            (dh_dd_num - dh_dd_ana).abs() < 1e-4,
            "∂H/∂δ*: numerical={}, analytical={}",
            dh_dd_num,
            dh_dd_ana
        );
    }

    // =========================================================================
    // Physical Sanity Tests
    // =========================================================================

    #[test]
    fn test_blvar_physical_bounds() {
        // Test that blvar produces physically reasonable values across conditions
        for theta in [0.0005, 0.001, 0.002, 0.005] {
            for h_ratio in [1.5, 2.0, 2.6, 3.5] {
                let mut station = BlStation::new();
                station.theta = theta;
                station.delta_star = theta * h_ratio;
                station.u = 1.0;

                blvar(&mut station, FlowType::Laminar, 0.0, 1e6);

                // Check physical bounds
                assert!(station.h > 1.0, "H should be > 1");
                assert!(station.hk > 1.0, "Hk should be > 1");
                assert!(station.hs > 1.0, "Hs should be > 1");
                assert!(station.cf > 0.0 || station.cf < 0.0, "Cf should be finite");
                assert!(station.cd > 0.0, "CD should be positive");
                assert!(station.r_theta > 0.0, "Rθ should be positive");
            }
        }
    }

    #[test]
    fn test_flow_type_affects_results() {
        let mut s_lam = BlStation::new();
        s_lam.theta = 0.002;
        s_lam.delta_star = 0.003;
        s_lam.u = 1.0;
        s_lam.ctau = 0.1;

        let mut s_turb = s_lam.clone();
        let mut s_wake = s_lam.clone();

        blvar(&mut s_lam, FlowType::Laminar, 0.0, 1e6);
        blvar(&mut s_turb, FlowType::Turbulent, 0.0, 1e6);
        blvar(&mut s_wake, FlowType::Wake, 0.0, 1e6);

        // Cf should be zero for wake
        assert_eq!(s_wake.cf, 0.0, "Wake Cf should be 0");

        // Cf for turbulent and laminar should differ
        assert!(
            (s_lam.cf - s_turb.cf).abs() > 1e-10,
            "Laminar and turbulent Cf should differ"
        );

        // Hs should differ between laminar and turbulent
        assert!(
            (s_lam.hs - s_turb.hs).abs() > 1e-6,
            "Laminar and turbulent Hs should differ"
        );
    }
}
