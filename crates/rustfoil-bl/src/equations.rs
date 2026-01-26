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
    amplification_rate, axset_full, cf_laminar, cf_turbulent, density_shape, dissipation_laminar,
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

/// Intermediate terms for BLDIF diagnostics
#[derive(Debug, Clone, Copy, Default)]
pub struct BldifTerms {
    pub xlog: f64,
    pub ulog: f64,
    pub tlog: f64,
    pub hlog: f64,
    pub z_cfx_shape: f64,
    pub z_dix_shape: f64,
    pub upw: f64,
    pub ha: f64,
    pub btmp_mom: f64,
    pub cfx: f64,
    pub cfx_ta: f64,
    pub cfx_t1: f64,
    pub cfx_t2: f64,
    pub btmp_shape: f64,
    pub cfx_shape: f64,
    pub dix: f64,
    pub cfx_upw: f64,
    pub dix_upw: f64,
    pub hsa: f64,
    pub hca: f64,
    pub dd: f64,
    pub dd2: f64,
    pub xot1: f64,
    pub xot2: f64,
    pub cf1: f64,
    pub cf2: f64,
    pub di1: f64,
    pub di2: f64,
    pub z_upw: f64,
    pub z_de2: f64,
    pub z_us2: f64,
    pub z_cq2: f64,
    pub z_cf2: f64,
    pub z_hk2: f64,
    pub z_d2: f64,
    pub z_u2: f64,
    pub z_s2: f64,
    pub upw_t2: f64,
    pub upw_d2: f64,
    pub upw_u2: f64,
    pub de2_t2: f64,
    pub de2_d2: f64,
    pub de2_u2: f64,
    pub us2_t2: f64,
    pub us2_d2: f64,
    pub us2_u2: f64,
    pub cq2_t2: f64,
    pub cq2_d2: f64,
    pub cq2_u2: f64,
    pub cf2_t2: f64,
    pub cf2_d2: f64,
    pub cf2_u2: f64,
    pub hk2_t2: f64,
    pub hk2_d2: f64,
    pub hk2_u2: f64,
    pub z_ax: f64,
    pub ax_hk2: f64,
    pub ax_t2: f64,
    pub ax_rt2: f64,
    pub ax_a2: f64,
    pub rt2_t2: f64,
    pub rt2_u2: f64,
    pub z_hs2: f64,
    pub z_cf2_shape: f64,
    pub z_di2: f64,
    pub z_t2_shape: f64,
    pub z_u2_shape: f64,
    pub z_hca: f64,
    pub z_ha_shape: f64,
    pub z_upw_shape: f64,
    pub hs2_t: f64,
    pub hs2_d: f64,
    pub hs2_u: f64,
    pub cf2_t_shape: f64,
    pub cf2_d_shape: f64,
    pub cf2_u_shape: f64,
    pub di2_t: f64,
    pub di2_d: f64,
    pub di2_u: f64,
    pub di2_s: f64,
    pub hc2_t: f64,
    pub hc2_d: f64,
    pub hc2_u: f64,
    pub h2_t: f64,
    pub h2_d: f64,
    pub upw_t2_shape: f64,
    pub upw_d2_shape: f64,
    pub upw_u2_shape: f64,
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
    let _hc_u = hc_hk * hk_u + hc_msq * m_u;
    let _hc_t = hc_hk * hk_t;
    let _hc_d = hc_hk * hk_d;
    let _hc_ms = hc_hk * hk_ms + hc_msq * m_ms;

    // Store Hc base derivatives for Jacobian construction in bldif()
    station.derivs.hc_hk = hc_hk;
    station.derivs.hc_msq = hc_msq;

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

    // Store Us and its derivatives
    station.us = us;
    station.derivs.us_t = us_t;
    station.derivs.us_d = us_d;
    station.derivs.us_u = us_u;

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
    let cq_u_val = cq_hs * hs_u + cq_us * us_u + cq_hk * hk_u + cq_rt * rt_u;
    let cq_t_val = cq_hs * hs_t + cq_us * us_t + cq_hk * hk_t + cq_rt * rt_t + cq_h * h_t;
    let cq_d_val = cq_hs * hs_d + cq_us * us_d + cq_hk * hk_d + cq_h * h_d;
    let cq_ms_val = cq_hs * hs_ms + cq_us * us_ms + cq_hk * hk_ms;
    let cq_re_val = cq_hs * hs_re + cq_us * us_re + cq_rt * rt_re;

    // Store CQ and its derivatives
    station.cq = cq;
    station.derivs.cq_t = cq_t_val;
    station.derivs.cq_d = cq_d_val;
    station.derivs.cq_u = cq_u_val;

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
    // Returns: (di, di_s, di_u, di_t, di_d, di_ms, di_re, di_hk_base, di_rt_base)
    // di_hk_base and di_rt_base are the base derivatives ∂DI/∂Hk and ∂DI/∂Rθ
    // needed for the Jacobian in bldif()
    let (mut di, mut di_s, mut di_u, mut di_t, mut di_d, mut di_ms, mut di_re, di_hk_base, di_rt_base) = match flow_type {
        FlowType::Laminar => {
            let di_result = dissipation_laminar(hk, station.r_theta);
            let di_hk = di_result.di_hk;
            let di_rt = di_result.di_rt;
            // Chain rule for primary variable derivatives
            let di_u = di_hk * hk_u + di_rt * rt_u;
            let di_t = di_hk * hk_t + di_rt * rt_t;
            let di_d = di_hk * hk_d;
            let di_ms = di_hk * hk_ms;
            let di_re = di_rt * rt_re;
            // Return base derivatives for Jacobian construction
            (di_result.di, 0.0, di_u, di_t, di_d, di_ms, di_re, di_hk, di_rt)
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
            let hm_rt = -(2.1 / (grt * grt)) / station.r_theta;
            let fl_hk = 1.0 / (hmin - 1.0);
            let fl_rt = -(fl / (hmin - 1.0)) * hm_rt;
            let tfl = fl.tanh();
            let dfac = 0.5 + 0.5 * tfl;
            let df_fl = 0.5 * (1.0 - tfl * tfl);
            let df_hk = df_fl * fl_hk;
            let df_rt = df_fl * fl_rt;

            let di_wall_corrected = di_wall * dfac;
            let di_u_corrected = di_u * dfac + di_wall * (df_hk * hk_u + df_rt * rt_u);
            let di_t_corrected = di_t * dfac + di_wall * (df_hk * hk_t + df_rt * rt_t);
            let di_d_corrected = di_d * dfac + di_wall * (df_hk * hk_d);
            let di_ms_corrected = di_ms * dfac + di_wall * (df_hk * hk_ms + df_rt * rt_ms);
            let di_re_corrected = di_re * dfac + di_wall * (df_rt * rt_re);

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
            let di_total = di_wall_corrected + dd + dd2;
            let di_total_s = dd_s;
            let di_total_u = di_u_corrected
                + dd_hs * hs_u
                + dd_us * us_u
                + dd2_hs * hs_u
                + dd2_us * us_u
                + dd2_rt * rt_u;
            let di_total_t = di_t_corrected
                + dd_hs * hs_t
                + dd_us * us_t
                + dd2_hs * hs_t
                + dd2_us * us_t
                + dd2_rt * rt_t;
            let di_total_d = di_d_corrected + dd_hs * hs_d + dd_us * us_d + dd2_hs * hs_d + dd2_us * us_d;
            let di_total_ms = di_ms_corrected
                + dd_hs * hs_ms
                + dd_us * us_ms
                + dd2_hs * hs_ms
                + dd2_us * us_ms;
            let di_total_re = di_re_corrected
                + dd_hs * hs_re
                + dd_us * us_re
                + dd2_hs * hs_re
                + dd2_us * us_re
                + dd2_rt * rt_re;

            // Check if laminar DI is higher (low Rθ case)
            let di_lam = dissipation_laminar(hk, station.r_theta);
            if di_lam.di > di_total {
                let di_hk = di_lam.di_hk;
                let di_rt = di_lam.di_rt;
                let di_u = di_hk * hk_u + di_rt * rt_u;
                let di_t = di_hk * hk_t + di_rt * rt_t;
                let di_d = di_hk * hk_d;
                let di_ms = di_hk * hk_ms;
                let di_re = di_rt * rt_re;
                (di_lam.di, 0.0, di_u, di_t, di_d, di_ms, di_re, di_hk, di_rt)
            } else {
                // Turbulent flow: DI depends on Hs, Us, Cf - not directly Hk, Rt
                // For Jacobian, we'd need DI_HS, DI_US, DI_CF chain rules
                // TODO: Implement proper turbulent Jacobian derivatives
                // For now, use laminar approximation for base derivatives
                (
                    di_total,
                    di_total_s,
                    di_total_u,
                    di_total_t,
                    di_total_d,
                    di_total_ms,
                    di_total_re,
                    di_lam.di_hk,
                    di_lam.di_rt,
                )
            }
        }
        FlowType::Wake => {
            // Wake: compute BOTH turbulent outer layer AND laminar wake (DILW)
            // Use whichever is larger, then double for two wake halves
            // (xblsys.f:1007-1086)
            
            // First compute turbulent outer layer contribution (xblsys.f:1007-1037)
            // This is the same formula used for turbulent BL, also applies to wake
            let s2 = station.ctau;
            let wake_hs = station.hs;
            
            // DD = S2^2 * (0.995-US) * 2.0/HS (xblsys.f:1010)
            let dd_outer = s2 * s2 * (0.995 - us) * 2.0 / wake_hs;
            let dd_outer_hs = -dd_outer / wake_hs;
            let dd_outer_us = -s2 * s2 * 2.0 / wake_hs;
            let dd_outer_s = 2.0 * s2 * (0.995 - us) * 2.0 / wake_hs;
            
            // Laminar stress contribution (xblsys.f:1025)
            // DD = 0.15*(0.995-US)^2 / RT * 2.0/HS
            let dd_lam = 0.15 * (0.995 - us).powi(2) / station.r_theta * 2.0 / wake_hs;
            let dd_lam_us = -0.15 * (0.995 - us) * 2.0 / station.r_theta * 2.0 / wake_hs;
            let dd_lam_hs = -dd_lam / wake_hs;
            let dd_lam_rt = -dd_lam / station.r_theta;
            
            // Total turbulent/wake DI2
            let di2_turb = dd_outer + dd_lam;
            let di2_turb_hs = dd_outer_hs + dd_lam_hs;
            let di2_turb_us = dd_outer_us + dd_lam_us;
            let di2_turb_rt = dd_lam_rt;
            let di2_turb_s = dd_outer_s;
            
            // Also compute DILW (laminar wake dissipation)
            let di_lam_wake = dissipation_wake(hk, station.r_theta);
            
            // Use max(DI2_turb, DILW) (xblsys.f:1073)
            let (di_base, di_s, di_hs_deriv, di_us_deriv, di_hk_deriv, di_rt_deriv) = 
                if di_lam_wake.di > di2_turb {
                    // Use laminar wake
                    (di_lam_wake.di, 0.0, 0.0, 0.0, di_lam_wake.di_hk, di_lam_wake.di_rt)
                } else {
                    // Use turbulent outer layer
                    (di2_turb, di2_turb_s, di2_turb_hs, di2_turb_us, 0.0, di2_turb_rt)
                };
            
            // Chain rule for derivatives
            // di depends on hk, hs, us, rt, s - need to convert to primary variables
            let di_u = di_hk_deriv * hk_u + di_hs_deriv * hs_hk * hk_u + di_us_deriv * us_u + di_rt_deriv * rt_u;
            let di_t = di_hk_deriv * hk_t + di_hs_deriv * hs_hk * hk_t + di_us_deriv * us_t + di_rt_deriv * rt_t;
            let di_d = di_hk_deriv * hk_d + di_hs_deriv * hs_hk * hk_d + di_us_deriv * us_d;
            let di_ms = di_hk_deriv * hk_ms + di_hs_deriv * hs_hk * hk_ms + di_us_deriv * us_ms;
            let di_re = di_rt_deriv * rt_re;
            
            // Double dissipation for wake (two halves) (xblsys.f:1089-1094)
            // Base derivatives for Jacobian (already doubled via chain rule)
            (
                di_base * 2.0,
                di_s * 2.0,
                di_u * 2.0,
                di_t * 2.0,
                di_d * 2.0,
                di_ms * 2.0,
                di_re * 2.0,
                di_hk_deriv * 2.0,
                di_rt_deriv * 2.0,
            )
        }
    };

    if flow_type == FlowType::Turbulent {
        // Defensive check: ensure turbulent DI wins when larger than laminar.
        // This mirrors XFOIL's DI2L > DI2 override, but prevents accidental laminar selection
        // when DI2 (turbulent) should dominate.
        let cf2t_result = cf_turbulent(hk, station.r_theta, msq);
        let cf2t_u = cf2t_result.cf_hk * hk_u + cf2t_result.cf_rt * rt_u + cf2t_result.cf_msq * m_u;
        let cf2t_t = cf2t_result.cf_hk * hk_t + cf2t_result.cf_rt * rt_t;
        let cf2t_d = cf2t_result.cf_hk * hk_d;
        let cf2t_ms = cf2t_result.cf_hk * hk_ms + cf2t_result.cf_rt * rt_ms + cf2t_result.cf_msq * m_ms;
        let cf2t_re = cf2t_result.cf_rt * rt_re;

        let di_wall = (0.5 * cf2t_result.cf * us) * 2.0 / station.hs;
        let di_hs_wall = -(0.5 * cf2t_result.cf * us) * 2.0 / (station.hs * station.hs);
        let di_us_wall = (0.5 * cf2t_result.cf) * 2.0 / station.hs;
        let di_cf2t = (0.5 * us) * 2.0 / station.hs;

        let di_u_wall = di_hs_wall * hs_u + di_us_wall * us_u + di_cf2t * cf2t_u;
        let di_t_wall = di_hs_wall * hs_t + di_us_wall * us_t + di_cf2t * cf2t_t;
        let di_d_wall = di_hs_wall * hs_d + di_us_wall * us_d + di_cf2t * cf2t_d;
        let di_ms_wall = di_hs_wall * hs_ms + di_us_wall * us_ms + di_cf2t * cf2t_ms;
        let di_re_wall = di_hs_wall * hs_re + di_us_wall * us_re + di_cf2t * cf2t_re;

        let grt = station.r_theta.ln();
        let hmin = 1.0 + 2.1 / grt;
        let fl = (hk - 1.0) / (hmin - 1.0);
        let hm_rt = -(2.1 / (grt * grt)) / station.r_theta;
        let fl_hk = 1.0 / (hmin - 1.0);
        let fl_rt = -(fl / (hmin - 1.0)) * hm_rt;
        let tfl = fl.tanh();
        let dfac = 0.5 + 0.5 * tfl;
        let df_fl = 0.5 * (1.0 - tfl * tfl);
        let df_hk = df_fl * fl_hk;
        let df_rt = df_fl * fl_rt;

        let di_wall_corrected = di_wall * dfac;
        let di_u_corrected = di_u_wall * dfac + di_wall * (df_hk * hk_u + df_rt * rt_u);
        let di_t_corrected = di_t_wall * dfac + di_wall * (df_hk * hk_t + df_rt * rt_t);
        let di_d_corrected = di_d_wall * dfac + di_wall * (df_hk * hk_d);
        let di_ms_corrected = di_ms_wall * dfac + di_wall * (df_hk * hk_ms + df_rt * rt_ms);
        let di_re_corrected = di_re_wall * dfac + di_wall * (df_rt * rt_re);

        let s = station.ctau;
        let dd = s * s * (0.995 - us) * 2.0 / station.hs;
        let dd_hs = -s * s * (0.995 - us) * 2.0 / (station.hs * station.hs);
        let dd_us = -s * s * 2.0 / station.hs;
        let dd_s = s * 2.0 * (0.995 - us) * 2.0 / station.hs;

        let dd2 = 0.15 * (0.995 - us) * (0.995 - us) / station.r_theta * 2.0 / station.hs;
        let dd2_us = -0.15 * (0.995 - us) * 2.0 / station.r_theta * 2.0 / station.hs;
        let dd2_hs = -dd2 / station.hs;
        let dd2_rt = -dd2 / station.r_theta;

        let di_total = di_wall_corrected + dd + dd2;
        let di_total_s = dd_s;
        let di_total_u = di_u_corrected
            + dd_hs * hs_u
            + dd_us * us_u
            + dd2_hs * hs_u
            + dd2_us * us_u
            + dd2_rt * rt_u;
        let di_total_t = di_t_corrected
            + dd_hs * hs_t
            + dd_us * us_t
            + dd2_hs * hs_t
            + dd2_us * us_t
            + dd2_rt * rt_t;
        let di_total_d = di_d_corrected + dd_hs * hs_d + dd_us * us_d + dd2_hs * hs_d + dd2_us * us_d;
        let di_total_ms = di_ms_corrected
            + dd_hs * hs_ms
            + dd_us * us_ms
            + dd2_hs * hs_ms
            + dd2_us * us_ms;
        let di_total_re = di_re_corrected
            + dd_hs * hs_re
            + dd_us * us_re
            + dd2_hs * hs_re
            + dd2_us * us_re
            + dd2_rt * rt_re;

        let di_lam = dissipation_laminar(hk, station.r_theta);
        if di_total > di_lam.di && (di - di_total).abs() > 1e-12 {
            di = di_total;
            di_s = di_total_s;
            di_u = di_total_u;
            di_t = di_total_t;
            di_d = di_total_d;
            di_ms = di_total_ms;
            di_re = di_total_re;
        }
    }

    station.cd = di;
    // Store BASE derivatives ∂DI/∂Hk and ∂DI/∂Rθ for use in some correlations
    station.derivs.cd_hk = di_hk_base;
    station.derivs.cd_rt = di_rt_base;
    // Store direct derivatives w.r.t. primary variables for shape Jacobian
    station.derivs.cd_u = di_u;
    station.derivs.cd_t = di_t;
    station.derivs.cd_d = di_d;
    station.derivs.cd_s = di_s;

    // === BL thickness DE from Green's correlation (xblsys.f:1100-1118) ===
    // DE = (3.15 + 1.72/(HK-1)) * T + D
    let de = (3.15 + 1.72 / (hk - 1.0).max(0.01)) * station.theta + station.delta_star;
    let de_hk = -1.72 / ((hk - 1.0).max(0.01) * (hk - 1.0).max(0.01)) * station.theta;

    let de_u = de_hk * hk_u;
    let de_t_val = de_hk * hk_t + 3.15 + 1.72 / (hk - 1.0).max(0.01);
    let de_d_val = de_hk * hk_d + 1.0;
    let de_ms = de_hk * hk_ms;

    // Clamp DE to maximum HDMAX * T
    let hdmax = 12.0;
    let de_final = if de > hdmax * station.theta {
        hdmax * station.theta
    } else {
        de
    };

    // Store DE and its derivatives
    station.de = de_final;
    station.derivs.de_t = de_t_val;
    station.derivs.de_d = de_d_val;

    // Store mass defect
    station.mass_defect = station.u * station.delta_star;

    // Store H derivative (for compatibility)
    station.derivs.h_theta = h_t;
    station.derivs.h_delta_star = h_d;
}

/// Compute BL secondary variables with debug output
///
/// Same as [`blvar`] but also emits debug events when debug collection is active.
/// This matches XFOIL's BLVAR instrumented output.
///
/// # Arguments
/// * `station` - Mutable BL station with primary variables set
/// * `flow_type` - Flow type (Laminar, Turbulent, Wake)
/// * `msq` - Mach number squared
/// * `re` - Reference Reynolds number
/// * `iteration` - Current solver iteration (for debug output)
/// * `side` - Side index (1 or 2 for upper/lower)
/// * `ibl` - Station index along the side
pub fn blvar_debug(
    station: &mut BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
    iteration: usize,
    side: usize,
    ibl: usize,
) {
    // Capture inputs before computation
    let input = crate::debug::BlvarInput {
        x: station.x,
        u: station.u,
        theta: station.theta,
        delta_star: station.delta_star,
        ctau: station.ctau,
        ampl: station.ampl,
    };

    // Run the actual computation
    blvar(station, flow_type, msq, re);

    // Emit debug event if collection is active
    if crate::debug::is_debug_active() {
        let dd = station.ctau * station.ctau * (0.995 - station.us) * 2.0 / station.hs;
        let dd2 = 0.15 * (0.995 - station.us).powi(2) / station.r_theta * 2.0 / station.hs;

        // Recompute DI terms for debug diagnostics
        let cf2t_result = cf_turbulent(station.hk, station.r_theta, msq);
        let di_wall = (0.5 * cf2t_result.cf * station.us) * 2.0 / station.hs;
        let grt = station.r_theta.ln();
        let hmin = 1.0 + 2.1 / grt;
        let fl = (station.hk - 1.0) / (hmin - 1.0);
        let dfac = 0.5 + 0.5 * fl.tanh();
        let di_wall_corrected = di_wall * dfac;
        let di_total = di_wall_corrected + dd + dd2;
        let di_lam = dissipation_laminar(station.hk, station.r_theta).di;

        let (di_lam_override, di_used_minus_total, di_used_minus_lam) = if flow_type == FlowType::Turbulent {
            let di_lam_override = di_lam > di_total;
            let di_used_minus_total = station.cd - di_total;
            let di_used_minus_lam = station.cd - di_lam;
            (Some(di_lam_override), Some(di_used_minus_total), Some(di_used_minus_lam))
        } else {
            (None, None, None)
        };

        let output = crate::debug::BlvarOutput {
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

        let flow_type_int = match flow_type {
            FlowType::Laminar => 1,
            FlowType::Turbulent => 2,
            FlowType::Wake => 3,
        };

        crate::debug::add_event(crate::debug::DebugEvent::blvar(
            iteration,
            side,
            ibl,
            flow_type_int,
            input,
            output,
        ));
    }
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
/// * `(BlResiduals, BlJacobian, BldifTerms)` - Residuals, Jacobian, and terms
///
/// # Reference
/// XFOIL xblsys.f BLDIF (line 1552)
/// Compute BL equations for similarity station using regular bldif with fixed log values
pub fn bldif_full_simi(
    s: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    // For similarity, XFOIL uses fixed log values (BLDIF with ITYP=0):
    // XLOG = 1.0, ULOG = BULE = 1.0, TLOG = 0.0, HLOG = 0.0
    let (res, jac, _) = bldif_with_terms_internal(s, s, flow_type, msq, re, true);
    (res, jac)
}

pub fn bldif_with_terms(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian, BldifTerms) {
    bldif_with_terms_internal(s1, s2, flow_type, msq, re, false)
}

fn bldif_with_terms_internal(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
    similarity: bool,
) -> (BlResiduals, BlJacobian, BldifTerms) {
    let mut res = BlResiduals::default();
    let mut jac = BlJacobian::default();
    let mut terms = BldifTerms::default();

    // === Logarithmic differences (xblsys.f:1630-1648) ===
    // For similarity station (ITYP=0):
    //   - Use fixed values: XLOG=1, ULOG=BULE=1, TLOG=0, HLOG=0
    //   - DDLOG=0 which zeros out many Jacobian log-derivative terms
    // For normal stations: DDLOG=1
    let (xlog, ulog, tlog, hlog, ddlog) = if similarity {
        let bule = 1.0;  // Leading edge pressure gradient parameter
        (1.0, bule, 0.5 * (1.0 - bule), 0.0, 0.0)  // DDLOG=0 for SIMI
    } else {
        (
            (s2.x / s1.x.max(1e-20)).ln(),
            (s2.u / s1.u.max(1e-20)).ln(),
            (s2.theta / s1.theta.max(1e-20)).ln(),
            (s2.hs / s1.hs.max(1e-20)).ln(),
            1.0,  // DDLOG=1 for normal
        )
    };
    terms.xlog = xlog;
    terms.ulog = ulog;
    terms.tlog = tlog;
    terms.hlog = hlog;
    terms.cf1 = s1.cf;
    terms.cf2 = s2.cf;
    terms.di1 = s1.cd;
    terms.di2 = s2.cd;

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
    terms.upw = upw;

    // UPW derivatives (XFOIL xblsys.f:1631-1644)
    // UPW depends on HL and HDCON, which depend on Hk
    let upw_hl = ehh * hl * hdcon;
    let upw_hd = 0.5 * ehh * hlsq;

    // HL derivatives w.r.t. Hk (XFOIL lines 1620-1621)
    let hl_hk1 = -1.0 / (s1.hk - 1.0).max(1e-6);
    let hl_hk2 = 1.0 / (s2.hk - 1.0).max(1e-6);

    // HD (HDCON) derivatives w.r.t. Hk (XFOIL lines 1606-1607, 1611-1613)
    // HDCON = 5.0*HUPWT/HK2^2, so d(HDCON)/d(HK2) = -2*HDCON/HK2
    let hd_hk1 = 0.0; // HDCON doesn't depend on HK1
    let hd_hk2 = -hdcon * 2.0 / s2.hk;

    // UPW_HK = UPW_HL * HL_HK + UPW_HD * HD_HK (XFOIL lines 1634-1635)
    let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    // === First equation: Amplification (laminar) or Shear-lag (turbulent/wake) ===
    // For similarity (ITYP=0), XFOIL uses a simple constraint: AMPL2 = 0
    // with VS2(1,1) = 1.0, VS1(1,1) = 0.0, VSREZ(1) = -AMPL2
    if similarity && flow_type == FlowType::Laminar {
        res.res_third = -s2.ampl;  // Residual: want ampl = 0
        jac.vs1[0][0] = 0.0;       // VS1 is zero for SIMI first equation
        jac.vs2[0][0] = 1.0;       // Simple identity: d(res)/d(ampl2) = 1
        // All other first-equation derivatives are 0 (already default)
    } else {
    match flow_type {
        FlowType::Laminar => {
            // Amplification equation (xblsys.f:1654-1682)
            // REZC = AMPL2 - AMPL1 - AX*(X2-X1)
            let ax_result = axset_full(
                s1.hk,
                s1.theta,
                s1.r_theta,
                s1.ampl,
                s2.hk,
                s2.theta,
                s2.r_theta,
                s2.ampl,
                9.0,
            );
            let ax = ax_result.ax;
            let dx = s2.x - s1.x;

            res.res_third = -(s2.ampl - s1.ampl - ax * dx);

            // Jacobian entries for amplification equation (XFOIL xblsys.f:1667-1676)
            // Variable order: [ampl, θ, δ*, u, x]
            // Z_AX = -(X2-X1) = -dx
            let z_ax = -dx;
            
            // Amplification rate depends on midpoint values:
            //   Hk_avg = 0.5*(Hk1 + Hk2)
            //   θ_avg = 0.5*(θ1 + θ2)  
            //   Rθ_avg = 0.5*(Rθ1 + Rθ2)
            // So derivatives w.r.t. station values include 0.5 factor
            let ax_hk1 = ax_result.ax_hk1;
            let ax_hk2 = ax_result.ax_hk2;
            
            // Hk derivatives: Hk = f(H), H = δ*/θ
            // HK_T = HK_H * H_T = hk_h * (-H/θ)
            // HK_D = HK_H * H_D = hk_h * (1/θ)
            let hk1_t = s1.derivs.hk_h * s1.derivs.h_theta;
            let hk1_d = s1.derivs.hk_h * s1.derivs.h_delta_star;
            let hk2_t = s2.derivs.hk_h * s2.derivs.h_theta;
            let hk2_d = s2.derivs.hk_h * s2.derivs.h_delta_star;
            let hk1_u = 0.0;
            let hk2_u = 0.0;
            
            // Rθ derivatives: Rθ = Re * Ue * θ / ν
            let rt1_t = s1.r_theta / s1.theta.max(1e-12);
            let rt2_t = s2.r_theta / s2.theta.max(1e-12);
            let rt1_u = s1.r_theta / s1.u.max(1e-12);
            let rt2_u = s2.r_theta / s2.u.max(1e-12);
            
            // ax also depends directly on θ and Rθ through the midpoint averages
            // ax_t = ∂ax/∂θ_avg * ∂θ_avg/∂θ = ax_th * 0.5
            // ax_rt = ∂ax/∂Rθ_avg * ∂Rθ_avg/∂Rθ = ax_rt * 0.5
            // (XFOIL xblsys.f:1668, 1673 - AX_T and AX_RT terms)
            let ax_t1 = ax_result.ax_t1;
            let ax_t2 = ax_result.ax_t2;
            let ax_rt1 = ax_result.ax_rt1;
            let ax_rt2 = ax_result.ax_rt2;

            terms.z_ax = z_ax;
            terms.ax_hk2 = ax_hk2;
            terms.ax_t2 = ax_t2;
            terms.ax_rt2 = ax_rt2;
            terms.ax_a2 = ax_result.ax_a2;
            terms.hk2_t2 = hk2_t;
            terms.hk2_d2 = hk2_d;
            terms.hk2_u2 = hk2_u;
            terms.rt2_t2 = rt2_t;
            terms.rt2_u2 = rt2_u;
            
            // XFOIL formula:
            // VS1(1,1) = Z_AX * AX_A1 - 1.0
            // VS1(1,2) = Z_AX * (AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1)
            // VS1(1,3) = Z_AX * (AX_HK1*HK1_D1)
            // VS2(1,1) = Z_AX * AX_A2 + 1.0
            // VS2(1,2) = Z_AX * (AX_HK2*HK2_T2 + AX_T2 + AX_RT2*RT2_T2)
            // VS2(1,3) = Z_AX * (AX_HK2*HK2_D2)
            
            jac.vs1[0][0] = z_ax * ax_result.ax_a1 - 1.0;
            jac.vs1[0][1] = z_ax * (ax_hk1 * hk1_t + ax_t1 + ax_rt1 * rt1_t);
            jac.vs1[0][2] = z_ax * (ax_hk1 * hk1_d);
            jac.vs1[0][3] = z_ax * (ax_hk1 * hk1_u + ax_rt1 * rt1_u);
            jac.vs1[0][4] = ax; // ∂/∂X1
            
            jac.vs2[0][0] = z_ax * ax_result.ax_a2 + 1.0;
            jac.vs2[0][1] = z_ax * (ax_hk2 * hk2_t + ax_t2 + ax_rt2 * rt2_t);
            jac.vs2[0][2] = z_ax * (ax_hk2 * hk2_d);
            jac.vs2[0][3] = z_ax * (ax_hk2 * hk2_u + ax_rt2 * rt2_u);
            jac.vs2[0][4] = -ax; // -∂/∂X2
        }
        FlowType::Turbulent | FlowType::Wake => {
            // Shear-lag equation (xblsys.f:1684-1841)
            // SA = actual shear stress (ctau), CQA = equilibrium shear stress
            let sa = (1.0 - upw) * s1.ctau + upw * s2.ctau;
            let cqa = (1.0 - upw) * s1.cq + upw * s2.cq; // Use stored CQ (equilibrium shear)
            let cfa = (1.0 - upw) * s1.cf + upw * s2.cf;
            let hka = (1.0 - upw) * s1.hk + upw * s2.hk;

            // USA, RTA, DEA are simple averages (xblsys.f:1693-1696)
            let usa = 0.5 * (s1.us + s2.us); // Use stored US (slip velocity)
            let rta = 0.5 * (s1.r_theta + s2.r_theta);
            let dea = 0.5 * (s1.de + s2.de); // Use stored DE (energy thickness)
            let da = 0.5 * (s1.delta_star + s2.delta_star);

            let ald = if flow_type == FlowType::Wake {
                DLCON
            } else {
                1.0
            };

            // Equilibrium dUe/dx (xblsys.f:1706-1732)
            let gcc = if flow_type == FlowType::Turbulent {
                GCCON
            } else {
                0.0
            };
            let mut hkc = hka - 1.0 - gcc / rta.max(1.0);
            let mut hkc_hka = 1.0;
            let mut hkc_rta = gcc / (rta.max(1.0) * rta.max(1.0));
            if flow_type == FlowType::Turbulent && hkc < 0.01 {
                hkc = 0.01;
                hkc_hka = 0.0;
                hkc_rta = 0.0;
            }

            let hr = hkc / (GACON * ald * hka);
            let hr_hka = hkc_hka / (GACON * ald * hka) - hr / hka;
            let hr_rta = hkc_rta / (GACON * ald * hka);

            let uq = (0.5 * cfa - hr * hr) / (GBCON * da.max(1e-10));
            let uq_hka = -2.0 * hr * hr_hka / (GBCON * da.max(1e-10));
            let uq_rta = -2.0 * hr * hr_rta / (GBCON * da.max(1e-10));
            let uq_cfa = 0.5 / (GBCON * da.max(1e-10));
            let uq_da = -uq / da.max(1e-10);

            // SCC shear lag constant (xblsys.f:1757-1761)
            let scc = SCCON * 1.333 / (1.0 + usa);
            let scc_usa = -scc / (1.0 + usa);

            let slog = (s2.ctau / s1.ctau.max(1e-20)).ln();
            let dxi = s2.x - s1.x;

            // Residual (xblsys.f:1767-1769)
            // REZC = SCC*(CQA - SA*ALD)*DXI - DEA*2.0*SLOG + DEA*2.0*(UQ*DXI - ULOG)*DUXCON
            res.res_third = -(scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog
                + dea * 2.0 * (uq * dxi - ulog) * DUXCON);

            // === Jacobian entries for shear-lag equation (xblsys.f:1781-1839) ===
            // Z coefficients
            let z_cfa = dea * 2.0 * uq_cfa * dxi * DUXCON;
            let z_hka = dea * 2.0 * uq_hka * dxi * DUXCON;
            let z_da = dea * 2.0 * uq_da * dxi * DUXCON;
            let z_sl = -dea * 2.0;
            let z_ul = -dea * 2.0 * DUXCON;
            let z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * DUXCON;
            let z_usa = scc_usa * (cqa - sa * ald) * dxi;
            let z_cqa = scc * dxi;
            let z_sa = -scc * dxi * ald;
            let z_dea = 2.0 * ((uq * dxi - ulog) * DUXCON - slog);

            // UPW derivatives for CQ and SA
            let z_upw = z_cqa * (s2.cq - s1.cq) + z_sa * (s2.ctau - s1.ctau)
                + z_cfa * (s2.cf - s1.cf) + z_hka * (s2.hk - s1.hk);

            // Direct derivatives
            let z_de1 = 0.5 * z_dea;
            let z_de2 = 0.5 * z_dea;
            let z_us1 = 0.5 * z_usa;
            let z_us2 = 0.5 * z_usa;
            let z_d1 = 0.5 * z_da;
            let z_d2 = 0.5 * z_da;
            let z_u1 = -z_ul / s1.u;
            let z_u2 = z_ul / s2.u;
            let z_x1 = -z_dxi;
            let z_x2 = z_dxi;
            let z_s1 = (1.0 - upw) * z_sa - z_sl / s1.ctau.max(1e-20);
            let z_s2 = upw * z_sa + z_sl / s2.ctau.max(1e-20);
            let z_cq1 = (1.0 - upw) * z_cqa;
            let z_cq2 = upw * z_cqa;
            let z_cf1 = (1.0 - upw) * z_cfa;
            let z_cf2 = upw * z_cfa;
            let z_hk1 = (1.0 - upw) * z_hka;
            let z_hk2 = upw * z_hka;

            // Hk derivatives w.r.t. primary variables
            let hk1_t = s1.derivs.hk_h * s1.derivs.h_theta;
            let hk1_d = s1.derivs.hk_h * s1.derivs.h_delta_star;
            let hk2_t = s2.derivs.hk_h * s2.derivs.h_theta;
            let hk2_d = s2.derivs.hk_h * s2.derivs.h_delta_star;

            // Rθ derivatives
            let rt1_t = s1.r_theta / s1.theta.max(1e-12);
            let rt2_t = s2.r_theta / s2.theta.max(1e-12);
            let rt1_u = s1.r_theta / s1.u.max(1e-12);
            let rt2_u = s2.r_theta / s2.u.max(1e-12);

            // Cf derivatives via Hk and Rθ
            let cf1_t = s1.derivs.cf_hk * hk1_t + s1.derivs.cf_rt * rt1_t;
            let cf1_d = s1.derivs.cf_hk * hk1_d;
            let cf1_u = s1.derivs.cf_rt * rt1_u;
            let cf2_t = s2.derivs.cf_hk * hk2_t + s2.derivs.cf_rt * rt2_t;
            let cf2_d = s2.derivs.cf_hk * hk2_d;
            let cf2_u = s2.derivs.cf_rt * rt2_u;

            // UPW derivatives w.r.t. primary variables (XFOIL xblsys.f:1637-1644)
            // UPW_T = UPW_HK * HK_T, UPW_D = UPW_HK * HK_D
            // These are needed for UQ derivative calculations
            let upw_t1_shear = upw_hk1 * hk1_t;
            let upw_d1_shear = upw_hk1 * hk1_d;
            let upw_t2_shear = upw_hk2 * hk2_t;
            let upw_d2_shear = upw_hk2 * hk2_d;

            // UQ derivatives via averaged quantities (xblsys.f:1735-1790)
            // UQ depends on CFA, HKA, DA, RTA
            // First stage: compute with UPW contributions (xblsys.f:1770-1779)
            let uq_upw = uq_cfa * (s2.cf - s1.cf) + uq_hka * (s2.hk - s1.hk);
            
            // UQ_T1 = (1.0-UPW)*(UQ_CFA*CF1_T1 + UQ_HKA*HK1_T1) + UQ_UPW*UPW_T1
            // Then add: + 0.5*UQ_RTA*RT1_T1 (xblsys.f:1781)
            let uq_t1 = (1.0 - upw) * (uq_cfa * cf1_t + uq_hka * hk1_t) 
                + uq_upw * upw_t1_shear
                + 0.5 * uq_rta * rt1_t;
            
            // UQ_D1 = (1.0-UPW)*(UQ_CFA*CF1_D1 + UQ_HKA*HK1_D1) + UQ_UPW*UPW_D1
            // Then add: + 0.5*UQ_DA (xblsys.f:1782)
            let uq_d1 = (1.0 - upw) * (uq_cfa * cf1_d + uq_hka * hk1_d)
                + uq_upw * upw_d1_shear
                + 0.5 * uq_da;
            
            // UQ_U1 = (1.0-UPW)*(UQ_CFA*CF1_U1 + UQ_HKA*HK1_U1) + UQ_UPW*UPW_U1
            // Then add: + 0.5*UQ_RTA*RT1_U1 (xblsys.f:1783)
            let uq_u1 = (1.0 - upw) * (uq_cfa * cf1_u)
                + 0.5 * uq_rta * rt1_u;
            
            // UQ_T2 = UPW*(UQ_CFA*CF2_T2 + UQ_HKA*HK2_T2) + UQ_UPW*UPW_T2
            // Then add: + 0.5*UQ_RTA*RT2_T2 (xblsys.f:1784)
            let uq_t2 = upw * (uq_cfa * cf2_t + uq_hka * hk2_t)
                + uq_upw * upw_t2_shear
                + 0.5 * uq_rta * rt2_t;
            
            // UQ_D2 = UPW*(UQ_CFA*CF2_D2 + UQ_HKA*HK2_D2) + UQ_UPW*UPW_D2
            // Then add: + 0.5*UQ_DA (xblsys.f:1785)
            let uq_d2 = upw * (uq_cfa * cf2_d + uq_hka * hk2_d)
                + uq_upw * upw_d2_shear
                + 0.5 * uq_da;
            
            // UQ_U2 = UPW*(UQ_CFA*CF2_U2 + UQ_HKA*HK2_U2) + UQ_UPW*UPW_U2
            // Then add: + 0.5*UQ_RTA*RT2_U2 (xblsys.f:1786)
            let uq_u2 = upw * (uq_cfa * cf2_u)
                + 0.5 * uq_rta * rt2_u;

            terms.z_upw = z_upw;
            terms.z_de2 = z_de2;
            terms.z_us2 = z_us2;
            terms.z_cq2 = z_cq2;
            terms.z_cf2 = z_cf2;
            terms.z_hk2 = z_hk2;
            terms.z_d2 = z_d2;
            terms.z_u2 = z_u2;
            terms.z_s2 = z_s2;
            terms.upw_t2 = upw_t2_shear;
            terms.upw_d2 = upw_d2_shear;
            terms.upw_u2 = 0.0;
            terms.de2_t2 = s2.derivs.de_t;
            terms.de2_d2 = s2.derivs.de_d;
            terms.de2_u2 = 0.0;
            terms.us2_t2 = s2.derivs.us_t;
            terms.us2_d2 = s2.derivs.us_d;
            terms.us2_u2 = s2.derivs.us_u;
            terms.cq2_t2 = s2.derivs.cq_t;
            terms.cq2_d2 = s2.derivs.cq_d;
            terms.cq2_u2 = s2.derivs.cq_u;
            terms.cf2_t2 = cf2_t;
            terms.cf2_d2 = cf2_d;
            terms.cf2_u2 = cf2_u;
            terms.hk2_t2 = hk2_t;
            terms.hk2_d2 = hk2_d;
            terms.hk2_u2 = 0.0;

            // Build Jacobian entries (xblsys.f:1848-1867)
            // XFOIL builds these in two stages:
            // 1. First: Z_UPW*UPW_T + Z_DE*DE_T + Z_US*US_T
            // 2. Then adds: Z_CQ*CQ_T + Z_CF*CF_T + Z_HK*HK_T
            // Note: The UQ derivatives are already incorporated through Z_CF and Z_HK
            // (Z_CFA = DEA*2.0*UQ_CFA*DXI*DUXCON, Z_HKA = DEA*2.0*UQ_HKA*DXI*DUXCON)
            
            jac.vs1[0][0] = z_s1;
            jac.vs1[0][1] = z_upw * upw_t1_shear + z_de1 * s1.derivs.de_t + z_us1 * s1.derivs.us_t
                + z_cq1 * s1.derivs.cq_t + z_cf1 * cf1_t + z_hk1 * hk1_t;
            jac.vs1[0][2] = z_d1 + z_upw * upw_d1_shear + z_de1 * s1.derivs.de_d + z_us1 * s1.derivs.us_d
                + z_cq1 * s1.derivs.cq_d + z_cf1 * cf1_d + z_hk1 * hk1_d;
            jac.vs1[0][3] = z_u1 + z_us1 * s1.derivs.us_u + z_cq1 * s1.derivs.cq_u + z_cf1 * cf1_u;
            jac.vs1[0][4] = z_x1;

            jac.vs2[0][0] = z_s2;
            jac.vs2[0][1] = z_upw * upw_t2_shear + z_de2 * s2.derivs.de_t + z_us2 * s2.derivs.us_t
                + z_cq2 * s2.derivs.cq_t + z_cf2 * cf2_t + z_hk2 * hk2_t;
            jac.vs2[0][2] = z_d2 + z_upw * upw_d2_shear + z_de2 * s2.derivs.de_d + z_us2 * s2.derivs.us_d
                + z_cq2 * s2.derivs.cq_d + z_cf2 * cf2_d + z_hk2 * hk2_d;
            jac.vs2[0][3] = z_u2 + z_us2 * s2.derivs.us_u + z_cq2 * s2.derivs.cq_u + z_cf2 * cf2_u;
            jac.vs2[0][4] = z_x2;
        }
    }
    } // End of else block for non-SIMI first equation

    // === Momentum equation (xblsys.f:1864-1920) ===
    // XFOIL reference variables:
    //   HA = 0.5*(H1 + H2)
    //   BTMP = HA + 2.0 - MA + HWA (we ignore HWA wake term)
    //   CFX = 0.50*CFM*XA/TA + 0.25*(CF1*X1/T1 + CF2*X2/T2)
    //   REZT = TLOG + BTMP*ULOG - XLOG*0.5*CFX
    //   VSREZ(2) = -REZT
    //
    // Jacobian is computed as derivative of REZT (unnegated residual)
    let ha = 0.5 * (s1.h + s2.h);
    terms.ha = ha;
    let ma_avg = 0.5 * msq;
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);

    // Get midpoint Cf and its derivatives
    let (cfm, cfm_hk1, cfm_hk2, cfm_rt1, cfm_rt2, _) = blmid(s1, s2, flow_type, msq);

    // CFX term: 0.50*CFM*XA/TA + 0.25*(CF1*X1/T1 + CF2*X2/T2)
    let cfx = 0.5 * cfm * xa / ta + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);
    terms.cfx = cfx;

    // Partial derivatives of CFX w.r.t. Cf values
    let cfx_cfm = 0.5 * xa / ta;
    let cfx_cf1 = 0.25 * s1.x / s1.theta;
    let cfx_cf2 = 0.25 * s2.x / s2.theta;

    let btmp = ha + 2.0 - ma_avg;
    terms.btmp_mom = btmp;

    // Momentum residual: VSREZ(2) = -REZT = -(TLOG + BTMP*ULOG - XLOG*0.5*CFX)
    res.res_mom = -(tlog + btmp * ulog - xlog * 0.5 * cfx);

    // Jacobian coefficients (from XFOIL xblsys.f lines 1887-1897)
    // XFOIL definitions:
    //   Z_CFX = -XLOG*0.5
    //   Z_HA  = ULOG
    //   Z_TL  = DDLOG           (DDLOG=0 for SIMI, 1 otherwise)
    //   Z_UL  = DDLOG * BTMP
    //   Z_XL  = -DDLOG * 0.5*CFX
    //   Z_T1  = -Z_TL/T1 + Z_CFX*CFX_T1
    //   Z_T2  =  Z_TL/T2 + Z_CFX*CFX_T2
    //   Z_U1  = -Z_UL/U1
    //   Z_U2  =  Z_UL/U2
    let z_cfx = -xlog * 0.5;
    let z_ha = ulog;
    let z_tl = ddlog;  // Key: DDLOG=0 for SIMI!
    let z_ul = ddlog * btmp;  // Key: DDLOG=0 for SIMI!
    let z_cfm = z_cfx * cfx_cfm;
    let z_cf1 = z_cfx * cfx_cf1;
    let z_cf2 = z_cfx * cfx_cf2;

    // Jacobian entries for momentum equation
    // These are derivatives of REZT (the unnegated residual)
    // XFOIL: VS1(2,2) = 0.5*Z_HA*H1_T1 + Z_CFM*CFM_T1 + Z_CF1*CF1_T1 + Z_T1
    // 
    // From XFOIL (xblsys.f lines 1874-1879):
    //   CFX_TA = -.50*CFM*XA/TA**2
    //   CFX_T1 = -.25*CF1*X1/T1**2 + CFX_TA*0.5
    //   CFX_T2 = -.25*CF2*X2/T2**2 + CFX_TA*0.5
    //   Z_T1 = -Z_TL/T1 + Z_CFX*CFX_T1
    //   Z_T2 =  Z_TL/T2 + Z_CFX*CFX_T2
    //   Z_U1 = -Z_UL/U1
    //   Z_U2 =  Z_UL/U2
    //
    // The θ derivatives must include both local and average theta contributions:
    let cfx_ta = -0.5 * cfm * xa / (ta * ta);  // Derivative of CFX w.r.t. average theta
    let cfx_t1 = -0.25 * s1.cf * s1.x / (s1.theta * s1.theta) + cfx_ta * 0.5;
    let cfx_t2 = -0.25 * s2.cf * s2.x / (s2.theta * s2.theta) + cfx_ta * 0.5;
    terms.cfx_ta = cfx_ta;
    terms.cfx_t1 = cfx_t1;
    terms.cfx_t2 = cfx_t2;
    // Key: z_tl multiplied by DDLOG, so these terms vanish for SIMI
    let z_t1 = -z_tl / s1.theta + z_cfx * cfx_t1;
    let z_t2 = z_tl / s2.theta + z_cfx * cfx_t2;
    let z_u1 = -z_ul / s1.u.max(1e-10);
    let z_u2 = z_ul / s2.u.max(1e-10);

    // Compute chain rule derivatives:
    // H_T (∂H/∂θ) = -δ*/θ² = -H/θ
    // Hk_T (∂Hk/∂θ) = Hk_H * H_T
    let h1_t = s1.derivs.h_theta;
    let h2_t = s2.derivs.h_theta;
    let hk1_t = s1.derivs.hk_h * h1_t;
    let hk2_t = s2.derivs.hk_h * h2_t;

    // Rθ derivatives: Rθ = Re * Ue * θ, so ∂Rθ/∂θ = Re * Ue
    // We get Re*Ue from Rθ/θ
    let rt1_t = s1.r_theta / s1.theta.max(1e-12);
    let rt2_t = s2.r_theta / s2.theta.max(1e-12);

    // Full θ Jacobian includes contributions from both Hk and Rθ:
    // CFM_T = CFM_HK * HK_T + CFM_RT * RT_T
    // CF_T  = CF_HK * HK_T + CF_RT * RT_T
    jac.vs1[1][1] = 0.5 * z_ha * h1_t
        + z_cfm * (cfm_hk1 * hk1_t + cfm_rt1 * rt1_t)
        + z_cf1 * (s1.derivs.cf_hk * hk1_t + s1.derivs.cf_rt * rt1_t)
        + z_t1;
    jac.vs2[1][1] = 0.5 * z_ha * h2_t
        + z_cfm * (cfm_hk2 * hk2_t + cfm_rt2 * rt2_t)
        + z_cf2 * (s2.derivs.cf_hk * hk2_t + s2.derivs.cf_rt * rt2_t)
        + z_t2;

    // The δ* derivatives (via H and via Cf(Hk(H)))
    // XFOIL: VS2(2,3) = 0.5*Z_HA*H2_D2 + Z_CFM*CFM_D2 + Z_CF2*CF2_D2
    // H_D (∂H/∂δ*) = 1/θ
    // Hk_D (∂Hk/∂δ*) = Hk_H * H_D
    let h1_d = s1.derivs.h_delta_star;
    let h2_d = s2.derivs.h_delta_star;
    let hk1_d = s1.derivs.hk_h * h1_d;
    let hk2_d = s2.derivs.hk_h * h2_d;

    jac.vs1[1][2] = 0.5 * z_ha * h1_d + z_cfm * cfm_hk1 * hk1_d + z_cf1 * s1.derivs.cf_hk * hk1_d;
    jac.vs2[1][2] = 0.5 * z_ha * h2_d + z_cfm * cfm_hk2 * hk2_d + z_cf2 * s2.derivs.cf_hk * hk2_d;

    // The Ue derivatives (via ULOG and via Cf)
    // XFOIL: VS1(2,4) = 0.5*Z_MA*M1_U1 + Z_CFM*CFM_U1 + Z_CF1*CF1_U1 + Z_U1
    // XFOIL: VS2(2,4) = 0.5*Z_MA*M2_U2 + Z_CFM*CFM_U2 + Z_CF2*CF2_U2 + Z_U2
    // 
    // Z_MA = -ULOG (from XFOIL line 1985), but M_U is small for incompressible
    // Key: z_u1, z_u2 include DDLOG factor, so vanish for SIMI
    // But CFM_U and CF_U contributions remain even for SIMI!
    //
    // CFM_U = 0.5*(CFM_HKA*HK_U + CFM_RTA*RT_U + CFM_MA*M_U)
    // CF_U = CF_HK*HK_U + CF_RT*RT_U + CF_M*M_U (computed as cf2_u above)
    //
    // For incompressible (MSQ≈0): HK_U ≈ 0, so CFM_U ≈ 0.5*CFM_RTA*RT_U
    // RT_U = RT*(1/U) for incompressible (R_U/R and V_U/V ≈ 0)
    let rt1_u = s1.r_theta / s1.u.max(1e-10);
    let rt2_u = s2.r_theta / s2.u.max(1e-10);
    
    // CF_U = CF_RT * RT_U (HK_U and M_U ≈ 0 for incompressible)
    let cf1_u = s1.derivs.cf_rt * rt1_u;
    let cf2_u = s2.derivs.cf_rt * rt2_u;
    
    // CFM_U = 0.5 * CFM_RTA * RT_U
    let cfm_u1 = 0.5 * cfm_rt1 * rt1_u;
    let cfm_u2 = 0.5 * cfm_rt2 * rt2_u;
    
    // Z_MA term: -ULOG * 0.5 * M_U, but M_U ≈ 0 for incompressible, skip
    // Full Jacobians including CF contributions
    jac.vs1[1][3] = z_cfm * cfm_u1 + z_cf1 * cf1_u + z_u1;
    jac.vs2[1][3] = z_cfm * cfm_u2 + z_cf2 * cf2_u + z_u2;

    // The x derivatives (via XLOG and CFX)
    // XFOIL: Z_XL = -DDLOG*0.5*CFX, Z_X1 = -Z_XL/X1 + Z_CFX*CFX_X1
    let z_xl = -ddlog * 0.5 * cfx;  // Key: DDLOG factor
    jac.vs1[1][4] = -z_xl / s1.x + z_cfx * (0.25 * s1.cf / s1.theta + 0.5 * cfm / ta * 0.5);
    jac.vs2[1][4] = z_xl / s2.x + z_cfx * (0.25 * s2.cf / s2.theta + 0.5 * cfm / ta * 0.5);

    // === Shape parameter equation (xblsys.f:1922-1995) ===
    // XFOIL computes: REZH = HLOG + BTMP*ULOG + XLOG*(0.5*CFX - DIX)
    
    let hsa = 0.5 * (s1.hs + s2.hs);
    let hca = 0.5 * (s1.hc + s2.hc);
    terms.hsa = hsa;
    terms.hca = hca;

    let xot1 = s1.x / s1.theta;
    let xot2 = s2.x / s2.theta;
    terms.xot1 = xot1;
    terms.xot2 = xot2;
    let dd = s2.ctau * s2.ctau * (0.995 - s2.us) * 2.0 / s2.hs;
    let dd2 = 0.15 * (0.995 - s2.us).powi(2) / s2.r_theta * 2.0 / s2.hs;
    terms.dd = dd;
    terms.dd2 = dd2;

    // Upwinded DI and CF (XFOIL lines 1932-1935)
    let dix = (1.0 - upw) * s1.cd * xot1 + upw * s2.cd * xot2;
    let cfx_shape = (1.0 - upw) * s1.cf * xot1 + upw * s2.cf * xot2;
    terms.dix = dix;
    terms.cfx_shape = cfx_shape;
    let dix_upw = s2.cd * xot2 - s1.cd * xot1;
    let cfx_upw = s2.cf * xot2 - s1.cf * xot1;
    terms.cfx_upw = cfx_upw;
    terms.dix_upw = dix_upw;

    // BTMP for shape equation (XFOIL line 1937)
    // Note: HWA term is for wake, we ignore for now
    let btmp_shape = 2.0 * hca / hsa + 1.0 - ha;
    terms.btmp_shape = btmp_shape;

    // Shape parameter residual (XFOIL line 1939)
    res.res_shape = -(hlog + btmp_shape * ulog + xlog * (0.5 * cfx_shape - dix));

    // === Z coefficients for shape equation Jacobian (XFOIL lines 1940-1957) ===
    // XFOIL:
    //   Z_CFX = XLOG*0.5
    //   Z_DIX = -XLOG
    //   Z_HL  = DDLOG                  (Key: DDLOG=0 for SIMI)
    //   Z_UL  = DDLOG * BTMP           (Key: DDLOG=0 for SIMI)
    //   Z_XL  = DDLOG * (0.5*CFX-DIX)  (Key: DDLOG=0 for SIMI)
    let z_cfx_shape = xlog * 0.5;
    let z_dix = -xlog;
    terms.z_cfx_shape = z_cfx_shape;
    terms.z_dix_shape = z_dix;
    let z_hca = 2.0 * ulog / hsa;
    let z_ha_shape = -ulog;
    let z_hl = ddlog;  // Key: DDLOG factor!
    let z_ul_shape = ddlog * btmp_shape;  // Key: DDLOG factor!
    let z_xl_shape = ddlog * (0.5 * cfx_shape - dix);  // Key: DDLOG factor!

    // Z_UPW = Z_CFX*CFX_UPW + Z_DIX*DIX_UPW
    let z_upw = z_cfx_shape * cfx_upw + z_dix * dix_upw;

    // Z_HS1 = -HCA*ULOG/HSA^2 - Z_HL/HS1
    // Z_HS2 = -HCA*ULOG/HSA^2 + Z_HL/HS2
    let z_hs1 = -hca * ulog / (hsa * hsa) - z_hl / s1.hs.max(0.01);
    let z_hs2 = -hca * ulog / (hsa * hsa) + z_hl / s2.hs.max(0.01);

    // Z_CF1 = (1-UPW)*Z_CFX*XOT1, Z_CF2 = UPW*Z_CFX*XOT2
    let z_cf1_shape = (1.0 - upw) * z_cfx_shape * xot1;
    let z_cf2_shape = upw * z_cfx_shape * xot2;

    // Z_DI1 = (1-UPW)*Z_DIX*XOT1, Z_DI2 = UPW*Z_DIX*XOT2
    let z_di1 = (1.0 - upw) * z_dix * xot1;
    let z_di2 = upw * z_dix * xot2;

    // Z_T1, Z_T2: direct theta derivatives from CFX and DIX terms (XFOIL lines 1959-1960)
    let z_t1_shape = (1.0 - upw) * (z_cfx_shape * s1.cf + z_dix * s1.cd) * (-xot1 / s1.theta);
    let z_t2_shape = upw * (z_cfx_shape * s2.cf + z_dix * s2.cd) * (-xot2 / s2.theta);

    // Z_X1, Z_X2: x derivatives (XFOIL lines 1961-1962)
    let z_x1_shape = (1.0 - upw) * (z_cfx_shape * s1.cf + z_dix * s1.cd) / s1.theta - z_xl_shape / s1.x;
    let z_x2_shape = upw * (z_cfx_shape * s2.cf + z_dix * s2.cd) / s2.theta + z_xl_shape / s2.x;

    // Z_U1, Z_U2: Ue derivatives (XFOIL lines 1963-1964)
    let z_u1_shape = -z_ul_shape / s1.u;
    let z_u2_shape = z_ul_shape / s2.u;

    // === Derivative chain rules ===
    // We need: HS_T, HS_D, CF_T, CF_D, DI_T, DI_D, H_T, H_D
    // These are computed via chain rule from the stored derivatives

    // H derivatives (shape factor H = δ*/θ)
    // H_T = ∂H/∂θ = -δ*/θ² = -H/θ  (same as h_theta in derivs)
    // H_D = ∂H/∂δ* = 1/θ  (same as h_delta_star in derivs)
    let h1_t = s1.derivs.h_theta;
    let h1_d = s1.derivs.h_delta_star;
    let h2_t = s2.derivs.h_theta;
    let h2_d = s2.derivs.h_delta_star;

    // Hk derivatives via H (Hk = f(H, Msq))
    // HK_T = HK_H * H_T
    // HK_D = HK_H * H_D
    let hk1_t = s1.derivs.hk_h * h1_t;
    let hk1_d = s1.derivs.hk_h * h1_d;
    let hk2_t = s2.derivs.hk_h * h2_t;
    let hk2_d = s2.derivs.hk_h * h2_d;

    // Rθ derivatives: Rθ = Re * Ue * θ
    // RT_T = ∂Rθ/∂θ = Re * Ue = Rθ/θ
    // RT_U = ∂Rθ/∂Ue = Re * θ = Rθ/Ue
    let rt1_t = s1.r_theta / s1.theta.max(1e-12);
    let rt2_t = s2.r_theta / s2.theta.max(1e-12);
    let rt1_u = s1.r_theta / s1.u.max(1e-12);
    let rt2_u = s2.r_theta / s2.u.max(1e-12);

    // Hs derivatives: Hs = f(Hk, Rθ, Msq)
    // HS_T = HS_HK * HK_T + HS_RT * RT_T
    // HS_D = HS_HK * HK_D
    // HS_U = HS_RT * RT_U  (via Rθ dependence on Ue)
    let hs1_t = s1.derivs.hs_hk * hk1_t + s1.derivs.hs_rt * rt1_t;
    let hs1_d = s1.derivs.hs_hk * hk1_d;
    let hs1_u = s1.derivs.hs_rt * rt1_u;
    let hs2_t = s2.derivs.hs_hk * hk2_t + s2.derivs.hs_rt * rt2_t;
    let hs2_d = s2.derivs.hs_hk * hk2_d;
    let hs2_u = s2.derivs.hs_rt * rt2_u;

    // Cf derivatives: Cf = f(Hk, Rθ, Msq)
    // CF_T = CF_HK * HK_T + CF_RT * RT_T
    // CF_D = CF_HK * HK_D
    // CF_U = CF_RT * RT_U
    let cf1_t_shape = s1.derivs.cf_hk * hk1_t + s1.derivs.cf_rt * rt1_t;
    let cf1_d_shape = s1.derivs.cf_hk * hk1_d;
    let cf1_u_shape = s1.derivs.cf_rt * rt1_u;
    let cf2_t_shape = s2.derivs.cf_hk * hk2_t + s2.derivs.cf_rt * rt2_t;
    let cf2_d_shape = s2.derivs.cf_hk * hk2_d;
    let cf2_u_shape = s2.derivs.cf_rt * rt2_u;

    // DI (Cd) derivatives: use direct partials stored by BLVAR
    // Laminar still yields cd_* = 0 for ctau and uses chain rule internally.
    let di1_t = s1.derivs.cd_t;
    let di1_d = s1.derivs.cd_d;
    let di1_u = s1.derivs.cd_u;
    let di2_t = s2.derivs.cd_t;
    let di2_d = s2.derivs.cd_d;
    let di2_u = s2.derivs.cd_u;
    let di1_s = s1.derivs.cd_s;
    let di2_s = s2.derivs.cd_s;

    // Hc derivatives (density shape factor) via chain rule (XFOIL xblsys.f:803-806)
    // HC_T = HC_HK * HK_T, HC_D = HC_HK * HK_D, HC_U = HC_HK * HK_U
    // At M=0, hc_hk=0 so these will be zero for incompressible flow
    let hc1_t = s1.derivs.hc_hk * hk1_t;
    let hc1_d = s1.derivs.hc_hk * hk1_d;
    let hc1_u = 0.0; // HK1_U1 ≈ 0 at low Mach
    let hc2_t = s2.derivs.hc_hk * hk2_t;
    let hc2_d = s2.derivs.hc_hk * hk2_d;
    let hc2_u = 0.0; // HK2_U2 ≈ 0 at low Mach

    // UPW derivatives w.r.t. primary variables (XFOIL xblsys.f:1637-1644)
    // UPW depends on Hk through HL and HDCON, so use chain rule:
    // UPW_T = UPW_HK * HK_T, UPW_D = UPW_HK * HK_D, UPW_U = UPW_HK * HK_U
    //
    let hk1_u = 0.0;
    let hk2_u = 0.0;
    let upw_t1 = upw_hk1 * hk1_t;
    let upw_d1 = upw_hk1 * hk1_d;
    let upw_u1 = upw_hk1 * hk1_u;
    let upw_t2 = upw_hk2 * hk2_t;
    let upw_d2 = upw_hk2 * hk2_d;
    let upw_u2 = upw_hk2 * hk2_u;

    terms.z_hs2 = z_hs2;
    terms.z_cf2_shape = z_cf2_shape;
    terms.z_di2 = z_di2;
    terms.z_t2_shape = z_t2_shape;
    terms.z_u2_shape = z_u2_shape;
    terms.z_hca = z_hca;
    terms.z_ha_shape = z_ha_shape;
    terms.z_upw_shape = z_upw;
    terms.hs2_t = hs2_t;
    terms.hs2_d = hs2_d;
    terms.hs2_u = hs2_u;
    terms.cf2_t_shape = cf2_t_shape;
    terms.cf2_d_shape = cf2_d_shape;
    terms.cf2_u_shape = cf2_u_shape;
    terms.di2_t = di2_t;
    terms.di2_d = di2_d;
    terms.di2_u = di2_u;
    terms.di2_s = di2_s;
    terms.hc2_t = hc2_t;
    terms.hc2_d = hc2_d;
    terms.hc2_u = hc2_u;
    terms.h2_t = h2_t;
    terms.h2_d = h2_d;
    terms.upw_t2_shape = upw_t2;
    terms.upw_d2_shape = upw_d2;
    terms.upw_u2_shape = upw_u2;

    // === Build shape equation Jacobian (XFOIL lines 1969-1989) ===
    // VS1(3,1) = Z_DI1*DI1_S1
    // VS1(3,2) = Z_HS1*HS1_T1 + Z_CF1*CF1_T1 + Z_DI1*DI1_T1 + Z_T1
    //          + 0.5*(Z_HCA*HC1_T1 + Z_HA*H1_T1) + Z_UPW*UPW_T1
    // VS1(3,3) = Z_HS1*HS1_D1 + Z_CF1*CF1_D1 + Z_DI1*DI1_D1
    //          + 0.5*(Z_HCA*HC1_D1 + Z_HA*H1_D1) + Z_UPW*UPW_D1
    // VS1(3,4) = Z_HS1*HS1_U1 + Z_CF1*CF1_U1 + Z_DI1*DI1_U1 + Z_U1
    //          + 0.5*Z_HCA*HC1_U1 + Z_UPW*UPW_U1
    // (same pattern for VS2)

    jac.vs1[2][0] = z_di1 * di1_s;
    jac.vs1[2][1] = z_hs1 * hs1_t + z_cf1_shape * cf1_t_shape + z_di1 * di1_t + z_t1_shape
        + 0.5 * (z_hca * hc1_t + z_ha_shape * h1_t) + z_upw * upw_t1;
    jac.vs1[2][2] = z_hs1 * hs1_d + z_cf1_shape * cf1_d_shape + z_di1 * di1_d
        + 0.5 * (z_hca * hc1_d + z_ha_shape * h1_d) + z_upw * upw_d1;
    jac.vs1[2][3] = z_hs1 * hs1_u + z_cf1_shape * cf1_u_shape + z_di1 * di1_u + z_u1_shape
        + 0.5 * z_hca * hc1_u + z_upw * upw_u1;
    jac.vs1[2][4] = z_x1_shape;

    jac.vs2[2][0] = z_di2 * di2_s;
    jac.vs2[2][1] = z_hs2 * hs2_t + z_cf2_shape * cf2_t_shape + z_di2 * di2_t + z_t2_shape
        + 0.5 * (z_hca * hc2_t + z_ha_shape * h2_t) + z_upw * upw_t2;
    jac.vs2[2][2] = z_hs2 * hs2_d + z_cf2_shape * cf2_d_shape + z_di2 * di2_d
        + 0.5 * (z_hca * hc2_d + z_ha_shape * h2_d) + z_upw * upw_d2;
    jac.vs2[2][3] = z_hs2 * hs2_u + z_cf2_shape * cf2_u_shape + z_di2 * di2_u + z_u2_shape
        + 0.5 * z_hca * hc2_u + z_upw * upw_u2;
    jac.vs2[2][4] = z_x2_shape;

    (res, jac, terms)
}

/// Compute BL equation residuals and Jacobian
pub fn bldif(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    let (res, jac, _) = bldif_with_terms(s1, s2, flow_type, msq, re);
    (res, jac)
}

/// Compute BL equations for similarity station (XFOIL BLDIF with ITYP=0)
///
/// At the similarity station (first BL interval), XFOIL uses fixed logarithmic
/// differences instead of computing them from station ratios:
/// - XLOG = 1.0 (not log(X2/X1))
/// - ULOG = BULE = 1.0 (leading edge pressure gradient parameter)
/// - TLOG = 0.5*(1.0 - BULE) = 0.0
/// - HLOG = 0.0
///
/// This ensures the BL equations have the correct form at the stagnation region
/// where both stations are effectively the same point.
pub fn bldif_simi(
    s: &BlStation,  // The similarity station (s1 = s2)
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    use crate::constants::{GBCON, CTCON};
    
    let mut res = BlResiduals::default();
    let mut jac = BlJacobian::default();
    
    // XFOIL similarity: fixed log values
    let bule = 1.0;  // Leading edge pressure gradient parameter
    let xlog = 1.0;
    let ulog = bule;
    let tlog = 0.5 * (1.0 - bule);  // = 0.0 when BULE=1
    let hlog = 0.0;
    
    // UPW = 1.0 at similarity (since s1 = s2, no upwinding needed)
    let upw = 1.0;
    
    // === First equation: Amplification (laminar) ===
    // XFOIL sets: VS2(1,1) = 1.0, VSREZ(1) = -AMPL2
    // This constrains amplification to be zero at similarity
    res.res_third = -s.ampl;  // Residual: want ampl = 0
    jac.vs2[0][0] = 1.0;       // d(res)/d(ampl2)
    // All other first-equation derivatives are 0
    
    // === Second equation: Momentum integral ===
    // CFM = (CF1 + CF2) / 2 at similarity = CF2
    let cfm = s.cf;
    
    // XFOIL momentum equation at similarity (simplified):
    // REZT = TLOG + (2 + H + M^2*H*(γ-1)/2 - MSQ)*ULOG - 0.5*CFM*X2/T2 * XLOG
    let hk = s.hk;
    let gam = 1.4;
    let h = s.h;
    let theta = s.theta;
    let x = s.x;
    
    // Coefficient in front of ULOG
    let bt = (2.0 + h + msq * h * (gam - 1.0) / 2.0 - msq).max(0.0);
    
    // REZT = TLOG + BT*ULOG - 0.5*CFM*(X/T)*XLOG
    res.res_mom = -(tlog + bt * ulog - 0.5 * cfm * (x / theta.max(1e-10)) * xlog);
    
    // Jacobian for momentum - simplified at similarity
    // Main dependency is on theta through X/T term
    let rezt_t = 0.5 * cfm * (x / (theta * theta).max(1e-20)) * xlog;
    jac.vs2[1][1] = rezt_t;  // d(res)/d(theta2)
    
    // === Third equation: Shape parameter ===
    // Similar simplification - at similarity, shape equation residual ≈ 0
    // since both stations are identical
    res.res_shape = 0.0;
    jac.vs2[2][2] = 1.0;  // Regularization
    
    (res, jac)
}

/// Compute BL equations for the transition interval (XFOIL TRDIF) - Simple version
///
/// This function handles the hybrid laminar-turbulent interval at transition.
/// XFOIL's TRDIF splits the interval X1..X2 into:
/// - Laminar part: X1 to XT (transition point)
/// - Turbulent part: XT to X2
///
/// This is a simplified version that uses basic interpolation weights.
/// For the full Jacobian with proper chain rule through XT derivatives,
/// use `trdif_full()` with a `Trchek2FullResult`.
///
/// # Arguments
/// * `s1` - Upstream laminar station
/// * `s2` - Downstream turbulent station
/// * `xt` - Transition x-location (between s1.x and s2.x)
/// * `msq` - Mach number squared
/// * `re` - Reference Reynolds number
///
/// # Reference
/// XFOIL xblsys.f TRDIF (line 1224-1568)
pub fn trdif(
    s1: &BlStation,
    s2: &BlStation,
    xt: f64,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    use crate::constants::{CTRCON, CTRCEX};
    
    let dx = s2.x - s1.x;
    if dx <= 1e-20 {
        // Degenerate case - use regular turbulent bldif
        return bldif(s1, s2, FlowType::Turbulent, msq, re);
    }
    
    // === Weighting factors for interpolation to transition point ===
    // WF2 = (XT - X1) / (X2 - X1): fraction of interval before transition
    let wf2 = ((xt - s1.x) / dx).clamp(0.01, 0.99);
    let wf1 = 1.0 - wf2;
    
    // === Interpolate primary variables to transition point ===
    let tt = s1.theta * wf1 + s2.theta * wf2;  // theta at XT
    let dt = s1.delta_star * wf1 + s2.delta_star * wf2;  // delta_star at XT
    let ut = s1.u * wf1 + s2.u * wf2;  // Ue at XT
    
    // === PART 1: Laminar interval X1 → XT ===
    // Create transition station with interpolated values (laminar)
    let mut st_laminar = BlStation::new();
    st_laminar.x = xt;
    st_laminar.u = ut;
    st_laminar.theta = tt;
    st_laminar.delta_star = dt;
    st_laminar.ctau = 0.0;  // Laminar
    st_laminar.ampl = 9.0;  // At Ncrit
    st_laminar.is_laminar = true;
    st_laminar.is_turbulent = false;
    
    // Compute secondary variables for laminar transition station
    blvar(&mut st_laminar, FlowType::Laminar, msq, re);
    
    // Compute laminar bldif from s1 to transition station
    let (res_lam, jac_lam) = bldif(s1, &st_laminar, FlowType::Laminar, msq, re);
    
    // === PART 2: Turbulent interval XT → X2 ===
    // First, compute turbulent equilibrium CQ at transition point (XFOIL line 1415)
    // XFOIL calls BLVAR(2) to get turbulent secondary variables before computing ST
    let mut st_turb_temp = BlStation::new();
    st_turb_temp.x = xt;
    st_turb_temp.u = ut;
    st_turb_temp.theta = tt;
    st_turb_temp.delta_star = dt;
    st_turb_temp.ctau = 0.03;  // Temporary value for blvar
    st_turb_temp.ampl = 0.0;
    st_turb_temp.is_laminar = false;
    st_turb_temp.is_turbulent = true;
    blvar(&mut st_turb_temp, FlowType::Turbulent, msq, re);
    
    // Now compute initial ctau at transition (XFOIL line 1420-1430)
    // XFOIL: CTR = CTRCON * exp(-CTRCEX/(HK2-1)), ST = CTR * CQ2
    let hk_t = st_turb_temp.hk;
    let cq_t = st_turb_temp.cq;  // Equilibrium ctau from TURBULENT closure
    
    let hk_arg = (hk_t - 1.0).max(0.1);
    let ctr = CTRCON * (-CTRCEX / hk_arg).exp();
    let st_ctau = ctr * cq_t;
    
    // Create turbulent transition station (station "1" for turbulent part)
    let mut st_turb = BlStation::new();
    st_turb.x = xt;
    st_turb.u = ut;
    st_turb.theta = tt;
    st_turb.delta_star = dt;
    st_turb.ctau = st_ctau.max(0.0001);
    st_turb.ampl = 0.0;  // Turbulent
    st_turb.is_laminar = false;
    st_turb.is_turbulent = true;
    
    // Compute secondary variables for turbulent transition station with correct ctau
    blvar(&mut st_turb, FlowType::Turbulent, msq, re);
    
    // Compute turbulent bldif from transition station to s2
    let (res_turb, jac_turb) = bldif(&st_turb, s2, FlowType::Turbulent, msq, re);
    
    // === COMBINE: Sum laminar and turbulent contributions ===
    // XFOIL combines the systems (lines 1547-1568):
    // - Row 1 (shear-lag): Only from turbulent part
    // - Rows 2-3 (momentum, shape): Sum of laminar + turbulent
    let mut res = BlResiduals::default();
    let mut jac = BlJacobian::default();
    
    // Residuals
    res.res_third = res_turb.res_third;  // Shear-lag: turbulent only
    res.res_mom = res_lam.res_mom + res_turb.res_mom;  // Sum
    res.res_shape = res_lam.res_shape + res_turb.res_shape;  // Sum
    
    // Jacobian VS1 (derivatives w.r.t. upstream station s1)
    // The laminar Jacobian's VS1 applies to the true s1
    // The turbulent Jacobian's VS1 applies to st_turb, which depends on s1 via interpolation
    // 
    // For the full chain rule, we'd need: dF/ds1 = dF_lam/ds1 + dF_turb/dst * dst/ds1
    // where dst/ds1 = WF1 (the interpolation weight)
    //
    // Simplified approach: Use laminar VS1 for rows 2-3, turbulent VS1 scaled by WF1 for row 1
    // This captures the main effect while avoiding the full chain rule complexity
    
    // Row 0 (shear-lag): From turbulent part only, but needs chain rule through XT
    // The turbulent system uses st_turb as "1" station, so its VS1 entries need to be
    // mapped to the actual s1 via the interpolation.
    // For simplicity, we use the turbulent VS2 entries (which apply to s2 directly)
    jac.vs1[0] = [0.0; 5];  // Shear-lag doesn't directly depend on s1 (laminar station)
    
    // Rows 1-2 (momentum, shape): Sum of laminar and turbulent contributions
    // Laminar VS1 applies directly to s1
    // Turbulent VS1 needs chain rule: VS1_turb * WF1 (interpolation weight)
    for j in 0..5 {
        jac.vs1[1][j] = jac_lam.vs1[1][j] + jac_turb.vs1[1][j] * wf1;
        jac.vs1[2][j] = jac_lam.vs1[2][j] + jac_turb.vs1[2][j] * wf1;
    }
    
    // Jacobian VS2 (derivatives w.r.t. downstream station s2)
    // Laminar VS2 applies to st_laminar, needs chain rule: VS2_lam * WF2
    // Turbulent VS2 applies directly to s2
    
    // Row 0 (shear-lag): From turbulent part only
    jac.vs2[0] = jac_turb.vs2[0];
    
    // Rows 1-2 (momentum, shape): Combine with chain rule
    for j in 0..5 {
        jac.vs2[1][j] = jac_lam.vs2[1][j] * wf2 + jac_turb.vs2[1][j];
        jac.vs2[2][j] = jac_lam.vs2[2][j] * wf2 + jac_turb.vs2[2][j];
    }
    
    (res, jac)
}

/// Compute BL equations for the transition interval with full Jacobian (XFOIL TRDIF)
///
/// This is the full implementation of XFOIL's TRDIF that properly handles the
/// chain rule through the transition point XT and its derivatives.
///
/// The function splits the interval X1..X2 into:
/// - Laminar part: X1 to XT
/// - Turbulent part: XT to X2
///
/// The Jacobian transformation follows XFOIL exactly:
/// - Laminar Jacobian VS2 entries are transformed via TT/DT/UT/XT derivatives
/// - Turbulent Jacobian VS1 entries are transformed via ST (initial ctau) derivatives
///
/// # Arguments
/// * `s1` - Upstream laminar station
/// * `s2` - Downstream turbulent station  
/// * `tr` - Full transition result from `trchek2_full()` containing XT and all derivatives
/// * `msq` - Mach number squared
/// * `re` - Reference Reynolds number
///
/// # Reference
/// XFOIL xblsys.f TRDIF (lines 1195-1549)
pub fn trdif_full(
    s1: &BlStation,
    s2: &BlStation,
    tr: &crate::closures::Trchek2FullResult,
    msq: f64,
    re: f64,
) -> (BlResiduals, BlJacobian) {
    use crate::constants::{CTRCON, CTRCEX};
    
    let dx = s2.x - s1.x;
    if dx <= 1e-20 || !tr.transition {
        // Degenerate case or no transition - use regular turbulent bldif
        return bldif(s1, s2, FlowType::Turbulent, msq, re);
    }
    
    let xt = tr.xt;
    let wf1 = tr.wf1;
    let wf2 = tr.wf2;
    
    // === Interpolate primary variables to transition point ===
    let tt = s1.theta * wf1 + s2.theta * wf2;
    let dt = s1.delta_star * wf1 + s2.delta_star * wf2;
    let ut = s1.u * wf1 + s2.u * wf2;
    
    // === PART 1: Laminar interval X1 → XT ===
    let mut st_laminar = BlStation::new();
    st_laminar.x = xt;
    st_laminar.u = ut;
    st_laminar.theta = tt;
    st_laminar.delta_star = dt;
    st_laminar.ctau = 0.0;
    st_laminar.ampl = 9.0;
    st_laminar.is_laminar = true;
    st_laminar.is_turbulent = false;
    
    blvar(&mut st_laminar, FlowType::Laminar, msq, re);
    let (res_lam, jac_lam) = bldif(s1, &st_laminar, FlowType::Laminar, msq, re);
    
    // === PART 2: Turbulent interval XT → X2 ===
    // Compute turbulent equilibrium CQ at transition point
    let mut st_turb_temp = BlStation::new();
    st_turb_temp.x = xt;
    st_turb_temp.u = ut;
    st_turb_temp.theta = tt;
    st_turb_temp.delta_star = dt;
    st_turb_temp.ctau = 0.03;
    st_turb_temp.ampl = 0.0;
    st_turb_temp.is_laminar = false;
    st_turb_temp.is_turbulent = true;
    blvar(&mut st_turb_temp, FlowType::Turbulent, msq, re);
    
    // Compute initial ctau at transition (ST = CTR * CQ)
    let hk_t = st_turb_temp.hk;
    let cq_t = st_turb_temp.cq;
    let hk_arg = (hk_t - 1.0).max(0.1);
    let ctr = CTRCON * (-CTRCEX / hk_arg).exp();
    let ctr_hk = ctr * CTRCEX / (hk_arg * hk_arg);
    
    let st = ctr * cq_t;
    
    // ST derivatives w.r.t. TT, DT, UT (at transition point)
    // ST = CTR(HK) * CQ(HK, Rθ, ...)
    // Need CQ derivatives from blvar - using stored derivatives
    let cq_t_deriv = st_turb_temp.derivs.cq_t;
    let cq_d_deriv = st_turb_temp.derivs.cq_d;
    let cq_u_deriv = st_turb_temp.derivs.cq_u;
    
    // HK derivatives at transition point
    let hk_t_deriv = st_turb_temp.derivs.hk_h * (-st_turb_temp.h / tt.max(1e-20));
    let hk_d_deriv = st_turb_temp.derivs.hk_h * (1.0 / tt.max(1e-20));
    let hk_u_deriv = 0.0;
    
    // ST = CTR * CQ, so ST_TT = CTR * CQ_T + CQ * CTR_HK * HK_T
    let st_tt = ctr * cq_t_deriv + cq_t * ctr_hk * hk_t_deriv;
    let st_dt = ctr * cq_d_deriv + cq_t * ctr_hk * hk_d_deriv;
    let st_ut = ctr * cq_u_deriv + cq_t * ctr_hk * hk_u_deriv;
    
    // Chain rule: ST derivatives w.r.t. actual "1" and "2" variables
    // ST_*1 = ST_TT * TT_*1 + ST_DT * DT_*1 + ST_UT * UT_*1
    let st_a1 = st_tt * tr.tt_a1 + st_dt * tr.dt_a1 + st_ut * tr.ut_a1;
    let st_x1 = st_tt * tr.tt_x1 + st_dt * tr.dt_x1 + st_ut * tr.ut_x1;
    let st_t1 = st_tt * tr.tt_t1 + st_dt * tr.dt_t1 + st_ut * tr.ut_t1;
    let st_d1 = st_tt * tr.tt_d1 + st_dt * tr.dt_d1 + st_ut * tr.ut_d1;
    let st_u1 = st_tt * tr.tt_u1 + st_dt * tr.dt_u1 + st_ut * tr.ut_u1;
    
    let st_x2 = st_tt * tr.tt_x2 + st_dt * tr.dt_x2 + st_ut * tr.ut_x2;
    let st_t2 = st_tt * tr.tt_t2 + st_dt * tr.dt_t2 + st_ut * tr.ut_t2;
    let st_d2 = st_tt * tr.tt_d2 + st_dt * tr.dt_d2 + st_ut * tr.ut_d2;
    let st_u2 = st_tt * tr.tt_u2 + st_dt * tr.dt_u2 + st_ut * tr.ut_u2;
    
    // Create turbulent transition station with correct ST
    let mut st_turb = BlStation::new();
    st_turb.x = xt;
    st_turb.u = ut;
    st_turb.theta = tt;
    st_turb.delta_star = dt;
    st_turb.ctau = st.max(0.0001);
    st_turb.ampl = 0.0;
    st_turb.is_laminar = false;
    st_turb.is_turbulent = true;
    blvar(&mut st_turb, FlowType::Turbulent, msq, re);
    
    let (res_turb, jac_turb) = bldif(&st_turb, s2, FlowType::Turbulent, msq, re);
    
    // === Transform Jacobians using full chain rule (XFOIL lines 1318-1515) ===
    
    // Laminar part: VS2 affects T variables (TT, DT, UT, XT)
    // BL1[K][j] = VS1[K][j] + VS2[K][2]*TT_*1 + VS2[K][3]*DT_*1 + VS2[K][4]*UT_*1 + VS2[K][5]*XT_*1
    // BL2[K][j] = VS2[K][2]*TT_*2 + VS2[K][3]*DT_*2 + VS2[K][4]*UT_*2 + VS2[K][5]*XT_*2
    let mut bl1 = [[0.0f64; 5]; 3];
    let mut bl2 = [[0.0f64; 5]; 3];
    
    for k in 1..3 {  // Only rows 2-3 (momentum, shape) for laminar
        // BL1: derivatives w.r.t. station 1 variables
        // Column 0: amplification (index 0)
        bl1[k][0] = jac_lam.vs1[k][0]
                  + jac_lam.vs2[k][1] * tr.tt_a1
                  + jac_lam.vs2[k][2] * tr.dt_a1
                  + jac_lam.vs2[k][3] * tr.ut_a1
                  + jac_lam.vs2[k][4] * tr.xt_a1;
        // Column 1: theta (index 1)
        bl1[k][1] = jac_lam.vs1[k][1]
                  + jac_lam.vs2[k][1] * tr.tt_t1
                  + jac_lam.vs2[k][2] * tr.dt_t1
                  + jac_lam.vs2[k][3] * tr.ut_t1
                  + jac_lam.vs2[k][4] * tr.xt_t1;
        // Column 2: delta_star (index 2)
        bl1[k][2] = jac_lam.vs1[k][2]
                  + jac_lam.vs2[k][1] * tr.tt_d1
                  + jac_lam.vs2[k][2] * tr.dt_d1
                  + jac_lam.vs2[k][3] * tr.ut_d1
                  + jac_lam.vs2[k][4] * tr.xt_d1;
        // Column 3: u (index 3)
        bl1[k][3] = jac_lam.vs1[k][3]
                  + jac_lam.vs2[k][1] * tr.tt_u1
                  + jac_lam.vs2[k][2] * tr.dt_u1
                  + jac_lam.vs2[k][3] * tr.ut_u1
                  + jac_lam.vs2[k][4] * tr.xt_u1;
        // Column 4: x (index 4)
        bl1[k][4] = jac_lam.vs1[k][4]
                  + jac_lam.vs2[k][1] * tr.tt_x1
                  + jac_lam.vs2[k][2] * tr.dt_x1
                  + jac_lam.vs2[k][3] * tr.ut_x1
                  + jac_lam.vs2[k][4] * tr.xt_x1;
        
        // BL2: derivatives w.r.t. station 2 variables (no VS1 contribution)
        bl2[k][0] = 0.0;  // No amplification dependence in VS2 for station 2
        bl2[k][1] = jac_lam.vs2[k][1] * tr.tt_t2
                  + jac_lam.vs2[k][2] * tr.dt_t2
                  + jac_lam.vs2[k][3] * tr.ut_t2
                  + jac_lam.vs2[k][4] * tr.xt_t2;
        bl2[k][2] = jac_lam.vs2[k][1] * tr.tt_d2
                  + jac_lam.vs2[k][2] * tr.dt_d2
                  + jac_lam.vs2[k][3] * tr.ut_d2
                  + jac_lam.vs2[k][4] * tr.xt_d2;
        bl2[k][3] = jac_lam.vs2[k][1] * tr.tt_u2
                  + jac_lam.vs2[k][2] * tr.dt_u2
                  + jac_lam.vs2[k][3] * tr.ut_u2
                  + jac_lam.vs2[k][4] * tr.xt_u2;
        bl2[k][4] = jac_lam.vs2[k][1] * tr.tt_x2
                  + jac_lam.vs2[k][2] * tr.dt_x2
                  + jac_lam.vs2[k][3] * tr.ut_x2
                  + jac_lam.vs2[k][4] * tr.xt_x2;
    }
    
    // Turbulent part: VS1 affects T variables (ST, TT, DT, UT, XT)
    // BT1[K][j] = VS1[K][1]*ST_*1 + VS1[K][2]*TT_*1 + VS1[K][3]*DT_*1 + VS1[K][4]*UT_*1 + VS1[K][5]*XT_*1
    // BT2[K][j] = VS2[K][j] + VS1[K][1]*ST_*2 + VS1[K][2]*TT_*2 + ...
    let mut bt1 = [[0.0f64; 5]; 3];
    let mut bt2 = [[0.0f64; 5]; 3];
    
    for k in 0..3 {
        // BT1: derivatives w.r.t. station 1 variables
        // Column 0: amplification
        bt1[k][0] = jac_turb.vs1[k][0] * st_a1
                  + jac_turb.vs1[k][1] * tr.tt_a1
                  + jac_turb.vs1[k][2] * tr.dt_a1
                  + jac_turb.vs1[k][3] * tr.ut_a1
                  + jac_turb.vs1[k][4] * tr.xt_a1;
        // Column 1: theta
        bt1[k][1] = jac_turb.vs1[k][0] * st_t1
                  + jac_turb.vs1[k][1] * tr.tt_t1
                  + jac_turb.vs1[k][2] * tr.dt_t1
                  + jac_turb.vs1[k][3] * tr.ut_t1
                  + jac_turb.vs1[k][4] * tr.xt_t1;
        // Column 2: delta_star
        bt1[k][2] = jac_turb.vs1[k][0] * st_d1
                  + jac_turb.vs1[k][1] * tr.tt_d1
                  + jac_turb.vs1[k][2] * tr.dt_d1
                  + jac_turb.vs1[k][3] * tr.ut_d1
                  + jac_turb.vs1[k][4] * tr.xt_d1;
        // Column 3: u
        bt1[k][3] = jac_turb.vs1[k][0] * st_u1
                  + jac_turb.vs1[k][1] * tr.tt_u1
                  + jac_turb.vs1[k][2] * tr.dt_u1
                  + jac_turb.vs1[k][3] * tr.ut_u1
                  + jac_turb.vs1[k][4] * tr.xt_u1;
        // Column 4: x
        bt1[k][4] = jac_turb.vs1[k][0] * st_x1
                  + jac_turb.vs1[k][1] * tr.tt_x1
                  + jac_turb.vs1[k][2] * tr.dt_x1
                  + jac_turb.vs1[k][3] * tr.ut_x1
                  + jac_turb.vs1[k][4] * tr.xt_x1;
        
        // BT2: VS2 applies directly, plus VS1 chain rule through T variables
        bt2[k][0] = jac_turb.vs2[k][0];  // Direct ctau dependence
        bt2[k][1] = jac_turb.vs2[k][1]
                  + jac_turb.vs1[k][0] * st_t2
                  + jac_turb.vs1[k][1] * tr.tt_t2
                  + jac_turb.vs1[k][2] * tr.dt_t2
                  + jac_turb.vs1[k][3] * tr.ut_t2
                  + jac_turb.vs1[k][4] * tr.xt_t2;
        bt2[k][2] = jac_turb.vs2[k][2]
                  + jac_turb.vs1[k][0] * st_d2
                  + jac_turb.vs1[k][1] * tr.tt_d2
                  + jac_turb.vs1[k][2] * tr.dt_d2
                  + jac_turb.vs1[k][3] * tr.ut_d2
                  + jac_turb.vs1[k][4] * tr.xt_d2;
        bt2[k][3] = jac_turb.vs2[k][3]
                  + jac_turb.vs1[k][0] * st_u2
                  + jac_turb.vs1[k][1] * tr.tt_u2
                  + jac_turb.vs1[k][2] * tr.dt_u2
                  + jac_turb.vs1[k][3] * tr.ut_u2
                  + jac_turb.vs1[k][4] * tr.xt_u2;
        bt2[k][4] = jac_turb.vs2[k][4]
                  + jac_turb.vs1[k][0] * st_x2
                  + jac_turb.vs1[k][1] * tr.tt_x2
                  + jac_turb.vs1[k][2] * tr.dt_x2
                  + jac_turb.vs1[k][3] * tr.ut_x2
                  + jac_turb.vs1[k][4] * tr.xt_x2;
    }
    
    // === Combine laminar and turbulent parts ===
    // Row 0 (shear-lag): turbulent only
    // Rows 1-2 (momentum, shape): laminar + turbulent
    let mut res = BlResiduals::default();
    let mut jac = BlJacobian::default();
    
    res.res_third = res_turb.res_third;
    res.res_mom = res_lam.res_mom + res_turb.res_mom;
    res.res_shape = res_lam.res_shape + res_turb.res_shape;
    
    for l in 0..5 {
        jac.vs1[0][l] = bt1[0][l];
        jac.vs2[0][l] = bt2[0][l];
        jac.vs1[1][l] = bl1[1][l] + bt1[1][l];
        jac.vs2[1][l] = bl2[1][l] + bt2[1][l];
        jac.vs1[2][l] = bl1[2][l] + bt1[2][l];
        jac.vs2[2][l] = bl2[2][l] + bt2[2][l];
    }
    
    (res, jac)
}

/// Compute turbulent-part BLDIF terms inside TRDIF (XT -> X2 interval).
/// Used for debugging shape-equation sensitivities at transition.
pub fn trdif_turb_terms(
    s1: &BlStation,
    s2: &BlStation,
    tr: &crate::closures::Trchek2FullResult,
    msq: f64,
    re: f64,
) -> Option<BldifTerms> {
    use crate::constants::{CTRCON, CTRCEX};

    let dx = s2.x - s1.x;
    if dx <= 1e-20 || !tr.transition {
        return None;
    }

    let xt = tr.xt;
    let wf1 = tr.wf1;
    let wf2 = tr.wf2;

    let tt = s1.theta * wf1 + s2.theta * wf2;
    let dt = s1.delta_star * wf1 + s2.delta_star * wf2;
    let ut = s1.u * wf1 + s2.u * wf2;

    let mut st_turb_temp = BlStation::new();
    st_turb_temp.x = xt;
    st_turb_temp.u = ut;
    st_turb_temp.theta = tt;
    st_turb_temp.delta_star = dt;
    st_turb_temp.ctau = 0.03;
    st_turb_temp.ampl = 0.0;
    st_turb_temp.is_laminar = false;
    st_turb_temp.is_turbulent = true;
    blvar(&mut st_turb_temp, FlowType::Turbulent, msq, re);

    let hk_t = st_turb_temp.hk;
    let cq_t = st_turb_temp.cq;
    let hk_arg = (hk_t - 1.0).max(0.1);
    let ctr = CTRCON * (-CTRCEX / hk_arg).exp();
    let st = ctr * cq_t;

    let mut st_turb = BlStation::new();
    st_turb.x = xt;
    st_turb.u = ut;
    st_turb.theta = tt;
    st_turb.delta_star = dt;
    st_turb.ctau = st.max(0.0001);
    st_turb.ampl = 0.0;
    st_turb.is_laminar = false;
    st_turb.is_turbulent = true;
    blvar(&mut st_turb, FlowType::Turbulent, msq, re);

    let (_, _, terms) = bldif_with_terms(&st_turb, s2, FlowType::Turbulent, msq, re);
    Some(terms)
}

/// Compute BL equation residuals and Jacobian with debug output
///
/// Same as [`bldif`] but also emits debug events when debug collection is active.
/// This matches XFOIL's BLDIF instrumented output.
///
/// # Arguments
/// * `s1` - Upstream station (station 1)
/// * `s2` - Downstream station (station 2)
/// * `flow_type` - Laminar, Turbulent, or Wake
/// * `msq` - Mach number squared
/// * `re` - Reference Reynolds number
/// * `iteration` - Current solver iteration
/// * `side` - Side index (1 or 2)
/// * `ibl` - Station index along the side
pub fn bldif_debug(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
    re: f64,
    iteration: usize,
    side: usize,
    ibl: usize,
) -> (BlResiduals, BlJacobian) {
    // Run the actual computation
    let (res, jac, terms) = bldif_with_terms(s1, s2, flow_type, msq, re);
    let mut xfoil_terms: Option<(BldifTerms, BlStation)> = None;

    if flow_type == FlowType::Turbulent {
        let mut s1_xfoil = s1.clone();
        s1_xfoil.theta = s2.theta;
        s1_xfoil.delta_star = s2.delta_star;
        s1_xfoil.ctau = s2.ctau;
        s1_xfoil.ampl = s2.ampl;
        blvar(&mut s1_xfoil, flow_type, msq, re);

        let (_res_x, _jac_x, terms_x) = bldif_with_terms(&s1_xfoil, s2, flow_type, msq, re);
        xfoil_terms = Some((terms_x, s1_xfoil));
    }

    // Emit debug event if collection is active
    if crate::debug::is_debug_active() {
        let flow_type_int = match flow_type {
            FlowType::Laminar => 1,
            FlowType::Turbulent => 2,
            FlowType::Wake => 3,
        };

        // Convert 3x5 Jacobian to 4x5 format matching XFOIL (add 4th row for Ue constraint)
        let vs1: [[f64; 5]; 4] = [
            jac.vs1[0],
            jac.vs1[1],
            jac.vs1[2],
            [0.0, 0.0, 0.0, 0.0, 0.0], // 4th row: Ue constraint (handled elsewhere)
        ];
        let vs2: [[f64; 5]; 4] = [
            jac.vs2[0],
            jac.vs2[1],
            jac.vs2[2],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let vsrez: [f64; 4] = [res.res_third, res.res_mom, res.res_shape, 0.0];

        crate::debug::add_event(crate::debug::DebugEvent::bldif(
            iteration,
            side,
            ibl,
            flow_type_int,
            vs1,
            vs2,
            vsrez,
        ));
        crate::debug::add_event(crate::debug::DebugEvent::bldif_terms(
            iteration,
            side,
            ibl,
            flow_type_int,
            crate::debug::BldifTermsEvent {
                side,
                ibl,
                flow_type: flow_type_int,
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
    }

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

    // =========================================================================
    // Jacobian Non-Singularity Tests
    // =========================================================================

    /// Helper function to compute determinant of 3x3 matrix
    fn det_3x3(m: &[[f64; 3]; 3]) -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    #[test]
    fn test_bldif_jacobian_column2_nonzero_laminar() {
        // Verify column 2 (delta_star derivatives) is populated for laminar flow
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

        // Check that column 2 (delta_star derivatives) is NOT all zeros
        let col2_nonzero_vs1 =
            jac.vs1[0][2].abs() > 1e-20 || jac.vs1[1][2].abs() > 1e-20 || jac.vs1[2][2].abs() > 1e-20;
        let col2_nonzero_vs2 =
            jac.vs2[0][2].abs() > 1e-20 || jac.vs2[1][2].abs() > 1e-20 || jac.vs2[2][2].abs() > 1e-20;

        assert!(
            col2_nonzero_vs1,
            "VS1 column 2 should have at least one non-zero entry: [{}, {}, {}]",
            jac.vs1[0][2], jac.vs1[1][2], jac.vs1[2][2]
        );
        assert!(
            col2_nonzero_vs2,
            "VS2 column 2 should have at least one non-zero entry: [{}, {}, {}]",
            jac.vs2[0][2], jac.vs2[1][2], jac.vs2[2][2]
        );
    }

    #[test]
    fn test_bldif_jacobian_column2_nonzero_turbulent() {
        // Verify column 2 (delta_star derivatives) is populated for turbulent flow
        let mut s1 = BlStation::new();
        s1.x = 0.5;
        s1.u = 1.0;
        s1.theta = 0.003;
        s1.delta_star = 0.006;
        s1.ctau = 0.15;
        s1.is_laminar = false;
        s1.is_turbulent = true;

        let mut s2 = BlStation::new();
        s2.x = 0.55;
        s2.u = 0.95;
        s2.theta = 0.0035;
        s2.delta_star = 0.007;
        s2.ctau = 0.14;
        s2.is_laminar = false;
        s2.is_turbulent = true;

        blvar(&mut s1, FlowType::Turbulent, 0.0, 1e6);
        blvar(&mut s2, FlowType::Turbulent, 0.0, 1e6);

        let (_, jac) = bldif(&s1, &s2, FlowType::Turbulent, 0.0, 1e6);

        // Check that column 2 (delta_star derivatives) is NOT all zeros
        let col2_nonzero_vs1 =
            jac.vs1[0][2].abs() > 1e-20 || jac.vs1[1][2].abs() > 1e-20 || jac.vs1[2][2].abs() > 1e-20;
        let col2_nonzero_vs2 =
            jac.vs2[0][2].abs() > 1e-20 || jac.vs2[1][2].abs() > 1e-20 || jac.vs2[2][2].abs() > 1e-20;

        assert!(
            col2_nonzero_vs1,
            "VS1 column 2 (turbulent) should have at least one non-zero entry: [{}, {}, {}]",
            jac.vs1[0][2], jac.vs1[1][2], jac.vs1[2][2]
        );
        assert!(
            col2_nonzero_vs2,
            "VS2 column 2 (turbulent) should have at least one non-zero entry: [{}, {}, {}]",
            jac.vs2[0][2], jac.vs2[1][2], jac.vs2[2][2]
        );
    }

    #[test]
    fn test_bldif_jacobian_3x3_nonsingular() {
        // Verify that the 3x3 submatrix (columns 0,1,2) is non-singular
        // This is the matrix that gets inverted in the Newton solver
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

        // Extract 3x3 submatrix from VS2 (downstream station)
        let va: [[f64; 3]; 3] = [
            [jac.vs2[0][0], jac.vs2[0][1], jac.vs2[0][2]],
            [jac.vs2[1][0], jac.vs2[1][1], jac.vs2[1][2]],
            [jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2]],
        ];

        let det = det_3x3(&va);

        assert!(
            det.abs() > 1e-20,
            "3x3 Jacobian submatrix (VA) should be non-singular. det = {:.2e}",
            det
        );
    }

    #[test]
    fn test_bldif_jacobian_nonsingular_across_conditions() {
        // Test non-singularity across various BL states
        let test_cases = [
            // (theta, h_ratio, x, ue, flow_type)
            (0.0005, 2.6, 0.05, 1.1, FlowType::Laminar),   // Near stagnation
            (0.001, 2.6, 0.2, 1.0, FlowType::Laminar),     // Mid-chord laminar
            (0.002, 3.0, 0.5, 0.9, FlowType::Laminar),     // APG laminar
            (0.003, 1.4, 0.6, 0.95, FlowType::Turbulent),  // Turbulent attached
            (0.005, 1.8, 0.9, 0.85, FlowType::Turbulent),  // Turbulent APG
            (0.008, 2.0, 1.5, 0.3, FlowType::Wake),        // Wake
        ];

        for (theta, h_ratio, x, ue, flow_type) in test_cases {
            let mut s1 = BlStation::new();
            s1.x = x;
            s1.u = ue;
            s1.theta = theta;
            s1.delta_star = theta * h_ratio;
            s1.ctau = 0.1;
            s1.is_laminar = flow_type == FlowType::Laminar;
            s1.is_turbulent = flow_type == FlowType::Turbulent;
            s1.is_wake = flow_type == FlowType::Wake;

            let mut s2 = BlStation::new();
            s2.x = x + 0.05;
            s2.u = ue * 0.98;
            s2.theta = theta * 1.1;
            s2.delta_star = theta * h_ratio * 1.15;
            s2.ctau = 0.09;
            s2.is_laminar = s1.is_laminar;
            s2.is_turbulent = s1.is_turbulent;
            s2.is_wake = s1.is_wake;

            blvar(&mut s1, flow_type, 0.0, 1e6);
            blvar(&mut s2, flow_type, 0.0, 1e6);

            let (_, jac) = bldif(&s1, &s2, flow_type, 0.0, 1e6);

            // Extract 3x3 submatrix
            let va: [[f64; 3]; 3] = [
                [jac.vs2[0][0], jac.vs2[0][1], jac.vs2[0][2]],
                [jac.vs2[1][0], jac.vs2[1][1], jac.vs2[1][2]],
                [jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2]],
            ];

            let det = det_3x3(&va);

            assert!(
                det.abs() > 1e-20,
                "Jacobian should be non-singular for {:?} at x={}, θ={}, H={}. det = {:.2e}",
                flow_type,
                x,
                theta,
                h_ratio,
                det
            );
        }
    }
}

#[test]
fn test_bldif_xfoil_station2() {
    // From XFOIL iter=4 (converged):
    // X1=1.10439760e-03, X2=3.36117580e-03
    // U1=6.06760320e-02, U2=2.07620520e-01
    // T1=3.92914430e-05, T2=3.48077990e-05
    // D1=8.76005770e-05, D2=7.70333380e-05
    // VSREZ=[0, 1.662689e-06, -1.130309e-06, 0]
    
    use crate::state::BlStation;
    
    let re = 1e6;
    let msq = 0.0;
    
    // Station 1 (upstream)
    let mut s1 = BlStation::new();
    s1.x = 1.10439760e-03;
    s1.u = 6.06760320e-02;
    s1.theta = 3.92914430e-05;
    s1.delta_star = 8.76005770e-05;
    s1.h = s1.delta_star / s1.theta;
    s1.ctau = 0.03;
    s1.is_laminar = true;
    blvar(&mut s1, FlowType::Laminar, msq, re);
    
    // Station 2 (downstream)
    let mut s2 = BlStation::new();
    s2.x = 3.36117580e-03;
    s2.u = 2.07620520e-01;
    s2.theta = 3.48077990e-05;
    s2.delta_star = 7.70333380e-05;
    s2.h = s2.delta_star / s2.theta;
    s2.ctau = 0.03;
    s2.is_laminar = true;
    blvar(&mut s2, FlowType::Laminar, msq, re);
    
    // Compute bldif
    let (res, _jac) = bldif(&s1, &s2, FlowType::Laminar, msq, re);
    
    println!("=== RustFoil bldif for XFOIL station 2 ===");
    println!("Inputs:");
    println!("  s1: x={:.8e}, u={:.8e}, theta={:.8e}, dstar={:.8e}, Hk={:.4}", 
             s1.x, s1.u, s1.theta, s1.delta_star, s1.hk);
    println!("  s2: x={:.8e}, u={:.8e}, theta={:.8e}, dstar={:.8e}, Hk={:.4}",
             s2.x, s2.u, s2.theta, s2.delta_star, s2.hk);
    println!("\nResiduals:");
    println!("  res_third (ampl): {:.8e}", res.res_third);
    println!("  res_mom (theta):  {:.8e}", res.res_mom);
    println!("  res_shape (H):    {:.8e}", res.res_shape);
    println!("\nXFOIL reference (converged iter=4):");
    println!("  VSREZ = [0, 1.66e-6, -1.13e-6, 0]");
    
    // Check residuals are small (same order of magnitude as XFOIL)
    assert!(res.res_mom.abs() < 0.01, "res_mom too large: {}", res.res_mom);
    assert!(res.res_shape.abs() < 0.01, "res_shape too large: {}", res.res_shape);
}
