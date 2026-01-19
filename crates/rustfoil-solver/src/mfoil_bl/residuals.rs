//! Residual equations for mfoil-style boundary layer solver.
//!
//! Each station has 3 residual equations:
//! 1. Momentum (von Kármán integral)
//! 2. Shape (H* evolution)
//! 3. Lag (turbulent shear lag) or Amplification (laminar N-factor)

use super::closures::*;
use super::state::MfoilParam;

/// Residual result from a station calculation.
#[derive(Debug, Clone)]
pub struct StationResidual {
    /// 3-element residual vector [R_mom, R_shape, R_lag]
    pub r: [f64; 3],
    /// 3x8 Jacobian w.r.t. [U1, U2] (flattened row-major)
    pub r_u: [[f64; 8]; 3],
    /// 3x2 Jacobian w.r.t. [x1, x2]
    pub r_x: [[f64; 2]; 3],
}

impl Default for StationResidual {
    fn default() -> Self {
        Self {
            r: [0.0; 3],
            r_u: [[0.0; 8]; 3],
            r_x: [[0.0; 2]; 3],
        }
    }
}

/// Calculate residual and Jacobians at a non-transition station.
///
/// # Arguments
/// * `param` - BL parameters
/// * `x` - Arc-length coordinates [x1, x2]
/// * `u` - States at nodes [[U1], [U2]], each is [θ, δ*, sa, Ue]
/// * `aux` - Auxiliary data [aux1, aux2] (wake gap)
///
/// # Returns
/// Station residual with R, R_U, R_x
pub fn residual_station(
    param: &MfoilParam,
    x: [f64; 2],
    u: [[f64; 4]; 2],
    aux: [f64; 2],
) -> StationResidual {
    let mut result = StationResidual::default();

    // Copy states, modifying delta_star to remove wake gap
    let mut u1 = u[0];
    let mut u2 = u[1];
    u1[1] -= aux[0];
    u2[1] -= aux[1];

    // Log changes
    let th1 = u1[0].max(1e-12);
    let th2 = u2[0].max(1e-12);
    let thlog = (th2 / th1).ln();
    let thlog_u = [-1.0 / th1, 0.0, 0.0, 0.0, 1.0 / th2, 0.0, 0.0, 0.0];

    // Compressibility-corrected edge velocity log
    let (uk1, uk1_u) = get_uk(u1[3], param);
    let (uk2, uk2_u) = get_uk(u2[3], param);
    let uelog = (uk2 / uk1).ln();
    let uelog_u = [
        0.0, 0.0, 0.0, -uk1_u / uk1,
        0.0, 0.0, 0.0, uk2_u / uk2,
    ];

    // x log and difference
    let x1 = x[0].max(1e-12);
    let x2 = x[1].max(1e-12);
    let xlog = (x2 / x1).ln();
    let xlog_x = [-1.0 / x1, 1.0 / x2];
    let dx = x2 - x1;
    let dx_x = [-1.0, 1.0];

    // Upwinding factor
    let (upw, upw_u) = get_upw(&u1, &u2, param);

    // Shape parameter H
    let (h1, h1_u1) = get_h(&u1);
    let (h2, h2_u2) = get_h(&u2);
    let h = 0.5 * (h1 + h2);
    let h_u = [
        0.5 * h1_u1[0], 0.5 * h1_u1[1], 0.5 * h1_u1[2], 0.5 * h1_u1[3],
        0.5 * h2_u2[0], 0.5 * h2_u2[1], 0.5 * h2_u2[2], 0.5 * h2_u2[3],
    ];

    // Kinetic energy shape parameter Hs (averaged)
    let (hs1, hs1_u1) = get_hs(&u1, param);
    let (hs2, hs2_u2) = get_hs(&u2, param);
    let (hs, hs_u) = upwind(0.5, &[0.0; 8], hs1, &hs1_u1, hs2, &hs2_u2);

    // Hs log change
    let hslog = (hs2 / hs1).ln();
    let mut hslog_u = [0.0; 8];
    for i in 0..4 {
        hslog_u[i] = -hs1_u1[i] / hs1;
        hslog_u[4 + i] = hs2_u2[i] / hs2;
    }

    // Wake shape parameter Hw
    let (hw1, hw1_u1) = get_hw(&u1, aux[0]);
    let (hw2, hw2_u2) = get_hw(&u2, aux[1]);
    let hw = 0.5 * (hw1 + hw2);
    let hw_u = [
        0.5 * hw1_u1[0], 0.5 * hw1_u1[1], 0.5 * hw1_u1[2], 0.5 * hw1_u1[3],
        0.5 * hw2_u2[0], 0.5 * hw2_u2[1], 0.5 * hw2_u2[2], 0.5 * hw2_u2[3],
    ];

    // Squared Mach number (averaged)
    let (ms1, ms1_u1) = get_mach2(&u1, param);
    let (ms2, ms2_u2) = get_mach2(&u2, param);
    let (ms, ms_u) = upwind(0.5, &[0.0; 8], ms1, &ms1_u1, ms2, &ms2_u2);

    // Handle similarity station
    let (thlog, thlog_u, hslog, hslog_u, uelog, uelog_u, xlog, xlog_x, dx, dx_x) = if param.simi {
        let thlog = 0.0;
        let thlog_u = [0.0; 8];
        let hslog = 0.0;
        let hslog_u = [0.0; 8];
        let uelog = 1.0;
        let uelog_u = [0.0; 8];
        let xlog = 1.0;
        let xlog_x = [0.0, 0.0];
        let dx = 0.5 * (x[0] + x[1]);
        let dx_x = [0.5, 0.5];
        (thlog, thlog_u, hslog, hslog_u, uelog, uelog_u, xlog, xlog_x, dx, dx_x)
    } else {
        (thlog, thlog_u, hslog, hslog_u, uelog, uelog_u, xlog, xlog_x, dx, dx_x)
    };

    // cf*x/θ (symmetrical average)
    let um = [
        0.5 * (u1[0] + u2[0]),
        0.5 * (u1[1] + u2[1]),
        0.5 * (u1[2] + u2[2]),
        0.5 * (u1[3] + u2[3]),
    ];
    let (cfxt1, cfxt1_u1, cfxt1_x1) = get_cfxt(&u1, x1, param);
    let (cfxt2, cfxt2_u2, cfxt2_x2) = get_cfxt(&u2, x2, param);
    let (cfxtm, cfxtm_um, cfxtm_xm) = get_cfxt(&um, 0.5 * (x1 + x2), param);
    
    let cfxt = 0.25 * cfxt1 + 0.5 * cfxtm + 0.25 * cfxt2;
    let mut cfxt_u = [0.0; 8];
    for i in 0..4 {
        cfxt_u[i] = 0.25 * cfxt1_u1[i] + 0.25 * cfxtm_um[i];
        cfxt_u[4 + i] = 0.25 * cfxtm_um[i] + 0.25 * cfxt2_u2[i];
    }
    let cfxt_x = [0.25 * cfxt1_x1 + 0.25 * cfxtm_xm, 0.25 * cfxtm_xm + 0.25 * cfxt2_x2];

    // Momentum equation: R_mom = thlog + (2+H+Hw-Ms)*uelog - 0.5*xlog*cfxt
    result.r[0] = thlog + (2.0 + h + hw - ms) * uelog - 0.5 * xlog * cfxt;
    for i in 0..8 {
        result.r_u[0][i] = thlog_u[i]
            + (h_u[i] + hw_u[i] - ms_u[i]) * uelog
            + (2.0 + h + hw - ms) * uelog_u[i]
            - 0.5 * xlog * cfxt_u[i];
    }
    result.r_x[0][0] = -0.5 * xlog_x[0] * cfxt - 0.5 * xlog * cfxt_x[0];
    result.r_x[0][1] = -0.5 * xlog_x[1] * cfxt - 0.5 * xlog * cfxt_x[1];

    // Dissipation function cDi*x/θ (upwinded)
    let (cdixt1, cdixt1_u1, cdixt1_x1) = get_cdixt(&u1, x1, param);
    let (cdixt2, cdixt2_u2, cdixt2_x2) = get_cdixt(&u2, x2, param);
    let (cdixt, cdixt_u) = upwind(upw, &upw_u, cdixt1, &cdixt1_u1, cdixt2, &cdixt2_u2);
    let cdixt_x = [(1.0 - upw) * cdixt1_x1, upw * cdixt2_x2];

    // cf*x/θ (upwinded for shape equation)
    let (cfxtu, cfxtu_u) = upwind(upw, &upw_u, cfxt1, &cfxt1_u1, cfxt2, &cfxt2_u2);
    let cfxtu_x = [(1.0 - upw) * cfxt1_x1, upw * cfxt2_x2];

    // Density shape parameter Hss (averaged)
    let (hss1, hss1_u1) = get_hss(&u1, param);
    let (hss2, hss2_u2) = get_hss(&u2, param);
    let (hss, hss_u) = upwind(0.5, &[0.0; 8], hss1, &hss1_u1, hss2, &hss2_u2);

    // Shape equation: R_shape = Hslog + (2*Hss/Hs + 1 - H - Hw)*uelog + xlog*(0.5*cfxtu - cDixt)
    result.r[1] = hslog + (2.0 * hss / hs + 1.0 - h - hw) * uelog + xlog * (0.5 * cfxtu - cdixt);
    for i in 0..8 {
        result.r_u[1][i] = hslog_u[i]
            + (2.0 * hss_u[i] / hs - 2.0 * hss / hs.powi(2) * hs_u[i] - h_u[i] - hw_u[i]) * uelog
            + (2.0 * hss / hs + 1.0 - h - hw) * uelog_u[i]
            + xlog * (0.5 * cfxtu_u[i] - cdixt_u[i]);
    }
    result.r_x[1][0] = xlog_x[0] * (0.5 * cfxtu - cdixt) + xlog * (0.5 * cfxtu_x[0] - cdixt_x[0]);
    result.r_x[1][1] = xlog_x[1] * (0.5 * cfxtu - cdixt) + xlog * (0.5 * cfxtu_x[1] - cdixt_x[1]);

    // Third equation: Shear lag (turbulent) or Amplification (laminar)
    if param.turb {
        // Turbulent shear lag equation
        let sa1 = u1[2].max(1e-12);
        let sa2 = u2[2].max(1e-12);
        let salog = (sa2 / sa1).ln();
        let salog_u = [0.0, 0.0, -1.0 / sa1, 0.0, 0.0, 0.0, 1.0 / sa2, 0.0];

        // BL thickness measure (averaged)
        let (de1, de1_u1) = get_de(&u1, param);
        let (de2, de2_u2) = get_de(&u2, param);
        let (de, de_u) = upwind(0.5, &[0.0; 8], de1, &de1_u1, de2, &de2_u2);

        // Slip velocity (averaged)
        let (us1, us1_u1) = get_us(&u1, param);
        let (us2, us2_u2) = get_us(&u2, param);
        let (us, us_u) = upwind(0.5, &[0.0; 8], us1, &us1_u1, us2, &us2_u2);

        // Kinematic shape parameter (upwinded)
        let (hk1, hk1_u1) = get_hk(&u1, param);
        let (hk2, hk2_u2) = get_hk(&u2, param);
        let (hk, hk_u) = upwind(upw, &upw_u, hk1, &hk1_u1, hk2, &hk2_u2);

        // Re_theta (averaged)
        let (ret1, ret1_u1) = get_ret(&u1, param);
        let (ret2, ret2_u2) = get_ret(&u2, param);
        let (ret, ret_u) = upwind(0.5, &[0.0; 8], ret1, &ret1_u1, ret2, &ret2_u2);

        // Skin friction (upwinded)
        let (cf1, cf1_u1) = get_cf(&u1, param);
        let (cf2, cf2_u2) = get_cf(&u2, param);
        let (cf, cf_u) = upwind(upw, &upw_u, cf1, &cf1_u1, cf2, &cf2_u2);

        // Delta star (averaged)
        let dsa = 0.5 * (u1[1] + u2[1]);
        let dsa_u = [0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0];

        // Equilibrium 1/Ue * dUe/dx
        let (uq, uq_u) = get_uq(dsa, &dsa_u, cf, &cf_u, hk, &hk_u, ret, &ret_u, param);

        // Equilibrium shear stress (upwinded)
        let (cteq1, cteq1_u1) = get_cteq(&u1, param);
        let (cteq2, cteq2_u2) = get_cteq(&u2, param);
        let (cteq, cteq_u) = upwind(upw, &upw_u, cteq1, &cteq1_u1, cteq2, &cteq2_u2);

        // Root of shear coefficient (upwinded)
        let (saa, saa_u) = upwind(
            upw,
            &upw_u,
            sa1,
            &[0.0, 0.0, 1.0, 0.0],
            sa2,
            &[0.0, 0.0, 1.0, 0.0],
        );

        // Lag coefficient
        let clag = param.slag_k / param.gb / (1.0 + us);
        let mut clag_u = [0.0; 8];
        for i in 0..8 {
            clag_u[i] = -clag / (1.0 + us) * us_u[i];
        }

        // Extra dissipation in wake
        let ald = if param.wake { param.dlr } else { 1.0 };

        // Shear lag equation
        result.r[2] = clag * (cteq - ald * saa) * dx - 2.0 * de * salog
            + 2.0 * de * (uq * dx - uelog) * param.cuq;
        for i in 0..8 {
            result.r_u[2][i] = clag_u[i] * (cteq - ald * saa) * dx
                + clag * (cteq_u[i] - ald * saa_u[i]) * dx
                - 2.0 * de_u[i] * salog
                - 2.0 * de * salog_u[i]
                + 2.0 * de_u[i] * (uq * dx - uelog) * param.cuq
                + 2.0 * de * (uq_u[i] * dx - uelog_u[i]) * param.cuq;
        }
        result.r_x[2][0] = clag * (cteq - ald * saa) * dx_x[0] + 2.0 * de * uq * dx_x[0];
        result.r_x[2][1] = clag * (cteq - ald * saa) * dx_x[1] + 2.0 * de * uq * dx_x[1];
    } else {
        // Laminar amplification equation
        if param.simi {
            // Similarity station: no amplification
            result.r[2] = u1[2] + u2[2];
            result.r_u[2] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
            result.r_x[2] = [0.0, 0.0];
        } else {
            // Amplification rate (averaged)
            let (damp1, damp1_u1) = get_damp(&u1, param);
            let (damp2, damp2_u2) = get_damp(&u2, param);
            let (damp, damp_u) = upwind(0.5, &[0.0; 8], damp1, &damp1_u1, damp2, &damp2_u2);

            // Amplification equation: sa2 - sa1 - damp*dx = 0
            result.r[2] = u2[2] - u1[2] - damp * dx;
            result.r_u[2] = [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
            for i in 0..8 {
                result.r_u[2][i] -= damp_u[i] * dx;
            }
            result.r_x[2][0] = -damp * dx_x[0];
            result.r_x[2][1] = -damp * dx_x[1];
        }
    }

    result
}

/// Calculate residual for a transition station.
///
/// The state U1 should be laminar; U2 should be turbulent.
/// This function calculates the transition location and linearizes it.
///
/// # Arguments
/// * `param` - BL parameters
/// * `x` - Arc-length coordinates [x1, x2]
/// * `u` - States at nodes [[U1], [U2]]
/// * `aux` - Auxiliary data [aux1, aux2]
///
/// # Returns
/// (StationResidual, transition_location_xt)
pub fn residual_transition(
    param: &MfoilParam,
    x: [f64; 2],
    u: [[f64; 4]; 2],
    aux: [f64; 2],
) -> (StationResidual, f64) {
    let mut result = StationResidual::default();

    let u1 = u[0];
    let u2 = u[1];
    let x1 = x[0];
    let x2 = x[1];
    let dx = x2 - x1;

    // Initial guess for transition location
    let mut xt = x1 + 0.5 * dx;

    // Newton iteration to find transition location where amplification = ncrit
    let ncrit = param.ncrit;
    let max_newton = 20;

    for iter in 0..max_newton {
        let w2 = (xt - x1) / dx;
        let w1 = 1.0 - w2;

        // State at transition
        let mut ut = [0.0; 4];
        let mut ut_xt = [0.0; 4];
        for i in 0..4 {
            ut[i] = w1 * u1[i] + w2 * u2[i];
            ut_xt[i] = (u2[i] - u1[i]) / dx;
        }
        ut[2] = ncrit; // Amplification at transition is ncrit
        ut_xt[2] = 0.0;

        // Amplification rates
        let (damp1, _) = get_damp(&u1, param);
        let (dampt, dampt_ut) = get_damp(&ut, param);

        // Residual for transition location: ncrit - sa1 - 0.5*(xt-x1)*(damp1 + dampt) = 0
        let rxt = ncrit - u1[2] - 0.5 * (xt - x1) * (damp1 + dampt);
        let mut rxt_xt = -0.5 * (damp1 + dampt);
        for i in 0..4 {
            rxt_xt -= 0.5 * (xt - x1) * dampt_ut[i] * ut_xt[i];
        }

        let dxt = -rxt / rxt_xt.max(1e-12);
        let dmax = 0.2 * dx * (1.1 - iter as f64 / max_newton as f64);
        let dxt_clamped = dxt.clamp(-dmax, dmax);

        if rxt.abs() < 1e-10 {
            break;
        }

        xt += dxt_clamped;
        xt = xt.clamp(x1 + 1e-6 * dx, x2 - 1e-6 * dx);
    }

    // Compute transition residual as sum of laminar and turbulent contributions
    let w2 = (xt - x1) / dx;
    let w1 = 1.0 - w2;

    // Laminar state at transition
    let mut utl = [0.0; 4];
    for i in 0..4 {
        utl[i] = w1 * u1[i] + w2 * u2[i];
    }
    utl[2] = ncrit;

    // Turbulent state at transition (with cttr)
    let mut param_turb = param.clone();
    param_turb.turb = true;

    let (cttr, _) = get_cttr(&utl, &param_turb);
    let mut utt = utl;
    utt[2] = cttr;

    // Laminar residual (x1 to xt)
    let mut param_lam = param.clone();
    param_lam.turb = false;
    let rl = residual_station(&param_lam, [x1, xt], [u1, utl], aux);

    // Turbulent residual (xt to x2)
    let rt = residual_station(&param_turb, [xt, x2], [utt, u2], aux);

    // Combined residual
    for i in 0..3 {
        result.r[i] = rl.r[i] + rt.r[i];
    }

    // Simplified Jacobian (finite difference for robustness)
    let eps = 1e-7;
    for j in 0..4 {
        // Perturb U1
        let mut u1p = u1;
        u1p[j] += eps;
        let (rp, _) = residual_transition(param, x, [u1p, u2], aux);
        for i in 0..3 {
            result.r_u[i][j] = (rp.r[i] - result.r[i]) / eps;
        }

        // Perturb U2
        let mut u2p = u2;
        u2p[j] += eps;
        let (rp, _) = residual_transition(param, x, [u1, u2p], aux);
        for i in 0..3 {
            result.r_u[i][4 + j] = (rp.r[i] - result.r[i]) / eps;
        }
    }

    // x derivatives (simplified)
    for j in 0..2 {
        let mut xp = x;
        xp[j] += eps;
        let (rp, _) = residual_transition(param, xp, u, aux);
        for i in 0..3 {
            result.r_x[i][j] = (rp.r[i] - result.r[i]) / eps;
        }
    }

    (result, xt)
}
