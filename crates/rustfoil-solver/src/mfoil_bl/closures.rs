//! Closure relations for mfoil-style boundary layer solver.
//!
//! These functions compute various boundary layer parameters and their
//! derivatives with respect to the state vector U = [θ, δ*, sa, Ue].

use super::state::MfoilParam;

/// Compute H = shape parameter = δ*/θ.
///
/// Returns (H, H_U) where H_U is the 4-element linearization.
#[inline]
pub fn get_h(u: &[f64; 4]) -> (f64, [f64; 4]) {
    let th = u[0].max(1e-12);
    let ds = u[1];
    let h = ds / th;
    let h_u = [-h / th, 1.0 / th, 0.0, 0.0];
    (h, h_u)
}

/// Compute Hw = wake gap shape parameter = wgap/θ.
///
/// Returns (Hw, Hw_U) where Hw_U is the 4-element linearization.
#[inline]
pub fn get_hw(u: &[f64; 4], wgap: f64) -> (f64, [f64; 4]) {
    let th = u[0].max(1e-12);
    let hw = wgap / th;
    let hw_u = [-hw / th, 0.0, 0.0, 0.0];
    (hw, hw_u)
}

/// Compute squared Mach number M².
///
/// Returns (M², M²_U) with compressibility correction.
pub fn get_mach2(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    if param.minf > 0.0 {
        let (uk, uk_u) = get_uk(u[3], param);
        let c2 = (param.gam - 1.0) * (param.h0 - 0.5 * uk * uk);
        let c2_uk = (param.gam - 1.0) * (-uk);
        let m2 = uk * uk / c2;
        let m2_uk = 2.0 * uk / c2 - m2 / c2 * c2_uk;
        let m2_u = [0.0, 0.0, 0.0, m2_uk * uk_u];
        (m2, m2_u)
    } else {
        (0.0, [0.0; 4])
    }
}

/// Compute Karman-Tsien corrected speed.
///
/// Returns (uk, uk_u) where uk_u is d(uk)/d(u).
pub fn get_uk(u: f64, param: &MfoilParam) -> (f64, f64) {
    if param.minf > 0.0 {
        let den = 1.0 - param.ktl * (u / param.vinf).powi(2);
        let den_u = -2.0 * param.ktl * u / param.vinf.powi(2);
        let uk = u * (1.0 - param.ktl) / den;
        let uk_u = (1.0 - param.ktl) / den - (uk / den) * den_u;
        (uk, uk_u)
    } else {
        (u, 1.0)
    }
}

/// Compute Hk = kinematic shape parameter.
///
/// Returns (Hk, Hk_U) with compressibility correction.
pub fn get_hk(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (h, h_u) = get_h(u);

    if param.minf > 0.0 {
        let (m2, m2_u) = get_mach2(u, param);
        let den = 1.0 + 0.113 * m2;
        let den_m2 = 0.113;
        let hk = (h - 0.29 * m2) / den;
        let mut hk_u = [0.0; 4];
        for i in 0..4 {
            hk_u[i] = (h_u[i] - 0.29 * m2_u[i]) / den - hk / den * den_m2 * m2_u[i];
        }
        (hk, hk_u)
    } else {
        (h, h_u)
    }
}

/// Compute Re_θ = Reynolds number based on momentum thickness.
///
/// Returns (Re_θ, Re_θ_U).
pub fn get_ret(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    if param.minf > 0.0 {
        let (m2, m2_u) = get_mach2(u, param);
        let (uk, uk_u) = get_uk(u[3], param);
        let gmi = param.gam - 1.0;
        let tr = 1.0 - 0.5 * uk * uk / param.h0;
        let tr_uk = -uk / param.h0;
        
        // Sutherland's law
        let f = tr.powf(1.5) * (1.0 + param.tsrat) / (tr + param.tsrat);
        let f_tr = 1.5 * f / tr - f / (tr + param.tsrat);
        let mu = param.mu0 * f;
        let mu_uk = param.mu0 * f_tr * tr_uk;
        
        // Density from isentropic relations
        let den = 1.0 + 0.5 * gmi * m2;
        let den_m2 = 0.5 * gmi;
        let rho = param.rho0 / den.powf(1.0 / gmi);
        let mut rho_u = [0.0; 4];
        for i in 0..4 {
            rho_u[i] = (-1.0 / gmi) * rho / den * den_m2 * m2_u[i];
        }
        
        let ret = rho * uk * u[0] / mu;
        let mut ret_u = [0.0; 4];
        for i in 0..4 {
            ret_u[i] = rho_u[i] * uk * u[0] / mu;
        }
        ret_u[0] += rho * uk / mu;
        ret_u[3] += (rho * u[0] / mu - ret / mu * mu_uk) * uk_u;
        
        (ret, ret_u)
    } else {
        let ret = param.rho0 * u[0] * u[3] / param.mu0;
        let ret_u = [u[3] / param.mu0, 0.0, 0.0, u[0] / param.mu0];
        (ret, ret_u)
    }
}

/// Compute Hs = H* = kinetic energy shape parameter.
///
/// Returns (Hs, Hs_U).
pub fn get_hs(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (mut hk, mut hk_u) = get_hk(u, param);

    // Limit Hk
    let hk_min = if param.wake { 1.00005 } else { 1.05 };
    if hk < hk_min {
        hk = hk_min;
        hk_u = [0.0; 4];
    }

    if param.turb {
        // Turbulent
        let hsmin = 1.5;
        let dhsinf = 0.015;
        let (mut ret, mut ret_u) = get_ret(u, param);

        // Limit Re_theta
        let mut ho = 4.0;
        let mut ho_u = [0.0; 4];
        if ret > 400.0 {
            ho = 3.0 + 400.0 / ret;
            for i in 0..4 {
                ho_u[i] = -400.0 / ret.powi(2) * ret_u[i];
            }
        }

        let mut reb = ret;
        let mut reb_u = ret_u;
        if ret < 200.0 {
            reb = 200.0;
            reb_u = [0.0; 4];
        }

        let (hs, hs_u) = if hk < ho {
            // Attached branch
            let hr = (ho - hk) / (ho - 1.0);
            let mut hr_u = [0.0; 4];
            for i in 0..4 {
                hr_u[i] = (ho_u[i] - hk_u[i]) / (ho - 1.0) - (ho - hk) / (ho - 1.0).powi(2) * ho_u[i];
            }
            let aa = (2.0 - hsmin - 4.0 / reb) * hr * hr;
            let mut aa_u = [0.0; 4];
            for i in 0..4 {
                aa_u[i] = (4.0 / reb.powi(2) * reb_u[i]) * hr * hr
                    + (2.0 - hsmin - 4.0 / reb) * 2.0 * hr * hr_u[i];
            }
            let hs = hsmin + 4.0 / reb + aa * 1.5 / (hk + 0.5);
            let mut hs_u = [0.0; 4];
            for i in 0..4 {
                hs_u[i] = -4.0 / reb.powi(2) * reb_u[i]
                    + aa_u[i] * 1.5 / (hk + 0.5)
                    - aa * 1.5 / (hk + 0.5).powi(2) * hk_u[i];
            }
            (hs, hs_u)
        } else {
            // Separated branch
            let lrb = reb.ln();
            let mut lrb_u = [0.0; 4];
            for i in 0..4 {
                lrb_u[i] = 1.0 / reb * reb_u[i];
            }
            let aa = hk - ho + 4.0 / lrb;
            let mut aa_u = [0.0; 4];
            for i in 0..4 {
                aa_u[i] = hk_u[i] - ho_u[i] - 4.0 / lrb.powi(2) * lrb_u[i];
            }
            let bb = 0.007 * lrb / aa.powi(2) + dhsinf / hk;
            let mut bb_u = [0.0; 4];
            for i in 0..4 {
                bb_u[i] = 0.007 * (lrb_u[i] / aa.powi(2) - 2.0 * lrb / aa.powi(3) * aa_u[i])
                    - dhsinf / hk.powi(2) * hk_u[i];
            }
            let hs = hsmin + 4.0 / reb + (hk - ho).powi(2) * bb;
            let mut hs_u = [0.0; 4];
            for i in 0..4 {
                hs_u[i] = -4.0 / reb.powi(2) * reb_u[i]
                    + 2.0 * (hk - ho) * (hk_u[i] - ho_u[i]) * bb
                    + (hk - ho).powi(2) * bb_u[i];
            }
            (hs, hs_u)
        };

        // Mach correction
        if param.minf > 0.0 {
            let (m2, m2_u) = get_mach2(u, param);
            let den = 1.0 + 0.014 * m2;
            let den_m2 = 0.014;
            let hs_corr = (hs + 0.028 * m2) / den;
            let mut hs_corr_u = [0.0; 4];
            for i in 0..4 {
                hs_corr_u[i] = (hs_u[i] + 0.028 * m2_u[i]) / den
                    - hs_corr / den * den_m2 * m2_u[i];
            }
            (hs_corr, hs_corr_u)
        } else {
            (hs, hs_u)
        }
    } else {
        // Laminar
        let a = hk - 4.35;
        let (hs, hs_hk) = if hk < 4.35 {
            let num = 0.0111 * a.powi(2) - 0.0278 * a.powi(3);
            let hs = num / (hk + 1.0) + 1.528 - 0.0002 * (a * hk).powi(2);
            let hs_hk = (0.0111 * 2.0 * a - 0.0278 * 3.0 * a.powi(2)) / (hk + 1.0)
                - num / (hk + 1.0).powi(2)
                - 0.0002 * 2.0 * a * hk * (hk + a);
            (hs, hs_hk)
        } else {
            let hs = 0.015 * a.powi(2) / hk + 1.528;
            let hs_hk = 0.015 * 2.0 * a / hk - 0.015 * a.powi(2) / hk.powi(2);
            (hs, hs_hk)
        };
        let mut hs_u = [0.0; 4];
        for i in 0..4 {
            hs_u[i] = hs_hk * hk_u[i];
        }
        (hs, hs_u)
    }
}

/// Compute Hss = density shape parameter.
///
/// Returns (Hss, Hss_U).
pub fn get_hss(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (m2, m2_u) = get_mach2(u, param);
    let (hk, hk_u) = get_hk(u, param);
    let num = 0.064 / (hk - 0.8) + 0.251;
    let mut num_u = [0.0; 4];
    for i in 0..4 {
        num_u[i] = -0.064 / (hk - 0.8).powi(2) * hk_u[i];
    }
    let hss = m2 * num;
    let mut hss_u = [0.0; 4];
    for i in 0..4 {
        hss_u[i] = m2_u[i] * num + m2 * num_u[i];
    }
    (hss, hss_u)
}

/// Compute cf = skin friction coefficient.
///
/// Returns (cf, cf_U).
pub fn get_cf(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    if param.wake {
        return (0.0, [0.0; 4]);
    }

    let (mut hk, mut hk_u) = get_hk(u, param);
    let (ret, ret_u) = get_ret(u, param);

    // Limit Hk
    if !param.turb && hk < 1.05 {
        hk = 1.05;
        hk_u = [0.0; 4];
    }

    if param.turb {
        // Turbulent cf
        let (m2, m2_u) = get_mach2(u, param);
        let fc = (1.0 + 0.5 * (param.gam - 1.0) * m2).sqrt();
        let mut fc_u = [0.0; 4];
        for i in 0..4 {
            fc_u[i] = 0.5 / fc * 0.5 * (param.gam - 1.0) * m2_u[i];
        }

        let mut aa = -1.33 * hk;
        let mut aa_u = [0.0; 4];
        for i in 0..4 {
            aa_u[i] = -1.33 * hk_u[i];
        }

        // Smooth limiting
        if aa < -17.0 {
            let aa_old = aa;
            aa = -20.0 + 3.0 * ((aa + 17.0) / 3.0).exp();
            for i in 0..4 {
                aa_u[i] = (aa + 20.0) / 3.0 * aa_u[i];
            }
        }

        let mut bb = (ret / fc).ln();
        let mut bb_u = [0.0; 4];
        for i in 0..4 {
            bb_u[i] = ret_u[i] / ret - fc_u[i] / fc;
        }
        if bb < 3.0 {
            bb = 3.0;
            bb_u = [0.0; 4];
        }
        bb /= 10.0_f64.ln();
        for i in 0..4 {
            bb_u[i] /= 10.0_f64.ln();
        }

        let cc = -1.74 - 0.31 * hk;
        let mut cc_u = [0.0; 4];
        for i in 0..4 {
            cc_u[i] = -0.31 * hk_u[i];
        }

        let dd = (4.0 - hk / 0.875).tanh();
        let mut dd_u = [0.0; 4];
        for i in 0..4 {
            dd_u[i] = (1.0 - dd * dd) * (-hk_u[i] / 0.875);
        }

        let cf0 = 0.3 * aa.exp() * bb.powf(cc);
        let mut cf0_u = [0.0; 4];
        for i in 0..4 {
            cf0_u[i] = cf0 * aa_u[i]
                + 0.3 * aa.exp() * cc * bb.powf(cc - 1.0) * bb_u[i]
                + cf0 * bb.ln() * cc_u[i];
        }

        let cf = (cf0 + 1.1e-4 * (dd - 1.0)) / fc;
        let mut cf_u = [0.0; 4];
        for i in 0..4 {
            cf_u[i] = (cf0_u[i] + 1.1e-4 * dd_u[i]) / fc - cf / fc * fc_u[i];
        }

        (cf, cf_u)
    } else {
        // Laminar cf
        let (num, num_hk) = if hk < 5.5 {
            let num = 0.0727 * (5.5 - hk).powi(3) / (hk + 1.0) - 0.07;
            let num_hk = 0.0727
                * (3.0 * (5.5 - hk).powi(2) / (hk + 1.0) * (-1.0)
                    - (5.5 - hk).powi(3) / (hk + 1.0).powi(2));
            (num, num_hk)
        } else {
            let num = 0.015 * (1.0 - 1.0 / (hk - 4.5)).powi(2) - 0.07;
            let num_hk = 0.015 * 2.0 * (1.0 - 1.0 / (hk - 4.5)) / (hk - 4.5).powi(2);
            (num, num_hk)
        };

        let cf = num / ret;
        let mut cf_u = [0.0; 4];
        for i in 0..4 {
            cf_u[i] = num_hk / ret * hk_u[i] - num / ret.powi(2) * ret_u[i];
        }

        (cf, cf_u)
    }
}

/// Compute de = BL thickness measure.
///
/// Returns (de, de_U).
pub fn get_de(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (mut hk, mut hk_u) = get_hk(u, param);
    
    // Limit Hk
    let hk_min = if param.wake { 1.00005 } else { 1.05 };
    if hk < hk_min {
        hk = hk_min;
        hk_u = [0.0; 4];
    }

    let aa = 3.15 + 1.72 / (hk - 1.0);
    let mut aa_u = [0.0; 4];
    for i in 0..4 {
        aa_u[i] = -1.72 / (hk - 1.0).powi(2) * hk_u[i];
    }

    let mut de = u[0] * aa + u[1];
    let mut de_u = [aa, 1.0, 0.0, 0.0];
    for i in 0..4 {
        de_u[i] += u[0] * aa_u[i];
    }

    // Cap at 12*theta
    let dmx = 12.0;
    if de > dmx * u[0] {
        de = dmx * u[0];
        de_u = [dmx, 0.0, 0.0, 0.0];
    }

    (de, de_u)
}

/// Compute cteq = √(equilibrium shear stress coefficient).
///
/// Returns (cteq, cteq_U).
pub fn get_cteq(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let cc = 0.5 / (param.ga.powi(2) * param.gb);
    let c = param.gc;
    let (mut hk, mut hk_u) = get_hk(u, param);
    let (hs, hs_u) = get_hs(u, param);
    let (h, h_u) = get_h(u);
    let (ret, ret_u) = get_ret(u, param);
    let (us, us_u) = get_us(u, param);

    // Limit Hk
    let hk_min = if param.wake { 1.00005 } else { 1.05 };
    if hk < hk_min {
        hk = hk_min;
        hk_u = [0.0; 4];
    }

    let (hkc, mut hkc_u) = if param.wake {
        let hkc = hk - 1.0;
        (hkc, hk_u)
    } else {
        let mut hkc = hk - 1.0 - c / ret;
        let mut hkc_u = [0.0; 4];
        for i in 0..4 {
            hkc_u[i] = hk_u[i] + c / ret.powi(2) * ret_u[i];
        }
        if hkc < 0.01 {
            hkc = 0.01;
            hkc_u = [0.0; 4];
        }
        (hkc, hkc_u)
    };

    let num = cc * hs * (hk - 1.0) * hkc.powi(2);
    let mut num_u = [0.0; 4];
    for i in 0..4 {
        num_u[i] = cc
            * (hs_u[i] * (hk - 1.0) * hkc.powi(2)
                + hs * hk_u[i] * hkc.powi(2)
                + hs * (hk - 1.0) * 2.0 * hkc * hkc_u[i]);
    }

    let den = (1.0 - us) * h * hk.powi(2);
    let mut den_u = [0.0; 4];
    for i in 0..4 {
        den_u[i] =
            (-us_u[i]) * h * hk.powi(2) + (1.0 - us) * h_u[i] * hk.powi(2) + (1.0 - us) * h * 2.0 * hk * hk_u[i];
    }

    let cteq = (num / den).sqrt();
    let mut cteq_u = [0.0; 4];
    for i in 0..4 {
        cteq_u[i] = 0.5 / cteq * (num_u[i] / den - num / den.powi(2) * den_u[i]);
    }

    (cteq, cteq_u)
}

/// Compute Us = normalized slip velocity.
///
/// Returns (Us, Us_U).
pub fn get_us(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (mut hk, mut hk_u) = get_hk(u, param);
    let (hs, hs_u) = get_hs(u, param);

    // Limit Hk
    let hk_min = if param.wake { 1.00005 } else { 1.05 };
    if hk < hk_min {
        hk = hk_min;
        hk_u = [0.0; 4];
    }

    let hk1 = hk - 1.0;
    let us = 0.5 * hs * hk1 / hk;
    let mut us_u = [0.0; 4];
    for i in 0..4 {
        us_u[i] = 0.5 * (hs_u[i] * hk1 / hk + hs * hk_u[i] / hk - hs * hk1 / hk.powi(2) * hk_u[i]);
    }

    (us, us_u)
}

/// Compute cttr = √(shear stress coefficient at transition).
///
/// Returns (cttr, cttr_U).
pub fn get_cttr(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let mut param_copy = param.clone();
    param_copy.wake = false; // Transition happens before wake

    let (cteq, cteq_u) = get_cteq(u, &param_copy);
    let (mut hk, mut hk_u) = get_hk(u, param);

    if hk < 1.05 {
        hk = 1.05;
        hk_u = [0.0; 4];
    }

    let c = param.ctau_c * ((-param.ctau_e) / (hk - 1.0)).exp();
    let mut c_u = [0.0; 4];
    for i in 0..4 {
        c_u[i] = c * param.ctau_e / (hk - 1.0).powi(2) * hk_u[i];
    }

    let cttr = c * cteq;
    let mut cttr_u = [0.0; 4];
    for i in 0..4 {
        cttr_u[i] = c_u[i] * cteq + c * cteq_u[i];
    }

    (cttr, cttr_u)
}

/// Compute amplification rate dN/ds.
///
/// Returns (damp, damp_U).
pub fn get_damp(u: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 4]) {
    let (mut hk, mut hk_u) = get_hk(u, param);
    let (ret, ret_u) = get_ret(u, param);

    // Limit Hk
    if hk < 1.05 {
        hk = 1.05;
        hk_u = [0.0; 4];
    }

    // Critical Re_theta for instability (Drela-Giles correlation)
    let hmi = 1.0 / (hk - 1.0);
    let mut hmi_u = [0.0; 4];
    for i in 0..4 {
        hmi_u[i] = -1.0 / (hk - 1.0).powi(2) * hk_u[i];
    }

    let aa = 2.492 * hmi.powf(0.43);
    let mut aa_u = [0.0; 4];
    for i in 0..4 {
        aa_u[i] = 2.492 * 0.43 * hmi.powf(0.43 - 1.0) * hmi_u[i];
    }

    let bb = (14.68 * hmi).tanh();
    let mut bb_u = [0.0; 4];
    for i in 0..4 {
        bb_u[i] = (1.0 - bb * bb) * 14.68 * hmi_u[i];
    }

    // Critical Re_theta
    let rth_crit = 10.0_f64.powf(aa + 0.7 * bb);
    let ln10 = 10.0_f64.ln();
    let mut rth_crit_u = [0.0; 4];
    for i in 0..4 {
        rth_crit_u[i] = rth_crit * ln10 * (aa_u[i] + 0.7 * bb_u[i]);
    }

    // Amplification rate
    if ret < rth_crit {
        // Below critical: no amplification
        (0.0, [0.0; 4])
    } else {
        // Above critical: amplification grows
        let rtr = ret / rth_crit;
        let mut rtr_u = [0.0; 4];
        for i in 0..4 {
            rtr_u[i] = ret_u[i] / rth_crit - ret / rth_crit.powi(2) * rth_crit_u[i];
        }

        // Rate coefficient (simplified from mfoil)
        let dn = 0.01 * ((rtr - 1.0).powi(2) + 0.25).sqrt();
        let dn_rtr = if rtr > 1.0 {
            0.01 * (rtr - 1.0) / ((rtr - 1.0).powi(2) + 0.25).sqrt()
        } else {
            0.0
        };

        let mut dn_u = [0.0; 4];
        for i in 0..4 {
            dn_u[i] = dn_rtr * rtr_u[i];
        }

        // Scale by 1/theta
        let damp = dn * ret / u[0].max(1e-12);
        let mut damp_u = [0.0; 4];
        for i in 0..4 {
            damp_u[i] = dn_u[i] * ret / u[0].max(1e-12) + dn * ret_u[i] / u[0].max(1e-12);
        }
        damp_u[0] -= damp / u[0].max(1e-12);

        (damp.min(10.0), damp_u) // Cap rate
    }
}

/// Compute upwinding factor based on shape parameter change.
///
/// Returns (upw, upw_U) where upw is between 0.5 (centered) and 1.0 (full backward).
pub fn get_upw(u1: &[f64; 4], u2: &[f64; 4], param: &MfoilParam) -> (f64, [f64; 8]) {
    let (hk1, hk1_u1) = get_hk(u1, param);
    let (hk2, hk2_u2) = get_hk(u2, param);

    let hut = 1.0; // Triggering constant
    let c = if param.wake { 1.0 } else { 5.0 };

    let huc = c * hut / hk2.powi(2);
    let mut huc_u = [0.0; 8];
    for i in 0..4 {
        huc_u[4 + i] = -2.0 * huc / hk2 * hk2_u2[i];
    }

    let aa = (hk2 - 1.0) / (hk1 - 1.0);
    let sga = aa.signum();
    let la = (sga * aa).ln();

    let mut la_u = [0.0; 8];
    for i in 0..4 {
        la_u[i] = -1.0 / (hk1 - 1.0) * hk1_u1[i];
        la_u[4 + i] = 1.0 / (hk2 - 1.0) * hk2_u2[i];
    }

    let mut hls = la.powi(2);
    let mut hls_u = [0.0; 8];
    for i in 0..8 {
        hls_u[i] = 2.0 * la * la_u[i];
    }

    if hls > 15.0 {
        hls = 15.0;
        hls_u = [0.0; 8];
    }

    let upw = 1.0 - 0.5 * (-hls * huc).exp();
    let mut upw_u = [0.0; 8];
    for i in 0..8 {
        upw_u[i] = -0.5 * (-hls * huc).exp() * (-hls_u[i] * huc - hls * huc_u[i]);
    }

    (upw, upw_u)
}

/// Compute upwind average of two scalars.
///
/// f = (1-upw)*f1 + upw*f2
///
/// Returns (f, f_U) where f_U is 8-element linearization [f_U1, f_U2].
pub fn upwind(
    upw: f64,
    upw_u: &[f64; 8],
    f1: f64,
    f1_u1: &[f64; 4],
    f2: f64,
    f2_u2: &[f64; 4],
) -> (f64, [f64; 8]) {
    let f = (1.0 - upw) * f1 + upw * f2;
    let mut f_u = [0.0; 8];
    for i in 0..4 {
        f_u[i] = (-upw_u[i]) * f1 + upw_u[i] * f2 + (1.0 - upw) * f1_u1[i];
        f_u[4 + i] = (-upw_u[4 + i]) * f1 + upw_u[4 + i] * f2 + upw * f2_u2[i];
    }
    (f, f_u)
}

/// Compute cf*x/θ for momentum equation.
///
/// Returns (cfxt, cfxt_U, cfxt_x).
pub fn get_cfxt(u: &[f64; 4], x: f64, param: &MfoilParam) -> (f64, [f64; 4], f64) {
    let (cf, cf_u) = get_cf(u, param);
    let th = u[0].max(1e-12);
    let cfxt = cf * x / th;
    let mut cfxt_u = [0.0; 4];
    for i in 0..4 {
        cfxt_u[i] = cf_u[i] * x / th;
    }
    cfxt_u[0] -= cfxt / th;
    let cfxt_x = cf / th;
    (cfxt, cfxt_u, cfxt_x)
}

/// Compute dissipation function cDi*x/θ for shape equation.
///
/// Returns (cdixt, cdixt_U, cdixt_x).
pub fn get_cdixt(u: &[f64; 4], x: f64, param: &MfoilParam) -> (f64, [f64; 4], f64) {
    let (hs, hs_u) = get_hs(u, param);
    let (cf, cf_u) = get_cf(u, param);
    let (us, us_u) = get_us(u, param);
    let th = u[0].max(1e-12);

    // Dissipation coefficient (simplified)
    let cdi = if param.turb {
        // Turbulent dissipation
        let sa = u[2].max(1e-10);
        let ctau = sa * sa;
        0.5 * ctau * us * (1.0 - us) + 0.5 * cf * (1.0 - us) * (1.0 - us)
    } else {
        // Laminar dissipation
        0.5 * cf * (1.0 - us) * (1.0 - us)
    };

    let cdixt = 2.0 * cdi / hs * x / th;

    // Approximate linearization
    let mut cdixt_u = [0.0; 4];
    let eps = 1e-7;
    for i in 0..4 {
        let mut u_pert = *u;
        u_pert[i] += eps;
        let (hs_p, _) = get_hs(&u_pert, param);
        let (cf_p, _) = get_cf(&u_pert, param);
        let (us_p, _) = get_us(&u_pert, param);
        let th_p = u_pert[0].max(1e-12);
        
        let cdi_p = if param.turb {
            let sa_p = u_pert[2].max(1e-10);
            let ctau_p = sa_p * sa_p;
            0.5 * ctau_p * us_p * (1.0 - us_p) + 0.5 * cf_p * (1.0 - us_p) * (1.0 - us_p)
        } else {
            0.5 * cf_p * (1.0 - us_p) * (1.0 - us_p)
        };
        
        let cdixt_p = 2.0 * cdi_p / hs_p * x / th_p;
        cdixt_u[i] = (cdixt_p - cdixt) / eps;
    }

    let cdixt_x = 2.0 * cdi / hs / th;

    (cdixt, cdixt_u, cdixt_x)
}

/// Get equilibrium 1/Ue * dUe/dx for shear lag.
///
/// Returns (uq, uq_U).
pub fn get_uq(
    ds: f64,
    ds_u: &[f64; 8],
    cf: f64,
    cf_u: &[f64; 8],
    hk: f64,
    hk_u: &[f64; 8],
    ret: f64,
    ret_u: &[f64; 8],
    param: &MfoilParam,
) -> (f64, [f64; 8]) {
    let (mut a, mut c) = (param.ga, param.gc);
    if param.wake {
        a *= param.dlr;
        c = 0.0;
    }

    let hk_lim = if param.wake { 1.00005 } else { 1.05 };
    let (hk, hk_u) = if hk < hk_lim {
        (hk_lim, [0.0; 8])
    } else {
        (hk, *hk_u)
    };

    let mut hkc = hk - 1.0 - c / ret;
    let mut hkc_u = [0.0; 8];
    for i in 0..8 {
        hkc_u[i] = hk_u[i] + c / ret.powi(2) * ret_u[i];
    }
    if hkc < 0.01 {
        hkc = 0.01;
        hkc_u = [0.0; 8];
    }

    let ut = 0.5 * cf - (hkc / (a * hk)).powi(2);
    let mut ut_u = [0.0; 8];
    for i in 0..8 {
        ut_u[i] = 0.5 * cf_u[i]
            - 2.0 * (hkc / (a * hk)) * (hkc_u[i] / (a * hk) - hkc / (a * hk.powi(2)) * hk_u[i]);
    }

    let uq = ut / (param.gb * ds);
    let mut uq_u = [0.0; 8];
    for i in 0..8 {
        uq_u[i] = ut_u[i] / (param.gb * ds) - uq / ds * ds_u[i];
    }

    (uq, uq_u)
}
