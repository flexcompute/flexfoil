//! XFOIL SETBL/BLPRV/BLKIN, with Qinf=1, gamma=1.4 and the pinned reference viscosity ratio.

/// Karman–Tsien edge speed and dU/dUI, where UI is the incompressible velocity.
pub fn edge_velocity(ui: f64, freestream_msq: f64) -> (f64, f64) {
    let beta = (1.0 - freestream_msq).sqrt();
    let tk = freestream_msq / (1.0 + beta).powi(2);
    let denominator = 1.0 - tk * ui * ui;
    let u = ui * (1.0 - tk) / denominator;
    let u_ui = (1.0 + tk * (2.0 * u * ui - 1.0)) / denominator;
    (u, u_ui)
}

pub fn local_mach_squared(u: f64, freestream_msq: f64) -> f64 {
    let hstinv = 0.4 * freestream_msq / (1.0 + 0.2 * freestream_msq);
    u * u * hstinv / (0.4 * (1.0 - 0.5 * u * u * hstinv))
}

#[derive(Debug, Clone, Copy)]
pub struct EdgeThermodynamics {
    pub msq: f64,
    pub m_u: f64,
    pub m_ms: f64,
    pub density: f64,
    pub r_theta: f64,
    pub rt_u: f64,
    pub rt_ms: f64,
}

/// Derivatives are at fixed physical edge speed U, as in BLKIN (not at fixed UI).
pub fn edge_thermodynamics(u: f64, theta: f64, freestream_msq: f64, re: f64) -> EdgeThermodynamics {
    // Official 6.99 SETBL leaves common-block HVRAT unset; instrumented full runs use zero.
    edge_thermodynamics_at_hvrat(u, theta, freestream_msq, re, 0.0)
}

fn edge_thermodynamics_at_hvrat(
    u: f64,
    theta: f64,
    freestream_msq: f64,
    re: f64,
    hvrat: f64,
) -> EdgeThermodynamics {
    let stagnation = 1.0 + 0.2 * freestream_msq;
    let hstinv = 0.4 * freestream_msq / stagnation;
    let hstinv_ms = 0.4 / stagnation.powi(2);
    let msq = local_mach_squared(u, freestream_msq);
    let tr = 1.0 + 0.2 * msq;
    let m_u = 2.0 * msq * tr / u;
    let m_ms = u * u * tr / (0.4 * (1.0 - 0.5 * u * u * hstinv)) * hstinv_ms;
    let density = stagnation.powf(2.5) * tr.powf(-2.5);
    let density_u = -density / tr * 0.5 * m_u;
    let density_ms = -density / tr * 0.5 * m_ms + 0.5 * density / stagnation;
    let he_inf = 1.0 - 0.5 * hstinv;
    let reybl = re * he_inf.powf(1.5) * (1.0 + hvrat) / (he_inf + hvrat);
    let reybl_ms = reybl * (1.5 / he_inf - 1.0 / (he_inf + hvrat)) * (-0.5 * hstinv_ms);
    let he = 1.0 - 0.5 * u * u * hstinv;
    let viscosity = he.powf(1.5) * (1.0 + hvrat) / (he + hvrat) / reybl;
    let viscosity_he = viscosity * (1.5 / he - 1.0 / (he + hvrat));
    let viscosity_u = viscosity_he * (-u * hstinv);
    let viscosity_ms = -viscosity / reybl * reybl_ms + viscosity_he * (-0.5 * u * u * hstinv_ms);
    let r_theta = density * u * theta / viscosity;
    EdgeThermodynamics {
        msq,
        m_u,
        m_ms,
        density,
        r_theta,
        rt_u: r_theta * (1.0 / u + density_u / density - viscosity_u / viscosity),
        rt_ms: r_theta * (density_ms / density - viscosity_ms / viscosity),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn matches_official_blprv_blkin_probe() {
        let (u, du) = edge_velocity(1.2, 0.25_f64.powi(2));
        let t = edge_thermodynamics_at_hvrat(u, 0.001, 0.25_f64.powi(2), 1e6, 0.35);
        assert!((u - 1.2087209494234783).abs() < 1e-14);
        assert!((du - 1.0551819772484472).abs() < 1e-14);
        assert!((t.msq - 0.09184214347131715).abs() < 1e-14);
        assert!((t.r_theta - 1196.6423833592621).abs() < 1e-9);
        assert!((t.rt_u * du - 977.9685407005923).abs() < 1e-9);
    }
    #[test]
    fn derivatives_match_finite_differences() {
        for msq in [0.0625, 0.25] {
            for u in [0.8, 1.0, 1.2] {
                let step = 1e-6;
                let t = edge_thermodynamics(u, 0.001, msq, 1e6);
                let lo = edge_thermodynamics(u - step, 0.001, msq, 1e6);
                let hi = edge_thermodynamics(u + step, 0.001, msq, 1e6);
                assert!(((hi.r_theta - lo.r_theta) / (2.0 * step) - t.rt_u).abs() < 1e-6);
                assert!(((hi.msq - lo.msq) / (2.0 * step) - t.m_u).abs() < 1e-8);
                let lo = edge_thermodynamics(u, 0.001, msq - step, 1e6);
                let hi = edge_thermodynamics(u, 0.001, msq + step, 1e6);
                assert!(((hi.r_theta - lo.r_theta) / (2.0 * step) - t.rt_ms).abs() < 1e-6);
            }
        }
    }
}

#[cfg(test)]
mod jacobian_tests {
    use crate::{bldif, blvar, BlStation, FlowType};

    #[test]
    fn compressible_velocity_jacobian_matches_residual_differences() {
        for flow in [FlowType::Laminar, FlowType::Turbulent, FlowType::Wake] {
            let mut first = BlStation::new();
            first.x = 0.2;
            first.u = 1.2;
            first.theta = 0.001;
            first.delta_star = if flow == FlowType::Laminar {
                0.0026
            } else {
                0.0016
            };
            first.ctau = 0.03;
            first.is_wake = flow == FlowType::Wake;
            let mut second = first.clone();
            second.x = 0.22;
            second.u *= 0.98;
            second.theta *= 1.1;
            second.delta_star *= 1.05;
            second.ctau *= 1.03;
            second.ampl = 0.5;
            blvar(&mut first, flow, 0.0625, 1e6);
            blvar(&mut second, flow, 0.0625, 1e6);
            let (_, jac) = bldif(&first, &second, flow, 0.0625, 1e6);
            for side in 0..2 {
                {
                    let variable = 3;
                    let step = 1e-6;
                    let evaluate = |offset: f64| {
                        let mut a = first.clone();
                        let mut b = second.clone();
                        let station = if side == 0 { &mut a } else { &mut b };
                        station.u += offset;
                        blvar(station, flow, 0.0625, 1e6);
                        let (res, _) = bldif(&a, &b, flow, 0.0625, 1e6);
                        [res.res_third, res.res_mom, res.res_shape]
                    };
                    let lo = evaluate(-step);
                    let hi = evaluate(step);
                    for equation in 0..3 {
                        let numerical = -(hi[equation] - lo[equation]) / (2.0 * step);
                        let analytic = if side == 0 {
                            jac.vs1[equation][variable]
                        } else {
                            jac.vs2[equation][variable]
                        };
                        assert!((numerical - analytic).abs() <= 2e-5 * (1.0 + analytic.abs()),
                            "{flow:?} side={side} variable={variable} equation={equation}: {analytic} vs {numerical}");
                    }
                }
            }
        }
    }
}
