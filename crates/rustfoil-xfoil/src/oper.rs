use rustfoil_bl::{
    add_event, is_debug_active, BlStation, DebugEvent, SurfaceBlState,
};
use rustfoil_core::Body;
use rustfoil_inviscid::{FactorizedSystem, InviscidSolver};
use serde::{Deserialize, Serialize};

use crate::{
    assembly::setbl,
    canonical_state::rows_to_stations,
    config::{OperatingMode, XfoilOptions},
    error::{Result, XfoilError},
    forces::{compute_panel_forces_from_gamma, update_force_state},
    march::blpini,
    result::XfoilViscousResult,
    solve::blsolv,
    state::XfoilState,
    state_ops::{compute_arc_lengths, gamqv, iblpan, qvfue, specal, stfind, stmove, uedginit, uicalc, xicalc},
    update::update,
    wake_panel::{qdcalc, qwcalc, xywake},
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum AlphaSpec {
    AlphaDeg(f64),
    TargetCl(f64),
}

pub fn coords_from_body(body: &Body) -> Vec<(f64, f64)> {
    let panels = body.panels();
    let mut node_x: Vec<f64> = panels.iter().map(|panel| panel.p1.x).collect();
    let mut node_y: Vec<f64> = panels.iter().map(|panel| panel.p1.y).collect();
    if let Some(last) = panels.last() {
        if (last.p2.x - panels[0].p1.x).abs() > 1.0e-10 || (last.p2.y - panels[0].p1.y).abs() > 1.0e-10 {
            node_x.push(last.p2.x);
            node_y.push(last.p2.y);
        }
    }
    node_x
        .iter()
        .zip(node_y.iter())
        .map(|(&x, &y)| (x, y))
        .collect()
}

pub fn solve_body_oper_point(
    body: &Body,
    spec: AlphaSpec,
    options: &XfoilOptions,
) -> Result<XfoilViscousResult> {
    let coords = coords_from_body(body);
    solve_coords_oper_point(&body.name, &coords, spec, options)
}

pub fn build_state_from_coords(
    name: &str,
    coords: &[(f64, f64)],
    spec: AlphaSpec,
    options: &XfoilOptions,
) -> Result<(XfoilState, FactorizedSystem)> {
    options
        .validate()
        .map_err(|msg| XfoilError::Message(msg.to_string()))?;

    let solver = InviscidSolver::new();
    let factorized = solver.factorize(coords)?;
    let alpha_deg = match spec {
        AlphaSpec::AlphaDeg(alpha) => alpha,
        AlphaSpec::TargetCl(_) => 0.0,
    };
    let alpha_rad = alpha_deg.to_radians();

    let panel_x: Vec<f64> = coords.iter().map(|(x, _)| *x).collect();
    let panel_y: Vec<f64> = coords.iter().map(|(_, y)| *y).collect();
    let panel_s = compute_arc_lengths(&panel_x, &panel_y);
    let (qinvu_0, qinvu_90) = factorized.surface_qinvu_basis();
    let gamu_0 = factorized.gamu_0.clone();
    let gamu_90 = factorized.gamu_90.clone();

    let operating_mode = match spec {
        AlphaSpec::AlphaDeg(_) => OperatingMode::PrescribedAlpha,
        AlphaSpec::TargetCl(target_cl) => OperatingMode::PrescribedCl { target_cl },
    };

    let mut state = XfoilState::new(
        name.to_string(),
        alpha_rad,
        operating_mode,
        panel_x,
        panel_y,
        panel_s,
        qinvu_0,
        qinvu_90,
        gamu_0,
        gamu_90,
        factorized.build_dij_with_default_wake()?,
    );
    state.panel_xp = factorized.geometry().xp.clone();
    state.panel_yp = factorized.geometry().yp.clone();
    state.xle = factorized.geometry().xle;
    state.yle = factorized.geometry().yle;
    state.xte = factorized.geometry().xte;
    state.yte = factorized.geometry().yte;
    state.chord = factorized.geometry().chord;
    state.sharp = factorized.geometry().sharp;
    state.ante = factorized.geometry().ante;

    Ok((state, factorized))
}

pub fn solve_coords_oper_point(
    name: &str,
    coords: &[(f64, f64)],
    spec: AlphaSpec,
    options: &XfoilOptions,
) -> Result<XfoilViscousResult> {
    let (mut state, factorized) = build_state_from_coords(name, coords, spec, options)?;
    solve_operating_point_from_state(&mut state, &factorized, options)
}

pub fn solve_operating_point_from_state(
    state: &mut XfoilState,
    factorized: &rustfoil_inviscid::FactorizedSystem,
    options: &XfoilOptions,
) -> Result<XfoilViscousResult> {
    specal(state, state.alpha_rad);
    if !state.lwake || (state.alpha_rad - state.awake).abs() > 1.0e-5 {
        xywake(state, factorized, options.wake_length_chords);
    }
    qwcalc(state, factorized);
    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        let (cl_inv, cm_inv) = compute_panel_forces_from_gamma(
            &state.panel_x,
            &state.panel_y,
            &state.qinv,
            state.alpha_rad,
        );
        eprintln!(
            "[CL DEBUG] initial inviscid cl={:.6e} cm={:.6e}",
            cl_inv, cm_inv
        );
    }
    stfind(state);
    iblpan(state);
    xicalc(state);
    uicalc(state);
    uedginit(state);
    if !state.lwdij || !state.ladij {
        qdcalc(state, factorized)?;
    }
    blpini(state, options.reynolds);

    if is_debug_active() {
        add_event(DebugEvent::viscal(
            0,
            state.alpha_rad,
            options.reynolds,
            options.mach,
            options.ncrit,
        ));
        emit_full_state(state, 0, options.mach, options.reynolds);
        add_event(DebugEvent::full_gamma_iter(0, state.gam.clone()));
    }

    for iter in 1..=options.max_iterations {
        let mut assembly = setbl(state, options.reynolds, options.ncrit, options.mach, iter);
        let solve = blsolv(state, &mut assembly, iter);
        update(
            state,
            &assembly,
            &solve,
            options.mach,
            options.reynolds,
            iter,
        );
        if let OperatingMode::PrescribedCl { .. } = state.operating_mode {
            specal(state, state.alpha_rad);
            uicalc(state);
        }
        if is_debug_active() {
            // Match XFOIL: dump BL state after UPDATE, before QVFUE/GAMQV/STMOVE.
            emit_full_state(state, iter, options.mach, options.reynolds);
        }
        qvfue(state);
        gamqv(state);
        if is_debug_active() {
            // Match XFOIL: dump panel gamma after GAMQV, before STMOVE.
            add_event(DebugEvent::full_gamma_iter(iter, state.gam.clone()));
        }
        if std::env::var("RUSTFOIL_DISABLE_STMOVE").is_err() {
            stmove(state);
        }
        update_force_state(state, options.mach, options.reynolds);
        state.iterations = iter;
        state.residual = solve.rms;
        state.converged = solve.rms <= options.tolerance;

        if is_debug_active() {
            add_event(DebugEvent::viscal_result(
                iter,
                solve.rms,
                solve.max,
                state.cl,
                state.cd,
                state.cm,
            ));
        }

        if state.converged {
            break;
        }
    }

    Ok(XfoilViscousResult {
        alpha_deg: state.alpha_rad.to_degrees(),
        cl: state.cl,
        cd: state.cd,
        cm: state.cm,
        x_tr_upper: transition_x(state, crate::state::XfoilSurface::Upper),
        x_tr_lower: transition_x(state, crate::state::XfoilSurface::Lower),
        converged: state.converged,
        iterations: state.iterations,
        residual: state.residual,
        cd_friction: state.cdf,
        cd_pressure: state.cdp,
        x_separation: separation_x(state, crate::state::XfoilSurface::Upper, options.mach, options.reynolds)
            .or_else(|| separation_x(state, crate::state::XfoilSurface::Lower, options.mach, options.reynolds)),
    })
}

fn emit_full_state(state: &XfoilState, iter: usize, mach: f64, reynolds: f64) {
    add_event(DebugEvent::full_bl_state(
        iter,
        surface_state(surface_stations(state, crate::state::XfoilSurface::Upper, mach, reynolds)),
        surface_state(surface_stations(state, crate::state::XfoilSurface::Lower, mach, reynolds)),
    ));
}

fn surface_state(stations: Vec<BlStation>) -> SurfaceBlState {
    SurfaceBlState {
        x: stations.iter().map(|station| station.x).collect(),
        theta: stations.iter().map(|station| station.theta).collect(),
        delta_star: stations.iter().map(|station| station.delta_star).collect(),
        ue: stations.iter().map(|station| station.u).collect(),
        hk: stations.iter().map(|station| station.hk).collect(),
        cf: stations.iter().map(|station| station.cf).collect(),
        mass_defect: stations.iter().map(|station| station.mass_defect).collect(),
    }
}

fn surface_stations(
    state: &XfoilState,
    surface: crate::state::XfoilSurface,
    mach: f64,
    reynolds: f64,
) -> Vec<BlStation> {
    let msq = mach * mach;
    match surface {
        crate::state::XfoilSurface::Upper => {
            let len = state.nbl_upper.min(state.upper_rows.len());
            rows_to_stations(
                &state.upper_rows[..len],
                state.iblte_upper,
                state.itran_upper,
                msq,
                reynolds,
            )
        }
        crate::state::XfoilSurface::Lower => {
            let len = state.nbl_lower.min(state.lower_rows.len());
            rows_to_stations(
                &state.lower_rows[..len],
                state.iblte_lower,
                state.itran_lower,
                msq,
                reynolds,
            )
        }
    }
}

fn transition_x(state: &crate::state::XfoilState, surface: crate::state::XfoilSurface) -> f64 {
    let (itran, iblte, xssitr, rows) = match surface {
        crate::state::XfoilSurface::Upper => (
            state.itran_upper,
            state.iblte_upper,
            state.xssitr_upper,
            &state.upper_rows,
        ),
        crate::state::XfoilSurface::Lower => (
            state.itran_lower,
            state.iblte_lower,
            state.xssitr_lower,
            &state.lower_rows,
        ),
    };
    if itran >= iblte || xssitr <= 0.0 || rows.len() < 2 || state.panel_x.is_empty() {
        return 1.0;
    }
    let surface_limit = iblte.min(rows.len().saturating_sub(1));
    let airfoil_rows = &rows[..=surface_limit];
    let (x, y) = interpolate_surface_point(airfoil_rows, xssitr);
    let chx = state.xte - state.xle;
    let chy = state.yte - state.yle;
    let chsq = chx * chx + chy * chy;
    if !x.is_finite() || !y.is_finite() || !chsq.is_finite() || chsq <= 1.0e-12 {
        1.0
    } else {
        (((x - state.xle) * chx + (y - state.yle) * chy) / chsq).clamp(0.0, 1.0)
    }
}

fn interpolate_surface_point(rows: &[crate::state::XfoilBlRow], target: f64) -> (f64, f64) {
    if rows.is_empty() {
        return (0.0, 0.0);
    }
    if target <= rows[0].x {
        return (rows[0].x_coord, rows[0].y_coord);
    }
    for i in 1..rows.len() {
        if target <= rows[i].x {
            let ds = (rows[i].x - rows[i - 1].x).abs();
            if ds <= 1.0e-12 {
                return (rows[i].x_coord, rows[i].y_coord);
            }
            let w = (target - rows[i - 1].x) / (rows[i].x - rows[i - 1].x);
            return (
                rows[i - 1].x_coord + w * (rows[i].x_coord - rows[i - 1].x_coord),
                rows[i - 1].y_coord + w * (rows[i].y_coord - rows[i - 1].y_coord),
            );
        }
    }
    let last = &rows[rows.len() - 1];
    (last.x_coord, last.y_coord)
}

fn separation_x(
    state: &XfoilState,
    surface: crate::state::XfoilSurface,
    mach: f64,
    reynolds: f64,
) -> Option<f64> {
    surface_stations(state, surface, mach, reynolds)
        .into_iter()
        .find(|station| !station.is_wake && station.cf <= 0.0)
        .map(|station| station.x_coord)
}
