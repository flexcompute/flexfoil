use rustfoil_bl::{
    add_event, is_debug_active, DebugEvent, SurfaceBlState,
};
use rustfoil_core::Body;
use rustfoil_inviscid::InviscidSolver;

use crate::{
    assembly::setbl,
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

#[derive(Debug, Clone, Copy)]
pub enum AlphaSpec {
    AlphaDeg(f64),
    TargetCl(f64),
}

pub fn solve_body_oper_point(
    body: &Body,
    spec: AlphaSpec,
    options: &XfoilOptions,
) -> Result<XfoilViscousResult> {
    let panels = body.panels();
    let mut node_x: Vec<f64> = panels.iter().map(|panel| panel.p1.x).collect();
    let mut node_y: Vec<f64> = panels.iter().map(|panel| panel.p1.y).collect();
    if let Some(last) = panels.last() {
        if (last.p2.x - panels[0].p1.x).abs() > 1.0e-10 || (last.p2.y - panels[0].p1.y).abs() > 1.0e-10 {
            node_x.push(last.p2.x);
            node_y.push(last.p2.y);
        }
    }
    let coords: Vec<(f64, f64)> = node_x.iter().zip(node_y.iter()).map(|(&x, &y)| (x, y)).collect();
    solve_coords_oper_point(&body.name, &coords, spec, options)
}

pub fn solve_coords_oper_point(
    name: &str,
    coords: &[(f64, f64)],
    spec: AlphaSpec,
    options: &XfoilOptions,
) -> Result<XfoilViscousResult> {
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
    state.sharp = factorized.geometry().sharp;
    state.ante = factorized.geometry().ante;

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
        emit_full_state(state, 0);
        add_event(DebugEvent::full_gamma_iter(0, state.gam.clone()));
    }

    for iter in 1..=options.max_iterations {
        let mut assembly = setbl(state, options.reynolds, options.ncrit, options.mach, iter);
        let solve = blsolv(state, &mut assembly, iter);
        update(state, &assembly, &solve, options.mach, options.reynolds);
        if let OperatingMode::PrescribedCl { .. } = state.operating_mode {
            specal(state, state.alpha_rad);
            uicalc(state);
        }
        qvfue(state);
        gamqv(state);
        if std::env::var("RUSTFOIL_DISABLE_STMOVE").is_err() {
            stmove(state);
            qvfue(state);
            gamqv(state);
        }
        update_force_state(state);
        state.iterations = iter;
        state.residual = solve.rms;
        state.converged = solve.rms <= options.tolerance;

        if is_debug_active() {
            emit_full_state(state, iter);
            add_event(DebugEvent::full_gamma_iter(iter, state.gam.clone()));
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
        x_separation: separation_x(&state.upper_rows).or_else(|| separation_x(&state.lower_rows)),
    })
}

fn emit_full_state(state: &XfoilState, iter: usize) {
    add_event(DebugEvent::full_bl_state(
        iter,
        surface_state(&state.upper_rows),
        surface_state(&state.lower_rows),
    ));
}

fn surface_state(rows: &[crate::state::XfoilBlRow]) -> SurfaceBlState {
    SurfaceBlState {
        x: rows.iter().map(|row| row.x).collect(),
        theta: rows.iter().map(|row| row.theta).collect(),
        delta_star: rows.iter().map(|row| row.dstr).collect(),
        ue: rows.iter().map(|row| row.uedg).collect(),
        hk: rows.iter().map(|row| row.hk).collect(),
        cf: rows.iter().map(|row| row.cf).collect(),
        mass_defect: rows.iter().map(|row| row.mass).collect(),
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
    let x = interpolate_surface_x(airfoil_rows, xssitr);
    let x_min = state.panel_x.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = state.panel_x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let chord = (x_max - x_min).abs();
    if !x.is_finite() || !chord.is_finite() || chord <= 1.0e-12 {
        1.0
    } else {
        ((x - x_min) / chord).clamp(0.0, 1.0)
    }
}

fn interpolate_surface_x(rows: &[crate::state::XfoilBlRow], target: f64) -> f64 {
    if rows.is_empty() {
        return 0.0;
    }
    if target <= rows[0].x {
        return rows[0].x_coord;
    }
    for i in 1..rows.len() {
        if target <= rows[i].x {
            let ds = (rows[i].x - rows[i - 1].x).abs();
            if ds <= 1.0e-12 {
                return rows[i].x_coord;
            }
            let w = (target - rows[i - 1].x) / (rows[i].x - rows[i - 1].x);
            return rows[i - 1].x_coord + w * (rows[i].x_coord - rows[i - 1].x_coord);
        }
    }
    rows[rows.len() - 1].x_coord
}

fn separation_x(rows: &[crate::state::XfoilBlRow]) -> Option<f64> {
    rows.iter()
        .find(|row| !row.is_wake && row.cf <= 0.0)
        .map(|row| row.x_coord)
}
