use rustfoil_bl::{
    add_event, is_debug_active, DebugEvent, SurfaceBlState,
};
use rustfoil_core::Body;
use rustfoil_inviscid::InviscidSolver;

use crate::{
    assembly::setbl,
    config::{OperatingMode, XfoilOptions},
    error::{Result, XfoilError},
    forces::update_force_state,
    march::{blpini, mrchdu, mrchue},
    result::XfoilViscousResult,
    solve::blsolv,
    state::XfoilState,
    state_ops::{compute_arc_lengths, gamqv, iblpan, qvfue, specal, stfind, stmove, uedginit, uicalc, ueset, xicalc},
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
    stfind(state);
    iblpan(state);
    xicalc(state);
    uicalc(state);
    uedginit(state);
    if !state.lwdij || !state.ladij {
        qdcalc(state, factorized)?;
    }
    blpini(state, options.reynolds);
    mrchue(state, options.reynolds);
    mrchdu(state, options.reynolds);
    ueset(state);
    qvfue(state);
    gamqv(state);
    update_force_state(state);

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
        let assembly = setbl(state);
        let solve = blsolv(&assembly);
        update(state, &solve);
        if let OperatingMode::PrescribedCl { .. } = state.operating_mode {
            specal(state, state.alpha_rad);
            uicalc(state);
        }
        ueset(state);
        qvfue(state);
        gamqv(state);
        stmove(state);
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
        x_tr_upper: transition_x(&state.upper_rows),
        x_tr_lower: transition_x(&state.lower_rows),
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

fn transition_x(rows: &[crate::state::XfoilBlRow]) -> f64 {
    rows.iter()
        .find(|row| row.is_turbulent && !row.is_wake)
        .map(|row| row.x_coord)
        .unwrap_or(1.0)
}

fn separation_x(rows: &[crate::state::XfoilBlRow]) -> Option<f64> {
    rows.iter()
        .find(|row| !row.is_wake && row.cf <= 0.0)
        .map(|row| row.x_coord)
}
