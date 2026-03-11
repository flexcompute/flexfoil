use rustfoil_bl::{add_event, is_debug_active, DebugEvent};

use crate::state::{XfoilBlRow, XfoilState};

pub fn blpini(state: &mut XfoilState, reynolds: f64) {
    initialize_surface(state.upper_rows.as_mut_slice(), reynolds, 1);
    initialize_surface(state.lower_rows.as_mut_slice(), reynolds, 2);
}

pub fn mrchue(state: &mut XfoilState, reynolds: f64) {
    march_surface(state.upper_rows.as_mut_slice(), reynolds, 1);
    march_surface(state.lower_rows.as_mut_slice(), reynolds, 2);
}

pub fn mrchdu(state: &mut XfoilState, reynolds: f64) {
    march_surface(state.upper_rows.as_mut_slice(), reynolds, 1);
    march_surface(state.lower_rows.as_mut_slice(), reynolds, 2);
}

fn initialize_surface(rows: &mut [XfoilBlRow], reynolds: f64, side: usize) {
    for row in rows.iter_mut().skip(1) {
        let x_eff = row.x.max(1.0e-6);
        let ue_eff = row.uedg.abs().max(0.02);
        row.theta = (0.45 * x_eff / (6.0 * ue_eff * reynolds)).sqrt().max(1.0e-8);
        row.dstr = 2.2 * row.theta;
        row.mass = row.uedg * row.dstr;
        row.h = 2.2;
        row.hk = 2.2;
        row.r_theta = reynolds * row.theta * ue_eff;
        row.cf = if row.is_wake { 0.0 } else { 0.664 / row.r_theta.sqrt().max(10.0) };
        row.cd = row.cf * 0.5;
        row.is_laminar = !row.is_wake;
        row.is_turbulent = row.is_wake;
    }

    if is_debug_active() {
        add_event(DebugEvent::bl_init(
            side,
            rows.len(),
            rows.iter().map(|row| row.x).collect(),
            rows.iter().map(|row| row.uedg).collect(),
            rows.iter().map(|row| row.theta).collect(),
            rows.iter().map(|row| row.dstr).collect(),
        ));
    }
}

fn march_surface(rows: &mut [XfoilBlRow], reynolds: f64, side: usize) {
    for i in 2..rows.len() {
        let prev = rows[i - 1].clone();
        let row = &mut rows[i];
        let dx = (row.x - prev.x).abs().max(1.0e-6);
        let ue = row.uedg.abs().max(0.02);
        row.theta = (prev.theta + 0.036 * dx / (reynolds * ue).sqrt()).max(prev.theta);
        row.dstr = (2.2 + 0.08 * i as f64) * row.theta;
        row.mass = row.uedg * row.dstr;
        row.h = row.dstr / row.theta.max(1.0e-12);
        row.hk = row.h;
        row.r_theta = reynolds * row.theta * ue;
        row.cf = if row.is_wake { 0.0 } else { 0.664 / row.r_theta.sqrt().max(10.0) };
        row.cd = row.cf * 0.5;
        if !row.is_wake && row.h > 2.8 {
            row.is_laminar = false;
            row.is_turbulent = true;
            row.ampl = 9.0;
            row.ctau = 0.03;
        }
    }

    if is_debug_active() {
        let _ = side;
    }
}
