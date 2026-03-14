mod support;

use rustfoil_xfoil::state_ops::{dsset, gamqv, iblpan, qvfue, stfind, stmove, uedginit, uicalc, ueset, xicalc};

use support::{
    assert_close_scalar, assert_close_slice, build_state_topology_fixture, prepare_state_topology_state,
    shifted_stmove_state, state_topology_fortran,
};

#[test]
fn stfind_matches_fortran() {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    let reference = &state_topology_fortran().stfind;

    assert_eq!(state.ist, reference.ist);
    assert_close_scalar("sst", state.sst, reference.sst, 1.0e-12);
    assert_close_scalar("sst_go", state.sst_go, reference.sst_go, 1.0e-12);
    assert_close_scalar("sst_gp", state.sst_gp, reference.sst_gp, 1.0e-12);
}

#[test]
fn iblpan_matches_fortran() {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    iblpan(&mut state);
    let reference = &state_topology_fortran().iblpan;

    assert_eq!(state.ipan_upper.len(), reference.nbl_upper);
    assert_eq!(state.nbl_lower, reference.nbl_lower);
    assert_eq!(state.iblte_upper + 1, reference.iblte_upper);
    assert_eq!(state.iblte_lower + 1, reference.iblte_lower);
    assert_eq!(&state.ipan_upper[1..], reference.ipan_upper.as_slice());
    assert_eq!(&state.ipan_lower[1..], reference.ipan_lower.as_slice());
}

#[test]
fn xicalc_matches_fortran() {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    iblpan(&mut state);
    xicalc(&mut state);
    let reference = &state_topology_fortran().xicalc;

    let upper_x: Vec<f64> = state.upper_rows.iter().map(|row| row.x).collect();
    let lower_x: Vec<f64> = state.lower_rows.iter().map(|row| row.x).collect();
    assert_close_slice(
        "upper_x_airfoil",
        &upper_x[..=state.iblte_upper],
        &reference.upper_x[..=state.iblte_upper],
        1.0e-10,
    );
    assert_close_slice("lower_x", &lower_x, &reference.lower_x, 1.0e-10);
    assert_close_slice("wgap", &state.wgap, &reference.wgap, 1.0e-12);
}

#[test]
fn uicalc_matches_fortran() {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    iblpan(&mut state);
    xicalc(&mut state);
    uicalc(&mut state);
    let reference = &state_topology_fortran().uicalc;

    let upper_uinv: Vec<f64> = state.upper_rows.iter().map(|row| row.uinv).collect();
    let lower_uinv: Vec<f64> = state.lower_rows.iter().map(|row| row.uinv).collect();
    let upper_uinv_a: Vec<f64> = state.upper_rows.iter().map(|row| row.uinv_a).collect();
    let lower_uinv_a: Vec<f64> = state.lower_rows.iter().map(|row| row.uinv_a).collect();
    assert_close_slice(
        "upper_uinv_airfoil",
        &upper_uinv[..=state.iblte_upper],
        &reference.upper_uinv[..=state.iblte_upper],
        1.0e-12,
    );
    assert_close_slice("lower_uinv", &lower_uinv, &reference.lower_uinv, 1.0e-12);
    assert_close_slice(
        "upper_uinv_a_airfoil",
        &upper_uinv_a[..=state.iblte_upper],
        &reference.upper_uinv_a[..=state.iblte_upper],
        1.0e-12,
    );
    assert_close_slice("lower_uinv_a", &lower_uinv_a, &reference.lower_uinv_a, 1.0e-12);
}

#[test]
fn uicalc_preserves_existing_uedg_until_explicit_init() {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    iblpan(&mut state);
    xicalc(&mut state);

    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.uedg = 123.456;
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.uedg = -234.567;
    }

    uicalc(&mut state);

    for row in state.upper_rows.iter().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        assert_close_scalar("upper_uedg_preserved", row.uedg, 123.456, 1.0e-12);
    }
    for row in state.lower_rows.iter().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        assert_close_scalar("lower_uedg_preserved", row.uedg, -234.567, 1.0e-12);
    }

    uedginit(&mut state);

    for row in state.upper_rows.iter().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        assert_close_scalar("upper_uedg_init", row.uedg, row.uinv, 1.0e-12);
    }
    for row in state.lower_rows.iter().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        assert_close_scalar("lower_uedg_init", row.uedg, row.uinv, 1.0e-12);
    }
}

#[test]
fn qvfue_matches_fortran() {
    let mut state = prepare_state_topology_state();
    qvfue(&mut state);
    let reference = &state_topology_fortran().qvfue;
    assert_close_slice("qvis", &state.qvis, &reference.qvis, 1.0e-12);
}

#[test]
fn gamqv_matches_fortran() {
    let mut state = prepare_state_topology_state();
    qvfue(&mut state);
    gamqv(&mut state);
    let reference = &state_topology_fortran().gamqv;
    assert_close_slice("gam", &state.gam, &reference.gam, 1.0e-12);
    assert_close_slice("gam_a", &state.gam_a, &reference.gam_a, 1.0e-12);
}

#[test]
fn ueset_matches_fortran() {
    let mut state = prepare_state_topology_state();
    ueset(&mut state);
    let reference = &state_topology_fortran().ueset;
    let upper_uedg: Vec<f64> = state.upper_rows.iter().map(|row| row.uedg).collect();
    let lower_uedg: Vec<f64> = state.lower_rows.iter().map(|row| row.uedg).collect();
    assert_close_slice(
        "upper_uedg_airfoil",
        &upper_uedg[1..=state.iblte_upper],
        &reference.upper_uedg[1..=state.iblte_upper],
        1.0e-12,
    );
    assert_close_slice("lower_uedg", &lower_uedg[1..], &reference.lower_uedg[1..], 1.0e-12);
}

#[test]
fn dsset_matches_fortran() {
    let mut state = prepare_state_topology_state();
    ueset(&mut state);
    dsset(&mut state);
    let reference = &state_topology_fortran().dsset;
    let upper_dstr: Vec<f64> = state.upper_rows.iter().map(|row| row.dstr).collect();
    let lower_dstr: Vec<f64> = state.lower_rows.iter().map(|row| row.dstr).collect();
    assert_close_slice(
        "upper_dstr_airfoil",
        &upper_dstr[1..=state.iblte_upper],
        &reference.upper_dstr[1..=state.iblte_upper],
        1.0e-12,
    );
    assert_close_slice("lower_dstr", &lower_dstr[1..], &reference.lower_dstr[1..], 1.0e-12);
}

#[test]
fn stmove_matches_fortran_where_applicable() {
    let mut state = shifted_stmove_state();
    stmove(&mut state);
    let reference = &state_topology_fortran().stmove;

    assert_eq!(state.ist, reference.ist);
    assert_eq!(state.nbl_upper, reference.nbl_upper);
    assert_eq!(state.nbl_lower, reference.nbl_lower);
    assert_eq!(
        &state.ipan_upper[1..state.nbl_upper],
        reference.ipan_upper.as_slice()
    );
    assert_eq!(&state.ipan_lower[1..], reference.ipan_lower.as_slice());

    let upper_x: Vec<f64> = state.upper_rows.iter().map(|row| row.x).collect();
    let lower_x: Vec<f64> = state.lower_rows.iter().map(|row| row.x).collect();
    let upper_theta: Vec<f64> = state.upper_rows.iter().map(|row| row.theta).collect();
    let lower_theta: Vec<f64> = state.lower_rows.iter().map(|row| row.theta).collect();
    let upper_dstr: Vec<f64> = state.upper_rows.iter().map(|row| row.dstr).collect();
    let lower_dstr: Vec<f64> = state.lower_rows.iter().map(|row| row.dstr).collect();
    let upper_uedg: Vec<f64> = state.upper_rows.iter().map(|row| row.uedg).collect();
    let lower_uedg: Vec<f64> = state.lower_rows.iter().map(|row| row.uedg).collect();
    let upper_mass: Vec<f64> = state.upper_rows.iter().map(|row| row.mass).collect();
    let lower_mass: Vec<f64> = state.lower_rows.iter().map(|row| row.mass).collect();
    assert_close_slice(
        "stmove.upper_x_airfoil",
        &upper_x[..state.nbl_upper],
        &reference.upper_x,
        1.0e-10,
    );
    assert_close_slice("stmove.lower_x", &lower_x, &reference.lower_x, 1.0e-10);
    assert_close_slice(
        "stmove.upper_theta",
        &upper_theta[..state.nbl_upper],
        &reference.upper_theta,
        1.0e-12,
    );
    assert_close_slice("stmove.lower_theta", &lower_theta, &reference.lower_theta, 1.0e-12);
    assert_close_slice(
        "stmove.upper_dstr",
        &upper_dstr[..state.nbl_upper],
        &reference.upper_dstr,
        1.0e-12,
    );
    assert_close_slice("stmove.lower_dstr", &lower_dstr, &reference.lower_dstr, 1.0e-12);
    assert_close_slice(
        "stmove.upper_uedg",
        &upper_uedg[..state.nbl_upper],
        &reference.upper_uedg,
        1.0e-12,
    );
    assert_close_slice("stmove.lower_uedg", &lower_uedg, &reference.lower_uedg, 1.0e-12);
    assert_close_slice(
        "stmove.upper_mass",
        &upper_mass[..state.nbl_upper],
        &reference.upper_mass,
        1.0e-12,
    );
    assert_close_slice("stmove.lower_mass", &lower_mass, &reference.lower_mass, 1.0e-12);
}

#[test]
fn stmove_keeps_projected_views_in_sync() {
    let mut state = shifted_stmove_state();
    stmove(&mut state);

    let qvis_after = state.qvis.clone();
    let gam_after = state.gam.clone();
    let gam_a_after = state.gam_a.clone();

    qvfue(&mut state);
    gamqv(&mut state);

    assert_close_slice("stmove.qvis_sync", &qvis_after, &state.qvis, 1.0e-12);
    assert_close_slice("stmove.gam_sync", &gam_after, &state.gam, 1.0e-12);
    assert_close_slice("stmove.gam_a_sync", &gam_a_after, &state.gam_a, 1.0e-12);
}
