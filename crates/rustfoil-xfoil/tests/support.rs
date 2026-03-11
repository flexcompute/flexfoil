use std::path::PathBuf;
use std::sync::OnceLock;

use nalgebra::DMatrix;
use rustfoil_testkit::fortran_runner::{
    compile_driver, run_fortran_json, run_fortran_json_with_args, FortranDriverSpec,
    XFOIL_STATE_OBJS, XFOIL_WORKFLOW_OBJS,
};
use rustfoil_xfoil::{
    config::{OperatingMode, XfoilOptions},
    oper::{solve_coords_oper_point, AlphaSpec},
    state::XfoilState,
    state_ops::{compute_arc_lengths, iblpan, qvfue, specal, stfind, uicalc, xicalc},
    wake_panel::{qdcalc, qwcalc, xywake},
};
use serde::Deserialize;

pub const MICROBENCH_MAX_RATIO: f64 = 1.15;
pub const WORKFLOW_MAX_RATIO: f64 = 1.30;

#[derive(Debug, Deserialize)]
pub struct FortranPerfCase {
    pub inner_loops: usize,
    pub samples_seconds: Vec<f64>,
    pub median_seconds: f64,
}

#[derive(Debug, Deserialize)]
pub struct FortranStateFind {
    pub ist: usize,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,
}

#[derive(Debug, Deserialize)]
pub struct FortranIblpan {
    pub nbl_upper: usize,
    pub nbl_lower: usize,
    pub iblte_upper: usize,
    pub iblte_lower: usize,
    pub ipan_upper: Vec<usize>,
    pub ipan_lower: Vec<usize>,
}

#[derive(Debug, Deserialize)]
pub struct FortranXiCalc {
    pub upper_x: Vec<f64>,
    pub lower_x: Vec<f64>,
    pub wgap: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranUiCalc {
    pub upper_uinv: Vec<f64>,
    pub lower_uinv: Vec<f64>,
    pub upper_uinv_a: Vec<f64>,
    pub lower_uinv_a: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranQvfue {
    pub qvis: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranGamqv {
    pub gam: Vec<f64>,
    pub gam_a: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranUeSet {
    pub upper_uedg: Vec<f64>,
    pub lower_uedg: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranDsSet {
    pub upper_dstr: Vec<f64>,
    pub lower_dstr: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranStMove {
    pub ist: usize,
    pub nbl_upper: usize,
    pub nbl_lower: usize,
    pub ipan_upper: Vec<usize>,
    pub ipan_lower: Vec<usize>,
    pub upper_x: Vec<f64>,
    pub lower_x: Vec<f64>,
    pub upper_theta: Vec<f64>,
    pub lower_theta: Vec<f64>,
    pub upper_dstr: Vec<f64>,
    pub lower_dstr: Vec<f64>,
    pub upper_uedg: Vec<f64>,
    pub lower_uedg: Vec<f64>,
    pub upper_mass: Vec<f64>,
    pub lower_mass: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranStatePerf {
    pub stfind: FortranPerfCase,
    pub iblpan: FortranPerfCase,
    pub xicalc: FortranPerfCase,
    pub uicalc: FortranPerfCase,
    pub qvfue: FortranPerfCase,
    pub gamqv: FortranPerfCase,
    pub ueset: FortranPerfCase,
    pub dsset: FortranPerfCase,
}

#[derive(Debug, Deserialize)]
pub struct FortranStateTopologyOutput {
    pub stfind: FortranStateFind,
    pub iblpan: FortranIblpan,
    pub xicalc: FortranXiCalc,
    pub uicalc: FortranUiCalc,
    pub qvfue: FortranQvfue,
    pub gamqv: FortranGamqv,
    pub ueset: FortranUeSet,
    pub dsset: FortranDsSet,
    pub stmove: FortranStMove,
    pub perf: FortranStatePerf,
}

#[derive(Debug, Deserialize)]
pub struct FortranQdcalcOutput {
    pub wake_x: Vec<f64>,
    pub wake_y: Vec<f64>,
    pub diag_sample: Vec<f64>,
    pub row0_sample: Vec<f64>,
    pub matrix_size: usize,
    pub perf: FortranPerfCase,
}

#[derive(Debug, Deserialize)]
pub struct WorkflowCase {
    pub alpha_deg: f64,
    pub converged: bool,
    pub iterations: usize,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
}

#[derive(Debug, Deserialize)]
pub struct WorkflowReference {
    pub cases: Vec<WorkflowCase>,
}

pub fn state_topology_fortran() -> &'static FortranStateTopologyOutput {
    static OUTPUT: OnceLock<FortranStateTopologyOutput> = OnceLock::new();
    OUTPUT.get_or_init(|| {
        let exe = compile_driver(&FortranDriverSpec {
            driver_source: &fortran_driver_path("state_topology_driver.f90"),
            executable_name: "state_topology_driver",
            object_names: XFOIL_STATE_OBJS,
        })
        .expect("compile state topology driver");
        run_fortran_json(exe).expect("run state topology driver")
    })
}

pub fn qdcalc_fortran() -> &'static FortranQdcalcOutput {
    static OUTPUT: OnceLock<FortranQdcalcOutput> = OnceLock::new();
    OUTPUT.get_or_init(|| {
        let exe = compile_driver(&FortranDriverSpec {
            driver_source: &fortran_driver_path("qdcalc_driver.f90"),
            executable_name: "qdcalc_driver",
            object_names: XFOIL_WORKFLOW_OBJS,
        })
        .expect("compile qdcalc driver");
        run_fortran_json(exe).expect("run qdcalc driver")
    })
}

pub fn workflow_fortran(alpha_deg: f64) -> WorkflowReference {
    let exe = workflow_driver_exe();
    run_fortran_json_with_args(exe, [format!("{alpha_deg:.3}")], None).expect("run VISCAL driver")
}

pub fn build_state_topology_fixture() -> XfoilState {
    let panel_x = vec![1.0, 0.75, 0.40, 0.0, 0.35, 0.85];
    let panel_y = vec![0.0, 0.08, 0.12, 0.0, -0.07, -0.02];
    let panel_s = compute_arc_lengths(&panel_x, &panel_y);
    let qinv = vec![0.72, 0.38, 0.12, -0.20, -0.55, -0.25];
    let qinv_a = vec![0.10, 0.08, 0.03, -0.04, -0.07, -0.02];
    let wake_qinv = vec![0.22, 0.18];
    let wake_qinv_a = vec![0.01, 0.01];
    let total_nodes = panel_x.len() + wake_qinv.len();
    let dij = DMatrix::from_fn(total_nodes, total_nodes, |row, col| {
        5.0e-4 * (row as f64 + 1.0) - 3.0e-4 * (col as f64 + 1.0)
    });

    let mut state = XfoilState::new(
        "synthetic".to_string(),
        0.0,
        OperatingMode::PrescribedAlpha,
        panel_x,
        panel_y,
        panel_s,
        qinv.clone(),
        qinv_a.clone(),
        dij,
    );
    state.qinv = qinv.clone();
    state.qinv_a = qinv_a;
    state.gam = qinv;
    state.gam_a = state.qinv_a.clone();
    state.panel_xp = vec![1.0, 0.95, 0.80, 0.0, 0.80, 1.0];
    state.panel_yp = vec![0.10, 0.09, 0.05, 0.0, -0.05, -0.10];
    state.sharp = true;
    state.ante = 0.02;
    state.wake_x = vec![1.05, 1.25];
    state.wake_y = vec![0.0, 0.0];
    state.wake_s = vec![0.0, 0.20];
    state.wake_qinv = wake_qinv;
    state.wake_qinv_a = wake_qinv_a;
    state.qvis = state
        .qinv
        .iter()
        .copied()
        .chain(state.wake_qinv.iter().copied())
        .collect();
    state.wgap = vec![0.0; state.wake_x.len()];
    state
}

pub fn prepare_state_topology_state() -> XfoilState {
    let mut state = build_state_topology_fixture();
    stfind(&mut state);
    iblpan(&mut state);
    xicalc(&mut state);
    uicalc(&mut state);
    seed_state_rows(&mut state);
    state
}

pub fn shifted_stmove_state() -> XfoilState {
    let mut state = prepare_state_topology_state();
    state.gam = vec![0.72, 0.28, -0.05, -0.25, -0.40, -0.20];
    state
}

pub fn build_qdcalc_state() -> (XfoilState, rustfoil_inviscid::FactorizedSystem) {
    let coords = naca0012_coords(40);
    let solver = rustfoil_inviscid::InviscidSolver::new();
    let factorized = solver.factorize(&coords).expect("factorize qdcalc fixture");
    let alpha_rad = 4.0_f64.to_radians();
    let panel_x: Vec<f64> = coords.iter().map(|(x, _)| *x).collect();
    let panel_y: Vec<f64> = coords.iter().map(|(_, y)| *y).collect();
    let panel_s = compute_arc_lengths(&panel_x, &panel_y);
    let mut state = XfoilState::new(
        "NACA0012".to_string(),
        alpha_rad,
        OperatingMode::PrescribedAlpha,
        panel_x,
        panel_y,
        panel_s,
        factorized.gamu_0.clone(),
        factorized.gamu_90.clone(),
        factorized
            .build_dij_with_default_wake()
            .expect("default wake dij"),
    );
    xywake(&mut state, &factorized, 1.0);
    qwcalc(&mut state, &factorized);
    specal(&mut state, alpha_rad);
    (state, factorized)
}

pub fn rust_qdcalc_sample() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, usize) {
    let (mut state, factorized) = build_qdcalc_state();
    qdcalc(&mut state, &factorized).expect("qdcalc");
    let diag_sample: Vec<f64> = (0..20.min(state.dij.nrows())).map(|i| state.dij[(i, i)]).collect();
    let row0_sample: Vec<f64> = (0..20.min(state.dij.ncols())).map(|j| state.dij[(0, j)]).collect();
    (
        state.wake_x.clone(),
        state.wake_y.clone(),
        diag_sample,
        row0_sample,
        state.dij.nrows(),
    )
}

pub fn run_workflow_case(alpha_deg: f64) -> rustfoil_xfoil::result::XfoilViscousResult {
    solve_coords_oper_point(
        "NACA0012",
        &naca0012_coords(160),
        AlphaSpec::AlphaDeg(alpha_deg),
        &XfoilOptions {
            reynolds: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 10,
            tolerance: 1.0e-4,
            wake_length_chords: 1.0,
            operating_mode: OperatingMode::PrescribedAlpha,
        },
    )
    .expect("run workflow case")
}

pub fn naca0012_coords(npanel: usize) -> Vec<(f64, f64)> {
    let nhalf = npanel / 2;
    let mut coords = Vec::with_capacity(npanel + 1);
    let pi = std::f64::consts::PI;
    let t = 0.12;

    for i in 0..=nhalf {
        let beta = pi * i as f64 / nhalf as f64;
        let xx = 0.5 * (1.0 - beta.cos());
        let yt = 5.0
            * t
            * (0.2969 * xx.sqrt()
                - 0.126 * xx
                - 0.3516 * xx * xx
                + 0.2843 * xx * xx * xx
                - 0.1036 * xx * xx * xx * xx);
        coords.push((xx, yt));
    }

    for i in 1..=nhalf {
        let beta = pi * (nhalf - i) as f64 / nhalf as f64;
        let xx = 0.5 * (1.0 - beta.cos());
        let yt = 5.0
            * t
            * (0.2969 * xx.sqrt()
                - 0.126 * xx
                - 0.3516 * xx * xx
                + 0.2843 * xx * xx * xx
                - 0.1036 * xx * xx * xx * xx);
        coords.push((xx, -yt));
    }

    coords
}

pub fn assert_close_scalar(label: &str, lhs: f64, rhs: f64, tol: f64) {
    let diff = (lhs - rhs).abs();
    assert!(
        diff <= tol,
        "{} mismatch: left {:.12e}, right {:.12e}, abs diff {:.12e}, tol {:.12e}",
        label,
        lhs,
        rhs,
        diff,
        tol
    );
}

pub fn assert_close_slice(label: &str, lhs: &[f64], rhs: &[f64], tol: f64) {
    assert_eq!(
        lhs.len(),
        rhs.len(),
        "{} length mismatch: {} vs {}",
        label,
        lhs.len(),
        rhs.len()
    );
    for (idx, (&left, &right)) in lhs.iter().zip(rhs.iter()).enumerate() {
        let diff = (left - right).abs();
        assert!(
            diff <= tol,
            "{}[{}] mismatch: left {:.12e}, right {:.12e}, abs diff {:.12e}, tol {:.12e}",
            label,
            idx,
            left,
            right,
            diff,
            tol
        );
    }
}

fn seed_state_rows(state: &mut XfoilState) {
    for (surface_idx, rows) in [&mut state.upper_rows, &mut state.lower_rows].into_iter().enumerate() {
        for (offset, row) in rows.iter_mut().enumerate().skip(1) {
            let ibl = offset as f64 + 1.0;
            let is = surface_idx as f64 + 1.0;
            row.theta = 0.008 + 0.002 * ibl + 0.0005 * is;
            row.dstr = 0.012 + 0.003 * ibl + 0.0008 * is;
            row.mass = 0.025 + 0.004 * ibl + 0.001 * is;
            row.uedg = row.uinv + 0.03 * ibl;
        }
    }
}

fn fortran_driver_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fortran")
        .join(name)
}

fn workflow_driver_exe() -> &'static std::path::PathBuf {
    static EXE: OnceLock<std::path::PathBuf> = OnceLock::new();
    EXE.get_or_init(|| {
        compile_driver(&FortranDriverSpec {
            driver_source: &fortran_driver_path("viscal_driver.f90"),
            executable_name: "viscal_driver",
            object_names: XFOIL_WORKFLOW_OBJS,
        })
        .expect("compile VISCAL driver")
    })
}
