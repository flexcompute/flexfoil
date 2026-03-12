use std::path::PathBuf;
use std::sync::OnceLock;

use nalgebra::DMatrix;
use rustfoil_bl::{finalize_debug, init_debug};
use rustfoil_testkit::fortran_runner::{
    compile_driver, run_fortran_json, run_fortran_json_with_args, run_fortran_test,
    run_xfoil_instrumented, FortranDriverSpec, XFOIL_STATE_OBJS, XFOIL_WORKFLOW_OBJS,
};
use rustfoil_xfoil::{
    assembly::setbl,
    config::{OperatingMode, XfoilOptions},
    march::{blpini, mrchdu, mrchue},
    oper::{solve_coords_oper_point, AlphaSpec},
    solve::blsolv,
    state::XfoilState,
    state_ops::{compute_arc_lengths, iblpan, specal, stfind, uedginit, uicalc, xicalc},
    update::update,
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
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    pub chord: f64,
    pub wake_x: Vec<f64>,
    pub wake_y: Vec<f64>,
    pub wake_s: Vec<f64>,
    pub wake_nx: Vec<f64>,
    pub wake_ny: Vec<f64>,
    pub wake_apanel: Vec<f64>,
    pub wake_qinvu_0: Vec<f64>,
    pub wake_qinvu_90: Vec<f64>,
    pub wake_qinv: Vec<f64>,
    pub wake_qinv_a: Vec<f64>,
    pub dij_flat: Vec<f64>,
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

#[derive(Debug, Deserialize, Clone)]
pub struct FortranMarchState {
    pub nbl_upper: usize,
    pub nbl_lower: usize,
    pub iblte_upper: usize,
    pub iblte_lower: usize,
    pub itran_upper: usize,
    pub itran_lower: usize,
    pub xssitr_upper: f64,
    pub xssitr_lower: f64,
    pub lblini: bool,
    pub upper_x: Vec<f64>,
    pub lower_x: Vec<f64>,
    pub upper_theta: Vec<f64>,
    pub lower_theta: Vec<f64>,
    pub upper_dstr: Vec<f64>,
    pub lower_dstr: Vec<f64>,
    pub upper_uedg: Vec<f64>,
    pub lower_uedg: Vec<f64>,
    pub upper_ctau: Vec<f64>,
    pub lower_ctau: Vec<f64>,
    pub upper_mass: Vec<f64>,
    pub lower_mass: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct FortranMarchOutput {
    pub setbl: FortranMarchState,
    pub perf: FortranMarchPerf,
}

#[derive(Debug, Deserialize)]
pub struct FortranMarchPerf {
    pub setbl: FortranPerfCase,
}

#[derive(Debug, Clone)]
pub struct XfoilBinaryReference {
    pub mrchue: FortranMarchState,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
    pub converged: bool,
    pub iterations: usize,
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
        let stdout = run_fortran_test(exe).expect("run qdcalc driver");
        let json_start = stdout.find('{').expect("qdcalc driver JSON start");
        serde_json::from_str(&stdout[json_start..]).expect("parse qdcalc driver JSON")
    })
}

pub fn workflow_fortran(alpha_deg: f64) -> WorkflowReference {
    let exe = workflow_driver_exe();
    run_fortran_json_with_args(exe, [format!("{alpha_deg:.3}")], None).expect("run VISCAL driver")
}

pub fn march_fortran() -> &'static FortranMarchOutput {
    static OUTPUT: OnceLock<FortranMarchOutput> = OnceLock::new();
    OUTPUT.get_or_init(|| {
        let exe = compile_driver(&FortranDriverSpec {
            driver_source: &fortran_driver_path("setbl_driver.f90"),
            executable_name: "setbl_driver",
            object_names: XFOIL_WORKFLOW_OBJS,
        })
        .expect("compile setbl driver");
        let stdout = run_fortran_test(exe).expect("run setbl driver");
        let json_start = stdout.find('{').expect("setbl driver JSON start");
        serde_json::from_str(&stdout[json_start..]).expect("parse setbl driver JSON")
    })
}

pub fn mrchue_fortran() -> &'static FortranMarchState {
    static OUTPUT: OnceLock<FortranMarchState> = OnceLock::new();
    OUTPUT.get_or_init(|| {
        let exe = compile_driver(&FortranDriverSpec {
            driver_source: &fortran_driver_path("mrchue_driver.f90"),
            executable_name: "mrchue_driver",
            object_names: XFOIL_WORKFLOW_OBJS,
        })
        .expect("compile mrchue driver");
        let stdout = run_fortran_test(exe).expect("run mrchue driver");
        let json_start = stdout.find('{').expect("mrchue driver JSON start");
        let parsed: serde_json::Value =
            serde_json::from_str(&stdout[json_start..]).expect("parse mrchue driver JSON");
        serde_json::from_value(parsed["mrchue"].clone()).expect("extract mrchue payload")
    })
}

pub fn xfoil_binary_reference() -> &'static XfoilBinaryReference {
    static OUTPUT: OnceLock<XfoilBinaryReference> = OnceLock::new();
    OUTPUT.get_or_init(|| build_xfoil_binary_reference(15.0))
}

pub fn xfoil_binary_oper_point(alpha_deg: f64) -> XfoilBinaryReference {
    build_xfoil_binary_reference(alpha_deg)
}

pub fn xfoil_iteration_debug(alpha_deg: f64) -> serde_json::Value {
    let source_coords_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("naca0012_xfoil_paneled.dat");
    let workdir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("target")
        .join(format!("xfoil-debug-{alpha_deg:.1}"));
    std::fs::create_dir_all(&workdir).expect("create xfoil debug workdir");
    let debug_path = workdir.join("xfoil_debug.json");
    let local_coords_path = workdir.join("naca0012_buffer_real.dat");
    let _ = std::fs::remove_file(&debug_path);
    std::fs::copy(&source_coords_path, &local_coords_path)
        .expect("copy xfoil buffer geometry into stage debug workdir");
    let commands = format!(
        "LOAD {}\n\nOPER\nVISC 1000000\nMACH 0\nITER 1\nALFA {alpha_deg:.1}\n\nQUIT\n",
        local_coords_path
            .file_name()
            .expect("local buffer airfoil filename")
            .to_string_lossy(),
    );
    run_xfoil_instrumented(&commands, &workdir).expect("run xfoil_instrumented");
    let debug_text = read_debug_text(&debug_path);
    serde_json::from_str(&debug_text).expect("parse xfoil debug JSON")
}

pub fn xfoil_debug_mrchue_state(alpha_deg: f64) -> FortranMarchState {
    let debug = xfoil_iteration_debug(alpha_deg);
    let events = debug["events"].as_array().expect("xfoil debug events");
    let mut upper = Vec::new();
    let mut lower = Vec::new();
    for event in events {
        if event["subroutine"].as_str() != Some("MRCHUE") {
            continue;
        }
        let side = event["side"].as_u64().expect("MRCHUE side") as usize;
        let ibl = event["ibl"].as_u64().expect("MRCHUE ibl") as usize;
        let tuple = (
            ibl,
            event["x"].as_f64().expect("MRCHUE x"),
            event["Ue"].as_f64().expect("MRCHUE Ue"),
            event["theta"].as_f64().expect("MRCHUE theta"),
            event["delta_star"].as_f64().expect("MRCHUE delta_star"),
            event["ctau"].as_f64().expect("MRCHUE ctau"),
        );
        if side == 1 {
            upper.push(tuple);
        } else if side == 2 {
            lower.push(tuple);
        }
    }
    mrchue_state_from_events(&upper, &lower)
}

pub fn rust_iteration_debug(alpha_deg: f64) -> serde_json::Value {
    let debug_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("target")
        .join(format!("rustfoil-debug-{alpha_deg:.1}.json"));
    let _ = std::fs::remove_file(&debug_path);

    let (mut state, _) = build_march_state(alpha_deg);
    blpini(&mut state, 1.0e6);
    init_debug(&debug_path);
    let assembly = setbl(&mut state, 1.0e6, 9.0, 0.0, 1);
    let mut assembly = assembly;
    let solve = blsolv(&mut state, &mut assembly, 1);
    update(&mut state, &assembly, &solve, 0.0, 1.0e6);
    finalize_debug();

    let debug_text = std::fs::read_to_string(&debug_path).expect("read rust debug JSON");
    serde_json::from_str(&debug_text).expect("parse rust debug JSON")
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
    state.wake_nx = vec![0.0, 0.0];
    state.wake_ny = vec![1.0, 1.0];
    state.wake_apanel = vec![0.0, 0.0];
    state.wake_qinv = wake_qinv;
    state.wake_qinv_a = wake_qinv_a;
    state.wake_qinvu_0 = state.wake_qinv.clone();
    state.wake_qinvu_90 = state.wake_qinv_a.clone();
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
    let (qinvu_0, qinvu_90) = factorized.surface_qinvu_basis();
    let gamu_0 = factorized.gamu_0.clone();
    let gamu_90 = factorized.gamu_90.clone();
    let mut state = XfoilState::new(
        "NACA0012".to_string(),
        alpha_rad,
        OperatingMode::PrescribedAlpha,
        panel_x,
        panel_y,
        panel_s,
        qinvu_0,
        qinvu_90,
        gamu_0,
        gamu_90,
        factorized
            .build_dij_with_default_wake()
            .expect("default wake dij"),
    );
    specal(&mut state, alpha_rad);
    xywake(&mut state, &factorized, 1.0);
    qwcalc(&mut state, &factorized);
    specal(&mut state, alpha_rad);
    (state, factorized)
}

pub fn rust_qdcalc_output() -> FortranQdcalcOutput {
    let (mut state, factorized) = build_qdcalc_state();
    qdcalc(&mut state, &factorized).expect("qdcalc");
    let geom = factorized.geometry();
    FortranQdcalcOutput {
        xle: geom.xle,
        yle: geom.yle,
        xte: geom.xte,
        yte: geom.yte,
        chord: geom.chord,
        wake_x: state.wake_x.clone(),
        wake_y: state.wake_y.clone(),
        wake_s: state.wake_s.clone(),
        wake_nx: state.wake_nx.clone(),
        wake_ny: state.wake_ny.clone(),
        wake_apanel: state.wake_apanel.clone(),
        wake_qinvu_0: state.wake_qinvu_0.clone(),
        wake_qinvu_90: state.wake_qinvu_90.clone(),
        wake_qinv: state.wake_qinv.clone(),
        wake_qinv_a: state.wake_qinv_a.clone(),
        dij_flat: flatten_matrix_row_major(&state.dij),
        matrix_size: state.dij.nrows(),
        perf: FortranPerfCase {
            inner_loops: 0,
            samples_seconds: Vec::new(),
            median_seconds: 0.0,
        },
    }
}

pub fn run_workflow_case(alpha_deg: f64) -> rustfoil_xfoil::result::XfoilViscousResult {
    let coords_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("naca0012.dat");
    let coords = load_dat_coords(&coords_path);
    solve_coords_oper_point(
        "NACA0012",
        &coords,
        AlphaSpec::AlphaDeg(alpha_deg),
        &XfoilOptions {
            reynolds: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 1,
            tolerance: 1.0e-4,
            wake_length_chords: 1.0,
            operating_mode: OperatingMode::PrescribedAlpha,
        },
    )
    .expect("run workflow case")
}

pub fn build_march_state(alpha_deg: f64) -> (XfoilState, rustfoil_inviscid::FactorizedSystem) {
    let paneled_coords_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("naca0012_xfoil_paneled.dat");
    let coords = load_dat_coords(&paneled_coords_path);
    let solver = rustfoil_inviscid::InviscidSolver::new();
    let factorized = solver.factorize(&coords).expect("factorize march fixture");
    let alpha_rad = alpha_deg.to_radians();
    let panel_x: Vec<f64> = coords.iter().map(|(x, _)| *x).collect();
    let panel_y: Vec<f64> = coords.iter().map(|(_, y)| *y).collect();
    let panel_s = compute_arc_lengths(&panel_x, &panel_y);
    let (qinvu_0, qinvu_90) = factorized.surface_qinvu_basis();
    let gamu_0 = factorized.gamu_0.clone();
    let gamu_90 = factorized.gamu_90.clone();

    let mut state = XfoilState::new(
        "NACA0012".to_string(),
        alpha_rad,
        OperatingMode::PrescribedAlpha,
        panel_x,
        panel_y,
        panel_s,
        qinvu_0,
        qinvu_90,
        gamu_0,
        gamu_90,
        factorized
            .build_dij_with_default_wake()
            .expect("default wake dij"),
    );
    state.panel_xp = factorized.geometry().xp.clone();
    state.panel_yp = factorized.geometry().yp.clone();
    state.sharp = factorized.geometry().sharp;
    state.ante = factorized.geometry().ante;

    specal(&mut state, alpha_rad);
    xywake(&mut state, &factorized, 1.0);
    qwcalc(&mut state, &factorized);
    specal(&mut state, alpha_rad);
    stfind(&mut state);
    iblpan(&mut state);
    xicalc(&mut state);
    uicalc(&mut state);
    uedginit(&mut state);
    qdcalc(&mut state, &factorized).expect("qdcalc");
    (state, factorized)
}

pub fn rust_setbl_state(alpha_deg: f64) -> FortranMarchState {
    let (mut state, _) = build_march_state(alpha_deg);
    blpini(&mut state, 1.0e6);
    let _ = setbl(&mut state, 1.0e6, 9.0, 0.0, 1);
    march_state_snapshot(&state)
}

pub fn rust_mrchue_state(alpha_deg: f64) -> FortranMarchState {
    let (mut state, _) = build_march_state(alpha_deg);
    blpini(&mut state, 1.0e6);
    mrchue(&mut state, 1.0e6, 9.0);
    march_state_snapshot(&state)
}

fn march_state_snapshot(state: &XfoilState) -> FortranMarchState {
    let upper_count = state.nbl_upper.min(state.upper_rows.len());
    let lower_count = state.nbl_lower.min(state.lower_rows.len());
    FortranMarchState {
        nbl_upper: upper_count,
        nbl_lower: lower_count,
        iblte_upper: state.iblte_upper + 1,
        iblte_lower: state.iblte_lower + 1,
        itran_upper: state.itran_upper + 1,
        itran_lower: state.itran_lower + 1,
        xssitr_upper: state.xssitr_upper,
        xssitr_lower: state.xssitr_lower,
        lblini: state.lblini,
        upper_x: state.upper_rows.iter().take(upper_count).map(|row| row.x).collect(),
        lower_x: state.lower_rows.iter().take(lower_count).map(|row| row.x).collect(),
        upper_theta: state
            .upper_rows
            .iter()
            .take(upper_count)
            .map(|row| row.theta)
            .collect(),
        lower_theta: state
            .lower_rows
            .iter()
            .take(lower_count)
            .map(|row| row.theta)
            .collect(),
        upper_dstr: state
            .upper_rows
            .iter()
            .take(upper_count)
            .map(|row| row.dstr)
            .collect(),
        lower_dstr: state
            .lower_rows
            .iter()
            .take(lower_count)
            .map(|row| row.dstr)
            .collect(),
        upper_uedg: state
            .upper_rows
            .iter()
            .take(upper_count)
            .map(|row| row.uedg)
            .collect(),
        lower_uedg: state
            .lower_rows
            .iter()
            .take(lower_count)
            .map(|row| row.uedg)
            .collect(),
        upper_ctau: state
            .upper_rows
            .iter()
            .take(upper_count)
            .map(|row| row.ctau)
            .collect(),
        lower_ctau: state
            .lower_rows
            .iter()
            .take(lower_count)
            .map(|row| row.ctau)
            .collect(),
        upper_mass: state
            .upper_rows
            .iter()
            .take(upper_count)
            .map(|row| row.mass)
            .collect(),
        lower_mass: state
            .lower_rows
            .iter()
            .take(lower_count)
            .map(|row| row.mass)
            .collect(),
    }
}

pub fn naca0012_coords(npanel: usize) -> Vec<(f64, f64)> {
    let nhalf = npanel / 2;
    let mut coords = Vec::with_capacity(npanel + 1);
    let pi = std::f64::consts::PI;
    let t = 0.12;

    for i in (0..=nhalf).rev() {
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
        let beta = pi * i as f64 / nhalf as f64;
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

fn flatten_matrix_row_major(matrix: &DMatrix<f64>) -> Vec<f64> {
    let mut flat = Vec::with_capacity(matrix.nrows() * matrix.ncols());
    for i in 0..matrix.nrows() {
        for j in 0..matrix.ncols() {
            flat.push(matrix[(i, j)]);
        }
    }
    flat
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

fn load_dat_coords(path: &std::path::Path) -> Vec<(f64, f64)> {
    let text = std::fs::read_to_string(path).expect("read paneled airfoil dat");
    text.lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            let x = parts.next()?.parse::<f64>().ok()?;
            let y = parts.next()?.parse::<f64>().ok()?;
            Some((x, y))
        })
        .collect()
}

fn build_xfoil_binary_reference(alpha_deg: f64) -> XfoilBinaryReference {
    let source_coords_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("naca0012.dat");
    let workdir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("target")
        .join(format!("xfoil-binary-{alpha_deg:.1}"));
    std::fs::create_dir_all(&workdir).expect("create xfoil binary workdir");
    let debug_path = workdir.join("xfoil_debug.json");
    let local_coords_path = workdir.join("naca0012.dat");
    let _ = std::fs::remove_file(&debug_path);
    std::fs::copy(&source_coords_path, &local_coords_path).expect("copy naca0012.dat into xfoil binary workdir");

    let commands = format!(
        "LOAD {}\n\nOPER\nVISC 1000000\nMACH 0\nITER 1\nALFA {alpha_deg:.1}\n\nQUIT\n",
        local_coords_path.file_name().expect("local airfoil filename").to_string_lossy(),
    );
    run_xfoil_instrumented(&commands, &workdir).expect("run xfoil_instrumented");
    let debug_text = read_debug_text(&debug_path);
    let debug: serde_json::Value = serde_json::from_str(&debug_text).expect("parse xfoil debug JSON");
    let events = debug["events"].as_array().expect("xfoil debug events");

    let mut upper = Vec::new();
    let mut lower = Vec::new();
    let mut cl = None;
    let mut cm = None;
    let mut cd = None;
    let mut x_tr_upper = None;
    let mut x_tr_lower = None;
    let mut converged = None;
    let mut iterations = None;
    for event in events {
        match event["subroutine"].as_str() {
            Some("MRCHUE") => {
                let side = event["side"].as_u64().expect("MRCHUE side") as usize;
                let ibl = event["ibl"].as_u64().expect("MRCHUE ibl") as usize;
                let tuple = (
                    ibl,
                    event["x"].as_f64().expect("MRCHUE x"),
                    event["Ue"].as_f64().expect("MRCHUE Ue"),
                    event["theta"].as_f64().expect("MRCHUE theta"),
                    event["delta_star"].as_f64().expect("MRCHUE delta_star"),
                    event["ctau"].as_f64().expect("MRCHUE ctau"),
                );
                if side == 1 {
                    upper.push(tuple);
                } else if side == 2 {
                    lower.push(tuple);
                }
            }
            Some("CL_DETAIL") => {
                cl = event["cl"].as_f64();
                cm = event["cm"].as_f64();
            }
            Some("CD_BREAKDOWN") => {
                cd = event["cd_total"].as_f64();
            }
            Some("VISCOUS_FINAL") => {
                x_tr_upper = event["x_tr_upper"].as_f64();
                x_tr_lower = event["x_tr_lower"].as_f64();
                converged = event["converged"].as_bool();
                iterations = event["iterations"].as_u64().map(|value| value as usize);
            }
            _ => {}
        }
    }

    XfoilBinaryReference {
        mrchue: mrchue_state_from_events(&upper, &lower),
        cl: cl.expect("xfoil cl"),
        cd: cd.expect("xfoil cd"),
        cm: cm.expect("xfoil cm"),
        x_tr_upper: x_tr_upper.expect("xfoil x_tr_upper"),
        x_tr_lower: x_tr_lower.expect("xfoil x_tr_lower"),
        converged: converged.expect("xfoil converged"),
        iterations: iterations.expect("xfoil iterations"),
    }
}

fn read_debug_text(debug_path: &std::path::Path) -> String {
    if let Ok(text) = std::fs::read_to_string(debug_path) {
        return text;
    }
    std::fs::read_to_string("xfoil_debug.json").expect("read xfoil_debug.json")
}

fn mrchue_state_from_events(
    upper_events: &[(usize, f64, f64, f64, f64, f64)],
    lower_events: &[(usize, f64, f64, f64, f64, f64)],
) -> FortranMarchState {
    fn fill(events: &[(usize, f64, f64, f64, f64, f64)]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let nbl = events.iter().map(|(ibl, ..)| *ibl).max().unwrap_or(1);
        let mut x = vec![0.0; nbl];
        let mut theta = vec![0.0; nbl];
        let mut dstr = vec![0.0; nbl];
        let mut uedg = vec![0.0; nbl];
        let mut ctau = vec![0.0; nbl];
        let mut mass = vec![0.0; nbl];
        for (ibl, xi, ue, th, ds, ct) in events {
            let idx = ibl - 1;
            x[idx] = *xi;
            theta[idx] = *th;
            dstr[idx] = *ds;
            uedg[idx] = *ue;
            ctau[idx] = *ct;
            mass[idx] = ds * ue;
        }
        (x, theta, dstr, uedg, ctau, mass)
    }

    let (upper_x, upper_theta, upper_dstr, upper_uedg, upper_ctau, upper_mass) = fill(upper_events);
    let (lower_x, lower_theta, lower_dstr, lower_uedg, lower_ctau, lower_mass) = fill(lower_events);
    FortranMarchState {
        nbl_upper: upper_x.len(),
        nbl_lower: lower_x.len(),
        iblte_upper: 0,
        iblte_lower: 0,
        itran_upper: 0,
        itran_lower: 0,
        xssitr_upper: 0.0,
        xssitr_lower: 0.0,
        lblini: true,
        upper_x,
        lower_x,
        upper_theta,
        lower_theta,
        upper_dstr,
        lower_dstr,
        upper_uedg,
        lower_uedg,
        upper_ctau,
        lower_ctau,
        upper_mass,
        lower_mass,
    }
}
