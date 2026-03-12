use rustfoil_bl::state::BlStation;
use rustfoil_coupling::stmove::find_stagnation_with_derivs;
use rustfoil_solver::viscous::{
    compute_n_wake_panels, extract_surface_xfoil, initialize_surface_stations_with_panel_idx,
    initialize_wake_bl_stations, setup_from_coords, solve_viscous_two_surfaces,
    ViscousSolverConfig,
};
use serde::Serialize;
use serde_json::Value;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize)]
struct StructuralMarkers {
    alpha_deg: i32,
    initial_ist: usize,
    final_ist: usize,
    initial_upper_len: usize,
    initial_lower_len: usize,
    final_upper_len: usize,
    final_lower_len: usize,
    initial_upper_ue: [f64; 2],
    initial_lower_ue: [f64; 2],
    final_upper_ue: [f64; 2],
    final_lower_ue: [f64; 2],
    cl: f64,
}

#[derive(Debug, Clone, Serialize)]
struct StructuralComparison {
    alpha_deg: i32,
    rust: StructuralMarkers,
    xfoil: StructuralMarkers,
    initial_ist_diff: i32,
    final_ist_diff: i32,
    initial_upper_len_diff: i32,
    initial_lower_len_diff: i32,
    final_upper_len_diff: i32,
    final_lower_len_diff: i32,
    cl_diff: f64,
}

const STRUCTURAL_SWEEP_ALPHAS: &[i32] = &[-15, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 15];

fn workspace_root() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir.parent().unwrap().parent().unwrap().to_path_buf()
}

fn load_naca2412_coords() -> Vec<(f64, f64)> {
    let path = workspace_root().join("naca2412_xfoil_paneled.dat");
    let content = fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("failed to read {:?}", path));

    content
        .lines()
        .filter_map(|line| {
            let parts: Vec<_> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                Some((parts[0].parse().ok()?, parts[1].parse().ok()?))
            } else {
                None
            }
        })
        .collect()
}

fn alpha_trace_name(alpha_deg: i32) -> String {
    format!("alpha_{:+06.1}.json", alpha_deg as f64)
}

fn xfoil_trace_path(alpha_deg: i32) -> PathBuf {
    workspace_root()
        .join("traces")
        .join("xfoil")
        .join("naca2412")
        .join("re3e06")
        .join(alpha_trace_name(alpha_deg))
}

fn find_matching_brace(content: &str, start_idx: usize) -> Option<usize> {
    let mut depth = 0usize;
    let mut in_string = false;
    let mut escaped = false;

    for (offset, ch) in content[start_idx..].char_indices() {
        if in_string {
            if escaped {
                escaped = false;
            } else if ch == '\\' {
                escaped = true;
            } else if ch == '"' {
                in_string = false;
            }
            continue;
        }

        match ch {
            '"' => in_string = true,
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    return Some(start_idx + offset + ch.len_utf8());
                }
            }
            _ => {}
        }
    }

    None
}

fn extract_event_blocks<'a>(content: &'a str, needle: &str) -> Vec<&'a str> {
    let mut blocks = Vec::new();
    let mut search_from = 0usize;

    while let Some(relative_pos) = content[search_from..].find(needle) {
        let pos = search_from + relative_pos;
        let start = content[..pos].rfind('{').expect("event start");
        let end = find_matching_brace(content, start).expect("event end");
        blocks.push(&content[start..end]);
        search_from = end;
    }

    blocks
}

fn parse_event(block: &str) -> Value {
    serde_json::from_str(block).unwrap_or_else(|_| panic!("failed to parse event block"))
}

fn parse_usize_from_line(line: &str) -> Option<usize> {
    line.split(':')
        .nth(1)?
        .trim()
        .trim_end_matches(',')
        .parse()
        .ok()
}

fn extract_full_bl_lengths(content: &str) -> Vec<(usize, usize)> {
    let lines: Vec<&str> = content.lines().collect();
    let mut lengths = Vec::new();

    for (idx, line) in lines.iter().enumerate() {
        if !line.contains("\"subroutine\": \"FULL_BL_STATE\"") {
            continue;
        }

        let mut nbl_upper = None;
        let mut nbl_lower = None;
        for next_line in lines.iter().skip(idx + 1).take(6) {
            if next_line.contains("\"nbl_upper\"") {
                nbl_upper = parse_usize_from_line(next_line);
            } else if next_line.contains("\"nbl_lower\"") {
                nbl_lower = parse_usize_from_line(next_line);
            }
        }

        if let (Some(upper), Some(lower)) = (nbl_upper, nbl_lower) {
            lengths.push((upper, lower));
        }
    }

    lengths
}

fn event_i64(event: &Value, key: &str) -> i64 {
    event
        .get(key)
        .and_then(Value::as_i64)
        .unwrap_or_else(|| panic!("missing i64 field {key}"))
}

fn blvar_ue(events: &[Value], side: i64, ibl: i64, from_end: bool) -> f64 {
    let iter = if from_end {
        Box::new(events.iter().rev()) as Box<dyn Iterator<Item = &Value>>
    } else {
        Box::new(events.iter()) as Box<dyn Iterator<Item = &Value>>
    };

    iter.filter(|event| {
        event_i64(event, "side") == side
            && event_i64(event, "ibl") == ibl
            && event.get("subroutine").and_then(Value::as_str) == Some("BLVAR")
    })
    .find_map(|event| event.get("input")?.get("u")?.as_f64())
    .unwrap_or(0.0)
}

fn event_usize(event: &Value, key: &str) -> usize {
    event
        .get(key)
        .and_then(Value::as_u64)
        .unwrap_or_else(|| panic!("missing usize field {key}"))
        as usize
}

fn event_f64(event: &Value, key: &str) -> f64 {
    event
        .get(key)
        .and_then(Value::as_f64)
        .unwrap_or_else(|| panic!("missing f64 field {key}"))
}

fn first_two_station_ue(stations: &[BlStation]) -> [f64; 2] {
    let values: Vec<f64> = stations.iter().skip(1).take(2).map(|station| station.u).collect();
    [
        values.first().copied().unwrap_or(0.0),
        values.get(1).copied().unwrap_or(0.0),
    ]
}

fn reconstruct_airfoil_gamma(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_airfoil_panels: usize,
) -> Vec<f64> {
    let mut gamma = vec![0.0; n_airfoil_panels];

    for station in upper_stations.iter().skip(1) {
        if !station.is_wake && station.panel_idx < n_airfoil_panels {
            gamma[station.panel_idx] = station.u;
        }
    }
    for station in lower_stations.iter().skip(1) {
        if !station.is_wake && station.panel_idx < n_airfoil_panels {
            gamma[station.panel_idx] = -station.u;
        }
    }

    gamma
}

fn load_xfoil_markers(alpha_deg: i32) -> StructuralMarkers {
    let path = xfoil_trace_path(alpha_deg);
    let content = fs::read_to_string(&path).unwrap_or_else(|_| panic!("failed to read {:?}", path));
    let stfind_events: Vec<Value> = extract_event_blocks(&content, "\"subroutine\": \"STFIND\"")
        .into_iter()
        .map(parse_event)
        .collect();
    let iblpan_events: Vec<Value> = extract_event_blocks(&content, "\"subroutine\": \"IBLPAN\"")
        .into_iter()
        .map(parse_event)
        .collect();
    let viscal_events: Vec<Value> =
        extract_event_blocks(&content, "\"subroutine\": \"VISCAL_RESULT\"")
            .into_iter()
            .map(parse_event)
            .collect();
    let ueset_events: Vec<Value> = extract_event_blocks(&content, "\"subroutine\": \"UESET\"")
        .into_iter()
        .map(parse_event)
        .collect();
    let blvar_events: Vec<Value> = extract_event_blocks(&content, "\"subroutine\": \"BLVAR\"")
        .into_iter()
        .map(parse_event)
        .collect();

    let first_stfind = stfind_events.first().expect("first STFIND");
    let last_stfind = stfind_events.last().expect("last STFIND");
    let last_viscal = viscal_events.last().expect("last VISCAL_RESULT");
    let full_bl_lengths = extract_full_bl_lengths(&content);
    let first_ueset = ueset_events.first();
    let last_ueset = ueset_events.last();
    let (initial_upper_len, initial_lower_len, final_upper_len, final_lower_len) =
        if let (Some(first_iblpan), Some(last_iblpan)) = (iblpan_events.first(), iblpan_events.last()) {
            (
                event_usize(first_iblpan, "NBL_upper"),
                event_usize(first_iblpan, "NBL_lower"),
                event_usize(last_iblpan, "NBL_upper"),
                event_usize(last_iblpan, "NBL_lower"),
            )
        } else if let (Some(first_ueset), Some(last_ueset)) = (first_ueset, last_ueset) {
            (
                event_usize(first_ueset, "nbl_upper"),
                event_usize(first_ueset, "nbl_lower"),
                event_usize(last_ueset, "nbl_upper"),
                event_usize(last_ueset, "nbl_lower"),
            )
        } else {
            let first = *full_bl_lengths.first().expect("first FULL_BL_STATE lengths");
            let last = *full_bl_lengths.last().expect("last FULL_BL_STATE lengths");
            (first.0, first.1, last.0, last.1)
        };

    StructuralMarkers {
        alpha_deg,
        initial_ist: event_usize(first_stfind, "IST").saturating_sub(1),
        final_ist: event_usize(last_stfind, "IST").saturating_sub(1),
        initial_upper_len,
        initial_lower_len,
        final_upper_len,
        final_lower_len,
        initial_upper_ue: [
            blvar_ue(&blvar_events, 1, 2, false),
            blvar_ue(&blvar_events, 1, 3, false),
        ],
        initial_lower_ue: [
            blvar_ue(&blvar_events, 2, 2, false),
            blvar_ue(&blvar_events, 2, 3, false),
        ],
        final_upper_ue: [
            blvar_ue(&blvar_events, 1, 2, true),
            blvar_ue(&blvar_events, 1, 3, true),
        ],
        final_lower_ue: [
            blvar_ue(&blvar_events, 2, 2, true),
            blvar_ue(&blvar_events, 2, 3, true),
        ],
        cl: event_f64(last_viscal, "CL"),
    }
}

fn run_rust_markers(alpha_deg: i32) -> StructuralMarkers {
    let coords = load_naca2412_coords();
    let setup_result = setup_from_coords(&coords, alpha_deg as f64).expect("setup");
    let full_arc = &setup_result.setup.arc_lengths;
    let ue_inviscid = &setup_result.setup.ue_inviscid;
    let ue_inviscid_alpha = setup_result.operating_sensitivity();
    let ist = setup_result.ist;
    let config = ViscousSolverConfig::with_reynolds(3.0e6).with_max_iterations(20);
    let n_airfoil_panels = setup_result.node_x.len();
    let (mut upper_stations, mut lower_stations) =
        setup_result.derive_station_views(config.reynolds, config.msq());
    let (upper_ue, lower_ue_extended) =
        setup_result.derive_uedg_views(config.reynolds, config.msq());

    let initial_upper_len = upper_stations.len();
    let initial_lower_len = lower_stations.len();
    let initial_upper_ue = first_two_station_ue(&upper_stations);
    let initial_lower_ue = first_two_station_ue(&lower_stations);

    let viscous = solve_viscous_two_surfaces(
        &mut upper_stations,
        &mut lower_stations,
        &upper_ue,
        &lower_ue_extended,
        &setup_result.setup.dij,
        &config,
        (alpha_deg as f64).to_radians(),
        &setup_result.node_x,
        &setup_result.node_y,
        ue_inviscid,
        &ue_inviscid_alpha,
    )
    .expect("viscous solve");

    let final_gamma = reconstruct_airfoil_gamma(&upper_stations, &lower_stations, n_airfoil_panels);
    let final_ist = find_stagnation_with_derivs(&final_gamma, full_arc)
        .map(|stag| stag.ist)
        .unwrap_or(ist);

    StructuralMarkers {
        alpha_deg,
        initial_ist: ist,
        final_ist,
        initial_upper_len,
        initial_lower_len,
        final_upper_len: upper_stations.len(),
        final_lower_len: lower_stations.len(),
        initial_upper_ue,
        initial_lower_ue,
        final_upper_ue: first_two_station_ue(&upper_stations),
        final_lower_ue: first_two_station_ue(&lower_stations),
        cl: viscous.cl,
    }
}

fn compare_markers(alpha_deg: i32) -> StructuralComparison {
    let rust = run_rust_markers(alpha_deg);
    let xfoil = load_xfoil_markers(alpha_deg);

    StructuralComparison {
        alpha_deg,
        initial_ist_diff: rust.initial_ist as i32 - xfoil.initial_ist as i32,
        final_ist_diff: rust.final_ist as i32 - xfoil.final_ist as i32,
        initial_upper_len_diff: rust.initial_upper_len as i32 - xfoil.initial_upper_len as i32,
        initial_lower_len_diff: rust.initial_lower_len as i32 - xfoil.initial_lower_len as i32,
        final_upper_len_diff: rust.final_upper_len as i32 - xfoil.final_upper_len as i32,
        final_lower_len_diff: rust.final_lower_len as i32 - xfoil.final_lower_len as i32,
        cl_diff: rust.cl - xfoil.cl,
        rust,
        xfoil,
    }
}

fn write_summary(name: &str, comparisons: &[StructuralComparison]) {
    let path = workspace_root().join("target").join(name);
    fs::write(&path, serde_json::to_string_pretty(comparisons).unwrap())
        .unwrap_or_else(|_| panic!("failed to write {:?}", path));
}

fn write_json<T: Serialize>(name: &str, data: &T) {
    let path = workspace_root().join("target").join(name);
    fs::write(&path, serde_json::to_string_pretty(data).unwrap())
        .unwrap_or_else(|_| panic!("failed to write {:?}", path));
}

// ---------------------------------------------------------------------------
// First-divergence report types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
enum DivergenceStage {
    Setup,
    Relocation,
    Transition,
    Wake,
    Converged,
}

#[derive(Debug, Clone, Serialize)]
struct IterationSnapshot {
    iteration: usize,
    cl: f64,
    cd: f64,
    cm: f64,
    rms_residual: f64,
    ist: usize,
}

#[derive(Debug, Clone, Serialize)]
struct DivergenceReport {
    alpha_deg: i32,
    first_divergence_iter: Option<usize>,
    divergence_stage: DivergenceStage,
    xfoil_iterations: Vec<IterationSnapshot>,
    rust_cl: f64,
    rust_cd: f64,
    rust_converged: bool,
    rust_iterations: usize,
    xfoil_cl: f64,
    xfoil_cd: f64,
    cl_error_pct: f64,
    cd_error_pct: f64,
    xfoil_x_tr_upper: f64,
    xfoil_x_tr_lower: f64,
    rust_x_tr_upper: f64,
    rust_x_tr_lower: f64,
}

fn load_xfoil_iteration_snapshots(alpha_deg: i32) -> Vec<IterationSnapshot> {
    let path = xfoil_trace_path(alpha_deg);
    let content = fs::read_to_string(&path).unwrap_or_else(|_| panic!("failed to read {:?}", path));

    let viscal_events: Vec<Value> =
        extract_event_blocks(&content, "\"subroutine\": \"VISCAL_RESULT\"")
            .into_iter()
            .map(parse_event)
            .collect();
    let stfind_events: Vec<Value> = extract_event_blocks(&content, "\"subroutine\": \"STFIND\"")
        .into_iter()
        .map(parse_event)
        .collect();

    viscal_events
        .iter()
        .enumerate()
        .map(|(i, event)| {
            let ist = stfind_events
                .get(i)
                .and_then(|s| s.get("IST").and_then(Value::as_u64))
                .map(|v| v as usize)
                .unwrap_or(0);
            IterationSnapshot {
                iteration: event_usize(event, "iteration"),
                cl: event_f64(event, "CL"),
                cd: event_f64(event, "CD"),
                cm: event_f64(event, "CM"),
                rms_residual: event_f64(event, "rms_residual"),
                ist,
            }
        })
        .collect()
}

fn xfoil_transition_locations(alpha_deg: i32) -> (f64, f64) {
    let path = xfoil_trace_path(alpha_deg);
    let content = fs::read_to_string(&path).unwrap_or_else(|_| panic!("failed to read {:?}", path));

    let tr_events: Vec<Value> =
        extract_event_blocks(&content, "\"subroutine\": \"TRCHEK2_FINAL\"")
            .into_iter()
            .map(parse_event)
            .filter(|e| e.get("transition").and_then(Value::as_bool).unwrap_or(false))
            .collect();

    let x_tr_upper = tr_events
        .iter()
        .find(|e| event_i64(e, "side") == 1)
        .and_then(|e| e.get("x1").and_then(Value::as_f64))
        .unwrap_or(1.0);
    let x_tr_lower = tr_events
        .iter()
        .find(|e| event_i64(e, "side") == 2)
        .and_then(|e| e.get("x1").and_then(Value::as_f64))
        .unwrap_or(1.0);

    (x_tr_upper, x_tr_lower)
}

fn classify_divergence(
    xfoil_snapshots: &[IterationSnapshot],
    xfoil_ist_initial: usize,
) -> (Option<usize>, DivergenceStage) {
    if xfoil_snapshots.is_empty() {
        return (Some(0), DivergenceStage::Setup);
    }

    for (i, snap) in xfoil_snapshots.iter().enumerate() {
        if i > 0 && snap.ist != xfoil_ist_initial && snap.ist != xfoil_snapshots[i - 1].ist {
            return (Some(snap.iteration), DivergenceStage::Relocation);
        }
    }

    (None, DivergenceStage::Converged)
}

fn build_divergence_report(alpha_deg: i32) -> DivergenceReport {
    let xfoil_snapshots = load_xfoil_iteration_snapshots(alpha_deg);
    let xfoil_markers = load_xfoil_markers(alpha_deg);
    let (xfoil_x_tr_upper, xfoil_x_tr_lower) = xfoil_transition_locations(alpha_deg);

    let rust_markers = run_rust_markers(alpha_deg);

    let xfoil_final = xfoil_snapshots.last().cloned().unwrap_or(IterationSnapshot {
        iteration: 0,
        cl: 0.0,
        cd: 0.0,
        cm: 0.0,
        rms_residual: 1.0,
        ist: 0,
    });

    let cl_error_pct = if xfoil_final.cl.abs() > 1e-6 {
        100.0 * (rust_markers.cl - xfoil_final.cl).abs() / xfoil_final.cl.abs()
    } else {
        0.0
    };
    let cd_error_pct = if xfoil_final.cd.abs() > 1e-6 {
        100.0 * (0.0_f64 - xfoil_final.cd).abs() / xfoil_final.cd.abs()
    } else {
        0.0
    };

    let (first_div, stage) =
        classify_divergence(&xfoil_snapshots, xfoil_markers.initial_ist);

    let final_stage = if rust_markers.final_ist as i32 != xfoil_markers.final_ist as i32
        && (rust_markers.final_ist as i32 - xfoil_markers.final_ist as i32).abs() > 1
    {
        DivergenceStage::Relocation
    } else if (xfoil_x_tr_upper - 1.0).abs() > 0.01 && cl_error_pct > 5.0 {
        DivergenceStage::Transition
    } else if cd_error_pct > 30.0 {
        DivergenceStage::Wake
    } else if cl_error_pct > 5.0 {
        stage
    } else {
        DivergenceStage::Converged
    };

    DivergenceReport {
        alpha_deg,
        first_divergence_iter: first_div,
        divergence_stage: final_stage,
        xfoil_iterations: xfoil_snapshots,
        rust_cl: rust_markers.cl,
        rust_cd: 0.0,
        rust_converged: rust_markers.cl.is_finite(),
        rust_iterations: 0,
        xfoil_cl: xfoil_final.cl,
        xfoil_cd: xfoil_final.cd,
        cl_error_pct,
        cd_error_pct,
        xfoil_x_tr_upper,
        xfoil_x_tr_lower,
        rust_x_tr_upper: 1.0,
        rust_x_tr_lower: 1.0,
    }
}

#[test]
fn test_naca2412_high_alpha_initial_setup_signature() {
    let comparisons: Vec<_> = [12, 15].into_iter().map(compare_markers).collect();
    write_summary("naca2412_high_alpha_setup.json", &comparisons);

    for comparison in comparisons {
        assert!(
            comparison.initial_ist_diff.abs() <= 1,
            "alpha {} initial IST drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(
            comparison.initial_upper_len_diff.abs() <= 1
                && comparison.initial_lower_len_diff.abs() <= 2,
            "alpha {} initial split-length drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(comparison.rust.initial_upper_ue[0] > 0.0);
        assert!(comparison.rust.initial_lower_ue[0] > 0.0);
    }
}

#[test]
fn test_naca2412_full_structural_sweep_gate() {
    let comparisons: Vec<_> = STRUCTURAL_SWEEP_ALPHAS
        .iter()
        .copied()
        .map(compare_markers)
        .collect();
    write_summary("naca2412_structural_sweep.json", &comparisons);

    for comparison in &comparisons {
        assert!(
            comparison.initial_ist_diff.abs() <= 1,
            "alpha {} initial IST drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(
            comparison.final_ist_diff.abs() <= 2,
            "alpha {} final IST drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(
            comparison.initial_upper_len_diff.abs() <= 1
                && comparison.initial_lower_len_diff.abs() <= 2,
            "alpha {} initial length drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(
            comparison.final_upper_len_diff.abs() <= 2
                && comparison.final_lower_len_diff.abs() <= 2,
            "alpha {} final length drift too large: {:?}",
            comparison.alpha_deg,
            comparison
        );
        assert!(comparison.rust.cl.is_finite(), "alpha {} CL is not finite", comparison.alpha_deg);
    }
}

#[test]
fn test_naca2412_first_divergence_sweep() {
    let reports: Vec<DivergenceReport> = STRUCTURAL_SWEEP_ALPHAS
        .iter()
        .copied()
        .map(build_divergence_report)
        .collect();

    write_json("divergence_report.json", &reports);

    eprintln!("\n=== First-Divergence All-Alpha Report ===");
    eprintln!(
        "{:>6} {:>12} {:>12} {:>10} {:>10} {:>10} {:>10}",
        "alpha", "stage", "first_div", "CL_err%", "CD_err%", "xtr_up_xf", "xtr_lo_xf"
    );

    for report in &reports {
        let div_str = report
            .first_divergence_iter
            .map(|i| format!("{}", i))
            .unwrap_or_else(|| "-".to_string());
        eprintln!(
            "{:>6} {:>12?} {:>12} {:>10.2} {:>10.2} {:>10.4} {:>10.4}",
            report.alpha_deg,
            report.divergence_stage,
            div_str,
            report.cl_error_pct,
            report.cd_error_pct,
            report.xfoil_x_tr_upper,
            report.xfoil_x_tr_lower,
        );
    }

    for report in &reports {
        assert!(
            report.rust_cl.is_finite(),
            "alpha {} rust CL is not finite",
            report.alpha_deg,
        );
    }
}
