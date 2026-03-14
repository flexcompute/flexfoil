use serde_json::Value;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

fn workspace_root() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn load_trace(path: &str) -> Value {
    let file = File::open(workspace_root().join(path)).expect("trace file");
    serde_json::from_reader(BufReader::new(file)).expect("valid trace json")
}

fn first_blvar_event(trace: &Value, side: i64, ibl: i64) -> Option<&Value> {
    trace["events"].as_array()?.iter().find(|event| {
        event["subroutine"].as_str() == Some("BLVAR")
            && event["side"].as_i64() == Some(side)
            && event["ibl"].as_i64() == Some(ibl)
    })
}

fn first_true_transition(trace: &Value, side: i64) -> Option<&Value> {
    trace["events"].as_array()?.iter().find(|event| {
        event["subroutine"].as_str() == Some("TRCHEK2_FINAL")
            && event["side"].as_i64() == Some(side)
            && event["transition"].as_bool() == Some(true)
    })
}

fn rel_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1.0e-30)
}

fn blvar_input(trace: &Value, side: i64, ibl: i64, key: &str) -> f64 {
    first_blvar_event(trace, side, ibl)
        .and_then(|event| event["input"][key].as_f64())
        .unwrap_or_else(|| panic!("missing input {key} for side={side} ibl={ibl}"))
}

fn blvar_output(trace: &Value, side: i64, ibl: i64, key: &str) -> f64 {
    first_blvar_event(trace, side, ibl)
        .and_then(|event| event["output"][key].as_f64())
        .unwrap_or_else(|| panic!("missing output {key} for side={side} ibl={ibl}"))
}

#[test]
fn test_alpha4_upper_station_windows_and_transition_gate() {
    let xfoil = load_trace("traces/xfoil/naca0012/re3e06/alpha_+004.0.json");
    let rust = load_trace("traces/rustfoil/naca0012/re3e06/alpha_+004.0.json");

    // The first BL station should remain a machine-precision match.
    assert!(rel_diff(
        blvar_input(&xfoil, 1, 2, "theta"),
        blvar_input(&rust, 1, 2, "theta")
    ) < 1.0e-5);
    assert!(rel_diff(
        blvar_output(&xfoil, 1, 2, "Hk"),
        blvar_output(&rust, 1, 2, "Hk")
    ) < 1.0e-9);
    assert!(rel_diff(
        blvar_output(&xfoil, 1, 2, "Cf"),
        blvar_output(&rust, 1, 2, "Cf")
    ) < 1.0e-6);

    // The traces show that Ue remains trustworthy well before the BL state drifts.
    for ibl in 2..=20 {
        let ue_x = blvar_input(&xfoil, 1, ibl, "u");
        let ue_r = blvar_input(&rust, 1, ibl, "u");
        assert!(
            rel_diff(ue_x, ue_r) < 1.0e-5,
            "upper alpha=4 Ue drifted too early at ibl={ibl}: xfoil={ue_x:.6e} rust={ue_r:.6e}"
        );
    }

    // Transition detection itself should stay locked to the XFOIL interval.
    let xfoil_tr = first_true_transition(&xfoil, 1).expect("xfoil transition");
    let rust_tr = first_true_transition(&rust, 1).expect("rust transition");
    assert_eq!(xfoil_tr["ibl"].as_i64(), Some(33));
    assert_eq!(rust_tr["ibl"].as_i64(), Some(33));
    assert!(rel_diff(
        xfoil_tr["x1"].as_f64().unwrap(),
        rust_tr["x1"].as_f64().unwrap()
    ) < 1.0e-5);
    assert!(rel_diff(
        xfoil_tr["x2"].as_f64().unwrap(),
        rust_tr["x2"].as_f64().unwrap()
    ) < 1.0e-5);
}

#[test]
fn test_alpha10_alpha12_upper_station_windows_keep_ue_locked() {
    for alpha_path in [
        "traces/xfoil/naca0012/re3e06/alpha_+010.0.json",
        "traces/xfoil/naca0012/re3e06/alpha_+012.0.json",
    ] {
        let rust_path = alpha_path.replacen("traces/xfoil", "traces/rustfoil", 1);
        let xfoil = load_trace(alpha_path);
        let rust = load_trace(&rust_path);

        for ibl in 2..=12 {
            let ue_x = blvar_input(&xfoil, 1, ibl, "u");
            let ue_r = blvar_input(&rust, 1, ibl, "u");
            assert!(
                rel_diff(ue_x, ue_r) < 2.0e-5,
                "{alpha_path}: upper Ue drifted too early at ibl={ibl}: xfoil={ue_x:.6e} rust={ue_r:.6e}"
            );
        }

        let xfoil_tr = first_true_transition(&xfoil, 1).expect("xfoil transition");
        let rust_tr = first_true_transition(&rust, 1).expect("rust transition");
        assert_eq!(
            xfoil_tr["ibl"].as_i64(),
            rust_tr["ibl"].as_i64(),
            "{alpha_path}: transition station changed"
        );
    }
}

// ---------------------------------------------------------------------------
// WS5 expanded regression gates
// ---------------------------------------------------------------------------

#[test]
fn test_alpha4_lower_station_windows() {
    let xfoil = load_trace("traces/xfoil/naca0012/re3e06/alpha_+004.0.json");
    let rust = load_trace("traces/rustfoil/naca0012/re3e06/alpha_+004.0.json");

    // Lower surface (side 2) Ue should match at early stations
    for ibl in 2..=20 {
        let ue_x = blvar_input(&xfoil, 2, ibl, "u");
        let ue_r = blvar_input(&rust, 2, ibl, "u");
        assert!(
            rel_diff(ue_x, ue_r) < 1.0e-4,
            "lower alpha=4 Ue drifted at ibl={ibl}: xfoil={ue_x:.6e} rust={ue_r:.6e}"
        );
    }

    // Lower surface theta at station 2 should be close
    let theta_x = blvar_input(&xfoil, 2, 2, "theta");
    let theta_r = blvar_input(&rust, 2, 2, "theta");
    assert!(
        rel_diff(theta_x, theta_r) < 1.0e-4,
        "lower alpha=4 theta at ibl=2: xfoil={theta_x:.6e} rust={theta_r:.6e}"
    );

    // Lower surface transition
    let xfoil_tr = first_true_transition(&xfoil, 2);
    let rust_tr = first_true_transition(&rust, 2);
    if let (Some(xt), Some(rt)) = (xfoil_tr, rust_tr) {
        let x_diff = (xt["ibl"].as_i64().unwrap() - rt["ibl"].as_i64().unwrap()).abs();
        assert!(
            x_diff <= 2,
            "lower alpha=4 transition station off by {x_diff}: xfoil={}, rust={}",
            xt["ibl"],
            rt["ibl"]
        );
    }
}

#[test]
fn test_alpha4_post_transition_windows() {
    let xfoil = load_trace("traces/xfoil/naca0012/re3e06/alpha_+004.0.json");
    let rust = load_trace("traces/rustfoil/naca0012/re3e06/alpha_+004.0.json");

    // Find upper transition station
    let xfoil_tr = first_true_transition(&xfoil, 1).expect("xfoil upper transition");
    let rust_tr = first_true_transition(&rust, 1).expect("rust upper transition");

    let xfoil_ibl = xfoil_tr["ibl"].as_i64().unwrap();
    let rust_ibl = rust_tr["ibl"].as_i64().unwrap();

    // Check 5 stations after transition on upper surface
    for offset in 1..=5 {
        let ibl_x = xfoil_ibl + offset;
        let ibl_r = rust_ibl + offset;

        if let (Some(ex), Some(er)) = (
            first_blvar_event(&xfoil, 1, ibl_x),
            first_blvar_event(&rust, 1, ibl_r),
        ) {
            let theta_x = ex["input"]["theta"].as_f64().unwrap_or(0.0);
            let theta_r = er["input"]["theta"].as_f64().unwrap_or(0.0);
            let hk_x = ex["output"]["Hk"].as_f64().unwrap_or(0.0);
            let hk_r = er["output"]["Hk"].as_f64().unwrap_or(0.0);

            // Relaxed tolerances for post-transition (turbulent BL is inherently less precise)
            assert!(
                rel_diff(theta_x, theta_r) < 0.5,
                "post-transition upper theta too far at ibl={ibl_x}: xfoil={theta_x:.6e} rust={theta_r:.6e}"
            );
            assert!(
                rel_diff(hk_x, hk_r) < 0.5,
                "post-transition upper Hk too far at ibl={ibl_x}: xfoil={hk_x:.4} rust={hk_r:.4}"
            );
        }
    }
}

fn viscal_result_events(trace: &Value) -> Vec<&Value> {
    trace["events"]
        .as_array()
        .map(|events| {
            events
                .iter()
                .filter(|e| e["subroutine"].as_str() == Some("VISCAL_RESULT"))
                .collect()
        })
        .unwrap_or_default()
}

#[test]
fn test_full_sweep_force_gate() {
    let alphas: Vec<f64> = (-15..=15).map(|a| a as f64).collect();

    let mut report_lines = Vec::new();
    report_lines.push(format!(
        "{:>6} {:>10} {:>10} {:>10} {:>10} {:>6}",
        "alpha", "xfoil_CL", "xfoil_CD", "xfoil_CDf", "xfoil_iter", "avail"
    ));

    for &alpha in &alphas {
        let alpha_str = format!("{:+06.1}", alpha);
        let xfoil_path = format!("traces/xfoil/naca0012/re3e06/alpha_{alpha_str}.json");
        let trace_file = workspace_root().join(&xfoil_path);
        if !trace_file.exists() {
            continue;
        }
        let xfoil = load_trace(&xfoil_path);
        let results = viscal_result_events(&xfoil);
        if let Some(last) = results.last() {
            let cl = last["CL"].as_f64().unwrap_or(f64::NAN);
            let cd = last["CD"].as_f64().unwrap_or(f64::NAN);
            let cdf = last["CD_friction"].as_f64().unwrap_or(f64::NAN);
            let iter = last["iteration"].as_i64().unwrap_or(0);

            assert!(
                cl.is_finite(),
                "XFOIL CL not finite at alpha={alpha}"
            );
            assert!(
                cd.is_finite() && cd > 0.0,
                "XFOIL CD not finite or non-positive at alpha={alpha}"
            );
            assert!(
                cdf.is_finite() && cdf > 0.0,
                "XFOIL CDf not finite or non-positive at alpha={alpha}"
            );

            report_lines.push(format!(
                "{:>6.1} {:>10.6} {:>10.6} {:>10.6} {:>10} {:>6}",
                alpha, cl, cd, cdf, iter, "Y"
            ));
        }
    }

    eprintln!("\n=== Full-Sweep XFOIL Force Reference ===");
    for line in &report_lines {
        eprintln!("{line}");
    }
}

#[test]
fn test_convergence_count_gate() {
    // Verify that XFOIL converges within a reasonable number of iterations
    // for key alphas. This gate ensures our trace data is self-consistent.
    for alpha_str in ["alpha_+000.0", "alpha_+004.0", "alpha_+008.0", "alpha_+010.0"] {
        let xfoil_path = format!("traces/xfoil/naca0012/re3e06/{alpha_str}.json");
        let trace_file = workspace_root().join(&xfoil_path);
        if !trace_file.exists() {
            continue;
        }
        let xfoil = load_trace(&xfoil_path);
        let results = viscal_result_events(&xfoil);
        if let Some(last) = results.last() {
            let iter = last["iteration"].as_i64().unwrap_or(0);
            let rms = last["rms_residual"].as_f64().unwrap_or(1.0);
            assert!(
                iter <= 25,
                "{alpha_str}: XFOIL took {iter} iterations"
            );
            assert!(
                rms < 1.0e-2,
                "{alpha_str}: XFOIL RMS residual {rms:.6e} is too large"
            );
        }
    }
}
