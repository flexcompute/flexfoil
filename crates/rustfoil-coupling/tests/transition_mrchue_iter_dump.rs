use rustfoil_bl::{finalize_debug, init_debug};
use rustfoil_coupling::march::{march_fixed_ue, MarchConfig};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct MrchueReference {
    metadata: Metadata,
    sides: std::collections::HashMap<String, SideData>,
}

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct SideData {
    stations: Vec<StationRef>,
}

#[derive(Debug, Deserialize)]
struct StationRef {
    ibl: usize,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
}

#[derive(Debug, Deserialize)]
struct DebugOutput {
    events: Vec<serde_json::Value>,
}

fn load_reference() -> Option<MrchueReference> {
    let paths = [
        "testdata/mrchue_iterations.json",
        "../testdata/mrchue_iterations.json",
        "../../testdata/mrchue_iterations.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(data) = serde_json::from_str(&content) {
                return Some(data);
            }
        }
    }
    None
}

fn load_xfoil_debug() -> Option<DebugOutput> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").ok();
    let mut paths = vec![
        "Xfoil-instrumented/bin/xfoil_debug.json".to_string(),
        "../Xfoil-instrumented/bin/xfoil_debug.json".to_string(),
        "../../Xfoil-instrumented/bin/xfoil_debug.json".to_string(),
    ];
    if let Some(dir) = manifest_dir {
        paths.insert(0, format!("{}/../../Xfoil-instrumented/bin/xfoil_debug.json", dir));
    }

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(data) = serde_json::from_str(&content) {
                println!("Loaded XFOIL debug: {}", path);
                return Some(data);
            }
        }
    }
    None
}

#[ignore]
#[test]
fn dump_mrchue_iter_for_transition() {
    let ref_data = match load_reference() {
        Some(data) => data,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let side_data = match ref_data.sides.get("1") {
        Some(data) => data,
        None => {
            eprintln!("Skipping: side 1 data not found");
            return;
        }
    };

    let xfoil_debug = match load_xfoil_debug() {
        Some(data) => data,
        None => {
            eprintln!("Skipping: xfoil_debug.json not found");
            return;
        }
    };

    let mut xfoil_u2_by_ibl = std::collections::HashMap::new();
    for event in &xfoil_debug.events {
        if event.get("subroutine") == Some(&serde_json::Value::String("BLPRV_STATE".to_string())) {
            if event.get("side") == Some(&serde_json::Value::Number(1.into()))
                && event.get("newton_iter") == Some(&serde_json::Value::Number(1.into()))
            {
                if let (Some(ibl), Some(u2)) = (event.get("ibl"), event.get("U2")) {
                    if let (Some(ibl), Some(u2)) = (ibl.as_u64(), u2.as_f64()) {
                        xfoil_u2_by_ibl.insert(ibl as usize, u2);
                    }
                }
            }
        }
    }

    let x: Vec<f64> = side_data.stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = side_data
        .stations
        .iter()
        .map(|s| xfoil_u2_by_ibl.get(&s.ibl).copied().unwrap_or(s.ue))
        .collect();

    let re = ref_data.metadata.reynolds;
    let msq = ref_data.metadata.mach * ref_data.metadata.mach;

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        ..Default::default()
    };

    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let debug_dir = format!("{}/target", manifest_dir);
    let debug_path = format!("{}/rust_debug_transition.json", debug_dir);
    let _ = fs::create_dir_all(&debug_dir);
    init_debug(&debug_path);
    let _result = march_fixed_ue(&x, &ue, re, msq, &config);
    finalize_debug();

    let rust_debug: DebugOutput =
        serde_json::from_str(&fs::read_to_string(&debug_path).unwrap()).unwrap();
    let xfoil_trdif_count = xfoil_debug
        .events
        .iter()
        .filter(|ev| ev.get("subroutine").and_then(|v| v.as_str()) == Some("TRDIF"))
        .count();
    println!("XFOIL TRDIF event count: {}", xfoil_trdif_count);

    let mut rust_iters: Vec<(usize, f64, f64)> = Vec::new();
    let mut rust_sys: std::collections::HashMap<usize, (Vec<Vec<f64>>, Vec<f64>)> =
        std::collections::HashMap::new();
    let mut rust_trdif: std::collections::HashMap<
        usize,
        (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>),
    > = std::collections::HashMap::new();
    let mut rust_trdif_derivs: Vec<serde_json::Map<String, serde_json::Value>> = Vec::new();
    for ev in rust_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("MRCHUE_ITER") {
            continue;
        }
        match ev.get("side").and_then(|v| v.as_u64()) {
            Some(0) | Some(1) => {}
            _ => continue,
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(30) {
            continue;
        }
        let iter = ev.get("iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let ctau = ev.get("ctau").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let theta = ev.get("theta").and_then(|v| v.as_f64()).unwrap_or(0.0);
        rust_iters.push((iter, ctau, theta));
    }
    rust_iters.sort_by_key(|(iter, _, _)| *iter);

    let mut xfoil_iters: Vec<(usize, f64, f64)> = Vec::new();
    let mut xfoil_sys: std::collections::HashMap<usize, (Vec<Vec<f64>>, Vec<f64>)> =
        std::collections::HashMap::new();
    let mut xfoil_trdif: std::collections::HashMap<
        usize,
        (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>),
    > = std::collections::HashMap::new();
    let mut xfoil_trdif_derivs: Vec<serde_json::Map<String, serde_json::Value>> = Vec::new();
    for ev in xfoil_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("MRCHUE_ITER") {
            continue;
        }
        if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
            continue;
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(32) {
            continue;
        }
        let iter = ev.get("newton_iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let output = ev.get("output").and_then(|v| v.as_object()).unwrap();
        let ctau = output.get("ctau").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let theta = output.get("theta").and_then(|v| v.as_f64()).unwrap_or(0.0);
        xfoil_iters.push((iter, ctau, theta));
    }
    xfoil_iters.sort_by_key(|(iter, _, _)| *iter);

    for ev in xfoil_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("VS2_BEFORE") {
            continue;
        }
        if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
            continue;
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(32) {
            continue;
        }
        let iter = ev.get("newton_iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let vs2 = ev
            .get("VS2_4x4")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vsrez = ev
            .get("VSREZ_rhs")
            .and_then(|v| v.as_array())
            .map(|vals| vals.iter().map(|v| v.as_f64().unwrap()).collect::<Vec<f64>>())
            .unwrap_or_default();
        xfoil_sys.insert(iter, (vs2, vsrez));
    }

    for ev in xfoil_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("TRDIF") {
            continue;
        }
        if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
            continue;
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(31) {
            continue;
        }
        let iter = ev.get("newton_iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let vs1 = ev
            .get("VS1")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vs2 = ev
            .get("VS2")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vsrez = ev
            .get("VSREZ")
            .and_then(|v| v.as_array())
            .map(|vals| vals.iter().map(|v| v.as_f64().unwrap()).collect::<Vec<f64>>())
            .unwrap_or_default();
        xfoil_trdif.insert(iter, (vs1, vs2, vsrez));
    }

    for ev in xfoil_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("TRDIF_DERIVS") {
            continue;
        }
        if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
            continue;
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(31) {
            continue;
        }
        if let Some(obj) = ev.as_object() {
            xfoil_trdif_derivs.push(obj.clone());
        }
    }

    println!("Rust MRCHUE_ITER (ibl=32):");
    for (iter, ctau, theta) in rust_iters {
        println!("  iter {:2}: ctau={:.6} theta={:.8}", iter, ctau, theta);
    }

    println!("XFOIL MRCHUE_ITER (ibl=32):");
    for (iter, ctau, theta) in xfoil_iters {
        println!("  iter {:2}: ctau={:.6} theta={:.8}", iter, ctau, theta);
    }

    for ev in rust_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("VS2_BEFORE") {
            continue;
        }
        match ev.get("side").and_then(|v| v.as_u64()) {
            Some(0) | Some(1) => {}
            _ => continue,
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(30) {
            continue;
        }
        let iter = ev.get("iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let vs2 = ev
            .get("VS2_4x4")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vsrez = ev
            .get("VSREZ_rhs")
            .and_then(|v| v.as_array())
            .map(|vals| vals.iter().map(|v| v.as_f64().unwrap()).collect::<Vec<f64>>())
            .unwrap_or_default();
        rust_sys.insert(iter, (vs2, vsrez));
    }

    for ev in rust_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("TRDIF") {
            continue;
        }
        match ev.get("side").and_then(|v| v.as_u64()) {
            Some(0) | Some(1) => {}
            _ => continue,
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(30) {
            continue;
        }
        let iter = ev.get("iter").and_then(|v| v.as_u64()).unwrap_or(0) as usize;
        let vs1 = ev
            .get("VS1")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vs2 = ev
            .get("VS2")
            .and_then(|v| v.as_array())
            .map(|rows| {
                rows.iter()
                    .map(|row| {
                        row.as_array()
                            .unwrap()
                            .iter()
                            .map(|v| v.as_f64().unwrap())
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>()
            })
            .unwrap_or_default();
        let vsrez = ev
            .get("VSREZ")
            .and_then(|v| v.as_array())
            .map(|vals| vals.iter().map(|v| v.as_f64().unwrap()).collect::<Vec<f64>>())
            .unwrap_or_default();
        rust_trdif.insert(iter, (vs1, vs2, vsrez));
    }

    for ev in rust_debug.events.iter() {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("TRDIF_DERIVS") {
            continue;
        }
        match ev.get("side").and_then(|v| v.as_u64()) {
            Some(0) | Some(1) => {}
            _ => continue,
        }
        if ev.get("ibl").and_then(|v| v.as_u64()) != Some(30) {
            continue;
        }
        if let Some(obj) = ev.as_object() {
            rust_trdif_derivs.push(obj.clone());
        }
    }

    println!("System compare (rust vs XFOIL 4x4):");
    let mut iters: Vec<usize> = rust_sys.keys().cloned().collect();
    iters.sort();
    for iter in iters {
        let (r_vs2, r_vsr) = rust_sys.get(&iter).unwrap();
        if let Some((x_vs2, x_vsr)) = xfoil_sys.get(&iter) {
            let mut max_vs2: f64 = 0.0;
            for i in 0..r_vs2.len().min(x_vs2.len()) {
                for j in 0..r_vs2[i].len().min(x_vs2[i].len()) {
                    let rv = r_vs2[i][j];
                    let xv = x_vs2[i][j];
                    max_vs2 = max_vs2.max((rv - xv).abs());
                }
            }
            let mut max_vsr: f64 = 0.0;
            for k in 0..r_vsr.len().min(x_vsr.len()) {
                max_vsr = max_vsr.max((r_vsr[k] - x_vsr[k]).abs());
            }
            println!("  iter {:2}: max_vs2={:.3e} max_vsr={:.3e}", iter, max_vs2, max_vsr);
        } else {
            println!("  iter {:2}: no XFOIL system", iter);
        }
    }

    let mut xfoil_trdif_iters: Vec<usize> = xfoil_trdif.keys().cloned().collect();
    xfoil_trdif_iters.sort();
    let mut rust_trdif_iters: Vec<usize> = rust_trdif.keys().cloned().collect();
    rust_trdif_iters.sort();
    println!("XFOIL TRDIF iters: {:?}", xfoil_trdif_iters);
    println!("Rust TRDIF iters: {:?}", rust_trdif_iters);

    println!("TRDIF compare (rust vs XFOIL):");
    let mut tr_iters: Vec<usize> = rust_trdif.keys().cloned().collect();
    tr_iters.sort();
    for iter in tr_iters {
        let (r_vs1, r_vs2, r_vsr) = rust_trdif.get(&iter).unwrap();
        if let Some((x_vs1, x_vs2, x_vsr)) = xfoil_trdif.get(&iter) {
            let mut max_vs1: f64 = 0.0;
            let mut max_vs2: f64 = 0.0;
            for i in 0..r_vs1.len().min(x_vs1.len()) {
                for j in 0..r_vs1[i].len().min(x_vs1[i].len()) {
                    max_vs1 = max_vs1.max((r_vs1[i][j] - x_vs1[i][j]).abs());
                }
            }
            for i in 0..r_vs2.len().min(x_vs2.len()) {
                for j in 0..r_vs2[i].len().min(x_vs2[i].len()) {
                    max_vs2 = max_vs2.max((r_vs2[i][j] - x_vs2[i][j]).abs());
                }
            }
            let mut max_vsr: f64 = 0.0;
            for k in 0..r_vsr.len().min(x_vsr.len()) {
                max_vsr = max_vsr.max((r_vsr[k] - x_vsr[k]).abs());
            }
            println!(
                "  iter {:2}: max_vs1={:.3e} max_vs2={:.3e} max_vsr={:.3e}",
                iter, max_vs1, max_vs2, max_vsr
            );
        } else {
            println!("  iter {:2}: no XFOIL TRDIF", iter);
        }
    }

    println!("TRDIF derivs compare (rust vs XFOIL):");
    let Some(r_deriv) = rust_trdif_derivs.first() else {
        println!("  no Rust TRDIF_DERIVS");
        return;
    };
    let Some(x_deriv) = xfoil_trdif_derivs.first() else {
        println!("  no XFOIL TRDIF_DERIVS");
        return;
    };
    let keys = [
        "WF1", "WF2", "XT",
        "XT_A1", "XT_X1", "XT_X2", "XT_T1", "XT_T2", "XT_D1", "XT_D2", "XT_U1", "XT_U2", "XT_MS", "XT_RE",
        "TT_A1", "TT_X1", "TT_X2", "TT_T1", "TT_T2", "TT_D1", "TT_D2", "TT_U1", "TT_U2", "TT_MS", "TT_RE",
        "DT_A1", "DT_X1", "DT_X2", "DT_T1", "DT_T2", "DT_D1", "DT_D2", "DT_U1", "DT_U2", "DT_MS", "DT_RE",
        "UT_A1", "UT_X1", "UT_X2", "UT_T1", "UT_T2", "UT_D1", "UT_D2", "UT_U1", "UT_U2", "UT_MS", "UT_RE",
    ];
    for key in keys {
        let rv = r_deriv.get(key).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let xv = x_deriv.get(key).and_then(|v| v.as_f64()).unwrap_or(0.0);
        let diff = (rv - xv).abs();
        println!("  {:<6} diff={:.3e} rust={:.6e} xfoil={:.6e}", key, diff, rv, xv);
    }

    println!("Shape terms compare (rust BLDIF_TERMS vs XFOIL SHAPE_JACOBIAN):");
    let rust_shape_terms = rust_debug.events.iter().find(|ev| {
        ev.get("subroutine").and_then(|v| v.as_str()) == Some("BLDIF_TERMS")
            && ev.get("ibl").and_then(|v| v.as_u64()) == Some(30)
            && ev
                .get("iter")
                .and_then(|v| v.as_u64())
                .or_else(|| ev.get("iteration").and_then(|v| v.as_u64()))
                == Some(1)
            && ev.get("flow_type").and_then(|v| v.as_u64()) == Some(2)
    });
    let xfoil_shape_terms = xfoil_debug.events.iter().find(|ev| {
        ev.get("subroutine").and_then(|v| v.as_str()) == Some("SHAPE_JACOBIAN")
            && ev.get("ibl").and_then(|v| v.as_u64()) == Some(32)
            && ev.get("newton_iter").and_then(|v| v.as_u64()) == Some(1)
            && ev.get("flow_type").and_then(|v| v.as_u64()) == Some(2)
    });

    if let (Some(r_terms), Some(x_terms)) = (rust_shape_terms, xfoil_shape_terms) {
        let items = [
            ("z_di2", "Z_DI2"),
            ("di2_s", "DI2_S2"),
            ("z_upw_shape", "Z_UPW"),
            ("upw_t2_shape", "UPW_T2"),
            ("upw_d2_shape", "UPW_D2"),
            ("upw_u2_shape", "UPW_U2"),
        ];
        for (r_key, x_key) in items {
            let rv = r_terms.get(r_key).and_then(|v| v.as_f64()).unwrap_or(0.0);
            let xv = x_terms.get(x_key).and_then(|v| v.as_f64()).unwrap_or(0.0);
            let diff = (rv - xv).abs();
            println!("  {:<12} diff={:.3e} rust={:.6e} xfoil={:.6e}", r_key, diff, rv, xv);
        }
    } else {
        println!("  missing rust or xfoil shape terms event");
    }

    println!("BLVAR output compare (rust vs XFOIL):");
    let rust_blvar = rust_debug
        .events
        .iter()
        .filter(|ev| {
            ev.get("subroutine").and_then(|v| v.as_str()) == Some("BLVAR")
                && ev.get("ibl").and_then(|v| v.as_u64()) == Some(30)
                && ev.get("flow_type").and_then(|v| v.as_u64()) == Some(2)
        })
        .min_by_key(|ev| ev.get("iteration").and_then(|v| v.as_u64()).unwrap_or(u64::MAX));
    let xfoil_blvar = rust_blvar.and_then(|r_blvar| {
        let r_input = r_blvar.get("input")?;
        let r_u = r_input.get("u")?.as_f64().unwrap_or(0.0);
        let r_theta = r_input.get("theta")?.as_f64().unwrap_or(0.0);
        let r_delta = r_input.get("delta_star")?.as_f64().unwrap_or(0.0);
        let r_ctau = r_input.get("ctau")?.as_f64().unwrap_or(0.0);
        xfoil_debug
            .events
            .iter()
            .filter(|ev| {
                ev.get("subroutine").and_then(|v| v.as_str()) == Some("BLVAR")
                    && ev.get("flow_type").and_then(|v| v.as_u64()) == Some(2)
            })
            .min_by(|a, b| {
                let score = |ev: &serde_json::Value| {
                    let input = ev.get("input").unwrap_or(&serde_json::Value::Null);
                    let u = input.get("u").and_then(|v| v.as_f64()).unwrap_or(0.0);
                    let theta = input.get("theta").and_then(|v| v.as_f64()).unwrap_or(0.0);
                    let delta = input.get("delta_star").and_then(|v| v.as_f64()).unwrap_or(0.0);
                    let ctau = input.get("ctau").and_then(|v| v.as_f64()).unwrap_or(0.0);
                    (u - r_u).abs() + (theta - r_theta).abs() + (delta - r_delta).abs() + (ctau - r_ctau).abs()
                };
                score(a).partial_cmp(&score(b)).unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    if let (Some(r_blvar), Some(x_blvar)) = (rust_blvar, xfoil_blvar) {
        for key in ["H", "Hk", "Hs", "Us", "Cd", "Cq"] {
            let rv = r_blvar
                .get("output")
                .and_then(|o| o.get(key))
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let xv = x_blvar
                .get("output")
                .and_then(|o| o.get(key))
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let diff = (rv - xv).abs();
            println!("  {:<4} diff={:.3e} rust={:.6e} xfoil={:.6e}", key, diff, rv, xv);
        }
    } else {
        println!("  missing rust or xfoil BLVAR event");
    }

    println!("Theta-at-transition compare (TRCHEK2_FINAL tt):");
    let rust_tr = rust_debug.events.iter().find(|ev| {
        ev.get("subroutine").and_then(|v| v.as_str()) == Some("TRCHEK2_FINAL")
            && ev.get("ibl").and_then(|v| v.as_u64()) == Some(30)
            && ev.get("transition").and_then(|v| v.as_bool()) == Some(true)
    });
    let xfoil_tr = xfoil_debug.events.iter().find(|ev| {
        ev.get("subroutine").and_then(|v| v.as_str()) == Some("TRCHEK2_FINAL")
            && ev.get("ibl").and_then(|v| v.as_u64()) == Some(32)
            && ev.get("transition").and_then(|v| v.as_bool()) == Some(true)
    });
    if let (Some(r_tr), Some(x_tr)) = (rust_tr, xfoil_tr) {
        let r_tt = r_tr.get("tt").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let x_tt = x_tr.get("tt").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let diff = (r_tt - x_tt).abs();
        println!("  tt diff={:.3e} rust={:.6e} xfoil={:.6e}", diff, r_tt, x_tt);
    } else {
        println!("  missing rust or xfoil TRCHEK2_FINAL transition event");
    }

    println!("Theta-at-transition compare (TRCHEK2_ITER):");
    let rust_tr_iter = rust_debug
        .events
        .iter()
        .filter(|ev| {
            ev.get("subroutine").and_then(|v| v.as_str()) == Some("TRCHEK2_ITER")
                && ev.get("ibl").and_then(|v| v.as_u64()) == Some(30)
                && ev.get("transition").and_then(|v| v.as_bool()) == Some(true)
        })
        .max_by_key(|ev| ev.get("trchek_iter").and_then(|v| v.as_u64()).unwrap_or(0));
    let xfoil_tr_iter = xfoil_debug
        .events
        .iter()
        .filter(|ev| {
            ev.get("subroutine").and_then(|v| v.as_str()) == Some("TRCHEK2_ITER")
                && ev.get("ibl").and_then(|v| v.as_u64()) == Some(32)
                && ev.get("transition").and_then(|v| v.as_bool()) == Some(true)
        })
        .max_by_key(|ev| ev.get("trchek_iter").and_then(|v| v.as_u64()).unwrap_or(0));
    if let (Some(r_tr), Some(x_tr)) = (rust_tr_iter, xfoil_tr_iter) {
        let r_wf1 = r_tr.get("wf1").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let r_wf2 = r_tr.get("wf2").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let r_t1 = r_tr.get("T1").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let r_t2 = r_tr.get("T2").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let r_tt = r_wf1 * r_t1 + r_wf2 * r_t2;
        let x_wf1 = x_tr.get("wf1").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let x_wf2 = x_tr.get("wf2").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let x_t1 = x_tr.get("T1").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let x_t2 = x_tr.get("T2").and_then(|v| v.as_f64()).unwrap_or(0.0);
        let x_tt = x_wf1 * x_t1 + x_wf2 * x_t2;
        let diff = (r_tt - x_tt).abs();
        println!("  tt diff={:.3e} rust={:.6e} xfoil={:.6e}", diff, r_tt, x_tt);
    } else {
        println!("  missing rust or xfoil TRCHEK2_ITER transition event");
    }
}
