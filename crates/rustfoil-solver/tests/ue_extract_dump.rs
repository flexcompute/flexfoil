use rustfoil_solver::viscous::{extract_surface_xfoil, setup_from_coords};
use serde::Serialize;
use std::fs;
use std::path::PathBuf;

#[derive(Serialize)]
struct SurfaceDump {
    arc: Vec<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    ue: Vec<f64>,
}

#[derive(Serialize)]
struct ExtractDump {
    alpha_deg: f64,
    ist: usize,
    sst: f64,
    ue_stag: f64,
    upper: SurfaceDump,
    lower: SurfaceDump,
}

fn load_xfoil_dat(filename: &str) -> Vec<(f64, f64)> {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop(); // crates/
    path.pop(); // project root
    path.push(filename);

    let content = fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Failed to read file: {:?}", path));

    let mut points = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                points.push((x, y));
            }
        }
    }

    points
}

#[test]
fn dump_extracted_surface_ue() {
    let coords = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let alpha_deg = 4.0;
    let setup_result = setup_from_coords(&coords, alpha_deg)
        .expect("setup_from_coords failed");

    let node_x = &setup_result.node_x;
    let node_y = &setup_result.node_y;
    let ue_inviscid = &setup_result.setup.ue_inviscid;
    let full_arc = &setup_result.setup.arc_lengths;
    let ist = setup_result.ist;
    let sst = setup_result.sst;

    // Interpolate Ue at stagnation (same logic as run_viscous_analysis)
    let ue_stag = if ist + 1 < ue_inviscid.len() && full_arc[ist + 1] != full_arc[ist] {
        let frac = (sst - full_arc[ist]) / (full_arc[ist + 1] - full_arc[ist]);
        ue_inviscid[ist] + frac * (ue_inviscid[ist + 1] - ue_inviscid[ist])
    } else {
        ue_inviscid[ist]
    };

    let (upper_arc, upper_x, upper_y, upper_ue) =
        extract_surface_xfoil(ist, sst, ue_stag, full_arc, node_x, node_y, ue_inviscid, true);
    let (lower_arc, lower_x, lower_y, lower_ue) =
        extract_surface_xfoil(ist, sst, ue_stag, full_arc, node_x, node_y, ue_inviscid, false);

    let dump = ExtractDump {
        alpha_deg,
        ist,
        sst,
        ue_stag,
        upper: SurfaceDump {
            arc: upper_arc,
            x: upper_x,
            y: upper_y,
            ue: upper_ue,
        },
        lower: SurfaceDump {
            arc: lower_arc,
            x: lower_x,
            y: lower_y,
            ue: lower_ue,
        },
    };

    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let out_dir = format!("{}/target", manifest_dir);
    let _ = fs::create_dir_all(&out_dir);
    let out_path = format!("{}/ue_extract.json", out_dir);
    fs::write(&out_path, serde_json::to_string_pretty(&dump).unwrap())
        .expect("Failed to write ue_extract.json");

    println!("Wrote extracted surface Ue to {}", out_path);
}
