use serde::Deserialize;
use std::fs;

use rustfoil_bl::equations::{bldif_with_terms_ncrit, blvar, FlowType};
use rustfoil_bl::state::BlStation;

#[derive(Debug, Deserialize)]
struct DebugOutput {
    events: Vec<serde_json::Value>,
}

#[derive(Debug, Deserialize, Clone)]
struct BldifEvent {
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "flow_type")]
    flow_type: usize,
    #[serde(rename = "X1")]
    x1: f64,
    #[serde(rename = "U1")]
    u1: f64,
    #[serde(rename = "T1")]
    t1: f64,
    #[serde(rename = "D1")]
    d1: f64,
    #[serde(rename = "S1")]
    s1: f64,
    #[serde(rename = "A1")]
    a1: f64,
    #[serde(rename = "X2")]
    x2: f64,
    #[serde(rename = "U2")]
    u2: f64,
    #[serde(rename = "T2")]
    t2: f64,
    #[serde(rename = "D2")]
    d2: f64,
    #[serde(rename = "S2")]
    s2: f64,
    #[serde(rename = "A2")]
    a2: f64,
    #[serde(rename = "VS2")]
    vs2: Vec<Vec<f64>>,
    #[serde(rename = "VSREZ")]
    vsrez: Vec<f64>,
}

fn load_xfoil_debug() -> Option<DebugOutput> {
    let paths = [
        "xfoil_debug.json",
        "Xfoil-instrumented/bin/xfoil_debug.json",
        "../Xfoil-instrumented/bin/xfoil_debug.json",
        "../../Xfoil-instrumented/bin/xfoil_debug.json",
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

fn flow_type_from_xfoil(ityp: usize) -> Option<FlowType> {
    match ityp {
        1 => Some(FlowType::Laminar),
        2 => Some(FlowType::Turbulent),
        3 => Some(FlowType::Wake),
        _ => None,
    }
}

fn station_from_event(
    x: f64,
    u: f64,
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> BlStation {
    let mut s = BlStation::new();
    s.x = x;
    s.u = u;
    s.theta = theta;
    s.delta_star = delta_star;
    s.ctau = ctau;
    s.ampl = ampl;
    s.is_laminar = matches!(flow_type, FlowType::Laminar);
    s.is_turbulent = matches!(flow_type, FlowType::Turbulent | FlowType::Wake);
    s.is_wake = matches!(flow_type, FlowType::Wake);
    blvar(&mut s, flow_type, msq, re);
    s
}

#[ignore]
#[test]
fn probe_n4412_upper_bubble_bldif_contract() {
    let data = match load_xfoil_debug() {
        Some(data) => data,
        None => {
            eprintln!("Skipping: xfoil_debug.json not found");
            return;
        }
    };

    let re = 3.0e6;
    let msq = 0.0;
    let ncrit = 9.0;

    let mut events = Vec::new();
    for ev in &data.events {
        if ev.get("subroutine").and_then(|v| v.as_str()) != Some("BLDIF") {
            continue;
        }
        let parsed: BldifEvent = serde_json::from_value(ev.clone()).unwrap();
        if parsed.side == 1 && (23..=30).contains(&parsed.ibl) {
            events.push(parsed);
        }
    }

    events.sort_by_key(|ev| (ev.ibl, ev.newton_iter));

    for ev in events {
        let Some(flow_type) = flow_type_from_xfoil(ev.flow_type) else {
            continue;
        };

        let s1 = station_from_event(
            ev.x1, ev.u1, ev.t1, ev.d1, ev.s1, ev.a1, flow_type, msq, re,
        );
        let s2 = station_from_event(
            ev.x2, ev.u2, ev.t2, ev.d2, ev.s2, ev.a2, flow_type, msq, re,
        );

        let (res, jac, _) = bldif_with_terms_ncrit(&s1, &s2, flow_type, msq, re, ncrit);

        let res_vec = [res.res_third, res.res_mom, res.res_shape];
        let max_res = (0..3)
            .map(|k| (res_vec[k] - ev.vsrez[k]).abs())
            .fold(0.0_f64, f64::max);

        let max_vs2 = (0..3)
            .flat_map(|row| (0..4).map(move |col| (row, col)))
            .map(|(row, col)| (jac.vs2[row][col] - ev.vs2[row][col]).abs())
            .fold(0.0_f64, f64::max);

        println!(
            "ibl={} iter={} flow={:?} max_res={:.3e} max_vs2={:.3e} x2={:.6} hk2={:.3} cf2={:+.3e}",
            ev.ibl,
            ev.newton_iter,
            flow_type,
            max_res,
            max_vs2,
            s2.x,
            s2.hk,
            s2.cf,
        );
    }
}
