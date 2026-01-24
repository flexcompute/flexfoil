use serde::Deserialize;
use std::fs;

use rustfoil_bl::equations::{bldif_with_terms, blvar, FlowType};
use rustfoil_bl::state::BlStation;

#[derive(Debug, Deserialize)]
struct DebugOutput {
    events: Vec<serde_json::Value>,
}

#[derive(Debug, Deserialize, Clone)]
struct BldifEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
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

#[derive(Debug, Deserialize, Clone)]
struct ShapeJacobianEvent {
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "flow_type")]
    flow_type: usize,
    #[serde(rename = "Z_HS2")]
    z_hs2: f64,
    #[serde(rename = "Z_CF2")]
    z_cf2: f64,
    #[serde(rename = "Z_DI2")]
    z_di2: f64,
    #[serde(rename = "Z_T2")]
    z_t2: f64,
    #[serde(rename = "Z_U2")]
    z_u2: f64,
    #[serde(rename = "Z_HCA")]
    z_hca: f64,
    #[serde(rename = "Z_HA")]
    z_ha: f64,
    #[serde(rename = "Z_UPW")]
    z_upw: f64,
    #[serde(rename = "HS2_T2")]
    hs2_t2: f64,
    #[serde(rename = "HS2_D2")]
    hs2_d2: f64,
    #[serde(rename = "HS2_U2")]
    hs2_u2: f64,
    #[serde(rename = "CF2_T2")]
    cf2_t2: f64,
    #[serde(rename = "CF2_D2")]
    cf2_d2: f64,
    #[serde(rename = "CF2_U2")]
    cf2_u2: f64,
    #[serde(rename = "DI2_T2")]
    di2_t2: f64,
    #[serde(rename = "DI2_D2")]
    di2_d2: f64,
    #[serde(rename = "DI2_U2")]
    di2_u2: f64,
    #[serde(rename = "DI2_S2")]
    di2_s2: f64,
    #[serde(rename = "HC2_T2")]
    hc2_t2: f64,
    #[serde(rename = "HC2_D2")]
    hc2_d2: f64,
    #[serde(rename = "HC2_U2")]
    hc2_u2: f64,
    #[serde(rename = "H2_T2")]
    h2_t2: f64,
    #[serde(rename = "H2_D2")]
    h2_d2: f64,
    #[serde(rename = "UPW_T2")]
    upw_t2: f64,
    #[serde(rename = "UPW_D2")]
    upw_d2: f64,
    #[serde(rename = "UPW_U2")]
    upw_u2: f64,
    #[serde(rename = "VS2_3_1")]
    vs2_3_1: f64,
    #[serde(rename = "VS2_3_2")]
    vs2_3_2: f64,
    #[serde(rename = "VS2_3_3")]
    vs2_3_3: f64,
    #[serde(rename = "VS2_3_4")]
    vs2_3_4: f64,
}

fn load_xfoil_debug() -> Option<DebugOutput> {
    let paths = [
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

#[ignore]
#[test]
fn compare_transition_bldif_upper_ibl32() {
    let data = match load_xfoil_debug() {
        Some(data) => data,
        None => {
            eprintln!("Skipping: xfoil_debug.json not found");
            return;
        }
    };

    let re = 1.0e6;
    let msq = 0.0;

    let mut events = Vec::new();
    let mut shape_events = std::collections::HashMap::new();
    for ev in data.events.iter() {
        match ev.get("subroutine").and_then(|v| v.as_str()) {
            Some("BLDIF") => {
                if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
                    continue;
                }
                if ev.get("ibl").and_then(|v| v.as_u64()) != Some(32) {
                    continue;
                }
                let parsed: BldifEvent = serde_json::from_value(ev.clone()).unwrap();
                if parsed.flow_type == 2 {
                    events.push(parsed);
                }
            }
            Some("SHAPE_JACOBIAN") => {
                if ev.get("side").and_then(|v| v.as_u64()) != Some(1) {
                    continue;
                }
                if ev.get("ibl").and_then(|v| v.as_u64()) != Some(32) {
                    continue;
                }
                let parsed: ShapeJacobianEvent = serde_json::from_value(ev.clone()).unwrap();
                if parsed.flow_type == 2 {
                    shape_events.insert(parsed.newton_iter, parsed);
                }
            }
            _ => {}
        }
    }

    if events.is_empty() {
        eprintln!("No BLDIF events for side=1, ibl=32, flow_type=2");
        return;
    }

    events.sort_by_key(|e| e.newton_iter);

    for ev in events {
        let mut st1 = BlStation::new();
        st1.x = ev.x1;
        st1.u = ev.u1;
        st1.theta = ev.t1;
        st1.delta_star = ev.d1;
        st1.ctau = ev.s1;
        st1.ampl = ev.a1;
        st1.is_laminar = false;
        st1.is_turbulent = true;
        blvar(&mut st1, FlowType::Turbulent, msq, re);

        let mut st2 = BlStation::new();
        st2.x = ev.x2;
        st2.u = ev.u2;
        st2.theta = ev.t2;
        st2.delta_star = ev.d2;
        st2.ctau = ev.s2;
        st2.ampl = ev.a2;
        st2.is_laminar = false;
        st2.is_turbulent = true;
        blvar(&mut st2, FlowType::Turbulent, msq, re);

        let (res, jac, terms) = bldif_with_terms(&st1, &st2, FlowType::Turbulent, msq, re);

        let mut max_vs2: f64 = 0.0;
        let mut row_max = [0.0f64; 3];
        for i in 0..3 {
            let row = ev.vs2.get(i).expect("VS2 row missing");
            for j in 0..4 {
                let xfoil_val = *row.get(j).expect("VS2 col missing");
                let diff = (jac.vs2[i][j] - xfoil_val).abs();
                max_vs2 = max_vs2.max(diff);
                row_max[i] = row_max[i].max(diff);
            }
        }

        let res_vec = [res.res_third, res.res_mom, res.res_shape];
        let mut max_res: f64 = 0.0;
        for k in 0..3 {
            let xfoil_val = *ev.vsrez.get(k).expect("VSREZ entry missing");
            max_res = max_res.max((res_vec[k] - xfoil_val).abs());
        }

        println!(
            "iter {}: max_vs2={:.3e} row_max=[{:.3e},{:.3e},{:.3e}] max_res={:.3e} ctau2={:.6}",
            ev.newton_iter, max_vs2, row_max[0], row_max[1], row_max[2], max_res, ev.s2
        );

        if ev.newton_iter == 1 || ev.newton_iter == 5 {
            if let Some(shape) = shape_events.get(&ev.newton_iter) {
                println!(
                    "  SHAPE terms diff: z_hs2={:.3e} z_cf2={:.3e} z_di2={:.3e} z_t2={:.3e} z_u2={:.3e}",
                    (terms.z_hs2 - shape.z_hs2).abs(),
                    (terms.z_cf2_shape - shape.z_cf2).abs(),
                    (terms.z_di2 - shape.z_di2).abs(),
                    (terms.z_t2_shape - shape.z_t2).abs(),
                    (terms.z_u2_shape - shape.z_u2).abs()
                );
                println!(
                    "  SHAPE derivs diff: hs2_t={:.3e} hs2_d={:.3e} hs2_u={:.3e} cf2_t={:.3e} cf2_d={:.3e} cf2_u={:.3e}",
                    (terms.hs2_t - shape.hs2_t2).abs(),
                    (terms.hs2_d - shape.hs2_d2).abs(),
                    (terms.hs2_u - shape.hs2_u2).abs(),
                    (terms.cf2_t_shape - shape.cf2_t2).abs(),
                    (terms.cf2_d_shape - shape.cf2_d2).abs(),
                    (terms.cf2_u_shape - shape.cf2_u2).abs()
                );
                println!(
                    "  SHAPE di diff: di2_t={:.3e} di2_d={:.3e} di2_u={:.3e} di2_s={:.3e}",
                    (terms.di2_t - shape.di2_t2).abs(),
                    (terms.di2_d - shape.di2_d2).abs(),
                    (terms.di2_u - shape.di2_u2).abs(),
                    (terms.di2_s - shape.di2_s2).abs()
                );
                println!(
                    "  SHAPE h/hc/upw diff: h2_t={:.3e} h2_d={:.3e} hc2_t={:.3e} hc2_d={:.3e} upw_t2={:.3e} upw_d2={:.3e}",
                    (terms.h2_t - shape.h2_t2).abs(),
                    (terms.h2_d - shape.h2_d2).abs(),
                    (terms.hc2_t - shape.hc2_t2).abs(),
                    (terms.hc2_d - shape.hc2_d2).abs(),
                    (terms.upw_t2_shape - shape.upw_t2).abs(),
                    (terms.upw_d2_shape - shape.upw_d2).abs()
                );
                println!(
                    "  SHAPE vs2 diff: v31={:.3e} v32={:.3e} v33={:.3e} v34={:.3e}",
                    (jac.vs2[2][0] - shape.vs2_3_1).abs(),
                    (jac.vs2[2][1] - shape.vs2_3_2).abs(),
                    (jac.vs2[2][2] - shape.vs2_3_3).abs(),
                    (jac.vs2[2][3] - shape.vs2_3_4).abs()
                );
                println!(
                    "  SHAPE vs2 values: rust=[{:.3e},{:.3e},{:.3e},{:.3e}] xfoil=[{:.3e},{:.3e},{:.3e},{:.3e}]",
                    jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2], jac.vs2[2][3],
                    shape.vs2_3_1, shape.vs2_3_2, shape.vs2_3_3, shape.vs2_3_4
                );
            } else {
                println!("  SHAPE event missing for iter {}", ev.newton_iter);
            }
        }
    }
}
