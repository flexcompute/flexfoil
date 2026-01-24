use serde::Deserialize;

use rustfoil_bl::equations::{bldif_with_terms, blvar, FlowType};
use rustfoil_bl::state::BlStation;

#[derive(Debug, Deserialize)]
struct MrchueData {
    metadata: Metadata,
    sides: std::collections::HashMap<String, SideData>,
}

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
}

#[derive(Debug, Deserialize)]
struct SideData {
    stations: Vec<StationData>,
}

#[derive(Debug, Deserialize)]
struct StationData {
    side: usize,
    ibl: usize,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    initial: StationState,
    #[serde(rename = "final")]
    final_state: StationFinal,
    iterations: Vec<IterationData>,
}

#[derive(Debug, Deserialize)]
struct StationState {
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
}

#[derive(Debug, Deserialize)]
struct StationFinal {
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
    #[allow(dead_code)]
    dmax: f64,
    #[allow(dead_code)]
    #[serde(rename = "Hk")]
    hk: f64,
    #[allow(dead_code)]
    #[serde(rename = "Rtheta")]
    rtheta: f64,
}

#[derive(Debug, Deserialize)]
struct IterationData {
    iter: usize,
    theta_in: f64,
    delta_star_in: f64,
    theta_out: f64,
    delta_star_out: f64,
    vsrez: Vec<f64>,
    #[serde(rename = "VS2")]
    vs2: Vec<Vec<f64>>,
    #[allow(dead_code)]
    dmax: f64,
    #[allow(dead_code)]
    relaxation: f64,
    #[allow(dead_code)]
    converged: bool,
}

fn build_station(
    x: f64,
    ue: f64,
    state: &StationState,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> BlStation {
    let mut st = BlStation::new();
    st.x = x;
    st.u = ue;
    st.theta = state.theta;
    st.delta_star = state.delta_star;
    st.ctau = state.ctau;
    st.ampl = state.ampl;
    st.is_laminar = flow_type == FlowType::Laminar;
    st.is_turbulent = flow_type == FlowType::Turbulent;
    st.is_wake = flow_type == FlowType::Wake;
    blvar(&mut st, flow_type, msq, re);
    st
}

fn build_station_iter(
    x: f64,
    ue: f64,
    theta: f64,
    delta_star: f64,
    base: &StationState,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> BlStation {
    let mut st = BlStation::new();
    st.x = x;
    st.u = ue;
    st.theta = theta;
    st.delta_star = delta_star;
    st.ctau = base.ctau;
    st.ampl = base.ampl;
    st.is_laminar = flow_type == FlowType::Laminar;
    st.is_turbulent = flow_type == FlowType::Turbulent;
    st.is_wake = flow_type == FlowType::Wake;
    blvar(&mut st, flow_type, msq, re);
    st
}

fn flatten_vs2(vs2: &[[f64; 5]; 3]) -> Vec<f64> {
    vs2.iter().flat_map(|row| row.iter().copied()).collect()
}

#[test]
fn compare_mrchue_vs2_iterations_lower_surface() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/mrchue_iterations.json");
    let raw = std::fs::read_to_string(&path).expect("testdata/mrchue_iterations.json missing");
    let data: MrchueData = serde_json::from_str(&raw).expect("json parse failed");

    let msq = data.metadata.mach * data.metadata.mach;
    let re = data.metadata.reynolds;
    let flow_type = FlowType::Laminar;

    let side = data
        .sides
        .get("2")
        .expect("side 2 missing in testdata");

    let stations_by_ibl: std::collections::HashMap<usize, &StationData> = side
        .stations
        .iter()
        .map(|s| (s.ibl, s))
        .collect();

    let focus_ibls = [64usize, 65, 66, 67];

    for ibl in focus_ibls {
        let s2 = stations_by_ibl.get(&ibl).expect("station missing");
        let s1 = stations_by_ibl
            .get(&(ibl - 1))
            .expect("prev station missing");

        let st1_state = StationState::from(&s1.final_state);
        let st1 = build_station(s1.x, s1.ue, &st1_state, flow_type, msq, re);

        println!("\n== Side 2 IBL {} ==", ibl);

        for iter in &s2.iterations {
            let st2 = build_station_iter(
                s2.x,
                s2.ue,
                iter.theta_in,
                iter.delta_star_in,
                &s2.initial,
                flow_type,
                msq,
                re,
            );

            let (_res, jac, _terms) = bldif_with_terms(&st1, &st2, flow_type, msq, re);
            let rust_vs2 = flatten_vs2(&jac.vs2);
            let xf_vs2: Vec<f64> = iter
                .vs2
                .iter()
                .take(3)
                .flat_map(|row| row.iter())
                .copied()
                .collect();

            let mut max_abs = 0.0f64;
            for (a, b) in rust_vs2.iter().zip(xf_vs2.iter()) {
                max_abs = max_abs.max((a - b).abs());
            }

            println!(
                "iter {:>2}: max |VS2 diff| = {:.6e} (rust row2 {:?}, xf row2 {:?})",
                iter.iter,
                max_abs,
                jac.vs2[2],
                iter.vs2.get(2).unwrap_or(&vec![]),
            );
        }
    }
}

// Helper to reuse final state struct as StationState where needed.
impl From<&StationFinal> for StationState {
    fn from(value: &StationFinal) -> Self {
        Self {
            theta: value.theta,
            delta_star: value.delta_star,
            ctau: value.ctau,
            ampl: value.ampl,
        }
    }
}
