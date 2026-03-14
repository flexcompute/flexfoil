use serde::Deserialize;

use rustfoil_bl::equations::{bldif_with_terms, blvar, FlowType};
use rustfoil_bl::state::BlStation;

#[derive(Debug, Deserialize)]
struct DebugOutput {
    events: Vec<serde_json::Value>,
}

#[derive(Debug, Deserialize, Clone)]
struct BlprvEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "T2")]
    t2: f64,
    #[serde(rename = "D2")]
    d2: f64,
}

#[derive(Debug, Deserialize, Clone)]
struct ShearJacobianEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "Z_UPW")]
    z_upw: f64,
    #[serde(rename = "Z_DE2")]
    z_de2: f64,
    #[serde(rename = "Z_US2")]
    z_us2: f64,
    #[serde(rename = "Z_CQ2")]
    z_cq2: f64,
    #[serde(rename = "Z_CF2")]
    z_cf2: f64,
    #[serde(rename = "Z_HK2")]
    z_hk2: f64,
    #[serde(rename = "Z_D2")]
    z_d2: f64,
    #[serde(rename = "Z_U2")]
    z_u2: f64,
    #[serde(rename = "Z_S2")]
    z_s2: f64,
    #[serde(rename = "UPW_T2")]
    upw_t2: f64,
    #[serde(rename = "UPW_D2")]
    upw_d2: f64,
    #[serde(rename = "UPW_U2")]
    upw_u2: f64,
    #[serde(rename = "DE2_T2")]
    de2_t2: f64,
    #[serde(rename = "DE2_D2")]
    de2_d2: f64,
    #[serde(rename = "DE2_U2")]
    de2_u2: f64,
    #[serde(rename = "US2_T2")]
    us2_t2: f64,
    #[serde(rename = "US2_D2")]
    us2_d2: f64,
    #[serde(rename = "US2_U2")]
    us2_u2: f64,
    #[serde(rename = "CQ2_T2")]
    cq2_t2: f64,
    #[serde(rename = "CQ2_D2")]
    cq2_d2: f64,
    #[serde(rename = "CQ2_U2")]
    cq2_u2: f64,
    #[serde(rename = "CF2_T2")]
    cf2_t2: f64,
    #[serde(rename = "CF2_D2")]
    cf2_d2: f64,
    #[serde(rename = "CF2_U2")]
    cf2_u2: f64,
    #[serde(rename = "HK2_T2")]
    hk2_t2: f64,
    #[serde(rename = "HK2_D2")]
    hk2_d2: f64,
    #[serde(rename = "HK2_U2")]
    hk2_u2: f64,
    #[serde(rename = "VS2_1_2")]
    vs2_1_2: f64,
    #[serde(rename = "VS2_1_3")]
    vs2_1_3: f64,
    #[serde(rename = "VS2_1_4")]
    vs2_1_4: f64,
}

#[derive(Debug, Deserialize, Clone)]
struct LaminarJacobianEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "AX")]
    ax: f64,
    #[serde(rename = "AX_HK2")]
    ax_hk2: f64,
    #[serde(rename = "AX_T2")]
    ax_t2: f64,
    #[serde(rename = "AX_RT2")]
    ax_rt2: f64,
    #[serde(rename = "AX_A2")]
    ax_a2: f64,
    #[serde(rename = "HK2_T2")]
    hk2_t2: f64,
    #[serde(rename = "HK2_D2")]
    hk2_d2: f64,
    #[serde(rename = "HK2_U2")]
    hk2_u2: f64,
    #[serde(rename = "RT2_T2")]
    rt2_t2: f64,
    #[serde(rename = "RT2_U2")]
    rt2_u2: f64,
    #[serde(rename = "Z_AX")]
    z_ax: f64,
    #[serde(rename = "VS2_1_2")]
    vs2_1_2: f64,
    #[serde(rename = "VS2_1_3")]
    vs2_1_3: f64,
    #[serde(rename = "VS2_1_4")]
    vs2_1_4: f64,
}

#[derive(Debug, Deserialize, Clone)]
struct ShapeJacobianEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
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

#[derive(Debug, Deserialize, Clone)]
struct BldifEvent {
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
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
}

#[derive(Debug, Deserialize, Clone)]
struct MomJacobianEvent {
    #[serde(rename = "call_id")]
    call_id: usize,
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    #[serde(default)]
    newton_iter: Option<usize>,
    #[serde(rename = "T2")]
    t2: f64,
    #[serde(rename = "D2")]
    d2: f64,
    #[serde(rename = "X1")]
    x1: f64,
    #[serde(rename = "U1")]
    u1: f64,
    #[serde(rename = "T1")]
    t1: f64,
    #[serde(rename = "D1")]
    d1: f64,
    #[serde(rename = "X2")]
    x2: f64,
    #[serde(rename = "U2")]
    u2: f64,
    #[serde(rename = "Z_HA")]
    z_ha: f64,
    #[serde(rename = "Z_CFM")]
    z_cfm: f64,
    #[serde(rename = "Z_CF2")]
    z_cf2: f64,
    #[serde(rename = "Z_T2")]
    z_t2: f64,
    #[serde(rename = "Z_U2")]
    z_u2: f64,
    #[serde(rename = "H2_T2")]
    h2_t2: f64,
    #[serde(rename = "H2_D2")]
    h2_d2: f64,
    #[serde(rename = "CFM_T2")]
    cfm_t2: f64,
    #[serde(rename = "CFM_D2")]
    cfm_d2: f64,
    #[serde(rename = "CFM_U2")]
    cfm_u2: f64,
    #[serde(rename = "CF2_T2")]
    cf2_t2: f64,
    #[serde(rename = "CF2_D2")]
    cf2_d2: f64,
    #[serde(rename = "CF2_U2")]
    cf2_u2: f64,
    #[serde(rename = "VS2_2_2")]
    vs2_2_2: f64,
    #[serde(rename = "VS2_2_3")]
    vs2_2_3: f64,
    #[serde(rename = "VS2_2_4")]
    vs2_2_4: f64,
}

#[derive(Debug, Deserialize, Clone)]
struct Vs2BeforeEvent {
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    #[serde(rename = "VS2_4x4")]
    vs2_4x4: Vec<Vec<f64>>,
}

#[derive(Debug, Deserialize, Clone)]
struct MrchueIterEvent {
    side: usize,
    ibl: usize,
    #[serde(rename = "newton_iter")]
    newton_iter: usize,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    #[serde(rename = "d_s_wake")]
    d_s_wake: f64,
    input: StationState,
    output: StationStateOut,
    #[serde(rename = "VS2")]
    #[allow(dead_code)]
    vs2: Vec<Vec<f64>>,
}

#[derive(Debug, Deserialize, Clone)]
struct StationState {
    theta: f64,
    #[serde(rename = "delta_star")]
    delta_star: f64,
    ctau: f64,
    ampl: f64,
}

#[derive(Debug, Deserialize, Clone)]
struct StationStateOut {
    theta: f64,
    #[serde(rename = "delta_star")]
    delta_star: f64,
    ctau: f64,
    ampl: f64,
    #[serde(rename = "Hk")]
    hk: f64,
    #[serde(rename = "Rtheta")]
    rtheta: f64,
}

fn build_station(
    x: f64,
    ue: f64,
    state: &StationState,
    delta_star_override: Option<f64>,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> BlStation {
    let mut st = BlStation::new();
    st.x = x;
    st.u = ue;
    st.theta = state.theta;
    st.delta_star = delta_star_override.unwrap_or(state.delta_star);
    st.ctau = state.ctau;
    st.ampl = state.ampl;
    st.is_laminar = flow_type == FlowType::Laminar;
    st.is_turbulent = flow_type == FlowType::Turbulent;
    st.is_wake = flow_type == FlowType::Wake;
    blvar(&mut st, flow_type, msq, re);
    st
}

fn build_station_out(
    x: f64,
    ue: f64,
    state: &StationStateOut,
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

fn flatten_vs2_4x4(vs2: &[[f64; 5]; 3]) -> Vec<f64> {
    let mut out = Vec::with_capacity(12);
    for row in vs2.iter() {
        out.extend_from_slice(&row[0..4]);
    }
    out
}

#[derive(Debug)]
struct RustMomTerms {
    z_ha: f64,
    z_cfm: f64,
    z_cf2: f64,
    z_t2: f64,
    z_u2: f64,
    h2_t2: f64,
    h2_d2: f64,
    cfm_t2: f64,
    cfm_d2: f64,
    cfm_u2: f64,
    cf2_t2: f64,
    cf2_d2: f64,
    cf2_u2: f64,
    vs2_2_2: f64,
    vs2_2_3: f64,
    vs2_2_4: f64,
}

fn compute_rust_mom_terms(
    s1: &BlStation,
    s2: &BlStation,
    flow_type: FlowType,
    msq: f64,
) -> RustMomTerms {
    use rustfoil_bl::closures::{cf_laminar, cf_turbulent};

    let xlog = (s2.x / s1.x.max(1e-20)).ln();
    let ulog = (s2.u / s1.u.max(1e-20)).ln();
    let ha = 0.5 * (s1.h + s2.h);
    let ma_avg = 0.5 * msq;
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);

    let hka = 0.5 * (s1.hk + s2.hk);
    let rta = 0.5 * (s1.r_theta + s2.r_theta);
    let ma = 0.5 * msq;

    let (cfm, cfm_hka, cfm_rta) = match flow_type {
        FlowType::Wake => (0.0, 0.0, 0.0),
        FlowType::Laminar => {
            let cf = cf_laminar(hka, rta, ma);
            (cf.cf, cf.cf_hk, cf.cf_rt)
        }
        FlowType::Turbulent => {
            let cf_t = cf_turbulent(hka, rta, ma);
            let cf_l = cf_laminar(hka, rta, ma);
            if cf_l.cf > cf_t.cf {
                (cf_l.cf, cf_l.cf_hk, cf_l.cf_rt)
            } else {
                (cf_t.cf, cf_t.cf_hk, cf_t.cf_rt)
            }
        }
    };

    let cfx = 0.5 * cfm * xa / ta + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);
    let cfx_cfm = 0.5 * xa / ta;
    let cfx_cf2 = 0.25 * s2.x / s2.theta;
    let btmp = ha + 2.0 - ma_avg;

    let z_cfx = -xlog * 0.5;
    let z_ha = ulog;
    let z_cfm = z_cfx * cfx_cfm;
    let z_cf2 = z_cfx * cfx_cf2;

    let cfx_ta = -0.5 * cfm * xa / (ta * ta);
    let cfx_t2 = -0.25 * s2.cf * s2.x / (s2.theta * s2.theta) + cfx_ta * 0.5;
    let z_t2 = 1.0 / s2.theta + z_cfx * cfx_t2;
    let z_u2 = btmp / s2.u;

    let h2_t2 = s2.derivs.h_theta;
    let h2_d2 = s2.derivs.h_delta_star;
    let hk2_t = s2.derivs.hk_h * h2_t2;
    let hk2_d = s2.derivs.hk_h * h2_d2;
    let rt2_t = s2.r_theta / s2.theta.max(1e-12);
    let rt2_u = s2.r_theta / s2.u.max(1e-12);

    let cfm_hk2 = 0.5 * cfm_hka;
    let cfm_rt2 = 0.5 * cfm_rta;
    let cfm_t2 = cfm_hk2 * hk2_t + cfm_rt2 * rt2_t;
    let cfm_d2 = cfm_hk2 * hk2_d;
    let cfm_u2 = cfm_rt2 * rt2_u;

    let cf2_t2 = s2.derivs.cf_hk * hk2_t + s2.derivs.cf_rt * rt2_t;
    let cf2_d2 = s2.derivs.cf_hk * hk2_d;
    let cf2_u2 = s2.derivs.cf_rt * rt2_u;

    let vs2_2_2 = 0.5 * z_ha * h2_t2 + z_cfm * cfm_t2 + z_cf2 * cf2_t2 + z_t2;
    let vs2_2_3 = 0.5 * z_ha * h2_d2 + z_cfm * cfm_d2 + z_cf2 * cf2_d2;
    let vs2_2_4 = z_cfm * cfm_u2 + z_cf2 * cf2_u2 + z_u2;

    RustMomTerms {
        z_ha,
        z_cfm,
        z_cf2,
        z_t2,
        z_u2,
        h2_t2,
        h2_d2,
        cfm_t2,
        cfm_d2,
        cfm_u2,
        cf2_t2,
        cf2_d2,
        cf2_u2,
        vs2_2_2,
        vs2_2_3,
        vs2_2_4,
    }
}

#[test]
fn compare_xfoil_mrchue_iter_vs2_lower_surface() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../Xfoil-instrumented/bin/xfoil_debug.json");
    let raw = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(_) => {
            eprintln!("Skipping test: xfoil_debug.json not found at {:?}", path);
            return;
        }
    };
    let data: DebugOutput = serde_json::from_str(&raw).expect("json parse failed");

    let msq = 0.0;
    let re = 1.0e6;
    let flow_type = FlowType::Laminar;

    let mut iter_events: Vec<MrchueIterEvent> = Vec::new();
    let mut vs2_before_events: Vec<Vs2BeforeEvent> = Vec::new();
    let mut mom_events: Vec<MomJacobianEvent> = Vec::new();
    let mut blprv_events: Vec<BlprvEvent> = Vec::new();
    let mut shear_events: Vec<ShearJacobianEvent> = Vec::new();
    let mut laminar_events: Vec<LaminarJacobianEvent> = Vec::new();
    let mut shape_events: Vec<ShapeJacobianEvent> = Vec::new();
    let mut bldif_events: Vec<BldifEvent> = Vec::new();
    for event in data.events {
        let subroutine = event
            .get("subroutine")
            .and_then(|v| v.as_str())
            .unwrap_or("");
        match subroutine {
            "MRCHUE_ITER" => {
                let parsed: MrchueIterEvent =
                    serde_json::from_value(event).expect("parse event");
                if parsed.side == 2 {
                    iter_events.push(parsed);
                }
            }
            "VS2_BEFORE" => {
                let parsed: Vs2BeforeEvent =
                    serde_json::from_value(event).expect("parse vs2_before");
                if parsed.side == 2 {
                    vs2_before_events.push(parsed);
                }
            }
            "MOM_JACOBIAN" => {
                let parsed: MomJacobianEvent =
                    serde_json::from_value(event).expect("parse mom_jacobian");
                if parsed.side == 2 {
                    mom_events.push(parsed);
                }
            }
            "BLPRV_STATE" => {
                let parsed: BlprvEvent =
                    serde_json::from_value(event).expect("parse blprv_state");
                if parsed.side == 2 {
                    blprv_events.push(parsed);
                }
            }
            "SHEAR_JACOBIAN" => {
                let parsed: ShearJacobianEvent =
                    serde_json::from_value(event).expect("parse shear_jacobian");
                if parsed.side == 2 {
                    shear_events.push(parsed);
                }
            }
            "LAMINAR_JACOBIAN" => {
                let parsed: LaminarJacobianEvent =
                    serde_json::from_value(event).expect("parse laminar_jacobian");
                if parsed.side == 2 {
                    laminar_events.push(parsed);
                }
            }
            "SHAPE_JACOBIAN" => {
                let parsed: ShapeJacobianEvent =
                    serde_json::from_value(event).expect("parse shape_jacobian");
                if parsed.side == 2 {
                    shape_events.push(parsed);
                }
            }
            "BLDIF" => {
                let parsed: BldifEvent = serde_json::from_value(event).expect("parse bldif");
                if parsed.side == 2 {
                    bldif_events.push(parsed);
                }
            }
            _ => {}
        }
    }

    // Build map of final (converged) outputs per station (max newton_iter)
    let mut final_by_ibl: std::collections::HashMap<usize, MrchueIterEvent> = std::collections::HashMap::new();
    for ev in &iter_events {
        final_by_ibl
            .entry(ev.ibl)
            .and_modify(|cur| {
                if ev.newton_iter > cur.newton_iter {
                    *cur = ev.clone();
                }
            })
            .or_insert_with(|| ev.clone());
    }

    let focus_ibls = [64usize, 65, 66, 67];

    let mut evs_by_ibl: std::collections::HashMap<usize, Vec<&MrchueIterEvent>> =
        std::collections::HashMap::new();
    for ev in &iter_events {
        evs_by_ibl.entry(ev.ibl).or_default().push(ev);
    }
    for evs in evs_by_ibl.values_mut() {
        evs.sort_by_key(|e| e.newton_iter);
    }

    let mut moms_by_ibl: std::collections::HashMap<usize, Vec<&MomJacobianEvent>> =
        std::collections::HashMap::new();
    for mom in &mom_events {
        moms_by_ibl.entry(mom.ibl).or_default().push(mom);
    }
    for moms in moms_by_ibl.values_mut() {
        moms.sort_by_key(|e| e.newton_iter.unwrap_or(0));
    }

    let mut blprv_by_key: std::collections::HashMap<(usize, usize), BlprvEvent> =
        std::collections::HashMap::new();
    for blp in &blprv_events {
        blprv_by_key.insert((blp.ibl, blp.newton_iter), blp.clone());
    }

    let mut shear_by_key: std::collections::HashMap<(usize, usize), ShearJacobianEvent> =
        std::collections::HashMap::new();
    for shear in &shear_events {
        shear_by_key.insert((shear.ibl, shear.newton_iter), shear.clone());
    }

    let mut laminar_by_key: std::collections::HashMap<(usize, usize), LaminarJacobianEvent> =
        std::collections::HashMap::new();
    for lam in &laminar_events {
        laminar_by_key.insert((lam.ibl, lam.newton_iter), lam.clone());
    }

    let mut shape_by_key: std::collections::HashMap<(usize, usize), ShapeJacobianEvent> =
        std::collections::HashMap::new();
    for shape in &shape_events {
        shape_by_key.insert((shape.ibl, shape.newton_iter), shape.clone());
    }

    let mut bldif_by_key: std::collections::HashMap<(usize, usize), BldifEvent> =
        std::collections::HashMap::new();
    for bld in &bldif_events {
        bldif_by_key.insert((bld.ibl, bld.newton_iter), bld.clone());
    }

    let mut d2_final_by_ibl: std::collections::HashMap<usize, f64> =
        std::collections::HashMap::new();
    for (ibl, evs) in &evs_by_ibl {
        if let Some(ev_last) = evs.last() {
            let d2_final = ev_last.input.delta_star - ev_last.d_s_wake;
            d2_final_by_ibl.insert(*ibl, d2_final);
        }
    }

    for ibl in focus_ibls {
        let prev = final_by_ibl
            .get(&(ibl - 1))
            .expect("prev station missing");

        println!("\n== Side 2 IBL {} ==", ibl);

        let evs = evs_by_ibl.get(&ibl).cloned().unwrap_or_default();
        let moms = moms_by_ibl.get(&ibl).cloned().unwrap_or_default();

        for ev in evs.iter() {
            let bldif = bldif_by_key
                .get(&(ibl - 1, ev.newton_iter))
                .expect("BLDIF missing for iteration");
            let mom = match moms
                .iter()
                .find(|m| m.newton_iter == Some(ev.newton_iter))
            {
                Some(m) => *m,
                None => {
                    println!("  MOM terms: missing for iter {}", ev.newton_iter);
                    continue;
                }
            };
            let st2 = build_station(
                mom.x2,
                mom.u2,
                &StationState {
                    theta: mom.t2,
                    delta_star: mom.d2,
                    ctau: ev.input.ctau,
                    ampl: ev.input.ampl,
                },
                Some(mom.d2),
                flow_type,
                msq,
                re,
            );
            let st1 = build_station(
                mom.x1,
                mom.u1,
                &StationState {
                    theta: mom.t1,
                    delta_star: mom.d1,
                    ctau: ev.input.ctau,
                    ampl: ev.input.ampl,
                },
                Some(mom.d1),
                flow_type,
                msq,
                re,
            );
            let (_res, jac, terms) = bldif_with_terms(&st1, &st2, flow_type, msq, re);

            let st2_bldif = build_station(
                bldif.x2,
                bldif.u2,
                &StationState {
                    theta: bldif.t2,
                    delta_star: bldif.d2,
                    ctau: bldif.s2,
                    ampl: bldif.a2,
                },
                Some(bldif.d2),
                flow_type,
                msq,
                re,
            );
            let st1_bldif = build_station(
                bldif.x1,
                bldif.u1,
                &StationState {
                    theta: bldif.t1,
                    delta_star: bldif.d1,
                    ctau: bldif.s1,
                    ampl: bldif.a1,
                },
                Some(bldif.d1),
                flow_type,
                msq,
                re,
            );
            let (_res_b, jac_bldif, terms_bldif) =
                bldif_with_terms(&st1_bldif, &st2_bldif, flow_type, msq, re);
            let rust_mom = compute_rust_mom_terms(&st1, &st2, flow_type, msq);

            let rust_vs2 = flatten_vs2_4x4(&jac_bldif.vs2);
            let xf_vs2: Vec<f64> = bldif
                .vs2
                .iter()
                .take(3)
                .flat_map(|row| row.iter().take(4))
                .copied()
                .collect();

            let mut max_abs = 0.0f64;
            let mut max_idx = 0usize;
            for (idx, (a, b)) in rust_vs2.iter().zip(xf_vs2.iter()).enumerate() {
                let diff = (a - b).abs();
                if diff > max_abs {
                    max_abs = diff;
                    max_idx = idx;
                }
            }
            let max_row = max_idx / 4;
            let max_col = max_idx % 4;

            println!(
                "iter {:>2}: max |VS2 diff| = {:.6e} at row {} col {} (rust row2 {:?}, xf row2 {:?})",
                ev.newton_iter,
                max_abs,
                max_row + 1,
                max_col + 1,
                jac_bldif.vs2[2],
                bldif.vs2.get(2).unwrap_or(&vec![]),
            );
            println!(
                "  MOM terms: Z_T2 {:.6e}/{:.6e} H2_T2 {:.6e}/{:.6e} CFM_T2 {:.6e}/{:.6e} CF2_T2 {:.6e}/{:.6e} VS2_2_2 {:.6e}/{:.6e}",
                rust_mom.z_t2,
                mom.z_t2,
                rust_mom.h2_t2,
                mom.h2_t2,
                rust_mom.cfm_t2,
                mom.cfm_t2,
                rust_mom.cf2_t2,
                mom.cf2_t2,
                rust_mom.vs2_2_2,
                mom.vs2_2_2,
            );
            if let Some(row0) = bldif.vs2.get(0) {
                if let Some(xf_12) = row0.get(1) {
                    println!(
                        "  VS2(1,2) rust/xfoil {:.6e}/{:.6e}",
                        jac_bldif.vs2[0][1], xf_12
                    );
                }
            }
            if let Some(lam) = laminar_by_key.get(&(ibl - 1, ev.newton_iter)) {
                let rust_vs2_12 = terms_bldif.z_ax
                    * (terms_bldif.ax_hk2 * terms_bldif.hk2_t2
                        + terms_bldif.ax_t2
                        + terms_bldif.ax_rt2 * terms_bldif.rt2_t2);
                let rust_vs2_13 = terms_bldif.z_ax * (terms_bldif.ax_hk2 * terms_bldif.hk2_d2);
                let rust_vs2_14 = terms_bldif.z_ax
                    * (terms_bldif.ax_hk2 * terms_bldif.hk2_u2
                        + terms_bldif.ax_rt2 * terms_bldif.rt2_u2);

                let xfoil_vs2_12 =
                    lam.z_ax * (lam.ax_hk2 * lam.hk2_t2 + lam.ax_t2 + lam.ax_rt2 * lam.rt2_t2);
                let xfoil_vs2_13 = lam.z_ax * (lam.ax_hk2 * lam.hk2_d2);
                let xfoil_vs2_14 =
                    lam.z_ax * (lam.ax_hk2 * lam.hk2_u2 + lam.ax_rt2 * lam.rt2_u2);

                println!(
                    "  LAM terms: VS2_1_2 {:.6e}/{:.6e} VS2_1_3 {:.6e}/{:.6e} VS2_1_4 {:.6e}/{:.6e}",
                    rust_vs2_12, xfoil_vs2_12, rust_vs2_13, xfoil_vs2_13, rust_vs2_14, xfoil_vs2_14
                );
                println!(
                    "  LAM ax: AX_HK2 {:.6e}/{:.6e} AX_T2 {:.6e}/{:.6e} AX_RT2 {:.6e}/{:.6e} HK2_T2 {:.6e}/{:.6e} RT2_T2 {:.6e}/{:.6e}",
                    terms_bldif.ax_hk2,
                    lam.ax_hk2,
                    terms_bldif.ax_t2,
                    lam.ax_t2,
                    terms_bldif.ax_rt2,
                    lam.ax_rt2,
                    terms_bldif.hk2_t2,
                    lam.hk2_t2,
                    terms_bldif.rt2_t2,
                    lam.rt2_t2
                );
            } else {
                println!("  LAM terms: missing for iter {}", ev.newton_iter);
            }
            if let Some(shape) = shape_by_key.get(&(ibl - 1, ev.newton_iter)) {
                let rust_vs2_31 = terms_bldif.z_di2 * terms_bldif.di2_s;
                let rust_vs2_32_hs = terms_bldif.z_hs2 * terms_bldif.hs2_t;
                let rust_vs2_32_cf = terms_bldif.z_cf2_shape * terms_bldif.cf2_t_shape;
                let rust_vs2_32_di = terms_bldif.z_di2 * terms_bldif.di2_t;
                let rust_vs2_32_t = terms_bldif.z_t2_shape;
                let rust_vs2_32_hc = 0.5
                    * (terms_bldif.z_hca * terms_bldif.hc2_t
                        + terms_bldif.z_ha_shape * terms_bldif.h2_t);
                let rust_vs2_32_upw = terms_bldif.z_upw_shape * terms_bldif.upw_t2_shape;
                let rust_vs2_32 =
                    rust_vs2_32_hs + rust_vs2_32_cf + rust_vs2_32_di + rust_vs2_32_t
                        + rust_vs2_32_hc + rust_vs2_32_upw;

                let rust_vs2_33_hs = terms_bldif.z_hs2 * terms_bldif.hs2_d;
                let rust_vs2_33_cf = terms_bldif.z_cf2_shape * terms_bldif.cf2_d_shape;
                let rust_vs2_33_di = terms_bldif.z_di2 * terms_bldif.di2_d;
                let rust_vs2_33_hc = 0.5
                    * (terms_bldif.z_hca * terms_bldif.hc2_d
                        + terms_bldif.z_ha_shape * terms_bldif.h2_d);
                let rust_vs2_33_upw = terms_bldif.z_upw_shape * terms_bldif.upw_d2_shape;
                let rust_vs2_33 =
                    rust_vs2_33_hs + rust_vs2_33_cf + rust_vs2_33_di + rust_vs2_33_hc
                        + rust_vs2_33_upw;
                let rust_vs2_34 = terms_bldif.z_hs2 * terms_bldif.hs2_u
                    + terms_bldif.z_cf2_shape * terms_bldif.cf2_u_shape
                    + terms_bldif.z_di2 * terms_bldif.di2_u
                    + terms_bldif.z_u2_shape
                    + 0.5 * terms_bldif.z_hca * terms_bldif.hc2_u
                    + terms_bldif.z_upw_shape * terms_bldif.upw_u2_shape;

                let xfoil_vs2_31 = shape.z_di2 * shape.di2_s2;
                let xfoil_vs2_32_hs = shape.z_hs2 * shape.hs2_t2;
                let xfoil_vs2_32_cf = shape.z_cf2 * shape.cf2_t2;
                let xfoil_vs2_32_di = shape.z_di2 * shape.di2_t2;
                let xfoil_vs2_32_t = shape.z_t2;
                let xfoil_vs2_32_hc =
                    0.5 * (shape.z_hca * shape.hc2_t2 + shape.z_ha * shape.h2_t2);
                let xfoil_vs2_32_upw = shape.z_upw * shape.upw_t2;
                let xfoil_vs2_32 =
                    xfoil_vs2_32_hs + xfoil_vs2_32_cf + xfoil_vs2_32_di + xfoil_vs2_32_t
                        + xfoil_vs2_32_hc + xfoil_vs2_32_upw;

                let xfoil_vs2_33_hs = shape.z_hs2 * shape.hs2_d2;
                let xfoil_vs2_33_cf = shape.z_cf2 * shape.cf2_d2;
                let xfoil_vs2_33_di = shape.z_di2 * shape.di2_d2;
                let xfoil_vs2_33_hc =
                    0.5 * (shape.z_hca * shape.hc2_d2 + shape.z_ha * shape.h2_d2);
                let xfoil_vs2_33_upw = shape.z_upw * shape.upw_d2;
                let xfoil_vs2_33 =
                    xfoil_vs2_33_hs + xfoil_vs2_33_cf + xfoil_vs2_33_di + xfoil_vs2_33_hc
                        + xfoil_vs2_33_upw;
                let xfoil_vs2_34 = shape.z_hs2 * shape.hs2_u2
                    + shape.z_cf2 * shape.cf2_u2
                    + shape.z_di2 * shape.di2_u2
                    + shape.z_u2
                    + 0.5 * shape.z_hca * shape.hc2_u2
                    + shape.z_upw * shape.upw_u2;

                println!(
                    "  SHAPE terms: VS2_3_1 {:.6e}/{:.6e} VS2_3_2 {:.6e}/{:.6e} VS2_3_3 {:.6e}/{:.6e} VS2_3_4 {:.6e}/{:.6e}",
                    rust_vs2_31,
                    xfoil_vs2_31,
                    rust_vs2_32,
                    xfoil_vs2_32,
                    rust_vs2_33,
                    xfoil_vs2_33,
                    rust_vs2_34,
                    xfoil_vs2_34
                );
                if (rust_vs2_32 - xfoil_vs2_32).abs() > 1e-6 {
                    println!(
                        "  SHAPE 3_2 parts: HS {:.6e}/{:.6e} CF {:.6e}/{:.6e} DI {:.6e}/{:.6e} T {:.6e}/{:.6e} HC {:.6e}/{:.6e} UPW {:.6e}/{:.6e}",
                        rust_vs2_32_hs,
                        xfoil_vs2_32_hs,
                        rust_vs2_32_cf,
                        xfoil_vs2_32_cf,
                        rust_vs2_32_di,
                        xfoil_vs2_32_di,
                        rust_vs2_32_t,
                        xfoil_vs2_32_t,
                        rust_vs2_32_hc,
                        xfoil_vs2_32_hc,
                        rust_vs2_32_upw,
                        xfoil_vs2_32_upw
                    );
                }
                if (rust_vs2_33 - xfoil_vs2_33).abs() > 1e-6 {
                    println!(
                        "  SHAPE 3_3 parts: HS {:.6e}/{:.6e} CF {:.6e}/{:.6e} DI {:.6e}/{:.6e} HC {:.6e}/{:.6e} UPW {:.6e}/{:.6e}",
                        rust_vs2_33_hs,
                        xfoil_vs2_33_hs,
                        rust_vs2_33_cf,
                        xfoil_vs2_33_cf,
                        rust_vs2_33_di,
                        xfoil_vs2_33_di,
                        rust_vs2_33_hc,
                        xfoil_vs2_33_hc,
                        rust_vs2_33_upw,
                        xfoil_vs2_33_upw
                    );
                }
            } else {
                println!("  SHAPE terms: missing for iter {}", ev.newton_iter);
            }
            if let Some(shear) = shear_by_key.get(&(ibl, ev.newton_iter)) {
                let rust_vs2_12 = terms_bldif.z_upw * terms_bldif.upw_t2
                    + terms_bldif.z_de2 * terms_bldif.de2_t2
                    + terms_bldif.z_us2 * terms_bldif.us2_t2
                    + terms_bldif.z_cq2 * terms_bldif.cq2_t2
                    + terms_bldif.z_cf2 * terms_bldif.cf2_t2
                    + terms_bldif.z_hk2 * terms_bldif.hk2_t2;
                let rust_vs2_13 = terms_bldif.z_d2
                    + terms_bldif.z_upw * terms_bldif.upw_d2
                    + terms_bldif.z_de2 * terms_bldif.de2_d2
                    + terms_bldif.z_us2 * terms_bldif.us2_d2
                    + terms_bldif.z_cq2 * terms_bldif.cq2_d2
                    + terms_bldif.z_cf2 * terms_bldif.cf2_d2
                    + terms_bldif.z_hk2 * terms_bldif.hk2_d2;
                let rust_vs2_14 = terms_bldif.z_u2
                    + terms_bldif.z_upw * terms_bldif.upw_u2
                    + terms_bldif.z_de2 * terms_bldif.de2_u2
                    + terms_bldif.z_us2 * terms_bldif.us2_u2
                    + terms_bldif.z_cq2 * terms_bldif.cq2_u2
                    + terms_bldif.z_cf2 * terms_bldif.cf2_u2
                    + terms_bldif.z_hk2 * terms_bldif.hk2_u2;

                let xfoil_vs2_12 = shear.z_upw * shear.upw_t2
                    + shear.z_de2 * shear.de2_t2
                    + shear.z_us2 * shear.us2_t2
                    + shear.z_cq2 * shear.cq2_t2
                    + shear.z_cf2 * shear.cf2_t2
                    + shear.z_hk2 * shear.hk2_t2;
                let xfoil_vs2_13 = shear.z_d2
                    + shear.z_upw * shear.upw_d2
                    + shear.z_de2 * shear.de2_d2
                    + shear.z_us2 * shear.us2_d2
                    + shear.z_cq2 * shear.cq2_d2
                    + shear.z_cf2 * shear.cf2_d2
                    + shear.z_hk2 * shear.hk2_d2;
                let xfoil_vs2_14 = shear.z_u2
                    + shear.z_upw * shear.upw_u2
                    + shear.z_de2 * shear.de2_u2
                    + shear.z_us2 * shear.us2_u2
                    + shear.z_cq2 * shear.cq2_u2
                    + shear.z_cf2 * shear.cf2_u2
                    + shear.z_hk2 * shear.hk2_u2;

                println!(
                    "  SHEAR terms: VS2_1_2 {:.6e}/{:.6e} VS2_1_3 {:.6e}/{:.6e} VS2_1_4 {:.6e}/{:.6e}",
                    rust_vs2_12, xfoil_vs2_12, rust_vs2_13, xfoil_vs2_13, rust_vs2_14, xfoil_vs2_14
                );
            } else {
                println!("  SHEAR terms: missing for iter {}", ev.newton_iter);
            }
        }
    }
}
