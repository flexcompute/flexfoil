//! Debug output infrastructure for XFOIL comparison
//!
//! This module provides JSON debug output matching XFOIL's instrumented format,
//! enabling direct comparison of internal solver state.

use serde::Serialize;
use std::cell::RefCell;
use std::fs::File;
use std::io::Write;
use std::path::Path;

thread_local! {
    static DEBUG_COLLECTOR: RefCell<Option<DebugCollector>> = RefCell::new(None);
}

/// Initialize debug collection to a file
pub fn init_debug(path: impl AsRef<Path>) {
    DEBUG_COLLECTOR.with(|c| {
        *c.borrow_mut() = Some(DebugCollector::new(path));
    });
}

/// Finalize and write debug output
pub fn finalize_debug() {
    DEBUG_COLLECTOR.with(|c| {
        if let Some(collector) = c.borrow_mut().take() {
            collector.write();
        }
    });
}

/// Check if debug collection is active
pub fn is_debug_active() -> bool {
    DEBUG_COLLECTOR.with(|c| c.borrow().is_some())
}

/// Add an event to the debug collector
pub fn add_event(event: DebugEvent) {
    DEBUG_COLLECTOR.with(|c| {
        if let Some(ref mut collector) = *c.borrow_mut() {
            collector.add_event(event);
        }
    });
}

/// Debug event collector
pub struct DebugCollector {
    events: Vec<DebugEvent>,
    path: String,
    call_id: usize,
}

impl DebugCollector {
    pub fn new(path: impl AsRef<Path>) -> Self {
        Self {
            events: Vec::new(),
            path: path.as_ref().to_string_lossy().to_string(),
            call_id: 0,
        }
    }

    pub fn add_event(&mut self, mut event: DebugEvent) {
        self.call_id += 1;
        event.call_id = self.call_id;
        self.events.push(event);
    }

    pub fn write(self) {
        let output = DebugOutput {
            events: self.events,
        };
        if let Ok(mut file) = File::create(&self.path) {
            if let Ok(json) = serde_json::to_string_pretty(&output) {
                let _ = file.write_all(json.as_bytes());
            }
        }
    }
}

/// Container for all debug events
#[derive(Serialize)]
pub struct DebugOutput {
    pub events: Vec<DebugEvent>,
}

/// A single debug event matching XFOIL's format
#[derive(Serialize, Clone)]
pub struct DebugEvent {
    pub call_id: usize,
    pub subroutine: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub iteration: Option<usize>,
    #[serde(flatten)]
    pub data: DebugData,
}

/// Type-specific debug data
#[derive(Serialize, Clone)]
#[serde(untagged)]
pub enum DebugData {
    Viscal(ViscalEvent),
    ViscalResult(ViscalResultEvent),
    Blvar(BlvarEvent),
    Bldif(BldifEvent),
    BldifTerms(BldifTermsEvent),
    MrchueIter(MrchueIterEvent),
    MrchueSys(MrchueSysEvent),
    MrchueMode(MrchueModeEvent),
    Vs2Before(Vs2BeforeEvent),
    Trdif(TrdifEvent),
    TrdifDerivs(TrdifDerivsEvent),
    Trchek2Iter(Trchek2IterEvent),
    Trchek2Final(Trchek2FinalEvent),
    Mrchue(MrchueEvent),
    Update(UpdateEvent),
    Qdcalc(QdcalcEvent),
    Blsolv(BlsolvEvent),
}

/// VISCAL iteration start event
#[derive(Serialize, Clone)]
pub struct ViscalEvent {
    pub alpha_rad: f64,
    pub reynolds: f64,
    pub mach: f64,
    pub ncrit: f64,
}

/// VISCAL iteration result event
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct ViscalResultEvent {
    pub rms_residual: f64,
    pub max_residual: f64,
    pub CL: f64,
    pub CD: f64,
    pub CM: f64,
}

/// BLVAR secondary variables event
#[derive(Serialize, Clone)]
pub struct BlvarEvent {
    pub side: usize,
    pub ibl: usize,
    pub flow_type: usize,
    pub input: BlvarInput,
    pub output: BlvarOutput,
}

#[derive(Serialize, Clone)]
pub struct BlvarInput {
    pub x: f64,
    pub u: f64,
    pub theta: f64,
    pub delta_star: f64,
    pub ctau: f64,
    pub ampl: f64,
}

#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct BlvarOutput {
    pub H: f64,
    pub Hk: f64,
    pub Hs: f64,
    pub Hc: f64,
    pub Rtheta: f64,
    pub Cf: f64,
    pub Cd: f64,
    pub Us: f64,
    pub Cq: f64,
    pub De: f64,
}

/// BLDIF Jacobian event
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct BldifEvent {
    pub side: usize,
    pub ibl: usize,
    pub flow_type: usize,
    pub VS1: [[f64; 5]; 4],
    pub VS2: [[f64; 5]; 4],
    pub VSREZ: [f64; 4],
}

/// BLDIF intermediate terms for diagnostics
#[derive(Serialize, Clone)]
pub struct BldifTermsEvent {
    pub side: usize,
    pub ibl: usize,
    pub flow_type: usize,
    pub xlog: f64,
    pub ulog: f64,
    pub tlog: f64,
    pub hlog: f64,
    pub upw: f64,
    pub ha: f64,
    pub btmp_mom: f64,
    pub cfx: f64,
    pub cfx_ta: f64,
    pub cfx_t1: f64,
    pub cfx_t2: f64,
    pub btmp_shape: f64,
    pub cfx_shape: f64,
    pub dix: f64,
    pub hsa: f64,
    pub hca: f64,
    pub xot1: f64,
    pub xot2: f64,
    pub z_di2: f64,
    pub di2_s: f64,
    pub z_upw_shape: f64,
    pub upw_t2_shape: f64,
    pub upw_d2_shape: f64,
    pub upw_u2_shape: f64,
}

/// TRCHEK2 iteration state (transition check)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct Trchek2IterEvent {
    pub global_iter: usize,
    pub side: usize,
    pub ibl: usize,
    pub trchek_iter: usize,
    pub x1: f64,
    pub x2: f64,
    pub ampl1: f64,
    pub ampl2: f64,
    pub ax: f64,
    pub residual: f64,
    pub wf1: f64,
    pub wf2: f64,
    pub xt: f64,
    pub Hk1: f64,
    pub Hk2: f64,
    pub Rt1: f64,
    pub Rt2: f64,
    pub T1: f64,
    pub T2: f64,
    pub U1: f64,
    pub U2: f64,
    pub Ncrit: f64,
    pub transition: bool,
}

/// TRCHEK2 final result (transition check)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct Trchek2FinalEvent {
    pub global_iter: usize,
    pub side: usize,
    pub ibl: usize,
    pub n_iterations: usize,
    pub converged: bool,
    pub x1: f64,
    pub x2: f64,
    pub ampl1: f64,
    pub ampl2_final: f64,
    pub ax_final: f64,
    pub xt_final: f64,
    pub Ncrit: f64,
    pub transition: bool,
    pub forced: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tt: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dt: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ut: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub Hk_t: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub Rt_t: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub St: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub Cq_t: Option<f64>,
}

/// MRCHUE initial march event
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct MrchueEvent {
    pub side: usize,
    pub ibl: usize,
    pub x: f64,
    pub Ue: f64,
    pub theta: f64,
    pub delta_star: f64,
    pub Hk: f64,
    pub Cf: f64,
}

/// MRCHUE Newton iteration event (per-iteration state)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct MrchueIterEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub direct: bool,
    pub dmax: f64,
    pub relaxation: f64,
    pub res_third: f64,
    pub res_mom: f64,
    pub res_shape: f64,
    pub delta_s: f64,
    pub delta_theta: f64,
    pub delta_delta_star: f64,
    pub delta_ue: f64,
    pub theta: f64,
    pub delta_star: f64,
    pub Ue: f64,
    pub H: f64,
    pub Hk: f64,
    pub Cf: f64,
    pub Cd: f64,
    pub ampl: f64,
    pub ctau: f64,
}

/// MRCHUE Newton system event (per-iteration Jacobian/residual)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct MrchueSysEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub direct: bool,
    pub flow_type: usize,
    pub x1: f64,
    pub x2: f64,
    pub t1: f64,
    pub t2: f64,
    pub u1: f64,
    pub u2: f64,
    pub hk1: f64,
    pub hk2: f64,
    pub ctau1: f64,
    pub ctau2: f64,
    pub VS2: [[f64; 5]; 3],
    pub VSREZ: [f64; 3],
}

/// MRCHUE mode decision event (direct vs inverse)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct MrchueModeEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub direct_before: bool,
    pub direct_after: bool,
    pub flow_type: usize,
    pub dmax: f64,
    pub rlx: f64,
    pub Htest: f64,
    pub Hk_test: f64,
    pub Hmax: f64,
    pub Ue: f64,
    pub theta: f64,
    pub delta_star: f64,
}

/// VS2 system dump (4x4 Newton system before solve)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct Vs2BeforeEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub direct: bool,
    pub flow_type: usize,
    pub VS2_4x4: [[f64; 4]; 4],
    pub VSREZ_rhs: [f64; 4],
}

/// TRDIF combined system (transition interval)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct TrdifEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub VS1: [[f64; 5]; 4],
    pub VS2: [[f64; 5]; 4],
    pub VSREZ: [f64; 4],
}

/// TRDIF chain-rule derivatives (XT/WF/TT/DT/UT)
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct TrdifDerivsEvent {
    pub side: usize,
    pub ibl: usize,
    pub iter: usize,
    pub WF1: f64,
    pub WF2: f64,
    pub XT: f64,
    pub XT_A1: f64,
    pub XT_X1: f64,
    pub XT_X2: f64,
    pub XT_T1: f64,
    pub XT_T2: f64,
    pub XT_D1: f64,
    pub XT_D2: f64,
    pub XT_U1: f64,
    pub XT_U2: f64,
    pub XT_MS: f64,
    pub XT_RE: f64,
    pub TT_A1: f64,
    pub TT_X1: f64,
    pub TT_X2: f64,
    pub TT_T1: f64,
    pub TT_T2: f64,
    pub TT_D1: f64,
    pub TT_D2: f64,
    pub TT_U1: f64,
    pub TT_U2: f64,
    pub TT_MS: f64,
    pub TT_RE: f64,
    pub DT_A1: f64,
    pub DT_X1: f64,
    pub DT_X2: f64,
    pub DT_T1: f64,
    pub DT_T2: f64,
    pub DT_D1: f64,
    pub DT_D2: f64,
    pub DT_U1: f64,
    pub DT_U2: f64,
    pub DT_MS: f64,
    pub DT_RE: f64,
    pub UT_A1: f64,
    pub UT_X1: f64,
    pub UT_X2: f64,
    pub UT_T1: f64,
    pub UT_T2: f64,
    pub UT_D1: f64,
    pub UT_D2: f64,
    pub UT_U1: f64,
    pub UT_U2: f64,
    pub UT_MS: f64,
    pub UT_RE: f64,
}

/// UPDATE Newton delta event
#[derive(Serialize, Clone)]
pub struct UpdateEvent {
    pub side: usize,
    pub ibl: usize,
    pub delta_ctau: f64,
    pub delta_theta: f64,
    pub delta_mass: f64,
    #[serde(rename = "delta_Ue")]
    pub delta_ue: f64,
    pub relaxation: f64,
}

/// QDCALC DIJ matrix event
#[derive(Serialize, Clone)]
#[allow(non_snake_case)]
pub struct QdcalcEvent {
    pub n_airfoil: usize,
    pub n_wake: usize,
    pub n_total: usize,
    pub DIJ_diagonal_sample: Vec<f64>,
    pub DIJ_row1_sample: Vec<f64>,
}

/// BLSOLV Newton system event
#[derive(Serialize, Clone)]
pub struct BlsolvEvent {
    pub system_size: usize,
}

// Helper functions to create events easily
impl DebugEvent {
    pub fn viscal(iteration: usize, alpha_rad: f64, reynolds: f64, mach: f64, ncrit: f64) -> Self {
        Self {
            call_id: 0,
            subroutine: "VISCAL".to_string(),
            iteration: Some(iteration),
            data: DebugData::Viscal(ViscalEvent {
                alpha_rad,
                reynolds,
                mach,
                ncrit,
            }),
        }
    }

    pub fn viscal_result(iteration: usize, rms: f64, max: f64, cl: f64, cd: f64, cm: f64) -> Self {
        Self {
            call_id: 0,
            subroutine: "VISCAL_RESULT".to_string(),
            iteration: Some(iteration),
            data: DebugData::ViscalResult(ViscalResultEvent {
                rms_residual: rms,
                max_residual: max,
                CL: cl,
                CD: cd,
                CM: cm,
            }),
        }
    }

    pub fn blvar(
        iteration: usize,
        side: usize,
        ibl: usize,
        flow_type: usize,
        input: BlvarInput,
        output: BlvarOutput,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "BLVAR".to_string(),
            iteration: Some(iteration),
            data: DebugData::Blvar(BlvarEvent {
                side,
                ibl,
                flow_type,
                input,
                output,
            }),
        }
    }

    pub fn bldif(
        iteration: usize,
        side: usize,
        ibl: usize,
        flow_type: usize,
        vs1: [[f64; 5]; 4],
        vs2: [[f64; 5]; 4],
        vsrez: [f64; 4],
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "BLDIF".to_string(),
            iteration: Some(iteration),
            data: DebugData::Bldif(BldifEvent {
                side,
                ibl,
                flow_type,
                VS1: vs1,
                VS2: vs2,
                VSREZ: vsrez,
            }),
        }
    }

    pub fn bldif_terms(
        iteration: usize,
        side: usize,
        ibl: usize,
        flow_type: usize,
        mut terms: BldifTermsEvent,
    ) -> Self {
        terms.side = side;
        terms.ibl = ibl;
        terms.flow_type = flow_type;
        Self {
            call_id: 0,
            subroutine: "BLDIF_TERMS".to_string(),
            iteration: Some(iteration),
            data: DebugData::BldifTerms(terms),
        }
    }

    pub fn mrchue_iter(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: MrchueIterEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "MRCHUE_ITER".to_string(),
            iteration: Some(iteration),
            data: DebugData::MrchueIter(event),
        }
    }

    pub fn mrchue_sys(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: MrchueSysEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "MRCHUE_SYS".to_string(),
            iteration: Some(iteration),
            data: DebugData::MrchueSys(event),
        }
    }

    pub fn mrchue_mode(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: MrchueModeEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "MRCHUE_MODE".to_string(),
            iteration: Some(iteration),
            data: DebugData::MrchueMode(event),
        }
    }

    pub fn vs2_before(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: Vs2BeforeEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "VS2_BEFORE".to_string(),
            iteration: Some(iteration),
            data: DebugData::Vs2Before(event),
        }
    }

    pub fn trdif(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: TrdifEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "TRDIF".to_string(),
            iteration: Some(iteration),
            data: DebugData::Trdif(event),
        }
    }

    pub fn trdif_derivs(
        iteration: usize,
        side: usize,
        ibl: usize,
        mut event: TrdifDerivsEvent,
    ) -> Self {
        event.side = side;
        event.ibl = ibl;
        Self {
            call_id: 0,
            subroutine: "TRDIF_DERIVS".to_string(),
            iteration: Some(iteration),
            data: DebugData::TrdifDerivs(event),
        }
    }

    pub fn trchek2_iter(
        iteration: usize,
        side: usize,
        ibl: usize,
        trchek_iter: usize,
        x1: f64,
        x2: f64,
        ampl1: f64,
        ampl2: f64,
        ax: f64,
        residual: f64,
        wf1: f64,
        wf2: f64,
        xt: f64,
        hk1: f64,
        hk2: f64,
        rt1: f64,
        rt2: f64,
        t1: f64,
        t2: f64,
        u1: f64,
        u2: f64,
        ncrit: f64,
        transition: bool,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "TRCHEK2_ITER".to_string(),
            iteration: Some(iteration),
            data: DebugData::Trchek2Iter(Trchek2IterEvent {
                global_iter: iteration,
                side,
                ibl,
                trchek_iter,
                x1,
                x2,
                ampl1,
                ampl2,
                ax,
                residual,
                wf1,
                wf2,
                xt,
                Hk1: hk1,
                Hk2: hk2,
                Rt1: rt1,
                Rt2: rt2,
                T1: t1,
                T2: t2,
                U1: u1,
                U2: u2,
                Ncrit: ncrit,
                transition,
            }),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn trchek2_final(
        iteration: usize,
        side: usize,
        ibl: usize,
        n_iterations: usize,
        converged: bool,
        x1: f64,
        x2: f64,
        ampl1: f64,
        ampl2_final: f64,
        ax_final: f64,
        xt_final: f64,
        ncrit: f64,
        transition: bool,
        forced: bool,
        tt: Option<f64>,
        dt: Option<f64>,
        ut: Option<f64>,
        hk_t: Option<f64>,
        rt_t: Option<f64>,
        st: Option<f64>,
        cq_t: Option<f64>,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "TRCHEK2_FINAL".to_string(),
            iteration: Some(iteration),
            data: DebugData::Trchek2Final(Trchek2FinalEvent {
                global_iter: iteration,
                side,
                ibl,
                n_iterations,
                converged,
                x1,
                x2,
                ampl1,
                ampl2_final,
                ax_final,
                xt_final,
                Ncrit: ncrit,
                transition,
                forced,
                tt,
                dt,
                ut,
                Hk_t: hk_t,
                Rt_t: rt_t,
                St: st,
                Cq_t: cq_t,
            }),
        }
    }

    pub fn mrchue(
        side: usize,
        ibl: usize,
        x: f64,
        ue: f64,
        theta: f64,
        delta_star: f64,
        hk: f64,
        cf: f64,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "MRCHUE".to_string(),
            iteration: None,
            data: DebugData::Mrchue(MrchueEvent {
                side,
                ibl,
                x,
                Ue: ue,
                theta,
                delta_star,
                Hk: hk,
                Cf: cf,
            }),
        }
    }

    pub fn update(
        iteration: usize,
        side: usize,
        ibl: usize,
        delta_ctau: f64,
        delta_theta: f64,
        delta_mass: f64,
        delta_ue: f64,
        relaxation: f64,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "UPDATE".to_string(),
            iteration: Some(iteration),
            data: DebugData::Update(UpdateEvent {
                side,
                ibl,
                delta_ctau,
                delta_theta,
                delta_mass,
                delta_ue,
                relaxation,
            }),
        }
    }

    pub fn qdcalc(
        n_airfoil: usize,
        n_wake: usize,
        n_total: usize,
        diag_sample: Vec<f64>,
        row1_sample: Vec<f64>,
    ) -> Self {
        Self {
            call_id: 0,
            subroutine: "QDCALC".to_string(),
            iteration: None,
            data: DebugData::Qdcalc(QdcalcEvent {
                n_airfoil,
                n_wake,
                n_total,
                DIJ_diagonal_sample: diag_sample,
                DIJ_row1_sample: row1_sample,
            }),
        }
    }

    pub fn blsolv(iteration: usize, system_size: usize) -> Self {
        Self {
            call_id: 0,
            subroutine: "BLSOLV".to_string(),
            iteration: Some(iteration),
            data: DebugData::Blsolv(BlsolvEvent { system_size }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debug_event_serialization() {
        let event = DebugEvent::viscal(1, 0.0698, 3_000_000.0, 0.0, 9.0);
        let json = serde_json::to_string_pretty(&event).unwrap();
        assert!(json.contains("VISCAL"));
        assert!(json.contains("alpha_rad"));
    }

    #[test]
    fn test_blvar_event() {
        let input = BlvarInput {
            x: 0.05,
            u: 1.02,
            theta: 0.00012,
            delta_star: 0.00031,
            ctau: 0.03,
            ampl: 0.0,
        };
        let output = BlvarOutput {
            H: 2.58,
            Hk: 2.58,
            Hs: 1.62,
            Hc: 0.0,
            Rtheta: 1234.5,
            Cf: 0.0045,
            Cd: 0.0012,
            Us: 0.22,
            Cq: 0.07,
            De: 0.00014,
        };
        let event = DebugEvent::blvar(1, 1, 3, 1, input, output);
        let json = serde_json::to_string_pretty(&event).unwrap();
        assert!(json.contains("BLVAR"));
        assert!(json.contains("\"H\": 2.58"));
    }
}
