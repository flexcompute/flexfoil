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
