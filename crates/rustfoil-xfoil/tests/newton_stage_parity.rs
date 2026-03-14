mod support;

use serde_json::Value;

use support::{assert_close_scalar, assert_close_slice, rust_iteration_debug, xfoil_iteration_debug};

fn event_by_name<'a>(debug: &'a Value, name: &str) -> &'a Value {
    debug["events"]
        .as_array()
        .expect("debug events")
        .iter()
        .find(|event| event["subroutine"].as_str() == Some(name))
        .unwrap_or_else(|| panic!("missing debug event {name}"))
}

fn events_by_name<'a>(debug: &'a Value, name: &str) -> Vec<&'a Value> {
    debug["events"]
        .as_array()
        .expect("debug events")
        .iter()
        .filter(|event| event["subroutine"].as_str() == Some(name))
        .collect()
}

fn flatten_numbers(value: &Value, out: &mut Vec<f64>) {
    match value {
        Value::Array(values) => {
            for value in values {
                flatten_numbers(value, out);
            }
        }
        Value::Number(number) => out.push(number.as_f64().expect("numeric value")),
        _ => {}
    }
}

fn flattened_numbers(value: &Value) -> Vec<f64> {
    let mut out = Vec::new();
    flatten_numbers(value, &mut out);
    out
}

#[test]
fn setbl_system_matches_xfoil_debug_samples() {
    let rust = rust_iteration_debug(15.0);
    let xfoil = xfoil_iteration_debug(15.0);
    let rust_event = event_by_name(&rust, "SETBL_SYSTEM");
    let xfoil_event = event_by_name(&xfoil, "SETBL_SYSTEM");

    assert_eq!(rust_event["nsys"], xfoil_event["nsys"], "setbl_system.nsys");
    for field in ["VA", "VB", "VDEL", "VM_diagonal", "VM_row1"] {
        let rust_values = flattened_numbers(&rust_event[field]);
        let xfoil_values = flattened_numbers(&xfoil_event[field]);
        let tol = match field {
            // XFOIL's debug fixture rounds these large VM entries to coarse decimals.
            "VM_diagonal" | "VM_row1" => 5.0e-1,
            _ => 2.0e-2,
        };
        assert_close_slice(&format!("setbl_system.{field}"), &rust_values, &xfoil_values, tol);
    }
}

#[test]
fn blsolv_solution_matches_xfoil_debug_samples() {
    let rust = rust_iteration_debug(15.0);
    let xfoil = xfoil_iteration_debug(15.0);
    let rust_event = event_by_name(&rust, "BLSOLV_SOLUTION");
    let xfoil_event = event_by_name(&xfoil, "BLSOLV_SOLUTION");

    assert_eq!(rust_event["nsys"], xfoil_event["nsys"], "blsolv_solution.nsys");
    let rust_values = flattened_numbers(&rust_event["deltas"]);
    let xfoil_values = flattened_numbers(&xfoil_event["deltas"]);
    assert_close_slice("blsolv_solution.deltas", &rust_values, &xfoil_values, 2.0e-4);
}

#[test]
fn update_detailed_matches_xfoil_first_surface_samples() {
    let rust = rust_iteration_debug(15.0);
    let xfoil = xfoil_iteration_debug(15.0);
    let rust_events = events_by_name(&rust, "UPDATE_DETAILED");
    let xfoil_events: Vec<&Value> = events_by_name(&xfoil, "UPDATE_DETAILED")
        .into_iter()
        .filter(|event| event["side"].as_u64() == Some(1))
        .collect();

    assert!(
        !rust_events.is_empty(),
        "rust debug should emit UPDATE_DETAILED events"
    );
    assert!(
        !xfoil_events.is_empty(),
        "xfoil debug should emit UPDATE_DETAILED events"
    );

    let n = rust_events.len().min(xfoil_events.len()).min(20);
    assert!(n > 0, "need comparable UPDATE_DETAILED events");
    for idx in 0..n {
        let rust_event = rust_events[idx];
        let xfoil_event = xfoil_events[idx];
        assert_close_scalar(
            &format!("update_detailed.theta_before[{idx}]"),
            rust_event["theta_before"].as_f64().expect("rust theta_before"),
            xfoil_event["theta_before"].as_f64().expect("xfoil theta_before"),
            2.0e-6,
        );
        assert_close_scalar(
            &format!("update_detailed.delta_star_before[{idx}]"),
            rust_event["delta_star_before"].as_f64().expect("rust delta_star_before"),
            xfoil_event["delta_star_before"].as_f64().expect("xfoil delta_star_before"),
            2.0e-6,
        );
        assert_close_scalar(
            &format!("update_detailed.ue_before[{idx}]"),
            rust_event["ue_before"].as_f64().expect("rust ue_before"),
            xfoil_event["ue_before"].as_f64().expect("xfoil ue_before"),
            2.0e-3,
        );
        assert_close_scalar(
            &format!("update_detailed.delta_theta[{idx}]"),
            rust_event["delta_theta"].as_f64().expect("rust delta_theta"),
            xfoil_event["delta_theta"].as_f64().expect("xfoil delta_theta"),
            2.0e-6,
        );
        assert_close_scalar(
            &format!("update_detailed.delta_mass[{idx}]"),
            rust_event["delta_mass"].as_f64().expect("rust delta_mass"),
            xfoil_event["delta_mass"].as_f64().expect("xfoil delta_mass"),
            2.0e-5,
        );
        assert_close_scalar(
            &format!("update_detailed.delta_ue[{idx}]"),
            rust_event["delta_ue"].as_f64().expect("rust delta_ue"),
            xfoil_event["delta_ue"].as_f64().expect("xfoil delta_ue"),
            2.0e-2,
        );
        assert_close_scalar(
            &format!("update_detailed.relaxation[{idx}]"),
            rust_event["relaxation"].as_f64().expect("rust relaxation"),
            xfoil_event["relaxation"].as_f64().expect("xfoil relaxation"),
            2.0e-3,
        );
        assert_close_scalar(
            &format!("update_detailed.theta_after[{idx}]"),
            rust_event["theta_after"].as_f64().expect("rust theta_after"),
            xfoil_event["theta_after"].as_f64().expect("xfoil theta_after"),
            2.0e-6,
        );
        assert_close_scalar(
            &format!("update_detailed.delta_star_after[{idx}]"),
            rust_event["delta_star_after"].as_f64().expect("rust delta_star_after"),
            xfoil_event["delta_star_after"].as_f64().expect("xfoil delta_star_after"),
            2.0e-6,
        );
        assert_close_scalar(
            &format!("update_detailed.ue_after[{idx}]"),
            rust_event["ue_after"].as_f64().expect("rust ue_after"),
            xfoil_event["ue_after"].as_f64().expect("xfoil ue_after"),
            2.0e-3,
        );
    }
}
