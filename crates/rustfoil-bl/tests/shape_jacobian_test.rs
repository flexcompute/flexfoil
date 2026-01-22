//! Unit tests for shape equation Jacobian coefficients
//!
//! Tests the intermediate calculations needed for the shape equation Jacobian
//! by comparing against XFOIL's SHAPE_JACOBIAN debug output.

use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;
use std::fs;

/// Load the XFOIL debug trace
fn load_debug_trace() -> Option<serde_json::Value> {
    let paths = [
        "testdata/newton_debug_trace.json",
        "../testdata/newton_debug_trace.json",
        "../../testdata/newton_debug_trace.json",
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

/// Load station reference data
fn load_station_reference() -> Option<serde_json::Value> {
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

/// Extract SHAPE_JACOBIAN event for IBL=3, iteration 1
fn get_shape_jacobian_event(trace: &serde_json::Value) -> Option<serde_json::Value> {
    let events = trace["events"].as_array()?;
    events
        .iter()
        .find(|e| {
            e["subroutine"].as_str() == Some("SHAPE_JACOBIAN")
                && e["side"].as_i64() == Some(1)
                && e["ibl"].as_i64() == Some(3)
                && e["iteration"].as_i64() == Some(1)
        })
        .cloned()
}

#[test]
fn test_shape_z_coefficients() {
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let shape_event = match get_shape_jacobian_event(&trace) {
        Some(e) => e,
        None => {
            eprintln!("Skipping: SHAPE_JACOBIAN event not found for IBL=3");
            return;
        }
    };

    println!("\n=== XFOIL Shape Equation Z Coefficients (IBL=3, iter 1) ===\n");
    
    // Extract XFOIL values
    let z_hs2 = shape_event["Z_HS2"].as_f64().unwrap();
    let z_cf2 = shape_event["Z_CF2"].as_f64().unwrap();
    let z_di2 = shape_event["Z_DI2"].as_f64().unwrap();
    let z_t2 = shape_event["Z_T2"].as_f64().unwrap();
    let z_u2 = shape_event["Z_U2"].as_f64().unwrap();
    let z_hca = shape_event["Z_HCA"].as_f64().unwrap();
    let z_ha = shape_event["Z_HA"].as_f64().unwrap();
    let z_upw = shape_event["Z_UPW"].as_f64().unwrap();

    println!("Z_HS2 = {:+.6e}", z_hs2);
    println!("Z_CF2 = {:+.6e}", z_cf2);
    println!("Z_DI2 = {:+.6e}", z_di2);
    println!("Z_T2  = {:+.6e}", z_t2);
    println!("Z_U2  = {:+.6e}", z_u2);
    println!("Z_HCA = {:+.6e}", z_hca);
    println!("Z_HA  = {:+.6e}", z_ha);
    println!("Z_UPW = {:+.6e}", z_upw);

    // These are the XFOIL formulas for the Z coefficients:
    // Z_CFX =  XLOG*0.5
    // Z_DIX = -XLOG
    // Z_HCA = 2.0*ULOG/HSA
    // Z_HA  = -ULOG
    // Z_HL  = DDLOG = 1.0
    // Z_UL  = DDLOG * BTMP = BTMP
    // 
    // Z_HS2 = -HCA*ULOG/HSA**2 + Z_HL/HS2 = -HCA*ULOG/HSA**2 + 1/HS2
    // Z_CF2 = UPW * Z_CFX * XOT2 = UPW * 0.5 * XLOG * X2/T2
    // Z_DI2 = UPW * Z_DIX * XOT2 = -UPW * XLOG * X2/T2
    // Z_T2  = UPW*(Z_CFX*CF2 + Z_DIX*DI2)*(-XOT2/T2) + Z_HWA*0.5*(-DW2/T2**2)
    //       = UPW*(0.5*XLOG*CF2 - XLOG*DI2)*(-X2/T2**2) + (wake term)
    // Z_U2  = Z_UL/U2 = BTMP/U2
    
    println!("\nNote: These coefficients determine the contribution of each derivative");
    println!("to the final shape equation Jacobian.\n");
}

#[test]
fn test_shape_derivatives() {
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let shape_event = match get_shape_jacobian_event(&trace) {
        Some(e) => e,
        None => {
            eprintln!("Skipping: SHAPE_JACOBIAN event not found for IBL=3");
            return;
        }
    };

    println!("\n=== XFOIL Shape Equation Derivatives (IBL=3, iter 1) ===\n");
    
    // Extract XFOIL derivative values
    let hs2_t2 = shape_event["HS2_T2"].as_f64().unwrap();
    let hs2_d2 = shape_event["HS2_D2"].as_f64().unwrap();
    let hs2_u2 = shape_event["HS2_U2"].as_f64().unwrap();
    
    let cf2_t2 = shape_event["CF2_T2"].as_f64().unwrap();
    let cf2_d2 = shape_event["CF2_D2"].as_f64().unwrap();
    let cf2_u2 = shape_event["CF2_U2"].as_f64().unwrap();
    
    let di2_t2 = shape_event["DI2_T2"].as_f64().unwrap();
    let di2_d2 = shape_event["DI2_D2"].as_f64().unwrap();
    let di2_u2 = shape_event["DI2_U2"].as_f64().unwrap();
    let di2_s2 = shape_event["DI2_S2"].as_f64().unwrap();
    
    let h2_t2 = shape_event["H2_T2"].as_f64().unwrap();
    let h2_d2 = shape_event["H2_D2"].as_f64().unwrap();
    
    let hc2_t2 = shape_event["HC2_T2"].as_f64().unwrap();
    let hc2_d2 = shape_event["HC2_D2"].as_f64().unwrap();
    let hc2_u2 = shape_event["HC2_U2"].as_f64().unwrap();

    println!("Hs derivatives:");
    println!("  HS2_T2 = {:+.6e}", hs2_t2);
    println!("  HS2_D2 = {:+.6e}", hs2_d2);
    println!("  HS2_U2 = {:+.6e}", hs2_u2);
    
    println!("\nCf derivatives:");
    println!("  CF2_T2 = {:+.6e}", cf2_t2);
    println!("  CF2_D2 = {:+.6e}", cf2_d2);
    println!("  CF2_U2 = {:+.6e}", cf2_u2);
    
    println!("\nDI (Cd) derivatives:");
    println!("  DI2_T2 = {:+.6e}", di2_t2);
    println!("  DI2_D2 = {:+.6e}", di2_d2);
    println!("  DI2_U2 = {:+.6e}", di2_u2);
    println!("  DI2_S2 = {:+.6e}", di2_s2);
    
    println!("\nH derivatives:");
    println!("  H2_T2  = {:+.6e} (should be -δ*/θ² = -H/θ)", h2_t2);
    println!("  H2_D2  = {:+.6e} (should be 1/θ)", h2_d2);
    
    println!("\nHc derivatives:");
    println!("  HC2_T2 = {:+.6e}", hc2_t2);
    println!("  HC2_D2 = {:+.6e}", hc2_d2);
    println!("  HC2_U2 = {:+.6e}", hc2_u2);
}

#[test]
fn test_shape_final_jacobian() {
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let ref_data = match load_station_reference() {
        Some(r) => r,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let shape_event = match get_shape_jacobian_event(&trace) {
        Some(e) => e,
        None => {
            eprintln!("Skipping: SHAPE_JACOBIAN event not found for IBL=3");
            return;
        }
    };

    println!("\n=== Shape Equation Jacobian Comparison (IBL=3, iter 1) ===\n");
    
    // XFOIL final values
    let xfoil_vs2_31 = shape_event["VS2_3_1"].as_f64().unwrap();
    let xfoil_vs2_32 = shape_event["VS2_3_2"].as_f64().unwrap();
    let xfoil_vs2_33 = shape_event["VS2_3_3"].as_f64().unwrap();
    let xfoil_vs2_34 = shape_event["VS2_3_4"].as_f64().unwrap();

    println!("XFOIL VS2 row 3 (shape equation Jacobian):");
    println!("  VS2[3,1] = {:+.6e} (∂/∂ctau)", xfoil_vs2_31);
    println!("  VS2[3,2] = {:+.6e} (∂/∂θ)", xfoil_vs2_32);
    println!("  VS2[3,3] = {:+.6e} (∂/∂δ*)", xfoil_vs2_33);
    println!("  VS2[3,4] = {:+.6e} (∂/∂Ue)", xfoil_vs2_34);

    // Compute RustFoil values
    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;
    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0]; // IBL=2
    let s1 = &stations[1]; // IBL=3

    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    curr.theta = s1["initial"]["theta"].as_f64().unwrap();
    curr.delta_star = s1["initial"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    // Debug: print intermediate values
    println!("\n=== Debug: Station 2 (curr) Derivatives ===");
    println!("  hs_hk = {:+.6e}", curr.derivs.hs_hk);
    println!("  hs_rt = {:+.6e}", curr.derivs.hs_rt);
    println!("  cf_hk = {:+.6e}", curr.derivs.cf_hk);
    println!("  cf_rt = {:+.6e}", curr.derivs.cf_rt);
    println!("  cd_hk = {:+.6e}", curr.derivs.cd_hk);
    println!("  cd_rt = {:+.6e}", curr.derivs.cd_rt);
    println!("  h_theta = {:+.6e}", curr.derivs.h_theta);
    println!("  h_delta_star = {:+.6e}", curr.derivs.h_delta_star);
    println!("  hk_h = {:+.6e}", curr.derivs.hk_h);
    println!("\n=== Debug: Station 2 (curr) Values ===");
    println!("  theta = {:+.6e}", curr.theta);
    println!("  delta_star = {:+.6e}", curr.delta_star);
    println!("  hs = {:+.6e}", curr.hs);
    println!("  hk = {:+.6e}", curr.hk);
    println!("  cf = {:+.6e}", curr.cf);
    println!("  cd = {:+.6e}", curr.cd);
    println!("  r_theta = {:+.6e}", curr.r_theta);

    let (_res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);

    // Row 2 is the shape equation (0-indexed)
    let rust_vs2_31 = jac.vs2[2][0];
    let rust_vs2_32 = jac.vs2[2][1];
    let rust_vs2_33 = jac.vs2[2][2];
    let rust_vs2_34 = jac.vs2[2][3];

    println!("\nRustFoil VS2 row 2 (shape equation Jacobian):");
    println!("  vs2[2][0] = {:+.6e} (∂/∂ctau)", rust_vs2_31);
    println!("  vs2[2][1] = {:+.6e} (∂/∂θ)", rust_vs2_32);
    println!("  vs2[2][2] = {:+.6e} (∂/∂δ*)", rust_vs2_33);
    println!("  vs2[2][3] = {:+.6e} (∂/∂Ue)", rust_vs2_34);

    println!("\n=== Comparison ===");
    let entries = [
        ("∂/∂ctau", xfoil_vs2_31, rust_vs2_31),
        ("∂/∂θ", xfoil_vs2_32, rust_vs2_32),
        ("∂/∂δ*", xfoil_vs2_33, rust_vs2_33),
        ("∂/∂Ue", xfoil_vs2_34, rust_vs2_34),
    ];

    let mut errors = 0;
    for (name, xfoil, rust) in entries {
        let rel_err = if xfoil.abs() > 1e-30 {
            ((rust - xfoil) / xfoil).abs() * 100.0
        } else {
            if rust.abs() > 1e-30 { 100.0 } else { 0.0 }
        };
        let sign_ok = (xfoil >= 0.0) == (rust >= 0.0);
        
        if rel_err > 10.0 || !sign_ok {
            errors += 1;
            println!("  {} XFOIL={:+.4e}, Rust={:+.4e}, err={:5.1}%{}",
                name, xfoil, rust, rel_err,
                if sign_ok { "" } else { " !! SIGN" });
        }
    }

    if errors > 0 {
        println!("\n*** {} entries have significant errors ***", errors);
        println!("This confirms the shape equation Jacobian needs fixing.");
    } else {
        println!("\nAll entries match within tolerance!");
    }
}

#[test]
fn test_shape_jacobian_formula() {
    // Test the XFOIL formula for computing the shape equation Jacobian
    // VS2(3,2) = Z_HS2*HS2_T2 + Z_CF2*CF2_T2 + Z_DI2*DI2_T2 + Z_T2
    //          + 0.5*(Z_HCA*HC2_T2 + Z_HA*H2_T2) + Z_UPW*UPW_T2
    
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let shape_event = match get_shape_jacobian_event(&trace) {
        Some(e) => e,
        None => {
            eprintln!("Skipping: SHAPE_JACOBIAN event not found for IBL=3");
            return;
        }
    };

    println!("\n=== Verify XFOIL Shape Jacobian Formula ===\n");

    // Extract all values
    let z_hs2 = shape_event["Z_HS2"].as_f64().unwrap();
    let z_cf2 = shape_event["Z_CF2"].as_f64().unwrap();
    let z_di2 = shape_event["Z_DI2"].as_f64().unwrap();
    let z_t2 = shape_event["Z_T2"].as_f64().unwrap();
    let z_hca = shape_event["Z_HCA"].as_f64().unwrap();
    let z_ha = shape_event["Z_HA"].as_f64().unwrap();
    let z_upw = shape_event["Z_UPW"].as_f64().unwrap();

    let hs2_t2 = shape_event["HS2_T2"].as_f64().unwrap();
    let cf2_t2 = shape_event["CF2_T2"].as_f64().unwrap();
    let di2_t2 = shape_event["DI2_T2"].as_f64().unwrap();
    let hc2_t2 = shape_event["HC2_T2"].as_f64().unwrap();
    let h2_t2 = shape_event["H2_T2"].as_f64().unwrap();
    let upw_t2 = shape_event["UPW_T2"].as_f64().unwrap();

    let hs2_d2 = shape_event["HS2_D2"].as_f64().unwrap();
    let cf2_d2 = shape_event["CF2_D2"].as_f64().unwrap();
    let di2_d2 = shape_event["DI2_D2"].as_f64().unwrap();
    let hc2_d2 = shape_event["HC2_D2"].as_f64().unwrap();
    let h2_d2 = shape_event["H2_D2"].as_f64().unwrap();
    let upw_d2 = shape_event["UPW_D2"].as_f64().unwrap();

    let z_u2 = shape_event["Z_U2"].as_f64().unwrap();
    let hs2_u2 = shape_event["HS2_U2"].as_f64().unwrap();
    let cf2_u2 = shape_event["CF2_U2"].as_f64().unwrap();
    let di2_u2 = shape_event["DI2_U2"].as_f64().unwrap();
    let hc2_u2 = shape_event["HC2_U2"].as_f64().unwrap();
    let upw_u2 = shape_event["UPW_U2"].as_f64().unwrap();

    // XFOIL's actual final values
    let xfoil_vs2_32 = shape_event["VS2_3_2"].as_f64().unwrap();
    let xfoil_vs2_33 = shape_event["VS2_3_3"].as_f64().unwrap();
    let xfoil_vs2_34 = shape_event["VS2_3_4"].as_f64().unwrap();

    // Compute VS2(3,2) using XFOIL formula
    let computed_vs2_32 = z_hs2 * hs2_t2 + z_cf2 * cf2_t2 + z_di2 * di2_t2 + z_t2
        + 0.5 * (z_hca * hc2_t2 + z_ha * h2_t2) + z_upw * upw_t2;

    // Compute VS2(3,3) using XFOIL formula  
    let computed_vs2_33 = z_hs2 * hs2_d2 + z_cf2 * cf2_d2 + z_di2 * di2_d2
        + 0.5 * (z_hca * hc2_d2 + z_ha * h2_d2) + z_upw * upw_d2;

    // Compute VS2(3,4) using XFOIL formula
    let computed_vs2_34 = z_hs2 * hs2_u2 + z_cf2 * cf2_u2 + z_di2 * di2_u2 + z_u2
        + 0.5 * z_hca * hc2_u2 + z_upw * upw_u2;

    println!("VS2(3,2) = ∂(shape residual)/∂θ:");
    println!("  Components:");
    println!("    Z_HS2*HS2_T2 = {:+.4e} * {:+.4e} = {:+.4e}", z_hs2, hs2_t2, z_hs2 * hs2_t2);
    println!("    Z_CF2*CF2_T2 = {:+.4e} * {:+.4e} = {:+.4e}", z_cf2, cf2_t2, z_cf2 * cf2_t2);
    println!("    Z_DI2*DI2_T2 = {:+.4e} * {:+.4e} = {:+.4e}", z_di2, di2_t2, z_di2 * di2_t2);
    println!("    Z_T2         = {:+.4e}", z_t2);
    println!("    0.5*Z_HCA*HC2_T2 = {:+.4e}", 0.5 * z_hca * hc2_t2);
    println!("    0.5*Z_HA*H2_T2   = {:+.4e}", 0.5 * z_ha * h2_t2);
    println!("    Z_UPW*UPW_T2     = {:+.4e}", z_upw * upw_t2);
    println!("  Sum = {:+.6e}", computed_vs2_32);
    println!("  XFOIL = {:+.6e}", xfoil_vs2_32);
    println!("  Match: {}", (computed_vs2_32 - xfoil_vs2_32).abs() / xfoil_vs2_32.abs() < 0.001);

    println!("\nVS2(3,3) = ∂(shape residual)/∂δ*:");
    println!("  Components:");
    println!("    Z_HS2*HS2_D2 = {:+.4e} * {:+.4e} = {:+.4e}", z_hs2, hs2_d2, z_hs2 * hs2_d2);
    println!("    Z_CF2*CF2_D2 = {:+.4e} * {:+.4e} = {:+.4e}", z_cf2, cf2_d2, z_cf2 * cf2_d2);
    println!("    Z_DI2*DI2_D2 = {:+.4e} * {:+.4e} = {:+.4e}", z_di2, di2_d2, z_di2 * di2_d2);
    println!("    0.5*Z_HCA*HC2_D2 = {:+.4e}", 0.5 * z_hca * hc2_d2);
    println!("    0.5*Z_HA*H2_D2   = {:+.4e}", 0.5 * z_ha * h2_d2);
    println!("    Z_UPW*UPW_D2     = {:+.4e}", z_upw * upw_d2);
    println!("  Sum = {:+.6e}", computed_vs2_33);
    println!("  XFOIL = {:+.6e}", xfoil_vs2_33);
    println!("  Match: {}", (computed_vs2_33 - xfoil_vs2_33).abs() / xfoil_vs2_33.abs() < 0.001);

    println!("\nVS2(3,4) = ∂(shape residual)/∂Ue:");
    println!("  Components:");
    println!("    Z_HS2*HS2_U2 = {:+.4e} * {:+.4e} = {:+.4e}", z_hs2, hs2_u2, z_hs2 * hs2_u2);
    println!("    Z_CF2*CF2_U2 = {:+.4e} * {:+.4e} = {:+.4e}", z_cf2, cf2_u2, z_cf2 * cf2_u2);
    println!("    Z_DI2*DI2_U2 = {:+.4e} * {:+.4e} = {:+.4e}", z_di2, di2_u2, z_di2 * di2_u2);
    println!("    Z_U2         = {:+.4e}", z_u2);
    println!("    0.5*Z_HCA*HC2_U2 = {:+.4e}", 0.5 * z_hca * hc2_u2);
    println!("    Z_UPW*UPW_U2     = {:+.4e}", z_upw * upw_u2);
    println!("  Sum = {:+.6e}", computed_vs2_34);
    println!("  XFOIL = {:+.6e}", xfoil_vs2_34);
    println!("  Match: {}", (computed_vs2_34 - xfoil_vs2_34).abs() / xfoil_vs2_34.abs().max(1e-10) < 0.001);
}
