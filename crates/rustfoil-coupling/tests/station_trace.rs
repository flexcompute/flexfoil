//! Trace exact station values from mrchue_iterations.json

use rustfoil_bl::closures::{amplification_rate, axset};

#[test]
fn test_exact_station_values() {
    // Exact values from mrchue_iterations.json (alpha=5 degrees)
    let stations = [
        // (ibl, x, Hk, theta, Rtheta, ampl)
        (17, 0.034705, 2.6225_f64, 6.125481e-5_f64, 118.60_f64, 1.771892e-9_f64),
        (18, 0.037666, 2.6815, 6.800512e-5, 130.49, 1.866317e-9),
        (19, 0.040885, 2.7424, 7.533508e-5, 142.93, 9.204327e-4),
    ];
    
    println!("\n=== DAMPL at each station (mrchue_iterations.json inputs) ===");
    for (ibl, x, hk, theta, rtheta, ampl_xfoil) in &stations {
        // Compute Rcrit
        let hmi = 1.0 / (hk - 1.0);
        let aa = 2.492 * hmi.powf(0.43);
        let bb = (14.0 * hmi - 9.24).tanh();
        let grcrit = aa + 0.7 * (bb + 1.0);
        let rcrit = 10.0_f64.powf(grcrit);
        
        let gr = rtheta.log10();
        let dgr = 0.08;
        let rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
        
        let result = amplification_rate(*hk, *theta, *rtheta);
        
        let status = if rnorm < 0.0 {
            "SUBCRIT"
        } else if rnorm < 1.0 {
            "RAMP"
        } else {
            "SUPER"
        };
        
        println!("IBL={}: Rcrit={:.1}, Rt={:.1}, RNORM={:.4}, status={}, AX={:.4e}", 
                 ibl, rcrit, rtheta, rnorm, status, result.ax);
    }
    
    println!("\n=== Check AXSET computation station 17->18 ===");
    let (_, _, hk1, t1, rt1, ampl1) = stations[0];
    let (_, _, hk2, t2, rt2, _) = stations[1];
    let ncrit = 9.0;
    
    // Test 1: Pass the actual XFOIL ampl values
    println!("\nWith XFOIL's ampl1 = {:.6e}:", ampl1);
    let ax1 = amplification_rate(hk1, t1, rt1);
    let ax2 = amplification_rate(hk2, t2, rt2);
    println!("  Individual AX1 = {:.6e} (from DAMPL)", ax1.ax);
    println!("  Individual AX2 = {:.6e} (from DAMPL)", ax2.ax);
    
    let axset_result = axset(hk1, t1, rt1, ampl1, hk2, t2, rt2, ampl1, ncrit);
    println!("  AXSET: AX1={:.6e}, AX2={:.6e}, RMS={:.6e}, DAX={:.6e}, Final={:.6e}",
             axset_result.ax1, axset_result.ax2, axset_result.axa_rms, 
             axset_result.dax, axset_result.ax);
    
    // Test 2: Pass ampl = 0
    println!("\nWith ampl1 = 0:");
    let axset_result2 = axset(hk1, t1, rt1, 0.0, hk2, t2, rt2, 0.0, ncrit);
    println!("  AXSET: AX1={:.6e}, AX2={:.6e}, RMS={:.6e}, DAX={:.6e}, Final={:.6e}",
             axset_result2.ax1, axset_result2.ax2, axset_result2.axa_rms, 
             axset_result2.dax, axset_result2.ax);
    
    // Explain where difference is
    println!("\n=== Analysis ===");
    println!("The RMS average of AX1 and AX2 is:");
    let ax1_val = ax1.ax;
    let ax2_val = ax2.ax;
    let rms = (0.5 * (ax1_val.powi(2) + ax2_val.powi(2))).sqrt();
    println!("  sqrt(0.5 * ({:.4e}^2 + {:.4e}^2)) = {:.4e}", ax1_val, ax2_val, rms);
    
    println!("\nThe DAX correction term is:");
    let avg_ampl = 0.5 * (ampl1 + ampl1);
    let arg = (20.0 * (ncrit - avg_ampl)).min(20.0);
    let exn = (-arg).exp();
    let dax = exn * 0.002 / (t1 + t2);
    println!("  With ampl={:.2e}: ARG={:.1}, EXN={:.2e}, DAX={:.4e}", ampl1, arg, exn, dax);
    
    let avg_ampl_zero = 0.0;
    let arg_zero = (20.0 * (ncrit - avg_ampl_zero)).min(20.0);
    let exn_zero = (-arg_zero).exp();
    let dax_zero = exn_zero * 0.002 / (t1 + t2);
    println!("  With ampl=0: ARG={:.1}, EXN={:.2e}, DAX={:.4e}", arg_zero, exn_zero, dax_zero);
    
    println!("\nFinal AX = RMS + DAX");
    println!("  With XFOIL ampl: {:.4e} + {:.4e} = {:.4e}", rms, dax, rms + dax);
    println!("  With ampl=0: {:.4e} + {:.4e} = {:.4e}", rms, dax_zero, rms + dax_zero);
    
    println!("\nXFOIL's actual dN from 17->18 = {:.4e}", stations[1].5 - stations[0].5);
    println!("Our dN = AX * dx = {:.4e} * {:.6} = {:.4e}", 
             rms + dax, stations[1].1 - stations[0].1, 
             (rms + dax) * (stations[1].1 - stations[0].1));
}
