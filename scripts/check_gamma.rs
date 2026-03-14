// Quick test to print gamma values
use rustfoil_inviscid::{InviscidSolver, FlowConditions};
use std::f64::consts::PI;

fn make_naca0012(n_panels: usize) -> Vec<(f64, f64)> {
    let t = 0.12;
    let n_half = n_panels / 2;
    let x_coords: Vec<f64> = (0..=n_half)
        .map(|i| {
            let beta = PI * (i as f64) / (n_half as f64);
            0.5 * (1.0 - beta.cos())
        })
        .collect();
    
    let thickness = |x: f64| -> f64 {
        5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2)
            + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
    };
    
    let mut points = Vec::with_capacity(n_panels);
    for i in (0..=n_half).rev() {
        let x = x_coords[i];
        let y = thickness(x);
        points.push((x, y));
    }
    for i in 1..n_half {
        let x = x_coords[i];
        let y = -thickness(x);
        points.push((x, y));
    }
    points
}

fn main() {
    let points = make_naca0012(160);
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&points).expect("Factorization failed");
    
    let flow = FlowConditions::with_alpha_deg(4.0);
    let solution = factorized.solve_alpha(&flow);
    
    println!("Gamma at alpha=4:");
    println!("  max: {:.4}", solution.gamma.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    println!("  min: {:.4}", solution.gamma.iter().cloned().fold(f64::INFINITY, f64::min));
    
    let max_idx = solution.gamma.iter().enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap().0;
    println!("  max at index: {}", max_idx);
    println!("  gamma[{}]: {:.4}", max_idx, solution.gamma[max_idx]);
    
    println!("\nFirst 10 gamma: {:?}", &solution.gamma[..10]);
    println!("Around max ({}..{}): {:?}", max_idx.saturating_sub(3), max_idx+4, 
             &solution.gamma[max_idx.saturating_sub(3)..max_idx+4]);
}
