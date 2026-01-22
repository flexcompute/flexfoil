//! Solver debug test

use rustfoil_coupling::solve::solve_4x4;

#[test]
fn test_solve_4x4_with_test_matrix() {
    // Matrix from newton_debug test iteration 0
    let a = [
        [1.0000e0, 0.0000e0, 0.0000e0, 0.0000e0],
        [0.0000e0, -8.2100e4, 5.5145e4, 1.6362e1],
        [4.0592e1, 0.0000e0, 1.2098e4, 4.7563e0],
        [0.0000e0, 0.0000e0, 0.0000e0, 1.0000e0],
    ];
    let b = [0.0, -5.6862e-1, 1.5765e-1, 0.0];

    println!("Input matrix A:");
    for (i, row) in a.iter().enumerate() {
        println!("  A[{}] = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", 
            i, row[0], row[1], row[2], row[3]);
    }
    println!("Input b = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", b[0], b[1], b[2], b[3]);

    let x = solve_4x4(&a, &b);
    
    println!("\nSolution x = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", x[0], x[1], x[2], x[3]);
    
    // Verify: compute A*x - b
    println!("\nVerification (A*x - b):");
    for i in 0..4 {
        let ax_i: f64 = (0..4).map(|j| a[i][j] * x[j]).sum();
        println!("  Row {}: A*x = {:.4e}, b = {:.4e}, residual = {:.4e}", 
            i, ax_i, b[i], ax_i - b[i]);
    }

    // Expected: x ≈ [0, 1.57e-5, 1.30e-5, 0]
    println!("\nExpected (from NumPy):");
    println!("  x = [0, 1.57e-5, 1.30e-5, 0]");
    
    // Check correctness
    assert!((x[0]).abs() < 1e-10, "x[0] should be ~0");
    assert!((x[1] - 1.57e-5).abs() < 1e-6, "x[1] should be ~1.57e-5, got {}", x[1]);
    assert!((x[2] - 1.30e-5).abs() < 1e-6, "x[2] should be ~1.30e-5, got {}", x[2]);
    assert!((x[3]).abs() < 1e-10, "x[3] should be ~0");
}
