//! Solver trace test - trace through Gaussian elimination step by step

#[test]
fn test_gauss_trace() {
    // Matrix from newton_debug test iteration 0
    let mut a = [
        [1.0000e0, 0.0000e0, 0.0000e0, 0.0000e0],
        [0.0000e0, -8.2100e4, 5.5145e4, 1.6362e1],
        [4.0592e1, 0.0000e0, 1.2098e4, 4.7563e0],
        [0.0000e0, 0.0000e0, 0.0000e0, 1.0000e0],
    ];
    let mut b = [0.0, -5.6862e-1, 1.5765e-1, 0.0];

    println!("Initial matrix:");
    print_matrix(&a, &b);

    // Forward elimination step by step
    for k in 0..3 {
        println!("\n=== Step k={} ===", k);
        
        // Find pivot
        let mut max_val = a[k][k].abs();
        let mut max_row = k;
        for i in (k + 1)..4 {
            if a[i][k].abs() > max_val {
                max_val = a[i][k].abs();
                max_row = i;
            }
        }
        println!("Pivot search: column {} max at row {} (value {:.4e})", k, max_row, max_val);

        // Swap rows
        if max_row != k {
            a.swap(k, max_row);
            b.swap(k, max_row);
            println!("Swapped rows {} and {}", k, max_row);
            print_matrix(&a, &b);
        }

        let pivot = a[k][k];
        println!("Pivot = {:.4e}", pivot);

        if pivot.abs() < 1e-20 {
            println!("SINGULAR!");
            return;
        }

        // Eliminate
        for i in (k + 1)..4 {
            if a[i][k].abs() > 1e-30 {
                let factor = a[i][k] / pivot;
                println!("Eliminating row {} with factor {:.4e}", i, factor);
                a[i][k] = 0.0;
                for j in (k + 1)..4 {
                    a[i][j] -= factor * a[k][j];
                }
                b[i] -= factor * b[k];
            }
        }

        println!("After elimination:");
        print_matrix(&a, &b);
    }

    println!("\n=== Back Substitution ===");
    
    let mut x = [0.0; 4];
    
    println!("a[3][3] = {:.4e}", a[3][3]);
    x[3] = b[3] / a[3][3];
    println!("x[3] = b[3]/a[3][3] = {:.4e} / {:.4e} = {:.4e}", b[3], a[3][3], x[3]);

    println!("\na[2][2] = {:.4e}", a[2][2]);
    x[2] = (b[2] - a[2][3] * x[3]) / a[2][2];
    println!("x[2] = (b[2] - a[2][3]*x[3]) / a[2][2]");
    println!("     = ({:.4e} - {:.4e}*{:.4e}) / {:.4e}", b[2], a[2][3], x[3], a[2][2]);
    println!("     = {:.4e}", x[2]);

    println!("\na[1][1] = {:.4e}", a[1][1]);
    x[1] = (b[1] - a[1][2] * x[2] - a[1][3] * x[3]) / a[1][1];
    println!("x[1] = (b[1] - a[1][2]*x[2] - a[1][3]*x[3]) / a[1][1]");
    println!("     = ({:.4e} - {:.4e}*{:.4e} - {:.4e}*{:.4e}) / {:.4e}", 
        b[1], a[1][2], x[2], a[1][3], x[3], a[1][1]);
    println!("     = {:.4e}", x[1]);

    println!("\na[0][0] = {:.4e}", a[0][0]);
    x[0] = (b[0] - a[0][1] * x[1] - a[0][2] * x[2] - a[0][3] * x[3]) / a[0][0];
    println!("x[0] = (b[0] - a[0][1]*x[1] - a[0][2]*x[2] - a[0][3]*x[3]) / a[0][0]");
    println!("     = ({:.4e} - {:.4e}*{:.4e} - {:.4e}*{:.4e} - {:.4e}*{:.4e}) / {:.4e}", 
        b[0], a[0][1], x[1], a[0][2], x[2], a[0][3], x[3], a[0][0]);
    println!("     = {:.4e}", x[0]);

    println!("\n=== Solution ===");
    println!("x = [{:.6e}, {:.6e}, {:.6e}, {:.6e}]", x[0], x[1], x[2], x[3]);
}

fn print_matrix(a: &[[f64; 4]; 4], b: &[f64; 4]) {
    for (i, row) in a.iter().enumerate() {
        println!("  [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}] | {:+.4e}", 
            row[0], row[1], row[2], row[3], b[i]);
    }
}
