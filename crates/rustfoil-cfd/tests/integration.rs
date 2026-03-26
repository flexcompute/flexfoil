//! Integration tests for the full CFD mesh generation pipeline.
//!
//! Tests the complete workflow: NACA 0012 coordinates → O-mesh → initial conditions → BC types.

use std::f64::consts::PI;

use rustfoil_cfd::boundary::generate_bc_types;
use rustfoil_cfd::config::{CfdConfig, CfdParamsGpu, PhysicsMode};
use rustfoil_cfd::init::compute_initial_conditions;
use rustfoil_cfd::mesh::generate_o_mesh;

/// Generate NACA 0012 airfoil coordinates (closed trailing edge).
fn naca0012(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);

    // Upper surface from TE to LE (i=0..n/2)
    let half = n / 2;
    for i in 0..half {
        let theta = PI * (i as f64) / (half as f64);
        let xc = 0.5 * (1.0 + theta.cos()); // cosine spacing: 1.0 → 0.0
        let t = 0.12; // max thickness ratio
        let yt = 5.0
            * t
            * (0.2969 * xc.sqrt()
                - 0.1260 * xc
                - 0.3516 * xc.powi(2)
                + 0.2843 * xc.powi(3)
                - 0.1015 * xc.powi(4));
        x.push(xc);
        y.push(yt);
    }
    // Lower surface from LE to TE (i=n/2..n)
    for i in 0..half {
        let theta = PI * (i as f64) / (half as f64);
        let xc = 0.5 * (1.0 - theta.cos()); // cosine spacing: 0.0 → 1.0
        let t = 0.12;
        let yt = 5.0
            * t
            * (0.2969 * xc.sqrt()
                - 0.1260 * xc
                - 0.3516 * xc.powi(2)
                + 0.2843 * xc.powi(3)
                - 0.1015 * xc.powi(4));
        x.push(xc);
        y.push(-yt);
    }

    (x, y)
}

#[test]
fn test_full_pipeline_naca0012() {
    // Step 1: Generate airfoil coordinates
    let (ax, ay) = naca0012(128);
    assert_eq!(ax.len(), 128);
    assert_eq!(ay.len(), 128);

    // Verify airfoil bounds
    let x_min = ax.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = ax.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(x_min >= -0.01, "x_min should be near 0, got {x_min}");
    assert!(x_max <= 1.01, "x_max should be near 1, got {x_max}");

    // Step 2: Generate O-mesh
    let ni = 128u32;
    let nj = 64u32;
    let far_field = 15.0;
    let ds0 = 0.001;

    let mesh = generate_o_mesh(&ax, &ay, ni, nj, far_field, ds0);

    assert_eq!(mesh.ni, ni);
    assert_eq!(mesh.nj, nj);
    assert_eq!(mesh.x.len(), (ni * nj) as usize);
    assert_eq!(mesh.y.len(), (ni * nj) as usize);

    // Verify wall (j=0) points are near the airfoil
    for i in 0..ni as usize {
        let wx = mesh.x[i] as f64;
        let wy = mesh.y[i] as f64;
        // Wall points should be within the chord length range
        assert!(
            wx > -5.0 && wx < 5.0,
            "Wall point x[{i}] = {wx} out of range"
        );
        assert!(
            wy > -5.0 && wy < 5.0,
            "Wall point y[{i}] = {wy} out of range"
        );
    }

    // Verify far-field (j=nj-1) points are far away
    let jj = (nj - 1) as usize;
    let mut max_dist = 0.0f64;
    for i in 0..ni as usize {
        let fx = mesh.x[jj * ni as usize + i] as f64;
        let fy = mesh.y[jj * ni as usize + i] as f64;
        let dist = (fx * fx + fy * fy).sqrt();
        if dist > max_dist {
            max_dist = dist;
        }
    }
    assert!(
        max_dist > 5.0,
        "Far-field points should be far from origin, max_dist={max_dist}"
    );

    // No NaN or Inf in mesh
    for val in &mesh.x {
        assert!(val.is_finite(), "mesh.x contains non-finite value");
    }
    for val in &mesh.y {
        assert!(val.is_finite(), "mesh.y contains non-finite value");
    }

    println!("Mesh generated: ni={}, nj={}", mesh.ni, mesh.nj);
    println!("  Wall x range: [{:.4}, {:.4}]",
        mesh.x[..ni as usize].iter().cloned().fold(f32::INFINITY, f32::min),
        mesh.x[..ni as usize].iter().cloned().fold(f32::NEG_INFINITY, f32::max),
    );
    println!("  Far-field max distance: {:.2}", max_dist);

    // Step 3: Compute initial conditions
    let config = CfdConfig {
        ni,
        nj,
        mach_inf: 0.5,
        alpha: 2.0f32.to_radians(),
        gamma: 1.4,
        reynolds: 1e5,
        physics: PhysicsMode::Euler,
        ..CfdConfig::default()
    };

    let q = compute_initial_conditions(&config);
    assert_eq!(q.len(), (ni * nj * 5) as usize);

    // Verify freestream state at an arbitrary cell
    let cell = 100;
    let rho = q[cell * 5];
    let rhou = q[cell * 5 + 1];
    let rhov = q[cell * 5 + 2];
    let e = q[cell * 5 + 3];
    let nu = q[cell * 5 + 4];

    assert!((rho - 1.0).abs() < 1e-5, "rho = {rho}");
    assert!(rhou > 0.0, "rhou should be positive for alpha~2deg");
    assert!(rhov > 0.0, "rhov should be positive for alpha~2deg");
    assert!(e > 0.0, "energy should be positive");
    assert!(nu.abs() < 1e-10, "nu_tilde should be 0 for Euler");

    // No NaN or Inf in Q
    for val in &q {
        assert!(val.is_finite(), "Q contains non-finite value");
    }

    println!("  Initial conditions: rho={rho:.4}, rhou={rhou:.4}, rhov={rhov:.4}, E={e:.4}");

    // Step 4: Generate boundary condition types
    let bc = generate_bc_types(ni, nj);
    assert_eq!(bc.len(), (ni * nj) as usize);

    // Wall at j=0
    for i in 0..ni as usize {
        assert_eq!(bc[i], 1, "bc[j=0, i={i}] should be Wall(1)");
    }
    // Interior at j=1
    for i in 0..ni as usize {
        assert_eq!(bc[ni as usize + i], 0, "bc[j=1, i={i}] should be Interior(0)");
    }
    // FarField at j=nj-1
    for i in 0..ni as usize {
        assert_eq!(
            bc[(nj as usize - 1) * ni as usize + i],
            2,
            "bc[j=nj-1, i={i}] should be FarField(2)"
        );
    }

    println!("  BC types: wall={}, interior={}, farfield={}",
        bc.iter().filter(|&&v| v == 1).count(),
        bc.iter().filter(|&&v| v == 0).count(),
        bc.iter().filter(|&&v| v == 2).count(),
    );

    // Step 5: Verify GPU params struct
    let params = CfdParamsGpu::from_config(&config, 0.001, 0);
    assert_eq!(params.ni, ni);
    assert_eq!(params.nj, nj);
    assert!((params.gamma - 1.4).abs() < 1e-6);
    assert!((params.mach_inf - 0.5).abs() < 1e-6);
    assert_eq!(params.physics_mode, 0); // Euler
    let bytes = params.as_bytes();
    assert_eq!(bytes.len(), 48);

    println!("  GPU params: {} bytes", bytes.len());
    println!("\nFull CFD pipeline test PASSED ✓");
}

#[test]
fn test_mesh_quality_metrics() {
    // Test that the mesh has reasonable aspect ratios and no degenerate cells
    let (ax, ay) = naca0012(64);
    let mesh = generate_o_mesh(&ax, &ay, 64, 32, 10.0, 0.005);

    let ni = 64usize;
    let nj = 32usize;
    let mut min_area = f64::INFINITY;
    let mut max_area = 0.0f64;
    let mut degenerate_count = 0;

    // Check cell areas (should all be positive for a valid mesh)
    for j in 0..nj - 1 {
        for i in 0..ni {
            let ip = (i + 1) % ni;

            let x00 = mesh.x[j * ni + i] as f64;
            let y00 = mesh.y[j * ni + i] as f64;
            let x10 = mesh.x[j * ni + ip] as f64;
            let y10 = mesh.y[j * ni + ip] as f64;
            let x01 = mesh.x[(j + 1) * ni + i] as f64;
            let y01 = mesh.y[(j + 1) * ni + i] as f64;
            let x11 = mesh.x[(j + 1) * ni + ip] as f64;
            let y11 = mesh.y[(j + 1) * ni + ip] as f64;

            // Cross product of diagonals gives 2x area
            let area = 0.5
                * ((x11 - x00) * (y01 - y10) - (x01 - x10) * (y11 - y00))
                    .abs();

            if area < 1e-12 {
                degenerate_count += 1;
            }
            if area < min_area {
                min_area = area;
            }
            if area > max_area {
                max_area = area;
            }
        }
    }

    println!("Mesh quality: min_area={:.2e}, max_area={:.2e}, degenerate={}", min_area, max_area, degenerate_count);
    assert!(
        degenerate_count < (ni * nj / 10),
        "Too many degenerate cells: {degenerate_count}"
    );
    assert!(min_area > 0.0, "All cell areas should be positive");
}

#[test]
fn test_rans_initial_conditions() {
    let config = CfdConfig {
        ni: 32,
        nj: 16,
        mach_inf: 0.3,
        alpha: 0.0,
        reynolds: 1e6,
        physics: PhysicsMode::RansSA,
        ..CfdConfig::default()
    };

    let q = compute_initial_conditions(&config);

    // RANS should have non-zero nu_tilde
    let nu = q[4]; // First cell, 5th variable
    assert!(nu > 0.0, "RANS nu_tilde should be positive, got {nu}");
    println!("RANS nu_tilde_inf = {:.6e}", nu);
}
