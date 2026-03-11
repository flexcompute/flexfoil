use rustfoil_bl::{bldif, blvar, BlStation, FlowType};

fn station(x: f64, u: f64, theta: f64, delta_star: f64, ctau: f64, ampl: f64) -> BlStation {
    let mut s = BlStation::new();
    s.x = x;
    s.u = u;
    s.theta = theta;
    s.delta_star = delta_star;
    s.ctau = ctau;
    s.ampl = ampl;
    s.is_laminar = false;
    s.is_turbulent = true;
    blvar(&mut s, FlowType::Turbulent, 0.0, 1.0e6);
    s
}

#[test]
fn test_bldif_matches_xfoil_turbulent_station_30_iter_2() {
    // XFOIL alpha=10, upper surface, ibl=30, Newton iter 2.
    // Reference values come from the raw XFOIL `BLDIF` event.
    let s1 = station(
        0.076690641,
        2.2560485,
        1.7974431e-4,
        5.1315745e-4,
        1.0436721e-1,
        1.1119815e1,
    );
    let s2 = station(
        0.081532629,
        2.1805929,
        1.9430622e-4,
        3.5921021e-4,
        8.7001313e-2,
        1.1119815e1,
    );

    let (res, jac) = bldif(&s1, &s2, FlowType::Turbulent, 0.0, 1.0e6);

    let expect_vsrez: [f64; 3] = [4.59314e-4, 9.074668e-2, 5.361225e-2];
    let expect_vs1: [[f64; 4]; 3] = [
        [1.892881e-2, 6.568537e-1, 3.310339e-1, 1.157667e-3],
        [0.0, -5.384414e3, -3.605573e1, -1.925901e0],
        [-1.289999e0, -9.309027e2, 5.090474e2, 6.108043e-1],
    ];
    let expect_vs2: [[f64; 4]; 3] = [
        [-5.371316e-2, -3.751163e0, 2.789894e0, -1.199010e-3],
        [0.0, 5.282089e3, 3.315815e0, 2.001993e0],
        [-2.863031e0, 3.662854e3, -1.823628e3, -6.448377e-1],
    ];

    let got_vsrez = [res.res_third, res.res_mom, res.res_shape];
    for i in 0..3 {
        let tol = expect_vsrez[i].abs().max(1.0) * 5.0e-2;
        assert!(
            (got_vsrez[i] - expect_vsrez[i]).abs() <= tol,
            "vsrez[{i}] expected {:.6e}, got {:.6e}",
            expect_vsrez[i],
            got_vsrez[i]
        );
    }

    for row in 0..3 {
        for col in 0..4 {
            let expected = expect_vs1[row][col];
            let got = jac.vs1[row][col];
            let tol = expected.abs().max(1.0) * 5.0e-2;
            assert!(
                (got - expected).abs() <= tol,
                "vs1[{row}][{col}] expected {:.6e}, got {:.6e}",
                expected,
                got
            );
        }
    }

    for row in 0..3 {
        for col in 0..4 {
            let expected = expect_vs2[row][col];
            let got = jac.vs2[row][col];
            let tol = expected.abs().max(1.0) * 5.0e-2;
            assert!(
                (got - expected).abs() <= tol,
                "vs2[{row}][{col}] expected {:.6e}, got {:.6e}",
                expected,
                got
            );
        }
    }
}
