use rustfoil_bl::{bldif, blvar, BlStation, FlowType};

fn laminar_station(x: f64, u: f64, theta: f64, delta_star: f64, ampl: f64) -> BlStation {
    let mut s = BlStation::new();
    s.x = x;
    s.u = u;
    s.theta = theta;
    s.delta_star = delta_star;
    s.ctau = 0.03;
    s.ampl = ampl;
    s.is_laminar = true;
    blvar(&mut s, FlowType::Laminar, 0.0, 1.0e6);
    s
}

fn turbulent_station(
    x: f64,
    u: f64,
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
) -> BlStation {
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

fn assert_close_matrix(label: &str, got: [[f64; 5]; 3], expected: [[f64; 4]; 3]) {
    for row in 0..3 {
        for col in 0..4 {
            let exp = expected[row][col];
            let val = got[row][col];
            let tol = exp.abs().max(1.0) * 5.0e-2;
            assert!(
                (val - exp).abs() <= tol,
                "{label}[{row}][{col}] expected {:.6e}, got {:.6e}",
                exp,
                val
            );
        }
    }
}

#[test]
fn test_laminar_half_bldif_matches_xfoil_transition29_iter2() {
    let s1 = laminar_station(0.072368748, 2.4713247, 8.599435e-05, 8.4371111e-04, 8.847934);
    let s2 = laminar_station(0.072646974, 2.4564394, 8.599435e-05, 8.4371111e-04, 9.0);

    let (res, jac) = bldif(&s1, &s2, FlowType::Laminar, 0.0, 1.0e6);

    let expect_vsrez: [f64; 3] = [1.690677e-16, 7.089791e-2, -5.010513e-2];
    let expect_vs1: [[f64; 4]; 3] = [
        [-1.007072e0, 9.515416e2, -6.867505e0, 0.0],
        [0.0, -1.128899e4, -3.516519e1, -4.779408e0],
        [0.0, 5.592007e2, -5.330745e1, 3.566021e0],
    ];
    let expect_vs2: [[f64; 4]; 3] = [
        [9.929285e-1, 9.515416e2, -6.867505e0, 0.0],
        [0.0, 1.196834e4, -3.516538e1, 4.808183e0],
        [0.0, -1.188036e3, 1.248144e2, -3.586357e0],
    ];

    let got_vsrez = [res.res_third, res.res_mom, res.res_shape];
    for i in 0..3 {
        let tol = expect_vsrez[i].abs().max(1.0) * 5.0e-2;
        assert!(
            (got_vsrez[i] - expect_vsrez[i]).abs() <= tol,
            "laminar vsrez[{i}] expected {:.6e}, got {:.6e}",
            expect_vsrez[i],
            got_vsrez[i]
        );
    }

    assert_close_matrix("laminar vs1", jac.vs1, expect_vs1);
    assert_close_matrix("laminar vs2", jac.vs2, expect_vs2);
}

#[test]
fn test_turbulent_half_bldif_matches_xfoil_transition29_iter2() {
    let s1 = turbulent_station(0.072646974, 2.4564394, 8.599435e-05, 8.4371111e-04, 1.4943847e-1, 0.0);
    let s2 = turbulent_station(0.076690641, 2.2400991, 8.599435e-05, 8.4371111e-04, 3.0e-2, 1.1210092e1);

    let (res, jac) = bldif(&s1, &s2, FlowType::Turbulent, 0.0, 1.0e6);

    let expect_vsrez: [f64; 3] = [-4.380881e-3, 1.083738e0, 1.937351e-2];
    let expect_vs1: [[f64; 4]; 3] = [
        [-4.113472e-3, 1.428533e1, 6.363093e-1, 8.455678e-4],
        [0.0, -6.399153e3, -5.360373e2, -4.808276e0],
        [-9.855429e0, 5.204126e3, 4.130541e2, 3.599056e0],
    ];
    let expect_vs2: [[f64; 4]; 3] = [
        [-8.671975e-2, 1.437849e1, 6.386148e-1, -9.109802e-4],
        [0.0, 1.685737e4, -5.360371e2, 5.272641e0],
        [-2.087836e0, -5.046076e3, 6.292629e2, -3.919468e0],
    ];

    let got_vsrez = [res.res_third, res.res_mom, res.res_shape];
    for i in 0..3 {
        let tol = expect_vsrez[i].abs().max(1.0) * 5.0e-2;
        assert!(
            (got_vsrez[i] - expect_vsrez[i]).abs() <= tol,
            "turbulent vsrez[{i}] expected {:.6e}, got {:.6e}",
            expect_vsrez[i],
            got_vsrez[i]
        );
    }

    assert_close_matrix("turbulent vs1", jac.vs1, expect_vs1);
    assert_close_matrix("turbulent vs2", jac.vs2, expect_vs2);
}
