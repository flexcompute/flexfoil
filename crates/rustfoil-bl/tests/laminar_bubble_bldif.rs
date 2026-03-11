use rustfoil_bl::equations::bldif;
use rustfoil_bl::{blvar, BlStation, FlowType};

fn station(x: f64, u: f64, theta: f64, delta_star: f64, ampl: f64) -> BlStation {
    let mut s = BlStation::new();
    s.x = x;
    s.u = u;
    s.theta = theta;
    s.delta_star = delta_star;
    s.ampl = ampl;
    s.ctau = 0.03;
    s.is_laminar = true;
    blvar(&mut s, FlowType::Laminar, 0.0, 1.0e6);
    s
}

#[test]
fn test_bldif_matches_xfoil_laminar_bubble_station_24_iter_2() {
    // XFOIL alpha=10, upper surface, ibl=24, Newton iter 2.
    // These values come from the raw XFOIL `BLDIF` / `VS2_BEFORE` event,
    // not the post-solve `MRCHUE_ITER` snapshot.
    let s1 = station(0.056034968, 2.6023193, 5.9971693e-05, 1.8303967e-04, 5.4085861e-01);
    let s2 = station(0.058776452, 2.5467323, 6.7633870e-05, 2.2670739e-04, 1.1273327);

    let (res, jac) = bldif(&s1, &s2, FlowType::Laminar, 0.0, 1.0e6);

    let expect_vsrez: [f64; 3] = [2.356789e-17, 1.102936e-2, -5.925047e-3];
    let expect_vs2 = [
        [1.0_f64, 2.779048e4, -6.578308e3, 0.0],
        [0.0, 1.462092e4, 1.171022e2, 2.045651e0],
        [0.0, 1.591092e3, -3.038229e2, -8.570490e-1],
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
