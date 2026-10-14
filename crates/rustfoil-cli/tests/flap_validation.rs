use std::process::Command;
use std::path::Path;

#[test]
fn invalid_flaps_fail_before_any_polar_or_file_access() {
    for mode in ["polar", "faithful-polar"] {
        for spec in ["0.72:0.5", "0.72:0.5:12deg", "nan:0.5:12", "0.72:inf:12", "0.72:0.5:nan", "0:0.5:12", "1:0.5:12", "-0.1:0.5:12", "0.72:0.5:inf", "0.72:1.5:12", "0.72:-0.1:12"] {
            let output = Command::new(env!("CARGO_BIN_EXE_rustfoil"))
                .args([mode, "/file-that-does-not-exist.dat", "--flap", spec]).output().unwrap();
            assert!(!output.status.success(), "{mode}: {spec}");
            assert!(output.stdout.is_empty(), "must not print an unflapped polar");
            let error = String::from_utf8_lossy(&output.stderr);
            assert!(error.contains("finite hinge_x"), "{mode}: {spec}: {error}");
        }
    }
}

#[test]
fn valid_flap_changes_polar_and_zero_deflection_preserves_it() {
    let foil = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/naca0012.dat");
    let run = |flap: Option<&str>| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rustfoil"));
        command.arg("polar").arg(&foil).args(["--alpha-start", "4", "--alpha-end", "4"]);
        if let Some(spec) = flap { command.args(["--flap", spec]); }
        let output = command.output().unwrap();
        assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
        output.stdout
    };
    assert_eq!(run(None), run(Some("0.72:0.5:0")));
    assert_ne!(run(None), run(Some("0.72:0.5:12")));
}

#[test]
fn nonzero_flap_on_a_short_contour_is_rejected() {
    let path = std::env::temp_dir().join(format!("rustfoil-short-flap-{}.dat", std::process::id()));
    std::fs::write(&path,"short\n1 .01\n.5 .1\n0 0\n.5 -.1\n1 -.01\n").unwrap();
    for mode in ["polar", "faithful-polar"] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rustfoil"));
        command.arg(mode).arg(&path).args(["--flap","0.72:0.5:12"]);
        if mode == "faithful-polar" { command.args(["--alpha","0:0:1","--re","1000000"]); }
        let output = command.output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("at least six contour points"));
    }
    std::fs::remove_file(path).unwrap();
}

#[test]
fn inviscid_cli_flap_matches_the_healed_sdk_geometry() {
    use rustfoil_core::{point, Body, CubicSpline, PanelingParams};
    use rustfoil_solver::inviscid::{FlowConditions, InviscidSolver};
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/naca0012.dat");
    let text=std::fs::read_to_string(&path).unwrap();
    let points:Vec<_>=text.lines().skip(1).filter(|line| !line.trim().is_empty()).map(|line| {
        let values:Vec<f64>=line.split_whitespace().map(|value|value.parse().unwrap()).collect();
        point(values[0],values[1])
    }).collect();
    let buffer=rustfoil_core::flap::xfoil_flap(&points,0.72,0.5,12.);
    let healed=CubicSpline::from_points(&buffer).unwrap().resample_xfoil(points.len(),&PanelingParams::default());
    let body=Body::from_points("SDK reference",&healed).unwrap();
    let expected=InviscidSolver::new().solve(&[body],&FlowConditions::with_alpha_deg(4.)).unwrap();
    let output=Command::new(env!("CARGO_BIN_EXE_rustfoil")).arg("polar").arg(path).args(["--alpha-start","4","--alpha-end","4","--flap","0.72:0.5:12"]).output().unwrap();
    assert!(output.status.success());
    let stdout=String::from_utf8(output.stdout).unwrap();
    let values:Vec<f64>=stdout.lines().last().unwrap().split_whitespace().map(|value|value.parse().unwrap()).collect();
    assert!((values[1]-expected.cl).abs()<5.1e-5);
    assert!((values[2]-expected.cm).abs()<5.1e-5);
}
