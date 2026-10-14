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
