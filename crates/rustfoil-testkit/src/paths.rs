use std::path::PathBuf;

fn env_path(name: &str) -> Option<PathBuf> {
    std::env::var_os(name).map(PathBuf::from)
}

pub fn workspace_root() -> PathBuf {
    if let Some(path) = env_path("RUSTFOIL_WORKSPACE_ROOT") {
        return path;
    }

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .map(PathBuf::from)
        .expect("rustfoil-testkit should live under crates/")
}

pub fn target_dir() -> PathBuf {
    workspace_root().join("target")
}

pub fn fortran_driver_dir() -> PathBuf {
    target_dir().join("fortran-drivers")
}

pub fn gfortran_path() -> PathBuf {
    env_path("RUSTFOIL_GFORTRAN").unwrap_or_else(|| PathBuf::from("/opt/homebrew/bin/gfortran"))
}

pub fn xfoil_src() -> PathBuf {
    env_path("RUSTFOIL_XFOIL_SRC")
        .unwrap_or_else(|| workspace_root().join("Xfoil").join("src"))
}

pub fn xfoil_instrumented_root() -> PathBuf {
    env_path("RUSTFOIL_XFOIL_INSTRUMENTED_ROOT")
        .unwrap_or_else(|| workspace_root().join("Xfoil-instrumented"))
}

pub fn xfoil_instrumented_src() -> PathBuf {
    env_path("RUSTFOIL_XFOIL_INSTRUMENTED_SRC")
        .unwrap_or_else(|| xfoil_instrumented_root().join("src"))
}

pub fn xfoil_instrumented_bin() -> PathBuf {
    env_path("RUSTFOIL_XFOIL_INSTRUMENTED_BIN")
        .unwrap_or_else(|| xfoil_instrumented_root().join("bin"))
}
