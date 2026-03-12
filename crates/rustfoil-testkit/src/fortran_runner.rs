//! FORTRAN compilation and execution utilities.
//!
//! This module provides helpers for compiling tiny FORTRAN driver programs and
//! linking them against selected XFOIL object files.

use serde::de::DeserializeOwned;
use std::ffi::OsStr;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use thiserror::Error;

use crate::paths::{
    fortran_driver_dir, gfortran_path, xfoil_instrumented_bin, xfoil_instrumented_src, xfoil_src,
};

/// Default XFOIL object set for inviscid/state-topology drivers.
pub const XFOIL_STATE_OBJS: &[&str] = &[
    "xpanel.o",
    "xfoil_debug.o",
];

/// Broader XFOIL object set for coupled/wake/workflow drivers.
pub const XFOIL_WORKFLOW_OBJS: &[&str] = &[
    "xpanel.o",
    "xoper.o",
    "xsolve.o",
    "xbl.o",
    "xblsys.o",
    "xfoil_debug.o",
    "naca.o",
    "xutils.o",
    "xgeom.o",
    "spline.o",
    "userio.o",
];

/// Error type for FORTRAN operations
#[derive(Debug, Error)]
pub enum FortranError {
    #[error("Compilation failed (exit {code:?}):\n{stderr}")]
    CompileFailed { stderr: String, code: Option<i32> },
    #[error("Execution failed (exit {code:?}):\n{stderr}")]
    RunFailed { stderr: String, code: Option<i32> },
    #[error("JSON parse failed: {0}")]
    Json(#[from] serde_json::Error),
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
}

#[derive(Debug, Clone)]
pub struct FortranDriverSpec<'a> {
    pub driver_source: &'a Path,
    pub executable_name: &'a str,
    pub object_names: &'a [&'a str],
}

pub fn compile_fortran(source_files: &[&str], output: &str) -> Result<(), FortranError> {
    let include_dir = xfoil_src();
    let mut cmd = Command::new(gfortran_path());
    cmd.arg(format!("-I{}", include_dir.display()));
    cmd.arg("-fallow-argument-mismatch");
    cmd.arg("-o").arg(output);
    for src in source_files {
        cmd.arg(src);
    }
    let output_result = cmd.output()?;
    if !output_result.status.success() {
        return Err(FortranError::CompileFailed {
            stderr: String::from_utf8_lossy(&output_result.stderr).to_string(),
            code: output_result.status.code(),
        });
    }
    Ok(())
}

pub fn ensure_xfoil_objects(object_names: &[&str]) -> Result<Vec<PathBuf>, FortranError> {
    let bin_dir = xfoil_instrumented_bin();
    std::fs::create_dir_all(fortran_driver_dir())?;

    let missing = object_names
        .iter()
        .any(|name| !bin_dir.join(name).exists());
    if missing {
        let mut cmd = Command::new("make");
        cmd.arg("-C").arg(&bin_dir);
        for name in object_names {
            cmd.arg(name);
        }
        let output = cmd.output()?;
        if !output.status.success() {
            return Err(FortranError::CompileFailed {
                stderr: String::from_utf8_lossy(&output.stderr).to_string(),
                code: output.status.code(),
            });
        }
    }

    Ok(object_names.iter().map(|name| bin_dir.join(name)).collect())
}

pub fn compile_driver(spec: &FortranDriverSpec<'_>) -> Result<PathBuf, FortranError> {
    let exe = fortran_driver_dir().join(spec.executable_name);
    let object_paths = ensure_xfoil_objects(spec.object_names)?;
    let include_dir = prepare_freeform_include_dir()?;

    let mut cmd = Command::new(gfortran_path());
    cmd.arg("-O2")
        .arg("-fdefault-real-8")
        .arg("-std=legacy")
        .arg("-fallow-argument-mismatch")
        .arg(format!("-I{}", include_dir.display()))
        .arg("-o")
        .arg(&exe)
        .arg(spec.driver_source);
    for object in &object_paths {
        cmd.arg(object);
    }
    #[cfg(target_os = "macos")]
    {
        cmd.arg("-Wl,-undefined,dynamic_lookup");
    }
    cmd.arg("-lm");

    let output = cmd.output()?;
    if !output.status.success() {
        return Err(FortranError::CompileFailed {
            stderr: format_output(&output),
            code: output.status.code(),
        });
    }

    Ok(exe)
}

pub fn run_fortran_test(executable: impl AsRef<Path>) -> Result<String, FortranError> {
    let output = Command::new(executable.as_ref()).output()?;
    if !output.status.success() {
        return Err(FortranError::RunFailed {
            stderr: format_output(&output),
            code: output.status.code(),
        });
    }
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

pub fn run_fortran_with_args<I, S>(
    executable: impl AsRef<Path>,
    args: I,
    working_directory: Option<&Path>,
) -> Result<String, FortranError>
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let mut cmd = Command::new(executable.as_ref());
    cmd.args(args);
    if let Some(cwd) = working_directory {
        cmd.current_dir(cwd);
    }
    let output = cmd.output()?;
    if !output.status.success() {
        return Err(FortranError::RunFailed {
            stderr: format_output(&output),
            code: output.status.code(),
        });
    }
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

pub fn run_fortran_json<T: DeserializeOwned>(executable: impl AsRef<Path>) -> Result<T, FortranError> {
    let stdout = run_fortran_test(executable)?;
    Ok(serde_json::from_str(&stdout)?)
}

pub fn run_fortran_json_with_args<T, I, S>(
    executable: impl AsRef<Path>,
    args: I,
    working_directory: Option<&Path>,
) -> Result<T, FortranError>
where
    T: DeserializeOwned,
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let stdout = run_fortran_with_args(executable, args, working_directory)?;
    Ok(serde_json::from_str(&stdout)?)
}

pub fn compile_and_run(source_files: &[&str], executable: &str) -> Result<String, FortranError> {
    compile_fortran(source_files, executable)?;
    run_fortran_test(executable)
}

pub fn run_xfoil_instrumented(commands: &str, working_directory: &Path) -> Result<String, FortranError> {
    let executable = xfoil_instrumented_bin().join("xfoil_instrumented");
    let debug_path = working_directory.join("xfoil_debug.json");
    let fallback_debug_path = std::env::current_dir()
        .ok()
        .map(|cwd| cwd.join("xfoil_debug.json"));
    let _ = std::fs::remove_file(&debug_path);
    if let Some(path) = &fallback_debug_path {
        if path != &debug_path {
            let _ = std::fs::remove_file(path);
        }
    }
    let output = Command::new(executable)
        .current_dir(working_directory)
        .env("XFOIL_DEBUG_JSON", &debug_path)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .and_then(|mut child| {
            use std::io::Write;
            if let Some(stdin) = child.stdin.as_mut() {
                stdin.write_all(commands.as_bytes())?;
            }
            child.wait_with_output()
        })?;
    if !output.status.success() {
        return Err(FortranError::RunFailed {
            stderr: format_output(&output),
            code: output.status.code(),
        });
    }
    if !debug_path.exists() {
        if let Some(path) = fallback_debug_path {
            if path.exists() {
                let _ = std::fs::rename(&path, &debug_path);
                if !debug_path.exists() {
                    let _ = std::fs::copy(&path, &debug_path);
                }
            }
        }
    }
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

fn format_output(output: &Output) -> String {
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    if stdout.trim().is_empty() {
        stderr.to_string()
    } else if stderr.trim().is_empty() {
        stdout.to_string()
    } else {
        format!("stdout:\n{}\n\nstderr:\n{}", stdout, stderr)
    }
}

fn prepare_freeform_include_dir() -> Result<PathBuf, FortranError> {
    let source_dir = xfoil_instrumented_src();
    let include_dir = fortran_driver_dir().join("freeform-includes");
    std::fs::create_dir_all(&include_dir)?;
    convert_fixed_include(&source_dir.join("PINDEX.INC"), &include_dir.join("PINDEX.INC"))?;
    convert_fixed_include(&source_dir.join("XFOIL.INC"), &include_dir.join("XFOIL.INC"))?;
    Ok(include_dir)
}

fn convert_fixed_include(source: &Path, destination: &Path) -> Result<(), FortranError> {
    let contents = std::fs::read_to_string(source)?;
    let converted = fixed_to_freeform(&contents);
    std::fs::write(destination, converted)?;
    Ok(())
}

fn fixed_to_freeform(contents: &str) -> String {
    let mut output = String::new();
    let mut current_stmt: Option<String> = None;

    let flush_stmt = |output: &mut String, current_stmt: &mut Option<String>| {
        if let Some(stmt) = current_stmt.take() {
            output.push_str(stmt.trim_end());
            output.push('\n');
        }
    };

    for raw_line in contents.lines() {
        if raw_line.trim().is_empty() {
            flush_stmt(&mut output, &mut current_stmt);
            output.push('\n');
            continue;
        }

        let first = raw_line.chars().next().unwrap_or(' ');
        if matches!(first, 'c' | 'C' | '*') {
            flush_stmt(&mut output, &mut current_stmt);
            output.push('!');
            output.push_str(raw_line.get(1..).unwrap_or("").trim_end());
            output.push('\n');
            continue;
        }

        let continuation = raw_line.chars().nth(5).unwrap_or(' ') != ' ';
        let text = if raw_line.len() > 6 { &raw_line[6..] } else { "" };
        let text = text.split('!').next().unwrap_or("").trim_end();
        if continuation {
            if let Some(stmt) = current_stmt.as_mut() {
                stmt.push(' ');
                stmt.push_str(text.trim());
            } else {
                current_stmt = Some(text.trim().to_string());
            }
        } else {
            flush_stmt(&mut output, &mut current_stmt);
            current_stmt = Some(text.to_string());
        }
    }

    if let Some(stmt) = current_stmt.take() {
        output.push_str(stmt.trim_end());
        output.push('\n');
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paths::{gfortran_path, xfoil_src};

    #[test]
    fn test_gfortran_exists() {
        let exists = gfortran_path().exists();
        if !exists {
            eprintln!("Warning: gfortran not found at {}", gfortran_path().display());
        }
    }

    #[test]
    fn test_xfoil_src_exists() {
        let path = xfoil_src();
        assert!(path.exists(), "XFOIL source directory not found at {}", path.display());
    }
}
