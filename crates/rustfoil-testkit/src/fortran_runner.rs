//! FORTRAN compilation and execution utilities.
//!
//! This module provides functions to compile and run FORTRAN test programs
//! against the XFOIL source code.

use std::io;
use std::path::Path;
use std::process::{Command, Output};

/// Path to gfortran compiler (Homebrew on macOS ARM)
pub const GFORTRAN: &str = "/opt/homebrew/bin/gfortran";

/// Path to XFOIL source directory for includes
pub const XFOIL_SRC: &str = "/Users/harry/flexfoil-boundary-layer/Xfoil/src";

/// Error type for FORTRAN operations
#[derive(Debug)]
pub enum FortranError {
    /// Compilation failed
    CompileFailed { stderr: String, code: Option<i32> },
    /// Execution failed
    RunFailed { stderr: String, code: Option<i32> },
    /// IO error
    Io(io::Error),
}

impl std::fmt::Display for FortranError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FortranError::CompileFailed { stderr, code } => {
                write!(f, "Compilation failed (exit {:?}):\n{}", code, stderr)
            }
            FortranError::RunFailed { stderr, code } => {
                write!(f, "Execution failed (exit {:?}):\n{}", code, stderr)
            }
            FortranError::Io(e) => write!(f, "IO error: {}", e),
        }
    }
}

impl std::error::Error for FortranError {}

impl From<io::Error> for FortranError {
    fn from(e: io::Error) -> Self {
        FortranError::Io(e)
    }
}

/// Compile FORTRAN source files into an executable.
///
/// # Arguments
/// * `source_files` - List of FORTRAN source file paths
/// * `output` - Path for the output executable
///
/// # Example
/// ```ignore
/// compile_fortran(&["test_closures.f", "xblsys.f"], "test_closures")?;
/// ```
pub fn compile_fortran(source_files: &[&str], output: &str) -> Result<(), FortranError> {
    let mut cmd = Command::new(GFORTRAN);
    
    // Add XFOIL include path
    cmd.arg(format!("-I{}", XFOIL_SRC));
    
    // Allow argument mismatch (needed for legacy FORTRAN)
    cmd.arg("-fallow-argument-mismatch");
    
    // Output file
    cmd.arg("-o").arg(output);
    
    // Source files
    for src in source_files {
        cmd.arg(src);
    }
    
    let output_result: Output = cmd.output()?;
    
    if !output_result.status.success() {
        return Err(FortranError::CompileFailed {
            stderr: String::from_utf8_lossy(&output_result.stderr).to_string(),
            code: output_result.status.code(),
        });
    }
    
    Ok(())
}

/// Run a FORTRAN executable and capture its stdout.
///
/// # Arguments
/// * `executable` - Path to the executable
///
/// # Returns
/// The stdout output as a string (typically JSON)
pub fn run_fortran_test(executable: impl AsRef<Path>) -> Result<String, FortranError> {
    let output = Command::new(executable.as_ref()).output()?;
    
    if !output.status.success() {
        return Err(FortranError::RunFailed {
            stderr: String::from_utf8_lossy(&output.stderr).to_string(),
            code: output.status.code(),
        });
    }
    
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

/// Compile and run a FORTRAN test program in one step.
///
/// # Arguments
/// * `source_files` - List of FORTRAN source file paths
/// * `executable` - Path for the temporary executable
///
/// # Returns
/// The stdout output as a string
pub fn compile_and_run(source_files: &[&str], executable: &str) -> Result<String, FortranError> {
    compile_fortran(source_files, executable)?;
    run_fortran_test(executable)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfortran_exists() {
        let exists = Path::new(GFORTRAN).exists();
        if !exists {
            eprintln!("Warning: gfortran not found at {}", GFORTRAN);
        }
    }

    #[test]
    fn test_xfoil_src_exists() {
        let exists = Path::new(XFOIL_SRC).exists();
        assert!(exists, "XFOIL source directory not found at {}", XFOIL_SRC);
    }
}
