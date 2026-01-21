//! Test infrastructure for comparing Rust BL implementations against XFOIL FORTRAN.
//!
//! This crate provides utilities for:
//! - Compiling and running FORTRAN test programs
//! - Loading JSON reference data
//! - Float comparison with tolerance

pub mod approx;
pub mod fortran_runner;

use serde::de::DeserializeOwned;
use std::fs;
use std::path::Path;

/// Load reference test data from a JSON file.
///
/// # Panics
/// Panics if the file cannot be read or parsed.
pub fn load_reference<T: DeserializeOwned>(path: impl AsRef<Path>) -> T {
    let path = path.as_ref();
    let contents = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read reference file {}: {}", path.display(), e));
    serde_json::from_str(&contents)
        .unwrap_or_else(|e| panic!("Failed to parse JSON from {}: {}", path.display(), e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

    #[derive(Deserialize)]
    struct TestData {
        value: f64,
    }

    #[test]
    fn test_load_reference_works() {
        // This test would need actual test data to run
        // For now it just verifies the module compiles
    }
}
