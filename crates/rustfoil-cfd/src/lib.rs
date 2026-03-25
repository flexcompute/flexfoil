//! RustFoil CFD - Structured mesh generation and setup for GPU-accelerated CFD.
//!
//! This crate provides CPU-side support for the WebGPU-based 2D CFD solver:
//! - O-type structured mesh generation around airfoils
//! - Solver configuration and parameter management
//! - Initial condition computation (freestream state)
//! - Boundary condition type mapping
//!
//! All heavy computation (flux evaluation, time stepping, turbulence model)
//! runs on WebGPU compute shaders. This crate prepares the data for GPU upload.

pub mod boundary;
pub mod config;
pub mod init;
pub mod mesh;
