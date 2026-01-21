//! Viscous-inviscid coupling for RustFoil
//!
//! This crate implements XFOIL's viscous-inviscid interaction scheme,
//! coupling the boundary layer equations with the inviscid panel method.
//!
//! # Modules
//!
//! - [`dij`] - Mass defect influence matrix (QDCALC)
//! - [`newton`] - Newton system construction (BLSYS)
//! - [`solve`] - Block tridiagonal solver (BLSOLV)
//! - [`march`] - Boundary layer marching (MRCHUE)
//! - [`update`] - Solution update procedures (UPDATE)

pub mod dij;
pub mod newton;
pub mod solve;
pub mod march;
pub mod update;
