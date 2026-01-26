//! Viscous-inviscid coupling for RustFoil
//!
//! This crate implements XFOIL's viscous-inviscid interaction scheme,
//! coupling the boundary layer equations with the inviscid panel method.
//!
//! # Modules
//!
//! - [`dij`] - Mass defect influence matrix (QDCALC)
//! - [`newton`] - Newton system construction (BLSYS)
//! - [`global_newton`] - Global Newton system for full VI coupling (SETBL)
//! - [`solve`] - Block tridiagonal solver (BLSOLV)
//! - [`march`] - Boundary layer marching (MRCHUE)
//! - [`update`] - Solution update procedures (UPDATE)
//! - [`stmove`] - Stagnation point relocation (STMOVE)
//! - [`wake`] - Wake marching and TE combination

pub mod dij;
pub mod global_newton;
pub mod march;
pub mod newton;
pub mod solve;
pub mod stmove;
pub mod update;
pub mod wake;
