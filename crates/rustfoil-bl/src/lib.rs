//! Boundary layer equations and closures for RustFoil
//!
//! This crate implements XFOIL's integral boundary layer formulation.
//!
//! # Modules
//!
//! - [`closures`] - Closure relations (HKIN, HS, CF, DI, HCT, DAMPL)
//! - [`constants`] - BLPAR constants from XFOIL
//! - [`state`] - BlStation struct for BL state at a single point
//! - [`equations`] - BLVAR/BLDIF for computing secondary variables and residuals

pub mod closures;
pub mod constants;
pub mod equations;
pub mod state;

pub use constants::*;
pub use equations::{bldif, blvar, BlJacobian, BlResiduals, FlowType};
pub use state::{BlDerivatives, BlStation};
