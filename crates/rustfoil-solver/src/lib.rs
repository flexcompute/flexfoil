//! RustFoil Solver - Inviscid and Viscous Flow Analysis
//!
//! This crate provides the complete flow analysis capability for RustFoil,
//! combining the inviscid panel method with viscous boundary layer analysis.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────────────┐
//! │                          rustfoil-solver                                │
//! ├─────────────────────────────────────────────────────────────────────────┤
//! │                                                                         │
//! │  ┌─────────────┐                    ┌─────────────────────────────────┐ │
//! │  │  inviscid/  │ ──── gamma ────►  │           viscous/              │ │
//! │  │             │                    │                                 │ │
//! │  │  Panel      │                    │  VISCAL iteration:             │ │
//! │  │  Method     │ ◄── mass defect ── │  • march BL                    │ │
//! │  │  (Stub)     │                    │  • Newton update               │ │
//! │  │             │                    │  • Ue coupling                 │ │
//! │  │  CL, CM     │                    │  • forces (CD)                 │ │
//! │  └─────────────┘                    └─────────────────────────────────┘ │
//! │         │                                        │                      │
//! │         └────────────────┬───────────────────────┘                      │
//! │                          │                                              │
//! │                          ▼                                              │
//! │                   ViscousResult                                         │
//! │                   (CL, CD, CM, transition, separation)                  │
//! └─────────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Modules
//!
//! - [`error`] - Error types for solver operations
//! - [`inviscid`] - Inviscid panel method (stub - to be replaced after merge)
//! - [`viscous`] - Viscous boundary layer solver (VISCAL)
//!
//! # MERGE NOTE
//!
//! This crate is designed to be merged with the `flexfoil` repository's
//! `rustfoil-solver`. After merging:
//!
//! 1. Replace the `inviscid/` stub with flexfoil's real implementation
//! 2. Merge the `error.rs` variants
//! 3. Keep the `viscous/` module as-is
//!
//! See the plan file for detailed merge instructions.
//!
//! # Example (Post-Merge)
//!
//! ```ignore
//! use rustfoil_core::Body;
//! use rustfoil_solver::{
//!     inviscid::{InviscidSolver, FlowConditions},
//!     viscous::{ViscousSolverConfig, solve_viscous, setup_from_inviscid},
//! };
//!
//! // 1. Create airfoil and solve inviscid
//! let body = Body::from_naca("0012", 160)?;
//! let solver = InviscidSolver::new();
//! let factorized = solver.factorize(&[body.clone()])?;
//! let inv_sol = factorized.solve_alpha(&FlowConditions::with_alpha_deg(4.0));
//!
//! // 2. Setup viscous from inviscid
//! let setup = setup_from_inviscid(&body, &inv_sol);
//! let config = ViscousSolverConfig::with_reynolds(1e6);
//!
//! // 3. Initialize BL and solve
//! let mut stations = initialize_bl_stations(...);
//! let result = solve_viscous(&mut stations, &setup.ue_inviscid, &setup.dij, &config)?;
//!
//! println!("CL = {:.4}, CD = {:.5}", result.cl, result.cd);
//! ```

pub mod error;
pub mod inviscid;
pub mod viscous;

// Re-export main types for convenience
pub use error::{SolverError, SolverResult};

// Re-export inviscid types (stub for now)
pub use inviscid::{FlowConditions, InviscidSolution, InviscidSolver};

// Re-export viscous types
pub use viscous::{
    compute_forces, solve_viscous, solve_viscous_polar_parallel, AeroForces, ViscousResult,
    ViscousSolverConfig, ViscousSetup,
};
