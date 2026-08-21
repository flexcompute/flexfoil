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
//! │  │             │                    │  • Newton update               │ │
//! │  │  CL, CM     │                    │  • Ue coupling                 │ │
//! │  │             │                    │  • forces (CD)                 │ │
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
//! - [`inviscid`] - Inviscid panel method (Linear Vorticity), plus the
//!   flowfield/streamline/smoke visualization helpers built on it
//! - [`viscous`] - Viscous boundary layer solver (VISCAL)
//!
//! # Which inviscid solver is XFOIL-faithful?
//!
//! Not this crate's. The XFOIL-faithful inviscid solver is the
//! **`rustfoil-inviscid`** crate, and that is what [`viscous`] itself calls.
//! [`inviscid`] is an older, self-contained implementation of the same panel
//! method. It is still live on two paths: the visualization code
//! ([`inviscid::velocity`], [`inviscid::smoke`]) is built on it, and the
//! inviscid-only entry points of the WASM and Python bindings report its `cl`,
//! `cm` and `cp`. It is single-body only, and new inviscid solver work —
//! multi-element included — belongs in `rustfoil-inviscid`. See the
//! [`inviscid`] module docs for the per-entry-point breakdown.
//!
//! # Example
//!
//! Note: this sketch drives the inviscid solver in [`inviscid`]. For the
//! XFOIL-faithful path, factorize with `rustfoil_inviscid::InviscidSolver`
//! instead (see [`viscous::setup_from_body`]).
//!
//! ```ignore
//! use rustfoil_core::Body;
//! use rustfoil_solver::{
//!     inviscid::{InviscidSolver, FlowConditions},
//!     viscous::{ViscousSolverConfig, solve_viscous},
//! };
//!
//! // 1. Create airfoil and solve inviscid
//! let body = Body::from_naca("0012", 160)?;
//! let solver = InviscidSolver::new();
//! let factorized = solver.factorize(&[body.clone()])?;
//! let inv_sol = factorized.solve_alpha(&FlowConditions::with_alpha_deg(4.0));
//!
//! // 2. Setup viscous from inviscid
//! let setup = ViscousSetup::from_inviscid(&body, &inv_sol);
//! let config = ViscousSolverConfig::with_reynolds(1e6);
//!
//! // 3. Initialize BL and solve  
//! let mut stations = setup.initialize_bl_stations();
//! let result = solve_viscous(&mut stations, &setup.ue_inviscid, &setup.dij, &config)?;
//!
//! println!("CL = {:.4}, CD = {:.5}", result.cl, result.cd);
//! ```

// Not the XFOIL-faithful inviscid solver, and single-body only. Live for the
// flowfield/streamline/smoke visualizations (`inviscid::velocity`,
// `inviscid::smoke`) built on it and for the bindings' inviscid-only entry
// points. The XFOIL-faithful solver is the `rustfoil-inviscid` crate — put new
// inviscid work (including multi-element) there.
pub mod inviscid;
pub mod viscous;

// Re-export main types for convenience
pub use inviscid::{SolverError, SolverResult};

// Re-export inviscid types (see the note on `mod inviscid`)
pub use inviscid::{FlowConditions, InviscidSolution, InviscidSolver, FactorizedSolution};

// Re-export viscous types
pub use viscous::{
    compute_forces, solve_viscous, solve_viscous_polar_parallel, AeroForces, ViscousResult,
    ViscousSolverConfig, ViscousSetup,
    // New integration with rustfoil-inviscid
    setup_from_body, setup_from_coords, SetupError, ViscousSetupResult,
};
