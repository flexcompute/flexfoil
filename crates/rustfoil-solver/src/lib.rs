//! RustFoil Solver - Aerodynamic analysis solvers.
//!
//! This crate provides the core aerodynamic solvers for RustFoil:
//!
//! # Solver Hierarchy
//!
//! ## Phase 2: Inviscid Solver (Linear Vorticity Panel Method)
//! - Fast linear system solve (O(N²) assembly, O(N³) solve for dense)
//! - Real-time capable for geometry manipulation
//! - Returns pressure coefficient (Cp) distribution
//! - Multi-body support built-in
//!
//! ## Phase 3: Boundary Layer Solver
//! - Integral boundary layer equations
//! - Thwaites (laminar) + Head/Green (turbulent)
//! - Transition prediction (eN method)
//!
//! ## Phase 4: Viscous-Inviscid Interaction
//! - Global Newton-Raphson coupling
//! - Transpiration velocity model
//! - Full polar generation
//!
//! # Mathematical Background
//!
//! ## Linear Vorticity Panel Method
//!
//! The inviscid solver uses a vortex panel method where:
//! - Each panel carries a linearly-varying vorticity distribution
//! - Boundary condition: no flow through surface (V·n = 0)
//! - Kutta condition: smooth flow departure at trailing edge
//!
//! The resulting linear system is:
//! ```text
//! [A]{γ} = {b}
//! ```
//! where:
//! - A_ij = influence of panel j's vorticity on panel i's normal velocity
//! - γ = vorticity strengths at panel nodes
//! - b = -V∞·n (freestream normal velocity at each panel)
//!
//! ## Kutta Condition
//!
//! For a sharp trailing edge, the vorticity must satisfy:
//! ```text
//! γ_upper(TE) + γ_lower(TE) = 0
//! ```
//! This ensures finite velocity at the trailing edge and determines
//! the circulation (and hence lift) around the airfoil.

#![warn(missing_docs)]
#![warn(clippy::all)]

pub mod inviscid;
pub mod boundary_layer;
pub mod viscous;
/// Mfoil-style global Newton boundary layer solver.
pub mod mfoil_bl;

// Re-exports for convenience
pub use boundary_layer::{
    BLSolution, BLConfig, BLSolver, TurbulentModel,
    XfoilConstants, BLState, BLStation, Surface,
};
pub use viscous::{
    ViscousSolver, ViscousSolution, ViscousConfig, CouplingMethod,
    NewtonConfig, NewtonGeometry, NewtonState, NewtonResult, NewtonVIISolver,
};
pub use mfoil_bl::{
    MfoilSolver, MfoilConfig, MfoilSolution,
};
