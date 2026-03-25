/**
 * CFD Solver - WebGPU-accelerated 2D compressible flow solver.
 *
 * Solves Euler / Navier-Stokes / RANS-SA equations on structured O-grids
 * using Roe approximate Riemann solver with MUSCL/WENO5 reconstruction.
 */

export { CfdSolver } from './CfdSolver';
export type { CfdSolverResult } from './CfdSolver';
export { DEFAULT_CFD_CONFIG, packCfdParams, NVAR } from './CfdConfig';
export type { CfdConfig, PhysicsMode, ReconstructionMode, TimeSteppingMode } from './CfdConfig';
