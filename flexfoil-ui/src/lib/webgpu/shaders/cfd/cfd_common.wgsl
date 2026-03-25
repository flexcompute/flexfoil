// Common struct definitions for CFD compute shaders.
// This file is prepended to each shader at pipeline creation time.

struct CfdParams {
    ni: u32,
    nj: u32,
    gamma: f32,
    cfl: f32,
    mach_inf: f32,
    alpha: f32,
    reynolds: f32,
    prandtl: f32,
    dt: f32,
    iteration: u32,
    physics_mode: u32,    // 0=Euler, 1=LaminarNS, 2=RANS_SA
    reconstruction: u32,  // 0=MUSCL, 1=WENO5
};

// Boundary condition type constants
const BC_INTERIOR: u32 = 0u;
const BC_WALL: u32 = 1u;
const BC_FARFIELD: u32 = 2u;
const BC_WAKE: u32 = 3u;

// Number of conservative variables per cell
const NVAR: u32 = 5u;

// Helper: flat index for cell (i, j) in a ni x nj grid
fn cell_idx(i: u32, j: u32, ni: u32) -> u32 {
    return j * ni + i;
}

// Helper: flat index for variable `v` at cell (i, j)
fn q_idx(i: u32, j: u32, v: u32, ni: u32) -> u32 {
    return (j * ni + i) * NVAR + v;
}

// Helper: wrap i-index for periodic O-grid
fn wrap_i(i: i32, ni: u32) -> u32 {
    let n = i32(ni);
    return u32(((i % n) + n) % n);
}
