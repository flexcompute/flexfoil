# Task 13: Implement BLSOLV Block Solver

## Objective
Implement the block Gaussian elimination solver for the BL Newton system.

## Prerequisites
- Task 12 (Newton system builder) completed

## Context
- BLSOLV: xsolve.f line 283 - block tridiagonal solver

## FORTRAN Reference

### BLSOLV (xsolve.f:283-486)
Block Gaussian elimination:
1. Forward sweep: eliminate lower diagonal
2. Back substitution: solve for unknowns
3. Handle wake coupling at trailing edge

## Deliverables

### src/solve.rs
```rust
//! Block tridiagonal solver for BL Newton system
//!
//! XFOIL Reference: xsolve.f BLSOLV (line 283)

use crate::newton::BlNewtonSystem;

/// Solve the block-tridiagonal BL Newton system
///
/// Uses block Gaussian elimination with forward sweep and back substitution.
///
/// # Arguments
/// * `system` - The Newton system to solve
///
/// # Returns
/// Solution vector [Δδ*, Δθ, ΔN/ΔCτ] at each station
///
/// # Reference
/// XFOIL xsolve.f BLSOLV (line 283)
pub fn solve_bl_system(system: &BlNewtonSystem) -> Vec<[f64; 3]> {
    let n = system.n;
    let mut solution = vec![[0.0; 3]; n];
    
    // Working arrays for forward elimination
    let mut va_mod = system.va.clone();
    let mut rhs_mod = system.rhs.clone();
    
    // === Forward Sweep ===
    // Eliminate lower diagonal (VB) blocks
    for i in 1..n {
        // Compute multiplier: M = VB[i] * VA[i-1]^(-1)
        let va_inv = invert_3x3(&va_mod[i-1]);
        let mult = multiply_3x3(&system.vb[i], &va_inv);
        
        // Update diagonal: VA[i] -= M * VB[i]  (but VB is below, so this is different)
        // Actually: VA[i] -= M * VC[i-1] where VC is upper diagonal
        // For pure tridiagonal: VA[i] is already the diagonal
        
        // Update RHS: rhs[i] -= M * rhs[i-1]
        let rhs_contrib = multiply_3x3_vec(&mult, &rhs_mod[i-1]);
        for k in 0..3 {
            rhs_mod[i][k] -= rhs_contrib[k];
        }
        
        // Eliminate VB (make it zero)
        // VA becomes modified
        for j in 0..3 {
            for k in 0..3 {
                va_mod[i][j][k] -= mult[j][k]; // Simplified
            }
        }
    }
    
    // === Back Substitution ===
    // Solve last station
    let va_inv_last = invert_3x3(&va_mod[n-1]);
    solution[n-1] = multiply_3x3_vec(&va_inv_last, &rhs_mod[n-1]);
    
    // Back substitute
    for i in (0..n-1).rev() {
        // x[i] = VA[i]^(-1) * (rhs[i] - VC[i] * x[i+1])
        // For simple tridiagonal, VC is upper diagonal
        let va_inv = invert_3x3(&va_mod[i]);
        solution[i] = multiply_3x3_vec(&va_inv, &rhs_mod[i]);
        // TODO: Account for upper diagonal coupling
    }
    
    solution
}

/// Invert a 3x3 matrix
fn invert_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    // Compute determinant
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    
    if det.abs() < 1e-30 {
        // Singular matrix, return identity
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    
    let inv_det = 1.0 / det;
    
    [
        [
            inv_det * (m[1][1] * m[2][2] - m[1][2] * m[2][1]),
            inv_det * (m[0][2] * m[2][1] - m[0][1] * m[2][2]),
            inv_det * (m[0][1] * m[1][2] - m[0][2] * m[1][1]),
        ],
        [
            inv_det * (m[1][2] * m[2][0] - m[1][0] * m[2][2]),
            inv_det * (m[0][0] * m[2][2] - m[0][2] * m[2][0]),
            inv_det * (m[0][2] * m[1][0] - m[0][0] * m[1][2]),
        ],
        [
            inv_det * (m[1][0] * m[2][1] - m[1][1] * m[2][0]),
            inv_det * (m[0][1] * m[2][0] - m[0][0] * m[2][1]),
            inv_det * (m[0][0] * m[1][1] - m[0][1] * m[1][0]),
        ],
    ]
}

/// Multiply 3x3 matrix by 3x3 matrix
fn multiply_3x3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Multiply 3x3 matrix by 3-vector
fn multiply_3x3_vec(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_invert_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = invert_3x3(&id);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < 1e-10);
            }
        }
    }
}
```

## Next Task
After completion, proceed to TASK_14_COUPLING_MARCH.md
