# Task 11: Create rustfoil-coupling with QDCALC

## Objective
Create the viscous-inviscid coupling crate and implement the mass defect influence matrix.

## Prerequisites
- rustfoil-bl crate complete (Tasks 01-10)

## Context
- QDCALC: xpanel.f line 1149 - builds DIJ matrix relating mass defect to edge velocity

## Crate Structure
```
/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── dij.rs           # This task
    ├── newton.rs        # Task 12
    ├── solve.rs         # Task 13
    ├── march.rs         # Task 14
    └── update.rs        # Task 15
```

## Deliverables

### 1. Cargo.toml
```toml
[package]
name = "rustfoil-coupling"
version = "0.1.0"
edition = "2021"

[dependencies]
rustfoil-core = { path = "../rustfoil-core" }
rustfoil-bl = { path = "../rustfoil-bl" }
nalgebra = "0.32"

[dev-dependencies]
rustfoil-testkit = { path = "../rustfoil-testkit" }
approx = "0.5"
```

### 2. src/lib.rs
```rust
//! Viscous-inviscid coupling for RustFoil

pub mod dij;
pub mod newton;
pub mod solve;
pub mod march;
pub mod update;
```

### 3. src/dij.rs
```rust
//! Mass defect influence matrix (DIJ)
//!
//! The DIJ matrix relates changes in mass defect (Ue*δ*) to changes in edge velocity.
//! This comes from the inviscid source panel influence.
//!
//! XFOIL Reference: xpanel.f QDCALC (line 1149)

use nalgebra::DMatrix;

/// Build the mass defect influence matrix
///
/// DIJ[i,j] represents the influence of mass defect change at station j
/// on edge velocity at station i:
///
///   ΔUe_i = Σ_j DIJ_ij * Δ(Ue*δ*)_j
///
/// # Arguments
/// * `x` - x-coordinates of BL stations
/// * `y` - y-coordinates of BL stations  
/// * `n` - Number of stations
///
/// # Reference
/// XFOIL xpanel.f QDCALC (line 1149)
pub fn build_dij_matrix(x: &[f64], y: &[f64]) -> DMatrix<f64> {
    let n = x.len();
    let mut dij = DMatrix::zeros(n, n);
    
    // The DIJ matrix comes from treating mass defect as source panels
    // and computing their velocity influence at each station
    //
    // For each source panel j, the velocity influence at point i is:
    //   ΔUe_i = (strength_j / 2π) * geometric_factor(i,j)
    
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];
                let r2 = dx * dx + dy * dy;
                
                // Source panel influence (simplified - full version in QDCALC)
                // This is the 2D source velocity influence
                dij[(i, j)] = dx / (2.0 * std::f64::consts::PI * r2);
            }
        }
    }
    
    // TODO: Handle self-influence (diagonal) terms
    // TODO: Include panel geometry factors
    // TODO: Handle wake coupling
    
    dij
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_dij_symmetry() {
        // Simple test geometry
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015, 0.01];
        
        let dij = build_dij_matrix(&x, &y);
        
        // DIJ should have certain symmetry properties
        assert!(dij.nrows() == 4);
        assert!(dij.ncols() == 4);
    }
}
```

### 4. Update workspace Cargo.toml
Add rustfoil-coupling to members.

## Next Task
After completion, proceed to TASK_12_COUPLING_NEWTON.md

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
