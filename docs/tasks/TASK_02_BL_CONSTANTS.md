# Task 02: Create rustfoil-bl Crate with Constants

## Objective
Create the boundary layer crate and implement BLPAR constants.

## Prerequisites
- Task 01 (rustfoil-testkit) completed

## Context
- FORTRAN constants in: `/Users/harry/flexfoil-boundary-layer/Xfoil/src/BLPAR.INC`

## Deliverables

### 1. Create crate structure

```
/Users/harry/flexfoil/crates/rustfoil-bl/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── constants.rs
    └── closures/
        └── mod.rs
```

### 2. Cargo.toml
```toml
[package]
name = "rustfoil-bl"
version = "0.1.0"
edition = "2021"

[dependencies]
rustfoil-core = { path = "../rustfoil-core" }

[dev-dependencies]
rustfoil-testkit = { path = "../rustfoil-testkit" }
approx = "0.5"
serde_json = "1"
```

### 3. src/constants.rs
Port BLPAR.INC exactly:
```rust
//! Boundary layer closure constants from XFOIL's BLPAR.INC
//! 
//! These constants are initialized in XFOIL's INIT subroutine (xfoil.f)
//! and remain fixed throughout execution.

/// Shear stress coefficient constant
pub const SCCON: f64 = 5.6;

/// G-beta locus  A  constant
pub const GACON: f64 = 6.70;

/// G-beta locus  B  constant  
pub const GBCON: f64 = 0.75;

/// G-beta locus  C  constant
pub const GCCON: f64 = 18.0;

/// Dissipation length constant
pub const DLCON: f64 = 0.9;

/// CtEQ constant (XFOIL default, can be modified)
pub const CTRCON: f64 = 1.8;

/// CtEQ exponent (XFOIL default, can be modified)
pub const CTRCEX: f64 = 3.3;

/// Ux constant
pub const DUXCON: f64 = 1.0;

/// Derived: Shear-lag coefficient = 0.5/(GACON² * GBCON)
pub const CTCON: f64 = 0.5 / (GACON * GACON * GBCON);

/// Cf scaling factor (normally 1.0)
pub const CFFAC: f64 = 1.0;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_ctcon_derivation() {
        // Verify CTCON matches XFOIL's calculation
        let expected = 0.5 / (6.70_f64.powi(2) * 0.75);
        assert!((CTCON - expected).abs() < 1e-15);
    }
}
```

### 4. src/lib.rs
```rust
//! Boundary layer equations and closures for RustFoil
//!
//! This crate implements XFOIL's integral boundary layer formulation.

pub mod constants;
pub mod closures;

pub use constants::*;
```

### 5. src/closures/mod.rs
```rust
//! Boundary layer closure relations
//!
//! Each closure function returns both the value and all partial derivatives,
//! which are required for the Newton-Raphson Jacobian.

// Will be populated by subsequent tasks:
// pub mod hkin;
// pub mod hs;
// pub mod cf;
// pub mod dissipation;
// pub mod density;
```

### 6. Update workspace Cargo.toml
Add to `/Users/harry/flexfoil/Cargo.toml`:
```toml
members = [
    # ... existing members ...
    "crates/rustfoil-testkit",
    "crates/rustfoil-bl",
]
```

## Verification
```bash
cd /Users/harry/flexfoil
cargo build -p rustfoil-bl
cargo test -p rustfoil-bl
```

## Next Task
After completion, proceed to TASK_03_HKIN.md
