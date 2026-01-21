//! Boundary layer closure constants from XFOIL's BLPAR.INC
//!
//! These constants are initialized in XFOIL's BLPINI subroutine (xbl.f)
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
