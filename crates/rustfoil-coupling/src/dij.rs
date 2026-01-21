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
///
/// # Returns
/// An n×n DMatrix where n is the number of stations
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
    fn test_dij_dimensions() {
        // Simple test geometry
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015, 0.01];
        
        let dij = build_dij_matrix(&x, &y);
        
        // DIJ should be n x n
        assert_eq!(dij.nrows(), 4);
        assert_eq!(dij.ncols(), 4);
    }
    
    #[test]
    fn test_dij_diagonal_zero() {
        // Self-influence terms should be zero (for now)
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015, 0.01];
        
        let dij = build_dij_matrix(&x, &y);
        
        // Diagonal should be zero (self-influence not yet implemented)
        for i in 0..4 {
            assert_eq!(dij[(i, i)], 0.0, "Diagonal element [{},{}] should be zero", i, i);
        }
    }
    
    #[test]
    fn test_dij_off_diagonal_nonzero() {
        // Off-diagonal elements should be non-zero for non-coincident points
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015, 0.01];
        
        let dij = build_dij_matrix(&x, &y);
        
        // Off-diagonal should be non-zero
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    assert!(dij[(i, j)] != 0.0, 
                        "Off-diagonal element [{},{}] should be non-zero", i, j);
                }
            }
        }
    }
    
    #[test]
    fn test_dij_source_influence_sign() {
        // For a point downstream (larger x), the influence from an upstream source
        // should contribute positive velocity (dx > 0)
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 0.0];
        
        let dij = build_dij_matrix(&x, &y);
        
        // Point 1 (x=1.0) influenced by source at point 0 (x=0.0)
        // dx = 1.0 - 0.0 = 1.0 > 0, so DIJ[1,0] > 0
        assert!(dij[(1, 0)] > 0.0, "Downstream point should have positive influence from upstream source");
        
        // Point 0 (x=0.0) influenced by source at point 1 (x=1.0)
        // dx = 0.0 - 1.0 = -1.0 < 0, so DIJ[0,1] < 0
        assert!(dij[(0, 1)] < 0.0, "Upstream point should have negative influence from downstream source");
    }
    
    #[test]
    fn test_dij_influence_decay() {
        // Influence should decay with distance (1/r² in the formula)
        let x = vec![0.0, 0.1, 0.5];
        let y = vec![0.0, 0.0, 0.0];
        
        let dij = build_dij_matrix(&x, &y);
        
        // Influence from station 0 on station 1 (close) vs station 2 (far)
        // |DIJ[1,0]| should be > |DIJ[2,0]| since station 1 is closer
        assert!(dij[(1, 0)].abs() > dij[(2, 0)].abs(), 
            "Closer stations should have stronger influence");
    }
    
    #[test]
    fn test_dij_empty_input() {
        // Empty input should produce empty matrix
        let x: Vec<f64> = vec![];
        let y: Vec<f64> = vec![];
        
        let dij = build_dij_matrix(&x, &y);
        
        assert_eq!(dij.nrows(), 0);
        assert_eq!(dij.ncols(), 0);
    }
    
    #[test]
    fn test_dij_single_station() {
        // Single station should produce 1x1 matrix with zero
        let x = vec![0.5];
        let y = vec![0.1];
        
        let dij = build_dij_matrix(&x, &y);
        
        assert_eq!(dij.nrows(), 1);
        assert_eq!(dij.ncols(), 1);
        assert_eq!(dij[(0, 0)], 0.0);
    }
}
