//! Boundary condition type mapping for the structured O-grid.
//!
//! In an O-grid topology:
//! - j=0: airfoil surface (wall)
//! - j=nj-1: far-field boundary
//! - i-direction: periodic (wraps around the airfoil)
//! - Wake cut: handled via periodicity in i

/// Boundary condition types (matches WGSL constants).
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcType {
    /// Interior cell (no special treatment)
    Interior = 0,
    /// Solid wall (slip for Euler, no-slip for NS/RANS)
    Wall = 1,
    /// Far-field (characteristic-based)
    FarField = 2,
    /// Wake cut (periodic connection)
    WakeCut = 3,
}

/// Generate boundary condition type array for the entire grid.
///
/// Returns a flat u32 array of length ni*nj.
pub fn generate_bc_types(ni: u32, nj: u32) -> Vec<u32> {
    let ni = ni as usize;
    let nj = nj as usize;
    let mut bc = vec![BcType::Interior as u32; ni * nj];

    // j=0: wall boundary
    for i in 0..ni {
        bc[i] = BcType::Wall as u32;
    }

    // j=nj-1: far-field boundary
    for i in 0..ni {
        bc[(nj - 1) * ni + i] = BcType::FarField as u32;
    }

    bc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bc_types() {
        let bc = generate_bc_types(8, 4);
        assert_eq!(bc.len(), 32);

        // Wall at j=0
        for i in 0..8 {
            assert_eq!(bc[i], BcType::Wall as u32);
        }
        // Interior at j=1,2
        for j in 1..3 {
            for i in 0..8 {
                assert_eq!(bc[j * 8 + i], BcType::Interior as u32);
            }
        }
        // FarField at j=3
        for i in 0..8 {
            assert_eq!(bc[3 * 8 + i], BcType::FarField as u32);
        }
    }
}
