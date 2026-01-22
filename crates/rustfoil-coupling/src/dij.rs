//! Mass defect influence matrix (DIJ)
//!
//! The DIJ matrix relates changes in mass defect (Ue*δ*) to changes in edge velocity.
//! This comes from the inviscid source panel influence.
//!
//! XFOIL Reference: xpanel.f QDCALC (line 1149)

use nalgebra::DMatrix;

/// Compute panel lengths for a set of BL stations
///
/// Returns the arc length of each panel segment.
fn compute_panel_lengths(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return vec![1.0; n.max(1)];
    }

    let mut ds = Vec::with_capacity(n);

    // First panel: use forward difference
    let dx = x[1] - x[0];
    let dy = y[1] - y[0];
    ds.push((dx * dx + dy * dy).sqrt().max(1e-10));

    // Interior panels: use central difference
    for i in 1..n - 1 {
        let dx = x[i + 1] - x[i - 1];
        let dy = y[i + 1] - y[i - 1];
        ds.push((dx * dx + dy * dy).sqrt().max(1e-10) / 2.0);
    }

    // Last panel: use backward difference
    let dx = x[n - 1] - x[n - 2];
    let dy = y[n - 1] - y[n - 2];
    ds.push((dx * dx + dy * dy).sqrt().max(1e-10));

    ds
}

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
    if n == 0 {
        return DMatrix::zeros(0, 0);
    }

    let mut dij = DMatrix::zeros(n, n);

    // Compute panel lengths for self-influence terms
    let ds = compute_panel_lengths(x, y);

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

                // Source panel influence with panel length factor
                // The influence is scaled by the panel length ds[j]
                let panel_factor = ds[j];
                dij[(i, j)] = dx * panel_factor / (2.0 * std::f64::consts::PI * r2);
            }
        }
    }

    // Self-influence (diagonal) terms
    //
    // For a source distribution, the self-induced velocity at a panel depends on:
    // 1. The local panel curvature (affects the velocity jump across the panel)
    // 2. The panel length
    //
    // In viscous-inviscid coupling, the diagonal term represents how a change
    // in mass defect at station i affects the edge velocity at that same station.
    //
    // From XFOIL's formulation, this is approximately:
    //   dij[i,i] ≈ -0.5 / ds[i]
    //
    // The negative sign indicates that increasing mass defect (which is like
    // adding sources) tends to slow down the local edge velocity.
    //
    // This is a simplification of XFOIL's QDCALC which uses panel geometry.
    for i in 0..n {
        // Self-influence coefficient
        // The factor accounts for the local source distribution effect
        // and is inversely proportional to panel spacing
        dij[(i, i)] = -0.5 / ds[i];
    }

    dij
}

/// Build the mass defect influence matrix with debug output
///
/// Same as [`build_dij_matrix`] but emits QDCALC debug events when
/// debug collection is active.
pub fn build_dij_matrix_debug(x: &[f64], y: &[f64], n_wake: usize) -> DMatrix<f64> {
    let dij = build_dij_matrix(x, y);

    // Emit debug event if collection is active
    if rustfoil_bl::is_debug_active() {
        let n = x.len();
        let n_airfoil = n.saturating_sub(n_wake);

        // Sample diagonal elements
        let diag_sample: Vec<f64> = (0..n.min(10)).map(|i| dij[(i, i)]).collect();

        // Sample first row
        let row1_sample: Vec<f64> = (0..n.min(10)).map(|j| dij[(0, j)]).collect();

        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::qdcalc(
            n_airfoil,
            n_wake,
            n,
            diag_sample,
            row1_sample,
        ));
    }

    dij
}

/// Build DIJ matrix with explicit panel geometry
///
/// This version takes pre-computed panel normals for more accurate
/// influence computation.
///
/// # Arguments
/// * `x` - x-coordinates of BL stations
/// * `y` - y-coordinates of BL stations
/// * `nx` - x-components of panel normals
/// * `ny` - y-components of panel normals
///
/// # Returns
/// An n×n DMatrix
pub fn build_dij_matrix_with_normals(
    x: &[f64],
    y: &[f64],
    nx: &[f64],
    ny: &[f64],
) -> DMatrix<f64> {
    let n = x.len();
    if n == 0 {
        return DMatrix::zeros(0, 0);
    }

    let mut dij = DMatrix::zeros(n, n);
    let ds = compute_panel_lengths(x, y);

    for i in 0..n {
        // Local tangent direction (perpendicular to normal)
        let tx = -ny[i];
        let ty = nx[i];

        for j in 0..n {
            if i != j {
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];
                let r2 = dx * dx + dy * dy;

                // Velocity induced by source at j in tangential direction at i
                // v_tangent = (σ / 2π) * (r · t) / r²
                let dot_rt = dx * tx + dy * ty;
                let panel_factor = ds[j];
                dij[(i, j)] = dot_rt * panel_factor / (2.0 * std::f64::consts::PI * r2);
            }
        }

        // Self-influence term
        dij[(i, i)] = -0.5 / ds[i];
    }

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
    fn test_dij_diagonal_nonzero() {
        // Self-influence terms should now be non-zero (negative)
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015, 0.01];

        let dij = build_dij_matrix(&x, &y);

        // Diagonal should be negative (self-influence slows down local velocity)
        for i in 0..4 {
            assert!(
                dij[(i, i)] < 0.0,
                "Diagonal element [{},{}] should be negative, got {}",
                i,
                i,
                dij[(i, i)]
            );
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
                    assert!(
                        dij[(i, j)] != 0.0,
                        "Off-diagonal element [{},{}] should be non-zero",
                        i,
                        j
                    );
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
        assert!(
            dij[(1, 0)] > 0.0,
            "Downstream point should have positive influence from upstream source"
        );

        // Point 0 (x=0.0) influenced by source at point 1 (x=1.0)
        // dx = 0.0 - 1.0 = -1.0 < 0, so DIJ[0,1] < 0
        assert!(
            dij[(0, 1)] < 0.0,
            "Upstream point should have negative influence from downstream source"
        );
    }

    #[test]
    fn test_dij_influence_decay() {
        // Influence should decay with distance (1/r² in the formula)
        let x = vec![0.0, 0.1, 0.5];
        let y = vec![0.0, 0.0, 0.0];

        let dij = build_dij_matrix(&x, &y);

        // Influence from station 0 on station 1 (close) vs station 2 (far)
        // |DIJ[1,0]| should be > |DIJ[2,0]| since station 1 is closer
        assert!(
            dij[(1, 0)].abs() > dij[(2, 0)].abs(),
            "Closer stations should have stronger influence"
        );
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
        // Single station should produce 1x1 matrix with non-zero diagonal
        let x = vec![0.5];
        let y = vec![0.1];

        let dij = build_dij_matrix(&x, &y);

        assert_eq!(dij.nrows(), 1);
        assert_eq!(dij.ncols(), 1);
        // Single station has self-influence (negative)
        assert!(dij[(0, 0)] < 0.0, "Single station should have negative self-influence");
    }

    #[test]
    fn test_panel_lengths() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 0.0, 0.0, 0.0];

        let ds = compute_panel_lengths(&x, &y);

        // All panels should have length 1.0 for uniform spacing
        assert_eq!(ds.len(), 4);
        for (i, &len) in ds.iter().enumerate() {
            assert!(
                (len - 1.0).abs() < 1e-10,
                "Panel {} length should be 1.0, got {}",
                i,
                len
            );
        }
    }

    #[test]
    fn test_dij_with_normals() {
        let x = vec![0.0, 0.1, 0.2];
        let y = vec![0.0, 0.01, 0.0];
        // Normals pointing upward (perpendicular to surface)
        let nx = vec![0.0, 0.0, 0.0];
        let ny = vec![1.0, 1.0, 1.0];

        let dij = build_dij_matrix_with_normals(&x, &y, &nx, &ny);

        assert_eq!(dij.nrows(), 3);
        assert_eq!(dij.ncols(), 3);

        // Diagonal should be negative
        for i in 0..3 {
            assert!(dij[(i, i)] < 0.0, "Diagonal should be negative");
        }
    }
}
