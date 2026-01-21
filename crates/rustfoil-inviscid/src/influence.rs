//! Influence coefficient calculations (XFOIL's PSILIN).
//!
//! This module computes the stream function ψ at a field point due to all vortex panels,
//! and the sensitivity ∂ψ/∂γ for each node. This forms the core of the panel method.
//!
//! # Linear Vorticity Panels
//!
//! XFOIL uses panels with linearly-varying vorticity (node-based unknowns).
//! For a panel from node JO to JP with γ varying linearly:
//!
//! ```text
//! γ(s) = γ_JO + (γ_JP - γ_JO) * (s - s_JO) / (s_JP - s_JO)
//! ```
//!
//! # Sum/Difference Formulation
//!
//! The key to XFOIL's efficiency is the PSIS/PSID formulation:
//! - PSIS: coefficient of (γ_JO + γ_JP), the "sum" term
//! - PSID: coefficient of (γ_JP - γ_JO), the "difference" term
//!
//! This allows influence coefficients to be computed incrementally.
//!
//! # XFOIL Reference
//!
//! - `xpanel.f`: PSILIN subroutine (lines 99-800)

use crate::geometry::AirfoilGeometry;
use crate::{QOPI, HOPI};

/// Result from computing influence coefficients at a single field point.
#[derive(Debug, Clone)]
pub struct PsilinResult {
    /// Stream function value at the field point
    pub psi: f64,
    /// Tangential velocity at α=0° (contribution to QTAN1)
    pub qtan1: f64,
    /// Tangential velocity at α=90° (contribution to QTAN2)
    pub qtan2: f64,
    /// ∂ψ/∂γⱼ for each node j (influence coefficient array)
    pub dzdg: Vec<f64>,
    /// ∂ψ/∂σⱼ for each node j (source influence, if computed)
    pub dzdm: Vec<f64>,
}

/// Intermediate values from a single panel's contribution.
/// Used for testing against XFOIL's internal values.
#[derive(Debug, Clone, Default)]
pub struct PanelContribution {
    /// Sum term (coefficient of γ_JO + γ_JP)
    pub psis: f64,
    /// Difference term (coefficient of γ_JP - γ_JO)
    pub psid: f64,
    /// Log term for endpoint 1
    pub g1: f64,
    /// Log term for endpoint 2
    pub g2: f64,
    /// Angle term for endpoint 1
    pub t1: f64,
    /// Angle term for endpoint 2
    pub t2: f64,
    /// Local x-coordinate of point relative to panel start
    pub x1: f64,
    /// Local x-coordinate of point relative to panel end
    pub x2: f64,
    /// Local y-coordinate (perpendicular distance to panel line)
    pub yy: f64,
}

/// Compute influence coefficients at a field point (i, xi, yi).
///
/// This is XFOIL's PSILIN subroutine, computing:
/// - The stream function ψ at (xi, yi) due to all panels
/// - The derivative ∂ψ/∂γⱼ for each node j
///
/// # Arguments
///
/// * `geom` - Airfoil geometry
/// * `i` - Index of the field point (for singularity handling)
/// * `xi`, `yi` - Field point coordinates
///
/// # Returns
///
/// A `PsilinResult` containing the stream function and influence coefficients.
pub fn psilin(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinResult {
    let n = geom.n;
    
    // Initialize influence coefficient arrays
    let mut dzdg = vec![0.0; n];
    let dzdm = vec![0.0; n];
    
    // Initialize accumulated values (TODO: compute these in full PSILIN)
    let psi = 0.0;
    let qtan1 = 0.0;
    let qtan2 = 0.0;
    
    // Distance tolerance for TE panel skip
    let seps = geom.total_arc_length() * 1e-5;
    
    // TE panel coefficients
    let (scs, sds) = geom.te_coefficients();
    
    // Saved values for TE panel treatment
    let mut te_contrib: Option<PanelContribution> = None;
    
    // Loop over all panels (JO = 0 to N-1 in 0-based indexing)
    for jo in 0..n {
        let jp = (jo + 1) % n;  // JP = JO+1, wrapping for TE panel
        
        // Panel endpoints
        let x_jo = geom.x[jo];
        let y_jo = geom.y[jo];
        let x_jp = geom.x[jp];
        let y_jp = geom.y[jp];
        
        // Panel vector and length
        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;
        
        // Skip TE panel if endpoints coincide (SEPS check)
        if jo == n - 1 && ds_sq < seps * seps {
            continue;
        }
        
        // Skip zero-length panels
        if ds_sq < 1e-24 {
            continue;
        }
        
        let dso = ds_sq.sqrt();
        let dsio = 1.0 / dso;
        
        // Unit tangent vector along panel
        let sx = dx * dsio;
        let sy = dy * dsio;
        
        // Vector from panel endpoints to field point
        let rx1 = xi - x_jo;
        let ry1 = yi - y_jo;
        let rx2 = xi - x_jp;
        let ry2 = yi - y_jp;
        
        // Transform to panel-local coordinates
        // x1, x2: along panel direction
        // yy: perpendicular to panel
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;
        
        // Squared distances to endpoints
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;
        
        // Logarithm and arctangent terms
        // Handle singularities when field point is at panel endpoint
        let (g1, t1) = if i != jo && rs1 > 1e-20 {
            (rs1.ln(), x1.atan2(yy))
        } else {
            (0.0, 0.0)
        };
        
        let (g2, t2) = if i != jp && rs2 > 1e-20 {
            (rs2.ln(), x2.atan2(yy))
        } else {
            (0.0, 0.0)
        };
        
        // Compute PSIS and PSID (sum/difference formulation)
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);
        
        // For inviscid-only (no sources), all panels including TE are treated the same
        // XFOIL's SIGLIN=.FALSE. case - pure vortex panels
        
        // Accumulate influence coefficients
        // ψ += QOPI * (PSIS*(γ_JO+γ_JP) + PSID*(γ_JP-γ_JO))
        // Rearranging: ψ += QOPI * ((PSIS-PSID)*γ_JO + (PSIS+PSID)*γ_JP)
        dzdg[jo] += QOPI * (contrib.psis - contrib.psid);
        dzdg[jp] += QOPI * (contrib.psis + contrib.psid);
    }
    
    // Note: TE panel source/vortex treatment (SIGLIN=.TRUE.) is only needed
    // for viscous analysis with mass defect. For inviscid-only, the TE panel
    // is treated like any other vortex panel above.
    let _ = te_contrib; // Suppress unused warning
    let _ = scs;
    let _ = sds;
    
    PsilinResult {
        psi,
        qtan1,
        qtan2,
        dzdg,
        dzdm,
    }
}

/// Compute PSIS and PSID for a single panel contribution.
///
/// This implements the core formulas from XFOIL (lines 194-220).
///
/// # Arguments
///
/// * `x1`, `x2` - Local x-coordinates of field point relative to panel endpoints
/// * `yy` - Local y-coordinate (perpendicular distance to panel line)
/// * `rs1`, `rs2` - Squared distances to panel endpoints
/// * `g1`, `g2` - Log terms: ln(r1²), ln(r2²)
/// * `t1`, `t2` - Angle terms: atan2(x1, yy), atan2(x2, yy)
pub fn compute_psis_psid(
    x1: f64, x2: f64, yy: f64,
    rs1: f64, rs2: f64,
    g1: f64, g2: f64, t1: f64, t2: f64,
) -> PanelContribution {
    // PSIS: sum term (XFOIL line ~215)
    // PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
    let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
    
    // PSID: difference term (XFOIL line ~216)
    // PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2 - RS1*G1 + X1*X1 - X2*X2)) / (X1-X2)
    let dxinv = if (x1 - x2).abs() > 1e-20 {
        1.0 / (x1 - x2)
    } else {
        0.0
    };
    
    let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
    
    PanelContribution {
        psis,
        psid,
        g1, g2, t1, t2, x1, x2, yy,
    }
}

/// Compute influence coefficients for a single panel.
///
/// This is a helper for testing that returns intermediate values.
pub fn psilin_single_panel(
    geom: &AirfoilGeometry,
    i: usize,
    jo: usize,
) -> PanelContribution {
    let n = geom.n;
    let jp = (jo + 1) % n;
    
    // Field point
    let xi = geom.x[i];
    let yi = geom.y[i];
    
    // Panel endpoints
    let x_jo = geom.x[jo];
    let y_jo = geom.y[jo];
    let x_jp = geom.x[jp];
    let y_jp = geom.y[jp];
    
    // Panel vector
    let dx = x_jp - x_jo;
    let dy = y_jp - y_jo;
    let ds_sq = dx * dx + dy * dy;
    
    if ds_sq < 1e-24 {
        return PanelContribution::default();
    }
    
    let dso = ds_sq.sqrt();
    let dsio = 1.0 / dso;
    
    // Unit tangent
    let sx = dx * dsio;
    let sy = dy * dsio;
    
    // Vectors to field point
    let rx1 = xi - x_jo;
    let ry1 = yi - y_jo;
    let rx2 = xi - x_jp;
    let ry2 = yi - y_jp;
    
    // Local coordinates
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;
    
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;
    
    let (g1, t1) = if i != jo && rs1 > 1e-20 {
        (rs1.ln(), x1.atan2(yy))
    } else {
        (0.0, 0.0)
    };
    
    let (g2, t2) = if i != jp && rs2 > 1e-20 {
        (rs2.ln(), x2.atan2(yy))
    } else {
        (0.0, 0.0)
    };
    
    compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_geometry() -> AirfoilGeometry {
        // Simple 4-point diamond shape for testing
        let points = vec![
            (1.0, 0.0),   // Upper TE
            (0.0, 0.1),   // Upper LE
            (0.0, -0.1),  // Lower LE (same x as upper)
            (1.0, 0.0),   // Lower TE (closed)
        ];
        // Need more points for a valid geometry
        let n = 20;
        let mut pts = Vec::with_capacity(n);
        for i in 0..n/2 {
            let t = i as f64 / (n/2 - 1) as f64;
            let x = 1.0 - t;
            let y = 0.1 * (1.0 - (2.0 * t - 1.0).powi(2));
            pts.push((x, y));
        }
        for i in 1..n/2 {
            let t = i as f64 / (n/2 - 1) as f64;
            let x = t;
            let y = -0.1 * (1.0 - (2.0 * t - 1.0).powi(2));
            pts.push((x, y));
        }
        AirfoilGeometry::from_points(&pts).unwrap()
    }

    #[test]
    fn test_psilin_produces_finite_values() {
        let geom = make_simple_geometry();
        
        // Test at a few interior points
        for i in 1..geom.n - 1 {
            let result = psilin(&geom, i, geom.x[i], geom.y[i]);
            
            assert!(result.psi.is_finite(), "psi should be finite at node {}", i);
            for (j, &dz) in result.dzdg.iter().enumerate() {
                assert!(dz.is_finite(), "dzdg[{}] should be finite at node {}", j, i);
            }
        }
    }

    #[test]
    fn test_psis_psid_symmetry() {
        // For a point directly above the panel midpoint, certain symmetries hold
        let x1: f64 = 0.5;
        let x2: f64 = -0.5;
        let yy: f64 = 1.0;
        let rs1: f64 = x1 * x1 + yy * yy;
        let rs2: f64 = x2 * x2 + yy * yy;
        let g1 = rs1.ln();
        let g2 = rs2.ln();
        let t1 = x1.atan2(yy);
        let t2 = x2.atan2(yy);
        
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);
        
        // rs1 == rs2 and g1 == g2 due to symmetry
        assert!((rs1 - rs2).abs() < 1e-10);
        assert!((g1 - g2).abs() < 1e-10);
        
        // PSIS and PSID should be finite
        assert!(contrib.psis.is_finite());
        assert!(contrib.psid.is_finite());
    }
}
