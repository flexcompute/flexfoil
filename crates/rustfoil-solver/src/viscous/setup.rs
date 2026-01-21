//! Setup functions for viscous solver - maps inviscid → viscous data.
//!
//! This module provides the "glue code" that connects the inviscid panel method
//! solution to the boundary layer solver. It handles:
//!
//! - Computing arc lengths along the airfoil surface
//! - Finding the stagnation point
//! - Initializing BL stations from inviscid edge velocities
//! - Building the mass defect influence matrix (DIJ)
//!
//! # MERGE NOTE
//! After merging with flexfoil, the `setup_from_inviscid()` function
//! should be uncommented and will accept `Body` and `InviscidSolution`
//! types from the real inviscid solver.

use nalgebra::DMatrix;
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::dij::build_dij_matrix;

/// Data extracted from inviscid solution for viscous solver.
///
/// This struct holds all the geometric and flow data needed to
/// initialize and run the viscous solver.
#[derive(Debug, Clone)]
pub struct ViscousSetup {
    /// Arc length coordinates along surface.
    /// Starts at stagnation, goes around the airfoil.
    pub arc_lengths: Vec<f64>,

    /// Edge velocities from inviscid solution (= gamma at each node).
    /// This is the boundary layer edge condition Ue.
    pub ue_inviscid: Vec<f64>,

    /// Mass defect influence matrix.
    /// DIJ[i,j] = dUe_i / d(mass_defect)_j
    pub dij: DMatrix<f64>,

    /// Panel x-coordinates (for force integration).
    pub panel_x: Vec<f64>,

    /// Panel y-coordinates.
    pub panel_y: Vec<f64>,

    /// Index of stagnation point in the arrays.
    pub stagnation_index: usize,

    /// Number of panels on upper surface (stag to TE).
    pub n_upper: usize,

    /// Number of panels on lower surface (stag to TE).
    pub n_lower: usize,
}

impl ViscousSetup {
    /// Create a ViscousSetup from raw data.
    ///
    /// This is used when you have already computed arc lengths and
    /// velocities from some source (e.g., external inviscid solver).
    pub fn from_raw(
        arc_lengths: Vec<f64>,
        ue_inviscid: Vec<f64>,
        panel_x: Vec<f64>,
        panel_y: Vec<f64>,
    ) -> Self {
        let n = arc_lengths.len();
        let stagnation_index = find_stagnation(&ue_inviscid);

        // Build DIJ matrix from panel coordinates
        let dij = build_dij_matrix(&panel_x, &panel_y);

        // Estimate upper/lower panel counts (assuming stag is near middle)
        let n_upper = n - stagnation_index;
        let n_lower = stagnation_index;

        Self {
            arc_lengths,
            ue_inviscid,
            dij,
            panel_x,
            panel_y,
            stagnation_index,
            n_upper,
            n_lower,
        }
    }

    /// Get upper surface arc lengths and edge velocities.
    ///
    /// Returns (arc_lengths, ue) for stations from stagnation to TE on upper surface.
    pub fn upper_surface(&self) -> (&[f64], &[f64]) {
        let i0 = self.stagnation_index;
        (&self.arc_lengths[i0..], &self.ue_inviscid[i0..])
    }

    /// Get lower surface arc lengths and edge velocities.
    ///
    /// Returns (arc_lengths, ue) for stations from stagnation to TE on lower surface.
    /// Note: Lower surface runs from stagnation backward in the array.
    pub fn lower_surface(&self) -> (&[f64], &[f64]) {
        let i0 = self.stagnation_index;
        (&self.arc_lengths[..=i0], &self.ue_inviscid[..=i0])
    }
}

/// Compute cumulative arc length from panel coordinates.
///
/// Returns arc length at each panel endpoint, starting from 0.
///
/// # Arguments
/// * `x` - x-coordinates of panel endpoints
/// * `y` - y-coordinates of panel endpoints
///
/// # Returns
/// Vector of arc lengths, same length as input coordinates.
pub fn compute_arc_lengths(x: &[f64], y: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), y.len(), "x and y must have same length");
    let n = x.len();

    if n == 0 {
        return vec![];
    }

    let mut s = vec![0.0; n];

    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }

    s
}

/// Compute arc lengths from panel lengths.
///
/// Takes panel lengths (ds) and returns cumulative arc length.
///
/// # Arguments
/// * `ds` - Array of panel lengths
///
/// # Returns
/// Cumulative arc length at each panel center (length = ds.len())
pub fn cumulative_arc_length(ds: &[f64]) -> Vec<f64> {
    let n = ds.len();
    let mut s = Vec::with_capacity(n);

    let mut cumulative = ds[0] / 2.0; // Start at center of first panel
    s.push(cumulative);

    for i in 1..n {
        cumulative += (ds[i - 1] + ds[i]) / 2.0;
        s.push(cumulative);
    }

    s
}

/// Compute arc lengths from stagnation point for one surface.
///
/// This computes cumulative arc length starting from 0 at the first point
/// (which should be the stagnation point) and increasing downstream.
///
/// # Arguments
/// * `x` - x-coordinates (stagnation first, then downstream)
/// * `y` - y-coordinates
///
/// # Returns
/// Arc lengths starting at 0 at stagnation, increasing downstream.
pub fn compute_arc_from_stagnation(x: &[f64], y: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), y.len(), "x and y must have same length");
    let n = x.len();

    if n == 0 {
        return vec![];
    }

    let mut s = vec![0.0; n];

    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }

    s
}

/// Extract one surface from stagnation to trailing edge.
///
/// Splits the airfoil at the stagnation point and extracts coordinates
/// and velocities for either the upper or lower surface, ordered from
/// stagnation downstream toward the trailing edge.
///
/// # Arguments
/// * `stag_idx` - Index of the stagnation point
/// * `x`, `y` - Full airfoil coordinates
/// * `ue` - Edge velocities from inviscid solution
/// * `is_upper` - true for upper surface, false for lower
///
/// # Returns
/// (x, y, ue) arrays for the surface, ordered from stagnation downstream.
/// Edge velocities are converted to positive (absolute) values.
pub fn extract_surface(
    stag_idx: usize,
    x: &[f64],
    y: &[f64],
    ue: &[f64],
    is_upper: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // For Selig format airfoils: TE(upper) → upper surface → LE → lower surface → TE(lower)
    // Indices 0..stag_idx: TE-upper to LE (upper surface, needs reversal for LE→TE direction)
    // Indices stag_idx..: LE to TE-lower (lower surface, already correct direction)
    if is_upper {
        // Upper surface: indices 0 to stag_idx, reversed to go LE → TE (downstream)
        (
            x[..=stag_idx].iter().rev().cloned().collect(),
            y[..=stag_idx].iter().rev().cloned().collect(),
            ue[..=stag_idx].iter().rev().map(|u| u.abs()).collect(),
        )
    } else {
        // Lower surface: indices stag_idx to end, LE → TE (already downstream direction)
        (
            x[stag_idx..].to_vec(),
            y[stag_idx..].to_vec(),
            ue[stag_idx..].iter().map(|u| u.abs()).collect(),
        )
    }
}

/// Find stagnation point index (where Ue is minimum or changes sign).
///
/// The stagnation point is where the surface velocity is approximately zero.
/// For a typical airfoil, this is near the leading edge.
///
/// # Arguments
/// * `ue` - Edge velocity array
///
/// # Returns
/// Index of the stagnation point (where |Ue| is minimum).
pub fn find_stagnation(ue: &[f64]) -> usize {
    if ue.is_empty() {
        return 0;
    }

    // Find where |Ue| is minimum
    ue.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Find stagnation by sign change in velocity.
///
/// Alternative stagnation finder that looks for where velocity changes sign,
/// which is more robust for some cases.
///
/// # Arguments
/// * `ue` - Edge velocity array
///
/// # Returns
/// Index of the stagnation point (where Ue changes sign).
pub fn find_stagnation_by_sign_change(ue: &[f64]) -> usize {
    if ue.len() < 2 {
        return 0;
    }

    for i in 1..ue.len() {
        if ue[i] * ue[i - 1] < 0.0 {
            // Sign change - stagnation is between i-1 and i
            // Return the one with smaller |Ue|
            if ue[i].abs() < ue[i - 1].abs() {
                return i;
            } else {
                return i - 1;
            }
        }
    }

    // No sign change found - use minimum |Ue|
    find_stagnation(ue)
}

/// Interpolate exact stagnation point location (XFOIL STFIND).
///
/// Finds where gamma changes from positive to negative and interpolates
/// the exact arc length where Ue = 0. This is critical for proper BL setup
/// because both surfaces must start at the true stagnation point with Ue ≈ 0
/// and have monotonically increasing velocity downstream.
///
/// # Arguments
/// * `gamma` - Surface velocities (signed, from inviscid solver)
/// * `s` - Arc lengths along full airfoil (from `compute_arc_lengths`)
///
/// # Returns
/// Tuple of (ist, sst, ue_stag):
/// - `ist`: Panel index where sign change occurs (stagnation is between ist and ist+1)
/// - `sst`: Interpolated arc length at exact stagnation point
/// - `ue_stag`: Interpolated velocity at stagnation (should be ≈ 0)
///
/// # XFOIL Reference
/// STFIND in xpanel.f (line 1357)
pub fn interpolate_stagnation(gamma: &[f64], s: &[f64]) -> (usize, f64, f64) {
    assert_eq!(gamma.len(), s.len(), "gamma and s must have same length");
    
    if gamma.len() < 2 {
        return (0, s.get(0).copied().unwrap_or(0.0), gamma.get(0).copied().unwrap_or(0.0));
    }
    
    // Find where gamma changes from positive to negative (XFOIL convention)
    for i in 0..gamma.len() - 1 {
        if gamma[i] >= 0.0 && gamma[i + 1] < 0.0 {
            let dgam = gamma[i + 1] - gamma[i];
            let ds = s[i + 1] - s[i];
            
            // Linear interpolation to find exact crossing point
            // Use the formula that minimizes roundoff (from XFOIL)
            let sst = if gamma[i] < -gamma[i + 1] {
                s[i] - ds * (gamma[i] / dgam)
            } else {
                s[i + 1] - ds * (gamma[i + 1] / dgam)
            };
            
            // Tweak if stagnation falls exactly on a node (very unlikely)
            let sst = if sst <= s[i] {
                s[i] + 1.0e-10
            } else if sst >= s[i + 1] {
                s[i + 1] - 1.0e-10
            } else {
                sst
            };
            
            // Interpolate Ue at stagnation (should be very close to 0)
            let frac = (sst - s[i]) / ds;
            let ue_stag = gamma[i] + frac * dgam;
            
            return (i, sst, ue_stag);
        }
    }
    
    // No sign change found - fallback to minimum |gamma|
    let ist = find_stagnation(gamma);
    (ist, s[ist], gamma[ist])
}

/// Extract surface with XFOIL-style arc lengths from interpolated stagnation.
///
/// This function implements XFOIL's IBLPAN + XICALC approach:
/// - Includes a virtual stagnation point at arc=0 (like XFOIL's IBL=1)
/// - Upper surface: virtual stagnation + panels from IST down to 0
/// - Lower surface: virtual stagnation + panels from IST+1 to N
///
/// Both surfaces start at arc=0 (stagnation) with Ue≈0 (interpolated) and have
/// monotonically increasing arc lengths.
///
/// # Arguments
/// * `ist` - Stagnation panel index (from `interpolate_stagnation`)
/// * `sst` - Interpolated stagnation arc length
/// * `ue_stag` - Interpolated Ue at stagnation (should be ≈0)
/// * `s` - Full airfoil arc lengths
/// * `x`, `y` - Full airfoil coordinates
/// * `ue` - Edge velocities (signed)
/// * `is_upper` - true for upper surface, false for lower
///
/// # Returns
/// Tuple of (arc, x, y, ue) for the surface:
/// - `arc`: Arc lengths from stagnation (starts at 0 for virtual stagnation point)
/// - `x`, `y`: Surface coordinates
/// - `ue`: Edge velocities (absolute values)
///
/// # XFOIL Reference
/// IBLPAN (line 1395) + XICALC (line 1455) in xpanel.f
pub fn extract_surface_xfoil(
    ist: usize,
    sst: f64,
    ue_stag: f64,
    s: &[f64],
    x: &[f64],
    y: &[f64],
    ue: &[f64],
    is_upper: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = s.len();
    assert_eq!(n, x.len());
    assert_eq!(n, y.len());
    assert_eq!(n, ue.len());
    
    // Minimum arc length near stagnation (XFOIL XFEPS)
    let xeps = 1.0e-10 * (s[n - 1] - s[0]).abs().max(1.0);
    
    // Interpolate position at stagnation point
    // sst is between s[ist] and s[ist+1]
    let ist_next = (ist + 1).min(n - 1);
    let frac = if s[ist_next] != s[ist] {
        ((sst - s[ist]) / (s[ist_next] - s[ist])).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let x_stag = x[ist] + frac * (x[ist_next] - x[ist]);
    let y_stag = y[ist] + frac * (y[ist_next] - y[ist]);
    
    if is_upper {
        // Upper surface: virtual stagnation + panels IST down to 0
        // Arc length = SST - S[i] (S decreases as we go upstream, so arc increases)
        // IST is the first real panel before stagnation, matching XFOIL's IBL=2
        let indices: Vec<usize> = (0..=ist).rev().collect();
        
        // Start with virtual stagnation point (arc=0, Ue≈0)
        let mut arc = vec![xeps]; // Very small arc, not exactly 0 to avoid division issues
        let mut x_surf = vec![x_stag];
        let mut y_surf = vec![y_stag];
        let mut ue_surf = vec![ue_stag.abs().max(0.001)]; // Small but positive Ue
        
        // Add real panels
        for &i in &indices {
            arc.push((sst - s[i]).max(xeps));
            x_surf.push(x[i]);
            y_surf.push(y[i]);
            ue_surf.push(ue[i].abs());
        }
        
        (arc, x_surf, y_surf, ue_surf)
    } else {
        // Lower surface: virtual stagnation + panels IST+1 to N-1
        // Arc length = S[i] - SST (S increases as we go downstream)
        let start = (ist + 1).min(n - 1);
        let indices: Vec<usize> = (start..n).collect();
        
        // Start with virtual stagnation point (arc=0, Ue≈0)
        let mut arc = vec![xeps];
        let mut x_surf = vec![x_stag];
        let mut y_surf = vec![y_stag];
        let mut ue_surf = vec![ue_stag.abs().max(0.001)];
        
        // Add real panels
        for &i in &indices {
            arc.push((s[i] - sst).max(xeps));
            x_surf.push(x[i]);
            y_surf.push(y[i]);
            ue_surf.push(ue[i].abs());
        }
        
        (arc, x_surf, y_surf, ue_surf)
    }
}

/// Initialize BL stations from inviscid solution.
///
/// Creates a vector of BlStation objects initialized with the inviscid
/// edge velocity and ready for the viscous iteration.
///
/// # Arguments
/// * `arc_lengths` - Arc length coordinates (from stagnation point)
/// * `ue_inviscid` - Edge velocities from inviscid solution
/// * `re` - Reynolds number
/// * `stag_idx` - Index of stagnation point
///
/// # Returns
/// Vector of initialized BlStation objects.
pub fn initialize_bl_stations(
    arc_lengths: &[f64],
    ue_inviscid: &[f64],
    re: f64,
    stag_idx: usize,
) -> Vec<BlStation> {
    let n = arc_lengths.len();
    assert_eq!(n, ue_inviscid.len(), "Arrays must have same length");

    let mut stations = Vec::with_capacity(n);

    for i in 0..n {
        let mut station = if i == stag_idx {
            // At stagnation point, use Hiemenz solution
            let ue_stag = ue_inviscid[i].abs().max(0.01);
            BlStation::stagnation(ue_stag, re)
        } else {
            // Away from stagnation, use default initialization
            BlStation::new()
        };

        // Set position and edge velocity
        station.x = arc_lengths[i];
        station.u = ue_inviscid[i];

        // Compute mass defect
        station.mass_defect = station.u * station.delta_star;

        // Initialize R_theta
        station.r_theta = re * station.u.abs() * station.theta;

        stations.push(station);
    }

    stations
}

/// Initialize BL stations for one surface (upper or lower).
///
/// This creates stations going downstream from the stagnation point
/// toward the trailing edge on either the upper or lower surface.
///
/// # Arguments
/// * `arc_lengths` - Arc lengths from stagnation (should increase downstream)
/// * `ue` - Edge velocities (all positive for attached flow)
/// * `re` - Reynolds number
///
/// # Returns
/// Vector of BlStation objects for one surface.
pub fn initialize_surface_stations(arc_lengths: &[f64], ue: &[f64], re: f64) -> Vec<BlStation> {
    let n = arc_lengths.len();
    assert_eq!(n, ue.len(), "Arrays must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut stations = Vec::with_capacity(n);

    // First station is at stagnation
    let ue_stag = ue[0].abs().max(0.01);
    let mut station = BlStation::stagnation(ue_stag, re);
    station.x = arc_lengths[0];
    station.u = ue[0].abs();
    station.mass_defect = station.u * station.delta_star;
    stations.push(station);

    // Remaining stations
    for i in 1..n {
        let mut station = BlStation::new();
        station.x = arc_lengths[i];
        station.u = ue[i].abs();
        station.is_laminar = true;
        station.is_turbulent = false;
        station.r_theta = re * station.u * station.theta;
        station.mass_defect = station.u * station.delta_star;
        stations.push(station);
    }

    stations
}

// =============================================================================
// Integration with rustfoil-inviscid
// =============================================================================

use rustfoil_core::Body;
use rustfoil_inviscid::{FlowConditions as InviscidFlowConditions, InviscidSolution, InviscidSolver};

/// Error type for setup operations
#[derive(Debug)]
pub enum SetupError {
    /// Inviscid solver error
    InviscidError(String),
    /// Geometry error
    GeometryError(String),
}

impl std::fmt::Display for SetupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SetupError::InviscidError(msg) => write!(f, "Inviscid solver error: {}", msg),
            SetupError::GeometryError(msg) => write!(f, "Geometry error: {}", msg),
        }
    }
}

impl std::error::Error for SetupError {}

impl From<rustfoil_inviscid::InviscidError> for SetupError {
    fn from(e: rustfoil_inviscid::InviscidError) -> Self {
        SetupError::InviscidError(e.to_string())
    }
}

/// Result of viscous setup from body.
#[derive(Debug, Clone)]
pub struct ViscousSetupResult {
    /// Inviscid solution
    pub inviscid: InviscidSolution,
    /// Viscous solver setup data
    pub setup: ViscousSetup,
    /// Node x-coordinates
    pub node_x: Vec<f64>,
    /// Node y-coordinates  
    pub node_y: Vec<f64>,
    /// Stagnation panel index (IST)
    pub ist: usize,
    /// Interpolated stagnation arc length (SST)
    pub sst: f64,
}

/// Setup viscous solver from Body and angle of attack.
///
/// This function:
/// 1. Converts Body coordinates to the format expected by rustfoil-inviscid
/// 2. Runs the inviscid solver at the given angle of attack
/// 3. Extracts stagnation point using XFOIL's STFIND algorithm
/// 4. Builds the ViscousSetup structure
///
/// # Arguments
/// * `body` - Airfoil body with panel geometry
/// * `alpha_deg` - Angle of attack in degrees
///
/// # Returns
/// `ViscousSetupResult` containing the inviscid solution and viscous setup data.
///
/// # Example
/// ```rust,ignore
/// use rustfoil_core::Body;
/// use rustfoil_solver::viscous::setup_from_body;
///
/// let body = Body::from_naca("0012", 160)?;
/// let result = setup_from_body(&body, 4.0)?;
/// println!("CL = {:.4}", result.inviscid.cl);
/// ```
pub fn setup_from_body(body: &Body, alpha_deg: f64) -> Result<ViscousSetupResult, SetupError> {
    let panels = body.panels();
    let n_panels = panels.len();
    
    if n_panels < 3 {
        return Err(SetupError::GeometryError(
            "Need at least 3 panels".to_string()
        ));
    }
    
    // Extract node coordinates (p1 of each panel)
    let mut node_x: Vec<f64> = panels.iter().map(|p| p.p1.x).collect();
    let mut node_y: Vec<f64> = panels.iter().map(|p| p.p1.y).collect();
    
    // Check if contour is closed (sharp TE) or open (blunt TE)
    let last_p2 = panels.last().unwrap().p2;
    let first_p1 = panels[0].p1;
    let is_closed = (last_p2.x - first_p1.x).abs() < 1e-10 
                 && (last_p2.y - first_p1.y).abs() < 1e-10;
    
    if !is_closed {
        // Blunt TE: add the lower TE node
        node_x.push(last_p2.x);
        node_y.push(last_p2.y);
    }
    
    // Convert to tuples for rustfoil-inviscid
    let coords: Vec<(f64, f64)> = node_x.iter()
        .zip(node_y.iter())
        .map(|(&x, &y)| (x, y))
        .collect();
    
    // Run inviscid solver
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&coords)?;
    let flow = InviscidFlowConditions::with_alpha_deg(alpha_deg);
    let inv_solution = factorized.solve_alpha(&flow);
    
    // gamma IS the surface velocity in XFOIL's formulation
    let ue_inviscid = inv_solution.gamma.clone();
    
    // Compute arc lengths
    let arc_lengths = compute_arc_lengths(&node_x, &node_y);
    
    // Find stagnation using XFOIL's interpolation method
    let (ist, sst, _ue_stag) = interpolate_stagnation(&ue_inviscid, &arc_lengths);
    
    // Build DIJ matrix
    let dij = build_dij_matrix(&node_x, &node_y);
    
    // Determine upper/lower counts based on stagnation
    let n = arc_lengths.len();
    let n_upper = n - ist;
    let n_lower = ist + 1;
    
    let setup = ViscousSetup {
        arc_lengths: arc_lengths.clone(),
        ue_inviscid,
        dij,
        panel_x: node_x.clone(),
        panel_y: node_y.clone(),
        stagnation_index: ist,
        n_upper,
        n_lower,
    };
    
    Ok(ViscousSetupResult {
        inviscid: inv_solution,
        setup,
        node_x,
        node_y,
        ist,
        sst,
    })
}

/// Setup viscous solver from raw coordinates.
///
/// Alternative to `setup_from_body()` when you have coordinates directly
/// (e.g., from XFOIL coordinate file).
///
/// # Arguments
/// * `coords` - Airfoil coordinates as (x, y) pairs, counter-clockwise from upper TE
/// * `alpha_deg` - Angle of attack in degrees
///
/// # Returns
/// `ViscousSetupResult` containing the inviscid solution and viscous setup data.
pub fn setup_from_coords(coords: &[(f64, f64)], alpha_deg: f64) -> Result<ViscousSetupResult, SetupError> {
    if coords.len() < 4 {
        return Err(SetupError::GeometryError(
            "Need at least 4 coordinate points".to_string()
        ));
    }
    
    let node_x: Vec<f64> = coords.iter().map(|(x, _)| *x).collect();
    let node_y: Vec<f64> = coords.iter().map(|(_, y)| *y).collect();
    
    // Run inviscid solver
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(coords)?;
    let flow = InviscidFlowConditions::with_alpha_deg(alpha_deg);
    let inv_solution = factorized.solve_alpha(&flow);
    
    let ue_inviscid = inv_solution.gamma.clone();
    let arc_lengths = compute_arc_lengths(&node_x, &node_y);
    let (ist, sst, _) = interpolate_stagnation(&ue_inviscid, &arc_lengths);
    let dij = build_dij_matrix(&node_x, &node_y);
    
    let n = arc_lengths.len();
    let n_upper = n - ist;
    let n_lower = ist + 1;
    
    let setup = ViscousSetup {
        arc_lengths: arc_lengths.clone(),
        ue_inviscid,
        dij,
        panel_x: node_x.clone(),
        panel_y: node_y.clone(),
        stagnation_index: ist,
        n_upper,
        n_lower,
    };
    
    Ok(ViscousSetupResult {
        inviscid: inv_solution,
        setup,
        node_x,
        node_y,
        ist,
        sst,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_arc_lengths() {
        // Simple case: 3 points on a line
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 0.0, 0.0];

        let s = compute_arc_lengths(&x, &y);

        assert_eq!(s.len(), 3);
        assert!((s[0] - 0.0).abs() < 1e-10);
        assert!((s[1] - 1.0).abs() < 1e-10);
        assert!((s[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_arc_lengths_diagonal() {
        // Diagonal line: sqrt(2) between points
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, 2.0];

        let s = compute_arc_lengths(&x, &y);

        let sqrt2 = 2.0_f64.sqrt();
        assert!((s[0] - 0.0).abs() < 1e-10);
        assert!((s[1] - sqrt2).abs() < 1e-10);
        assert!((s[2] - 2.0 * sqrt2).abs() < 1e-10);
    }

    #[test]
    fn test_find_stagnation() {
        // Velocity distribution with minimum near middle
        let ue = vec![0.5, 0.3, 0.1, -0.05, -0.2, -0.4];

        let idx = find_stagnation(&ue);
        assert_eq!(idx, 3); // -0.05 has smallest |Ue|
    }

    #[test]
    fn test_find_stagnation_by_sign_change() {
        let ue = vec![0.5, 0.3, 0.1, -0.05, -0.2, -0.4];

        let idx = find_stagnation_by_sign_change(&ue);
        assert!(idx == 2 || idx == 3); // Sign change is between 2 and 3
    }

    #[test]
    fn test_find_stagnation_empty() {
        let ue: Vec<f64> = vec![];
        assert_eq!(find_stagnation(&ue), 0);
    }

    #[test]
    fn test_initialize_bl_stations() {
        let arc_lengths = vec![0.0, 0.1, 0.2, 0.3];
        let ue = vec![0.1, 0.5, 0.8, 1.0];
        let re = 1e6;

        let stations = initialize_bl_stations(&arc_lengths, &ue, re, 0);

        assert_eq!(stations.len(), 4);

        // First station should have stagnation initialization
        assert!((stations[0].h - 2.216).abs() < 0.01);

        // All stations should have correct x and u
        for (i, station) in stations.iter().enumerate() {
            assert!((station.x - arc_lengths[i]).abs() < 1e-10);
            assert!((station.u - ue[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_initialize_surface_stations() {
        let arc_lengths = vec![0.0, 0.1, 0.2, 0.3, 0.4];
        let ue = vec![0.1, 0.5, 0.8, 0.9, 1.0];
        let re = 1e6;

        let stations = initialize_surface_stations(&arc_lengths, &ue, re);

        assert_eq!(stations.len(), 5);

        // First station should be stagnation-like
        assert!((stations[0].h - 2.216).abs() < 0.01);

        // All should be laminar initially
        for station in &stations {
            assert!(station.is_laminar);
            assert!(!station.is_turbulent);
        }
    }

    #[test]
    fn test_viscous_setup_from_raw() {
        let arc = vec![0.0, 0.1, 0.2, 0.3, 0.4];
        // Design ue so that index 2 has the smallest |Ue|
        let ue = vec![0.5, 0.2, 0.05, -0.2, -0.5];
        let x = vec![0.0, 0.1, 0.2, 0.3, 0.4];
        let y = vec![0.05, 0.03, 0.0, -0.03, -0.05];

        let setup = ViscousSetup::from_raw(arc.clone(), ue.clone(), x, y);

        assert_eq!(setup.arc_lengths.len(), 5);
        assert_eq!(setup.stagnation_index, 2); // 0.05 has smallest |Ue|
        assert_eq!(setup.dij.nrows(), 5);
        assert_eq!(setup.dij.ncols(), 5);
    }

    #[test]
    fn test_cumulative_arc_length() {
        let ds = vec![0.1, 0.1, 0.1, 0.1];

        let s = cumulative_arc_length(&ds);

        assert_eq!(s.len(), 4);
        // First station at center of first panel: 0.05
        assert!((s[0] - 0.05).abs() < 1e-10);
        // Subsequent stations at centers
        assert!((s[1] - 0.15).abs() < 1e-10);
        assert!((s[2] - 0.25).abs() < 1e-10);
        assert!((s[3] - 0.35).abs() < 1e-10);
    }

    #[test]
    fn test_compute_arc_from_stagnation() {
        // Simple horizontal line from stagnation
        let x = vec![0.0, 0.1, 0.2, 0.3];
        let y = vec![0.0, 0.0, 0.0, 0.0];

        let s = compute_arc_from_stagnation(&x, &y);

        assert_eq!(s.len(), 4);
        assert!((s[0] - 0.0).abs() < 1e-10); // Starts at 0 (stagnation)
        assert!((s[1] - 0.1).abs() < 1e-10);
        assert!((s[2] - 0.2).abs() < 1e-10);
        assert!((s[3] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_extract_surface_upper() {
        // Airfoil-like geometry: TE -> LE (lower) -> LE (stag) -> TE (upper)
        // Indices:                 0       1           2           3      4
        let x = vec![1.0, 0.5, 0.0, 0.5, 1.0];
        let y = vec![-0.05, -0.08, 0.0, 0.08, 0.05];
        let ue = vec![-0.5, -0.3, 0.05, 0.3, 0.5]; // Signed velocities

        let stag_idx = 2;

        let (ux, uy, u_ue) = extract_surface(stag_idx, &x, &y, &ue, true);

        // Upper surface: indices 2, 3, 4
        assert_eq!(ux.len(), 3);
        assert!((ux[0] - 0.0).abs() < 1e-10); // Starts at stagnation
        assert!((ux[1] - 0.5).abs() < 1e-10);
        assert!((ux[2] - 1.0).abs() < 1e-10);

        // Ue should be absolute values
        assert!((u_ue[0] - 0.05).abs() < 1e-10);
        assert!((u_ue[1] - 0.3).abs() < 1e-10);
        assert!((u_ue[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_extract_surface_lower() {
        // Same geometry as above
        let x = vec![1.0, 0.5, 0.0, 0.5, 1.0];
        let y = vec![-0.05, -0.08, 0.0, 0.08, 0.05];
        let ue = vec![-0.5, -0.3, 0.05, 0.3, 0.5];

        let stag_idx = 2;

        let (lx, ly, l_ue) = extract_surface(stag_idx, &x, &y, &ue, false);

        // Lower surface: indices 2, 1, 0 (reversed to go downstream)
        assert_eq!(lx.len(), 3);
        assert!((lx[0] - 0.0).abs() < 1e-10); // Starts at stagnation
        assert!((lx[1] - 0.5).abs() < 1e-10);
        assert!((lx[2] - 1.0).abs() < 1e-10);

        // Ue should be absolute values (reversed order)
        assert!((l_ue[0] - 0.05).abs() < 1e-10);
        assert!((l_ue[1] - 0.3).abs() < 1e-10);
        assert!((l_ue[2] - 0.5).abs() < 1e-10);
    }
}
