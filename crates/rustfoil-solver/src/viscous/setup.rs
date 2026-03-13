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
use rustfoil_bl::debug;
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::stmove::find_stagnation_with_derivs;
use rustfoil_coupling::dij::build_dij_matrix;
use super::state::XfoilLikeViscousState;

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

    /// Inviscid wake x-coordinates used to build the wake-aware DIJ matrix.
    pub wake_x: Vec<f64>,

    /// Inviscid wake y-coordinates used to build the wake-aware DIJ matrix.
    pub wake_y: Vec<f64>,

    /// Wake arc-length coordinates measured from the airfoil trailing edge.
    pub wake_s: Vec<f64>,

    /// Inviscid wake-edge velocities used to seed the wake BL march.
    pub wake_ue: Vec<f64>,
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
            wake_x: Vec::new(),
            wake_y: Vec::new(),
            wake_s: Vec::new(),
            wake_ue: Vec::new(),
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
/// - `ue`: Edge velocities in the split XFOIL surface convention
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
        
        // Start with a geometric stagnation anchor. Like XFOIL's IBL=1, this
        // is not a real panel-owned velocity station.
        let mut arc = vec![0.0];
        let mut x_surf = vec![x_stag];
        let mut y_surf = vec![y_stag];
        let mut ue_surf = vec![ue_stag.signum() * ue_stag.abs().min(1.0e-12)];
        
        // Add real panels
        for &i in &indices {
            arc.push((sst - s[i]).max(xeps));
            x_surf.push(x[i]);
            y_surf.push(y[i]);
            ue_surf.push(ue[i]);
        }
        
        (arc, x_surf, y_surf, ue_surf)
    } else {
        // Lower surface: virtual stagnation + panels IST+1 to N-1
        // Arc length = S[i] - SST (S increases as we go downstream)
        let start = (ist + 1).min(n - 1);
        let indices: Vec<usize> = (start..n).collect();
        
        // Start with a geometric stagnation anchor. Like XFOIL's IBL=1, this
        // is not a real panel-owned velocity station.
        let mut arc = vec![0.0];
        let mut x_surf = vec![x_stag];
        let mut y_surf = vec![y_stag];
        let mut ue_surf = vec![ue_stag.signum() * ue_stag.abs().min(1.0e-12)];
        
        // Add real panels
        for &i in &indices {
            arc.push((s[i] - sst).max(xeps));
            x_surf.push(x[i]);
            y_surf.push(y[i]);
            // Full-contour QINV/GAM is negative on the lower surface, while
            // XFOIL's split UINV(IBL,2) is the positive edge-speed magnitude.
            ue_surf.push(-ue[i]);
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
            // At stagnation point, use Thwaites formula (XFOIL-compatible)
            let x_stag = arc_lengths[i].max(1e-12);
            let ue_stag = ue_inviscid[i].abs().max(0.01);
            BlStation::stagnation(x_stag, ue_stag, re)
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
///
/// # Note
/// This version does not populate `x_coord`. Use `initialize_surface_stations_with_coords`
/// if you need proper x/c transition location reporting.
pub fn initialize_surface_stations(arc_lengths: &[f64], ue: &[f64], re: f64) -> Vec<BlStation> {
    // Fallback: use arc lengths as x_coord (will give > 1.0 values for x_tr)
    initialize_surface_stations_with_coords(arc_lengths, ue, arc_lengths, re)
}

/// Initialize BL stations for one surface with panel x-coordinates for x/c reporting.
///
/// This creates stations going downstream from the stagnation point
/// toward the trailing edge on either the upper or lower surface.
///
/// # Arguments
/// * `arc_lengths` - Arc lengths from stagnation (should increase downstream)
/// * `ue` - Edge velocities (all positive for attached flow)
/// * `panel_x` - Panel x-coordinates (for transition x/c reporting)
/// * `re` - Reynolds number
///
/// # Returns
/// Vector of BlStation objects for one surface.
pub fn initialize_surface_stations_with_coords(
    arc_lengths: &[f64],
    ue: &[f64],
    panel_x: &[f64],
    re: f64,
) -> Vec<BlStation> {
    let n = arc_lengths.len();
    assert_eq!(n, ue.len(), "arc_lengths and ue must have same length");
    assert_eq!(n, panel_x.len(), "arc_lengths and panel_x must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut stations = Vec::with_capacity(n);

    // First station is a virtual stagnation anchor. Seed its thickness from
    // the first real downstream station, then zero the panel-owned quantities.
    let seed_x = arc_lengths.get(1).copied().unwrap_or(1.0e-6).abs().max(1.0e-12);
    let seed_ue = ue.get(1).copied().unwrap_or(0.01).abs().max(0.01);
    let mut station = BlStation::stagnation(seed_x, seed_ue, re);
    station.x = arc_lengths[0];
    station.x_coord = panel_x[0];
    station.u = 0.0;
    station.mass_defect = 0.0;
    stations.push(station);

    // Remaining stations
    for i in 1..n {
        let mut station = BlStation::new();
        station.x = arc_lengths[i];
        station.x_coord = panel_x[i];
        station.u = ue[i].abs();
        station.is_laminar = true;
        station.is_turbulent = false;
        station.r_theta = re * station.u * station.theta;
        station.mass_defect = station.u * station.delta_star;
        stations.push(station);
    }

    stations
}

/// Initialize BL stations with panel indices for VI coupling.
///
/// This version properly maps each BL station to its corresponding panel
/// in the global DIJ matrix, enabling proper viscous-inviscid coupling.
///
/// # Arguments
/// * `arc_lengths` - Arc lengths from stagnation (should increase downstream)
/// * `ue` - Edge velocities (all positive for attached flow)
/// * `panel_x` - Panel x-coordinates (for transition x/c reporting)
/// * `stagnation_idx` - Index of stagnation point in full panel array
/// * `n_airfoil_panels` - Total number of airfoil panels (for wake panel offset)
/// * `is_upper` - True for upper surface, false for lower
/// * `re` - Reynolds number
///
/// # Panel Index Mapping (Selig format)
/// - Upper surface: panel_idx = stagnation_idx - local_idx (decreasing toward panel 0 = upper TE)
///   - When local_idx > stagnation_idx, maps to wake panels N, N+1, ...
/// - Lower surface: panel_idx = stagnation_idx + local_idx (increasing toward panel N-1 = lower TE)
///   - When stagnation_idx + local_idx >= N, maps to wake panels N, N+1, ...
pub fn initialize_surface_stations_with_panel_idx(
    arc_lengths: &[f64],
    ue: &[f64],
    panel_x: &[f64],
    stagnation_idx: usize,
    n_airfoil_panels: usize,
    is_upper: bool,
    re: f64,
) -> Vec<BlStation> {
    let n = arc_lengths.len();
    assert_eq!(n, ue.len(), "arc_lengths and ue must have same length");
    assert_eq!(n, panel_x.len(), "arc_lengths and panel_x must have same length");
    assert!(n >= 2, "Need at least 2 stations");

    let mut stations = Vec::with_capacity(n);

    // First station is a virtual stagnation anchor. Seed its thickness from
    // the first real downstream station, then zero the panel-owned quantities.
    let seed_x = arc_lengths.get(1).copied().unwrap_or(1.0e-6).abs().max(1.0e-12);
    let seed_ue = ue.get(1).copied().unwrap_or(0.01).abs().max(0.01);
    let mut station = BlStation::stagnation(seed_x, seed_ue, re);
    station.x = arc_lengths[0];
    station.x_coord = panel_x[0];
    station.panel_idx = usize::MAX;
    station.u = 0.0;
    station.mass_defect = 0.0;
    stations.push(station);

    // Remaining stations
    for i in 1..n {
        let mut station = BlStation::new();
        station.x = arc_lengths[i];
        station.x_coord = panel_x[i];
        // Map BL station to panel index
        //
        // IMPORTANT: extract_surface_xfoil() prepends a virtual stagnation point,
        // so the station indices are offset by 1 from the panel data:
        //   - Station 0 = virtual stagnation (interpolated between panels)
        //   - Station 1 = data from panel ist
        //   - Station 2 = data from panel ist-1 (for upper) or ist+1 (for lower)
        //   - etc.
        //
        // Upper surface: stations go from virtual stag (0) toward upper TE (panel 0)
        //   - Station i (for i >= 1) corresponds to panel ist - (i - 1)
        //   - When ist - (i - 1) < 0, we've passed TE and entered wake
        //
        // Lower surface: stations go from virtual stag (0) toward lower TE (panel N-1)
        //   - Station i (for i >= 1) corresponds to panel ist + i
        //   - When ist + i >= N, we've entered wake
        station.panel_idx = if is_upper {
            // Upper surface: station i maps to panel ist - (i - 1)
            // The -1 accounts for the virtual stagnation point at station 0
            let panel = stagnation_idx as i64 - (i as i64 - 1);
            if panel >= 0 {
                panel as usize  // Normal airfoil panel
            } else {
                // Wake panel: station is past upper TE (panel 0)
                // XFOIL assigns wake panels indices N, N+1, N+2, ...
                n_airfoil_panels + (-panel - 1) as usize
            }
        } else {
            // Lower surface: station i maps to panel ist + i
            // This naturally gives N, N+1, ... for wake when large enough
            stagnation_idx + i
        };
        station.u = ue[i].abs();
        station.is_laminar = true;
        station.is_turbulent = false;
        station.r_theta = re * station.u * station.theta;
        station.mass_defect = station.u * station.delta_star;
        stations.push(station);
    }

    // Debug: print panel indices for first few and last few stations
    if rustfoil_bl::is_debug_active() {
        let surface_name = if is_upper { "Upper" } else { "Lower" };
        eprintln!("[DEBUG setup] {} surface: n={}, stagnation_idx={}, n_airfoil={}", 
            surface_name, n, stagnation_idx, n_airfoil_panels);
        eprintln!("[DEBUG setup]   First 5 panel_idx: {:?}", 
            stations.iter().take(5).map(|s| s.panel_idx).collect::<Vec<_>>());
        eprintln!("[DEBUG setup]   Last 5 panel_idx: {:?}", 
            stations.iter().rev().take(5).rev().map(|s| s.panel_idx).collect::<Vec<_>>());
    }

    stations
}

// =============================================================================
// Integration with rustfoil-inviscid
// =============================================================================

use rustfoil_core::Body;
use rustfoil_inviscid::influence::psilin_with_dqdm;
use rustfoil_inviscid::system::{compute_wake_count, setexp};
use rustfoil_inviscid::{
    FactorizedSystem, FlowConditions as InviscidFlowConditions, InviscidSolution, InviscidSolver,
};

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

impl ViscousSetupResult {
    /// Derive split `BlStation` views from the canonical setup state.
    pub fn derive_station_views(
        &self,
        reynolds: f64,
        msq: f64,
    ) -> (Vec<BlStation>, Vec<BlStation>) {
        let canonical =
            initialize_canonical_state_from_setup_result(self, reynolds, msq);
        (
            canonical.upper_station_view(),
            canonical.lower_station_view(),
        )
    }

    /// Derive split `Ue` views from the canonical setup state.
    pub fn derive_uedg_views(
        &self,
        reynolds: f64,
        msq: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let canonical =
            initialize_canonical_state_from_setup_result(self, reynolds, msq);
        (
            canonical.upper_uedg_view(),
            canonical.lower_uedg_view(),
        )
    }

    /// Expose the inviscid alpha-sensitivity column used by the operating-variable path.
    pub fn operating_sensitivity(&self) -> Vec<f64> {
        if self.inviscid.gamma_a.len() == self.setup.ue_inviscid.len() {
            self.inviscid.gamma_a.clone()
        } else {
            vec![0.0; self.setup.ue_inviscid.len()]
        }
    }
}

/// Build the canonical XFOIL-style viscous setup state from an inviscid setup result.
///
/// Phase 2 makes this the primary setup authoring path. Existing split station
/// vectors remain available only as transitional views derived from the returned
/// canonical state.
pub fn initialize_canonical_state_from_setup_result(
    setup_result: &ViscousSetupResult,
    reynolds: f64,
    msq: f64,
) -> XfoilLikeViscousState {
    initialize_canonical_state(
        &setup_result.setup,
        &setup_result.node_x,
        &setup_result.node_y,
        &setup_result.inviscid.gamma_a,
        setup_result.ist,
        setup_result.sst,
        reynolds,
        msq,
    )
}

fn initialize_canonical_state(
    setup: &ViscousSetup,
    node_x: &[f64],
    node_y: &[f64],
    ue_inviscid_alpha: &[f64],
    ist: usize,
    sst: f64,
    reynolds: f64,
    msq: f64,
) -> XfoilLikeViscousState {
    let full_arc = &setup.arc_lengths;
    let ue_inviscid = &setup.ue_inviscid;
    let ue_stag = if ist + 1 < ue_inviscid.len() && full_arc[ist + 1] != full_arc[ist] {
        let frac = (sst - full_arc[ist]) / (full_arc[ist + 1] - full_arc[ist]);
        ue_inviscid[ist] + frac * (ue_inviscid[ist + 1] - ue_inviscid[ist])
    } else {
        ue_inviscid[ist]
    };

    let (upper_arc, upper_x, _, upper_ue) = extract_surface_xfoil(
        ist,
        sst,
        ue_stag,
        full_arc,
        node_x,
        node_y,
        ue_inviscid,
        true,
    );
    let (lower_arc, lower_x, _, lower_ue) = extract_surface_xfoil(
        ist,
        sst,
        ue_stag,
        full_arc,
        node_x,
        node_y,
        ue_inviscid,
        false,
    );

    let n_airfoil_panels = node_x.len();
    let upper_stations = initialize_surface_stations_with_panel_idx(
        &upper_arc,
        &upper_ue,
        &upper_x,
        ist,
        n_airfoil_panels,
        true,
        reynolds,
    );
    let mut lower_stations = initialize_surface_stations_with_panel_idx(
        &lower_arc,
        &lower_ue,
        &lower_x,
        ist,
        n_airfoil_panels,
        false,
        reynolds,
    );

    let n_wake_panels = compute_n_wake_panels(setup.dij.nrows(), n_airfoil_panels);
    if n_wake_panels > 0 {
        let upper_te = upper_stations.last().cloned().unwrap_or_default();
        let lower_te = lower_stations.last().cloned().unwrap_or_default();
        let ante = if upper_te.panel_idx < node_x.len() && lower_te.panel_idx < node_x.len() {
            let dx = node_x[upper_te.panel_idx] - node_x[lower_te.panel_idx];
            let dy = node_y[upper_te.panel_idx] - node_y[lower_te.panel_idx];
            (dx * dx + dy * dy).sqrt()
        } else {
            0.0
        };
        let wake_stations = initialize_wake_bl_stations(
            &upper_te,
            &lower_te,
            ante,
            n_airfoil_panels,
            &setup.wake_x,
            &setup.wake_s,
            &setup.wake_ue,
            reynolds,
            msq,
        );
        lower_stations.extend(wake_stations);
    }

    let mut state =
        XfoilLikeViscousState::from_station_views(&upper_stations, &lower_stations, n_airfoil_panels);
    if let Some(stag) = find_stagnation_with_derivs(ue_inviscid, full_arc) {
        state.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
    } else {
        state.set_stagnation_metadata(ist, sst, 0.0, 0.0);
    }
    state.set_panel_inviscid_arrays(ue_inviscid, ue_inviscid_alpha);
    state.set_wake_geometry(&setup.wake_x, &setup.wake_y, &setup.wake_s);
    state
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
    
    // Build wake geometry and wake-edge velocity from the inviscid field,
    // matching XFOIL's XYWAKE/QWCALC ownership more closely than the older
    // straight-wake + heuristic-Ue setup path.
    let (wake_x, wake_y, wake_s, wake_ue) =
        build_xfoil_wake_geometry_and_ue(&factorized, &inv_solution.gamma, flow.alpha);
    let dij = factorized.build_dij_with_wake(&wake_x, &wake_y)?;
    
    // Emit FULLDIJ debug event if debug collection is active
    if debug::is_debug_active() {
        let nsys = dij.nrows();
        if nsys > 0 {
            // Flatten to row-major order
            let mut flattened = Vec::with_capacity(nsys * nsys);
            for i in 0..nsys {
                for j in 0..nsys {
                    flattened.push(dij[(i, j)]);
                }
            }
            
            // Extract diagonal sample (first 20 values)
            let diag_sample: Vec<f64> = (0..nsys.min(20)).map(|i| dij[(i, i)]).collect();
            
            // Extract row 1 sample (first 20 values)
            let row1_sample: Vec<f64> = (0..nsys.min(20)).map(|j| dij[(0, j)]).collect();
            
            debug::add_event(debug::DebugEvent::full_dij(nsys, flattened, diag_sample, row1_sample));
        }
    }
    
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
        wake_x,
        wake_y,
        wake_s,
        wake_ue,
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
    
    // Build wake geometry and wake-edge velocity from the inviscid field,
    // matching XFOIL's XYWAKE/QWCALC ownership more closely than the older
    // straight-wake + heuristic-Ue setup path.
    let (wake_x, wake_y, wake_s, wake_ue) =
        build_xfoil_wake_geometry_and_ue(&factorized, &inv_solution.gamma, flow.alpha);
    let dij = factorized.build_dij_with_wake(&wake_x, &wake_y)?;
    
    // Emit FULLDIJ debug event if debug collection is active
    if debug::is_debug_active() {
        let nsys = dij.nrows();
        if nsys > 0 {
            // Flatten to row-major order
            let mut flattened = Vec::with_capacity(nsys * nsys);
            for i in 0..nsys {
                for j in 0..nsys {
                    flattened.push(dij[(i, j)]);
                }
            }
            
            // Extract diagonal sample (first 20 values)
            let diag_sample: Vec<f64> = (0..nsys.min(20)).map(|i| dij[(i, i)]).collect();
            
            // Extract row 1 sample (first 20 values)
            let row1_sample: Vec<f64> = (0..nsys.min(20)).map(|j| dij[(0, j)]).collect();
            
            debug::add_event(debug::DebugEvent::full_dij(nsys, flattened, diag_sample, row1_sample));
        }
    }
    
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
        wake_x,
        wake_y,
        wake_s,
        wake_ue,
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

/// Initialize wake BL stations with proper panel indices for DIJ coupling.
///
/// XFOIL puts wake stations on the lower surface (IS=2) after the TE.
/// Wake stations are indexed IBLTE+1, IBLTE+2, ..., NBL(2) where
/// the panel index for wake station i is N_airfoil + (i - IBLTE - 1).
///
/// # Arguments
/// * `upper_te` - Upper surface trailing edge station (for wake initial condition)
/// * `lower_te` - Lower surface trailing edge station (for wake initial condition)
/// * `n_airfoil_panels` - Number of airfoil panels (N)
/// * `n_wake_panels` - Number of wake panels (from DIJ matrix)
/// * `wake_length` - Wake length in chord lengths
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
///
/// # Returns
/// Vector of wake BlStations with proper panel_idx for DIJ matrix lookup
///
/// # XFOIL Reference
/// - Wake panels: indices N, N+1, N+2, ... (0-indexed)
/// - Wake is placed on lower surface (IS=2) after IBLTE(2)
/// - IPAN(IBL,2) = N + (IBL - IBLTE(2) - 1) for wake stations
pub fn initialize_wake_bl_stations(
    upper_te: &BlStation,
    lower_te: &BlStation,
    ante: f64,
    n_airfoil_panels: usize,
    inviscid_wake_x: &[f64],
    inviscid_wake_s: &[f64],
    inviscid_wake_ue: &[f64],
    re: f64,
    msq: f64,
) -> Vec<BlStation> {
    let n_wake_panels = inviscid_wake_x
        .len()
        .min(inviscid_wake_s.len())
        .min(inviscid_wake_ue.len());
    if n_wake_panels == 0 {
        return Vec::new();
    }
    
    // Combine upper and lower TE for wake initial condition (TESYS) and place
    // the first wake row on the same inviscid wake node used by the DIJ matrix.
    use rustfoil_coupling::wake::{combine_te_for_wake, march_wake};
    
    let mut wake_initial = combine_te_for_wake(upper_te, lower_te, ante);
    let wake_arc_offset = lower_te.x;
    wake_initial.x = wake_arc_offset + inviscid_wake_s[0];
    wake_initial.x_coord = inviscid_wake_x[0];
    wake_initial.u = inviscid_wake_ue[0].abs().max(0.01);
    wake_initial.dw = wake_gap_from_te_distance(ante, inviscid_wake_s[0]);
    wake_initial.mass_defect = wake_initial.u * (wake_initial.delta_star + wake_initial.dw);
    
    let mut wake_stations = if n_wake_panels > 1 {
        let wake_dw: Vec<f64> = inviscid_wake_s[1..n_wake_panels]
            .iter()
            .map(|s| wake_gap_from_te_distance(ante, *s))
            .collect();
        march_wake(
            &wake_initial,
            &inviscid_wake_s[1..n_wake_panels]
                .iter()
                .map(|s| wake_arc_offset + *s)
                .collect::<Vec<_>>(),
            &inviscid_wake_ue[1..n_wake_panels],
            &wake_dw,
            re,
            msq,
        )
    } else {
        vec![wake_initial]
    };
    
    // Set proper panel indices for DIJ matrix lookup
    // Wake panel indices are N, N+1, N+2, ... (0-indexed)
    for (i, station) in wake_stations.iter_mut().enumerate() {
        station.panel_idx = n_airfoil_panels + i;
        station.is_wake = true;
        station.is_turbulent = true;
        station.is_laminar = false;
        station.dw = inviscid_wake_s
            .get(i)
            .map(|s| wake_gap_from_te_distance(ante, *s))
            .unwrap_or(station.dw);
        station.x = inviscid_wake_s
            .get(i)
            .map(|s| wake_arc_offset + *s)
            .unwrap_or(station.x);
        station.x_coord = inviscid_wake_x.get(i).copied().unwrap_or(station.x_coord);
    }
    
    if debug::is_debug_active() {
        eprintln!("[DEBUG wake] Created {} wake stations with panel_idx {} to {}",
            wake_stations.len(),
            n_airfoil_panels,
            n_airfoil_panels + wake_stations.len().saturating_sub(1));
    }
    
    wake_stations
}

fn wake_gap_from_te_distance(ante: f64, wake_s_from_te: f64) -> f64 {
    if ante <= 1.0e-10 {
        return 0.0;
    }

    let telrat = 2.5_f64;
    let dwdxte = 0.0_f64;
    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;
    let zn = 1.0 - wake_s_from_te / (telrat * ante);

    if zn >= 0.0 {
        ante * (aa + bb * zn) * zn * zn
    } else {
        0.0
    }
}

fn build_xfoil_wake_geometry_and_ue(
    factorized: &FactorizedSystem,
    _gamma: &[f64],
    alpha_rad: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let geom = factorized.geometry();
    let n = geom.n;
    if n < 2 {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    }

    let wake_length = 1.0;
    let n_wake = compute_wake_count(n, wake_length);
    if n_wake == 0 {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    }

    let ds1 = 0.5 * ((geom.s[1] - geom.s[0]) + (geom.s[n - 1] - geom.s[n - 2]));
    let snew = setexp(ds1, wake_length * geom.chord, n_wake);
    let xte = 0.5 * (geom.x[0] + geom.x[n - 1]);
    let yte = 0.5 * (geom.y[0] + geom.y[n - 1]);

    let sx = 0.5 * (geom.yp[n - 1] - geom.yp[0]);
    let sy = 0.5 * (geom.xp[0] - geom.xp[n - 1]);
    let smod = (sx * sx + sy * sy).sqrt().max(1.0e-12);
    let nx_init = sx / smod;
    let ny_init = sy / smod;

    let cosa = alpha_rad.cos();
    let sina = alpha_rad.sin();
    let (surface_qinvu_0, surface_qinvu_90) = factorized.surface_qinvu_basis();
    let mut gam = (0..n)
        .map(|i| cosa * factorized.gamu_0[i] + sina * factorized.gamu_90[i])
        .collect::<Vec<_>>();
    if let Some(last) = gam.last_mut() {
        *last = 0.0;
    }

    let mut wake_x = Vec::with_capacity(n_wake);
    let mut wake_y = Vec::with_capacity(n_wake);
    let mut wake_nx = vec![0.0; n_wake];
    let mut wake_ny = vec![0.0; n_wake];

    wake_x.push(xte - 0.0001 * ny_init);
    wake_y.push(yte + 0.0001 * nx_init);
    wake_nx[0] = nx_init;
    wake_ny[0] = ny_init;

    let (psi_x0, psi_y0) =
        compute_psi_gradient_at_wake_point(geom, n, wake_x[0], wake_y[0], &gam, cosa, sina);
    let grad_mag0 = (psi_x0 * psi_x0 + psi_y0 * psi_y0).sqrt().max(1.0e-12);
    let mut nx_next = -psi_x0 / grad_mag0;
    let mut ny_next = -psi_y0 / grad_mag0;
    if n_wake > 1 {
        wake_nx[1] = nx_next;
        wake_ny[1] = ny_next;
    }

    for i in 1..n_wake {
        let ds = snew[i] - snew[i - 1];
        let xi = wake_x[i - 1] - ds * ny_next;
        let yi = wake_y[i - 1] + ds * nx_next;
        wake_x.push(xi);
        wake_y.push(yi);

        if i < n_wake - 1 {
            let (psi_x, psi_y) =
                compute_psi_gradient_at_wake_point(geom, n + i, xi, yi, &gam, cosa, sina);
            let grad_mag = (psi_x * psi_x + psi_y * psi_y).sqrt().max(1.0e-12);
            nx_next = -psi_x / grad_mag;
            ny_next = -psi_y / grad_mag;
            wake_nx[i + 1] = nx_next;
            wake_ny[i + 1] = ny_next;
        }
    }

    let wake_s = snew.clone();
    let mut wake_qinvu_0 = vec![0.0; n_wake];
    let mut wake_qinvu_90 = vec![0.0; n_wake];
    wake_qinvu_0[0] = surface_qinvu_0[n - 1];
    wake_qinvu_90[0] = surface_qinvu_90[n - 1];
    for i in 1..n_wake {
        let (psi_x_0, psi_y_0) = compute_psi_gradient_at_wake_point(
            geom,
            n + i,
            wake_x[i],
            wake_y[i],
            &factorized.gamu_0,
            1.0,
            0.0,
        );
        let (psi_x_90, psi_y_90) = compute_psi_gradient_at_wake_point(
            geom,
            n + i,
            wake_x[i],
            wake_y[i],
            &factorized.gamu_90,
            0.0,
            1.0,
        );
        wake_qinvu_0[i] = psi_x_0 * wake_nx[i] + psi_y_0 * wake_ny[i];
        wake_qinvu_90[i] = psi_x_90 * wake_nx[i] + psi_y_90 * wake_ny[i];
    }
    let wake_ue = wake_qinvu_0
        .iter()
        .zip(wake_qinvu_90.iter())
        .map(|(&q0, &q90)| (cosa * q0 + sina * q90).abs())
        .collect();

    (wake_x, wake_y, wake_s, wake_ue)
}

fn compute_psi_gradient_at_wake_point(
    geom: &rustfoil_inviscid::geometry::AirfoilGeometry,
    point_index: usize,
    x: f64,
    y: f64,
    gamma: &[f64],
    cosa: f64,
    sina: f64,
) -> (f64, f64) {
    let result_x = psilin_with_dqdm(geom, point_index, x, y, 1.0, 0.0);
    let result_y = psilin_with_dqdm(geom, point_index, x, y, 0.0, 1.0);

    let mut psi_x = -sina;
    let mut psi_y = cosa;
    for (&dqdg_x, &g) in result_x.dqdg.iter().zip(gamma.iter()) {
        psi_x += dqdg_x * g;
    }
    for (&dqdg_y, &g) in result_y.dqdg.iter().zip(gamma.iter()) {
        psi_y += dqdg_y * g;
    }

    (psi_x, psi_y)
}

fn compute_wake_arc_lengths(te_x: f64, te_y: f64, wake_x: &[f64], wake_y: &[f64]) -> Vec<f64> {
    let n = wake_x.len().min(wake_y.len());
    let mut wake_s = Vec::with_capacity(n);
    let mut prev_x = te_x;
    let mut prev_y = te_y;
    let mut s = 0.0;

    for i in 0..n {
        let dx = wake_x[i] - prev_x;
        let dy = wake_y[i] - prev_y;
        s += (dx * dx + dy * dy).sqrt();
        wake_s.push(s);
        prev_x = wake_x[i];
        prev_y = wake_y[i];
    }

    wake_s
}

/// Compute the number of wake panels from DIJ matrix dimensions.
///
/// # Arguments
/// * `dij_size` - Size of the DIJ matrix (rows or columns)
/// * `n_airfoil_panels` - Number of airfoil panels
///
/// # Returns
/// Number of wake panels
pub fn compute_n_wake_panels(dij_size: usize, n_airfoil_panels: usize) -> usize {
    dij_size.saturating_sub(n_airfoil_panels)
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

    #[test]
    fn test_extract_surface_xfoil_uses_virtual_stagnation_anchor() {
        let x = vec![1.0, 0.5, 0.0, 0.5, 1.0];
        let y = vec![-0.05, -0.08, 0.0, 0.08, 0.05];
        let s = compute_arc_lengths(&x, &y);
        let gamma = vec![0.6, 0.2, -0.1, -0.4, -0.7];
        let (ist, sst, ue_stag) = interpolate_stagnation(&gamma, &s);

        let (upper_arc, _, _, upper_ue) =
            extract_surface_xfoil(ist, sst, ue_stag, &s, &x, &y, &gamma, true);
        let (lower_arc, _, _, lower_ue) =
            extract_surface_xfoil(ist, sst, ue_stag, &s, &x, &y, &gamma, false);

        assert_eq!(upper_arc[0], 0.0);
        assert_eq!(lower_arc[0], 0.0);
        assert!(upper_ue[0].abs() <= 1.0e-12);
        assert!(lower_ue[0].abs() <= 1.0e-12);
        assert!(upper_arc[1] > 0.0);
        assert!(lower_arc[1] > 0.0);
    }

    #[test]
    fn test_initialize_surface_stations_with_panel_idx_virtual_anchor_zero_velocity() {
        let arc = vec![0.0, 0.01, 0.04, 0.08];
        let ue = vec![0.0, 0.2, 0.4, 0.6];
        let panel_x = vec![0.0, 0.05, 0.2, 0.5];

        let stations = initialize_surface_stations_with_panel_idx(
            &arc,
            &ue,
            &panel_x,
            7,
            20,
            true,
            1.0e6,
        );

        assert_eq!(stations[0].x, 0.0);
        assert_eq!(stations[0].u, 0.0);
        assert_eq!(stations[0].mass_defect, 0.0);
        assert_eq!(stations[0].panel_idx, usize::MAX);
        assert!(stations[0].theta > 0.0);
        assert!(stations[1].u > 0.0);
        assert_eq!(stations[1].panel_idx, 7);
        assert_eq!(stations[2].panel_idx, 6);
    }

    #[test]
    fn test_initialize_canonical_state_appends_wake_on_lower_side_only() {
        let node_x = vec![1.0, 0.5, 0.0, 0.5, 1.0];
        let node_y = vec![0.05, 0.08, 0.0, -0.08, -0.05];
        let arc_lengths = compute_arc_lengths(&node_x, &node_y);
        let ue_inviscid = vec![0.6, 0.2, -0.1, -0.4, -0.7];
        let (ist, sst, _) = interpolate_stagnation(&ue_inviscid, &arc_lengths);
        let setup = ViscousSetup {
            arc_lengths: arc_lengths.clone(),
            ue_inviscid: ue_inviscid.clone(),
            dij: DMatrix::zeros(7, 7),
            panel_x: node_x.clone(),
            panel_y: node_y.clone(),
            stagnation_index: ist,
            n_upper: arc_lengths.len() - ist,
            n_lower: ist + 1,
            wake_x: vec![1.05, 1.25],
            wake_y: vec![-0.01, -0.02],
            wake_s: compute_wake_arc_lengths(1.0, 0.0, &[1.05, 1.25], &[-0.01, -0.02]),
            wake_ue: vec![0.8, 0.9],
        };

        let state = initialize_canonical_state(
            &setup,
            &node_x,
            &node_y,
            &vec![0.0; ue_inviscid.len()],
            ist,
            sst,
            1.0e6,
            0.0,
        );

        assert_eq!(state.ist, ist);
        assert_eq!(state.qinv, ue_inviscid);
        assert!(state.upper_rows.iter().all(|row| !row.is_wake));
        assert!(state.lower_rows.iter().any(|row| row.is_wake));
        assert_eq!(state.wake_x, setup.wake_x);
        assert_eq!(state.wake_y, setup.wake_y);
        assert_eq!(state.wake_s, setup.wake_s);
    }
}
