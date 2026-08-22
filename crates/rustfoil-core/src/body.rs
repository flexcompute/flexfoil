//! Aerodynamic body representation.
//!
//! A `Body` represents a single aerodynamic element (airfoil, slat, flap, etc.)
//! as a closed contour discretized into panels.
//!
//! # Multi-Element Groundwork
//! `Body` is a single-element type: each instance owns its own contour, panels,
//! and per-body trailing- and leading-edge indices, so a multi-element
//! configuration can be held as a `Vec<Body>` without special-casing. That
//! representation is groundwork rather than a working multi-body solver —
//! paneling a configuration as a whole, inviscid interaction between elements,
//! and viscous treatment of the resulting wakes and gaps are not implemented,
//! and the current solve path operates on a single body. See the development
//! phases in `README.md` for status.
//!
//! # Panel Ordering Convention
//! Points run in Selig/XFOIL order, starting from the trailing edge:
//! 1. Upper surface (TE → LE)
//! 2. Lower surface (LE → TE)
//!
//! This is the order [`crate::naca::naca4`] emits and the order the surface
//! extraction in rustfoil-solver (`viscous::setup`) expects, so node 0 is the
//! upper-surface trailing edge and the last node is the lower-surface trailing
//! edge. [`crate::layout::ElementSpan`] follows the same convention.
//!
//! This convention ensures:
//! - Normal vectors point outward (into the flow)
//! - The trailing edge is easily identified as the first/last point
//! - Boundary layer marching proceeds naturally from stagnation point

use crate::error::GeometryError;
use crate::panel::Panel;
use crate::point::Point;

/// Tolerance for deciding whether a contour's first and last points meet.
///
/// A gap smaller than this in **either** coordinate counts as closed; see
/// [`contour_is_closed`] for why the comparison is componentwise.
///
/// # Why this value
/// `1e-10` is the tolerance the solve path already applies when it decides
/// whether to synthesize a blunt trailing-edge node
/// (`rustfoil-solver/src/viscous/setup.rs` and
/// `rustfoil-solver/src/inviscid/mod.rs`). Adopting this helper there therefore
/// classifies every geometry exactly as it is classified today, which matters
/// because the classification changes the node count and so the solver result.
///
/// It is deliberately looser than [`crate::point::GEOMETRY_TOLERANCE`] (`1e-12`),
/// which answers a different question: whether two points are the *same point*.
/// A contour whose ends are 1e-11 apart is closed for paneling purposes without
/// its endpoints being duplicates.
pub const CONTOUR_CLOSURE_TOLERANCE: f64 = 1e-10;

/// Whether a contour's first and last points meet, i.e. whether the contour is
/// closed.
///
/// This is the single canonical closure test. It answers the question that
/// decides how a contour is paneled: a **closed** contour (sharp trailing edge)
/// has a duplicate endpoint, so its last panel already runs back to its first
/// node; an **open** contour (blunt trailing edge) has a real gap across the
/// trailing-edge base.
///
/// # Comparison
/// Componentwise (L∞): closed if `|Δx| < tolerance` *and* `|Δy| < tolerance`,
/// with `tolerance` = [`CONTOUR_CLOSURE_TOLERANCE`]. Use
/// [`contour_is_closed_within`] to state a different tolerance explicitly.
///
/// Fewer than two points cannot form a contour and report as not closed.
///
/// # Example
/// ```
/// use rustfoil_core::body::contour_is_closed;
/// use rustfoil_core::point::point;
///
/// let sharp = [point(1.0, 0.0), point(0.0, 0.0), point(1.0, 0.0)];
/// assert!(contour_is_closed(&sharp));
///
/// let blunt = [point(1.0, -0.005), point(0.0, 0.0), point(1.0, 0.005)];
/// assert!(!contour_is_closed(&blunt));
/// ```
#[inline]
pub fn contour_is_closed(points: &[Point]) -> bool {
    contour_is_closed_within(points, CONTOUR_CLOSURE_TOLERANCE)
}

/// [`contour_is_closed`] with an explicit tolerance.
///
/// Provided so a caller that has to reproduce an existing classification can
/// state the tolerance it depends on at the call site instead of inheriting a
/// default that may be revised.
#[inline]
pub fn contour_is_closed_within(points: &[Point], tolerance: f64) -> bool {
    if points.len() < 2 {
        return false;
    }
    let first = &points[0];
    let last = &points[points.len() - 1];
    (last.x - first.x).abs() < tolerance && (last.y - first.y).abs() < tolerance
}

/// A single aerodynamic body discretized into panels.
///
/// # Invariants
/// - The body forms a closed contour (first and last points are coincident,
///   or close enough that they form the trailing edge)
/// - All panels have non-zero length
/// - Panels are ordered counter-clockwise
#[derive(Debug, Clone)]
pub struct Body {
    /// Human-readable identifier (e.g., "main", "slat", "flap")
    pub name: String,

    /// Ordered panels forming the closed contour
    panels: Vec<Panel>,

    /// Index of the trailing edge panel (upper surface side).
    ///
    /// For the Kutta condition, we need to know which panel is at the TE.
    /// This is typically the last panel (index = n_panels - 1), but may
    /// differ for blunt trailing edges or special geometries.
    te_panel_upper: usize,

    /// Index of the trailing edge panel (lower surface side).
    ///
    /// For sharp trailing edges, this is panel 0. For blunt TE, there may
    /// be a "base" panel connecting the two TE points.
    te_panel_lower: usize,

    /// Index of the leading edge point (approximate, for reference).
    ///
    /// The LE is identified as the point with minimum x-coordinate.
    le_point_idx: usize,

    /// Whether the input contour's first and last points met — see
    /// [`Body::is_closed`].
    is_closed: bool,
}

impl Body {
    /// Construct a body from raw coordinate points.
    ///
    /// # Arguments
    /// * `name` - Identifier for the body (e.g., "main", "slat")
    /// * `points` - Ordered points forming the airfoil contour. Should run in
    ///   Selig/XFOIL order: from the upper-surface trailing edge, around the
    ///   leading edge, back to the lower-surface trailing edge.
    ///
    /// # Point Closure
    /// The points should form a closed contour:
    /// - Either the first and last points are coincident (sharp TE)
    /// - Or they are distinct (blunt TE), in which case no closing panel is added
    ///
    /// # Errors
    /// - `InsufficientPoints` if fewer than 3 points are provided
    /// - `DegeneratePanel` if any resulting panel has zero length
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::body::Body;
    /// use rustfoil_core::point::point;
    ///
    /// // Simple diamond airfoil
    /// let points = vec![
    ///     point(1.0, 0.0),   // TE
    ///     point(0.5, -0.1),  // Lower
    ///     point(0.0, 0.0),   // LE
    ///     point(0.5, 0.1),   // Upper
    ///     point(1.0, 0.0),   // Back to TE (closed)
    /// ];
    /// let body = Body::from_points("diamond", &points).unwrap();
    /// assert_eq!(body.n_panels(), 4);
    /// ```
    pub fn from_points(name: &str, points: &[Point]) -> Result<Self, GeometryError> {
        const MIN_POINTS: usize = 3;

        if points.len() < MIN_POINTS {
            return Err(GeometryError::InsufficientPoints {
                required: MIN_POINTS,
                provided: points.len(),
            });
        }

        // Check if the contour is closed (first ≈ last point). Recorded on the
        // body and reported by `is_closed()`; it does not affect the panel
        // count (see below).
        let is_closed = contour_is_closed(points);

        // Number of panels, the same either way:
        // - Closed contour: n_points - 1, the last point being a duplicate of
        //   the first, so the last panel already runs back to node 0.
        // - Open contour: n_points - 1, with no base panel added across the
        //   trailing-edge gap.
        let n_panels = points.len() - 1;

        // Build panels
        let mut panels = Vec::with_capacity(n_panels);
        for i in 0..n_panels {
            let panel = Panel::new_with_index(points[i], points[i + 1], i)?;
            panels.push(panel);
        }

        // No closing panel is added for an open contour: the body is paneled as
        // given. Callers that need to know which case they have read
        // `is_closed()`.

        // Find trailing edge panels (first and last by convention)
        let te_panel_lower = 0;
        let te_panel_upper = panels.len() - 1;

        // Find leading edge (minimum x-coordinate)
        let le_point_idx = points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        Ok(Self {
            name: name.to_string(),
            panels,
            te_panel_upper,
            te_panel_lower,
            le_point_idx,
            is_closed,
        })
    }

    /// Whether the contour this body was built from was closed, i.e. whether
    /// its first and last points met within [`CONTOUR_CLOSURE_TOLERANCE`].
    ///
    /// - `true` — sharp trailing edge. The last panel runs back to node 0, so
    ///   the body has `n_panels()` distinct nodes.
    /// - `false` — blunt trailing edge. There is a gap between the last panel's
    ///   end and node 0, so the body has `n_panels() + 1` distinct nodes.
    ///
    /// # This does not change the paneling
    /// A body is paneled as given either way: `n_panels()` is
    /// `points.len() - 1` in both cases and no base panel is synthesized across
    /// a blunt trailing edge. The flag reports which case the caller has; it is
    /// the caller that decides what to do about it. The solve path, for
    /// instance, appends the missing lower trailing-edge node itself before
    /// building the influence matrix.
    #[inline]
    pub fn is_closed(&self) -> bool {
        self.is_closed
    }

    /// Number of distinct nodes in the body's contour.
    ///
    /// `n_panels()` for a closed contour, `n_panels() + 1` for an open one —
    /// the count the panel method works in, as distinct from the number of
    /// input points.
    #[inline]
    pub fn n_nodes(&self) -> usize {
        if self.is_closed {
            self.panels.len()
        } else {
            self.panels.len() + 1
        }
    }

    /// Number of panels in the body.
    #[inline]
    pub fn n_panels(&self) -> usize {
        self.panels.len()
    }

    /// Access the panels as a slice.
    #[inline]
    pub fn panels(&self) -> &[Panel] {
        &self.panels
    }

    /// Mutable access to panels (for repaneling operations).
    #[inline]
    pub fn panels_mut(&mut self) -> &mut [Panel] {
        &mut self.panels
    }

    /// Index of the upper-surface trailing edge panel.
    ///
    /// This is used for enforcing the Kutta condition, which requires
    /// the vorticity at the trailing edge to satisfy γ_upper + γ_lower = 0
    /// (for a sharp TE) to ensure finite velocity.
    #[inline]
    pub fn te_upper_index(&self) -> usize {
        self.te_panel_upper
    }

    /// Index of the lower-surface trailing edge panel.
    #[inline]
    pub fn te_lower_index(&self) -> usize {
        self.te_panel_lower
    }

    /// Index of the leading edge point.
    #[inline]
    pub fn le_index(&self) -> usize {
        self.le_point_idx
    }

    /// Total surface arc length of the body.
    pub fn arc_length(&self) -> f64 {
        self.panels.iter().map(|p| p.length()).sum()
    }

    /// Get all panel midpoints (control points) as a vector.
    ///
    /// Useful for setting up the influence coefficient matrix.
    pub fn control_points(&self) -> Vec<Point> {
        self.panels.iter().map(|p| p.midpoint()).collect()
    }

    /// Compute the approximate chord length.
    ///
    /// Defined as the distance from leading edge to trailing edge.
    pub fn chord(&self) -> f64 {
        if self.panels.is_empty() {
            return 0.0;
        }

        let le = &self.panels[self.le_point_idx.min(self.panels.len() - 1)].p1;
        let te = &self.panels[self.te_panel_lower].p1;

        (te - le).norm()
    }

    /// Iterator over panel indices for the lower surface (TE → LE).
    pub fn lower_surface_indices(&self) -> impl Iterator<Item = usize> {
        0..=self.le_point_idx.min(self.panels.len() - 1)
    }

    /// Iterator over panel indices for the upper surface (LE → TE).
    pub fn upper_surface_indices(&self) -> impl Iterator<Item = usize> {
        self.le_point_idx.min(self.panels.len() - 1)..self.panels.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;
    use approx::assert_relative_eq;

    fn make_diamond_airfoil() -> Body {
        // Simple closed diamond shape
        let points = vec![
            point(1.0, 0.0),  // TE
            point(0.5, -0.1), // Lower
            point(0.0, 0.0),  // LE
            point(0.5, 0.1),  // Upper
            point(1.0, 0.0),  // Back to TE
        ];
        Body::from_points("diamond", &points).unwrap()
    }

    #[test]
    fn test_body_construction() {
        let body = make_diamond_airfoil();

        assert_eq!(body.n_panels(), 4);
        assert_eq!(body.name, "diamond");
        assert_eq!(body.te_lower_index(), 0);
        assert_eq!(body.te_upper_index(), 3);
    }

    #[test]
    fn test_leading_edge_detection() {
        let body = make_diamond_airfoil();

        // LE should be at index 2 (the point at x=0)
        assert_eq!(body.le_index(), 2);
    }

    #[test]
    fn test_chord_length() {
        let body = make_diamond_airfoil();

        // Chord from LE (0,0) to TE (1,0) should be 1.0
        assert_relative_eq!(body.chord(), 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_insufficient_points() {
        let points = vec![point(0.0, 0.0), point(1.0, 0.0)];
        let result = Body::from_points("test", &points);

        assert!(matches!(
            result,
            Err(GeometryError::InsufficientPoints {
                required: 3,
                provided: 2
            })
        ));
    }

    #[test]
    fn test_arc_length() {
        // Square: 4 panels of length 1
        let points = vec![
            point(1.0, 0.0),
            point(0.0, 0.0),
            point(0.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 0.0),
        ];
        let body = Body::from_points("square", &points).unwrap();

        assert_relative_eq!(body.arc_length(), 4.0);
    }

    #[test]
    fn test_control_points() {
        let body = make_diamond_airfoil();
        let cps = body.control_points();

        assert_eq!(cps.len(), 4);
        // First panel midpoint should be between (1,0) and (0.5,-0.1)
        assert_relative_eq!(cps[0].x, 0.75);
        assert_relative_eq!(cps[0].y, -0.05);
    }

    // --- Closed-contour test ---------------------------------------------

    #[test]
    fn contour_closure_needs_two_points() {
        assert!(!contour_is_closed(&[]));
        assert!(!contour_is_closed(&[point(1.0, 0.0)]));
    }

    #[test]
    fn contour_closure_detects_a_sharp_trailing_edge() {
        let sharp = [
            point(1.0, 0.0),
            point(0.5, -0.05),
            point(0.0, 0.0),
            point(0.5, 0.05),
            point(1.0, 0.0),
        ];
        assert!(contour_is_closed(&sharp));
    }

    #[test]
    fn contour_closure_detects_a_blunt_trailing_edge() {
        let blunt = [
            point(1.0, -0.005),
            point(0.5, -0.05),
            point(0.0, 0.0),
            point(0.5, 0.05),
            point(1.0, 0.005),
        ];
        assert!(!contour_is_closed(&blunt));
    }

    #[test]
    fn contour_closure_is_componentwise_at_the_tolerance() {
        let inside = CONTOUR_CLOSURE_TOLERANCE * 0.5;
        let outside = CONTOUR_CLOSURE_TOLERANCE * 2.0;

        // Just inside in both components.
        let closed = [
            point(1.0, 0.0),
            point(0.0, 0.0),
            point(1.0 + inside, inside),
        ];
        assert!(contour_is_closed(&closed));

        // Outside in x alone, and in y alone.
        let open_x = [point(1.0, 0.0), point(0.0, 0.0), point(1.0 + outside, 0.0)];
        let open_y = [point(1.0, 0.0), point(0.0, 0.0), point(1.0, outside)];
        assert!(!contour_is_closed(&open_x));
        assert!(!contour_is_closed(&open_y));
    }

    #[test]
    fn contour_closure_honours_an_explicit_tolerance() {
        // A gap of 1e-11 is closed at the default tolerance (1e-10) and open at
        // the point-coincidence tolerance (1e-12) — the two answer different
        // questions, so the caller can state which one it means.
        let pts = [point(1.0, 0.0), point(0.0, 0.0), point(1.0 + 1e-11, 0.0)];
        assert!(contour_is_closed_within(&pts, CONTOUR_CLOSURE_TOLERANCE));
        assert!(!contour_is_closed_within(
            &pts,
            crate::point::GEOMETRY_TOLERANCE
        ));
    }

    #[test]
    fn body_reports_its_closure_state() {
        let sharp = make_diamond_airfoil();
        assert!(sharp.is_closed());

        let blunt = Body::from_points(
            "blunt",
            &[
                point(1.0, -0.005),
                point(0.5, -0.05),
                point(0.0, 0.0),
                point(0.5, 0.05),
                point(1.0, 0.005),
            ],
        )
        .unwrap();
        assert!(!blunt.is_closed());
    }

    #[test]
    fn closure_state_does_not_change_the_panel_count() {
        // Both contours have 5 points and must produce 4 panels: the closure
        // flag is reported, not acted on. If this ever changes, every solver
        // number downstream changes with it.
        let sharp = Body::from_points(
            "sharp",
            &[
                point(1.0, 0.0),
                point(0.5, -0.05),
                point(0.0, 0.0),
                point(0.5, 0.05),
                point(1.0, 0.0),
            ],
        )
        .unwrap();
        let blunt = Body::from_points(
            "blunt",
            &[
                point(1.0, -0.005),
                point(0.5, -0.05),
                point(0.0, 0.0),
                point(0.5, 0.05),
                point(1.0, 0.005),
            ],
        )
        .unwrap();

        assert_eq!(sharp.n_panels(), 4);
        assert_eq!(blunt.n_panels(), 4);
        assert_eq!(sharp.te_upper_index(), 3);
        assert_eq!(blunt.te_upper_index(), 3);
    }

    #[test]
    fn node_count_follows_the_closure_state() {
        let sharp = make_diamond_airfoil();
        // Closed: the last panel ends on node 0, so 4 panels are 4 nodes.
        assert_eq!(sharp.n_nodes(), 4);

        let blunt = Body::from_points(
            "blunt",
            &[
                point(1.0, -0.005),
                point(0.5, -0.05),
                point(0.0, 0.0),
                point(0.5, 0.05),
                point(1.0, 0.005),
            ],
        )
        .unwrap();
        // Open: the last panel's end is a fifth distinct node.
        assert_eq!(blunt.n_nodes(), 5);
    }
}
