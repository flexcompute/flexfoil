//! Per-element paneling of a multi-element configuration.
//!
//! [`Configuration::panel_all`] turns a [`Configuration`] into a
//! [`PaneledConfiguration`]: every element re-paneled with its own
//! [`PanelingParams`](crate::spline::PanelingParams) and its own node count,
//! placed by its own [`Placement`](crate::placement::Placement), and the results
//! concatenated into the single global node array described by a [`Layout`].
//!
//! # One PANGEN call per element
//! [`CubicSpline::resample_xfoil`] is a port of XFOIL's PANGEN, which is written
//! for one contour. Three of its steps read the whole point set as that one
//! contour:
//!
//! - it locates a single leading edge, by Newton iteration over the whole array;
//! - it applies trailing-edge bunching at the two *ends* of the array, which for
//!   a concatenation of elements are the first element's lower trailing edge and
//!   the last element's upper trailing edge, and nothing in between;
//! - it normalises the curvature array by one global maximum, so all elements
//!   share a curvature scale. Whichever element has the tightest leading edge
//!   sets that scale for the rest — usually the smallest element, which leaves
//!   the *largest* element's own curvature peak the flattest of the three.
//!
//! Calling it once per element is what makes those three per-element: each call
//! sees exactly one contour, so it finds that element's leading edge, bunches at
//! that element's trailing edge, and normalises by that element's own peak
//! curvature. Nothing inside PANGEN is changed or reimplemented here; this
//! module is the loop around it.
//!
//! # Placement is applied after paneling
//! Each element is paneled in its **own** coordinates, and the resulting nodes
//! are then mapped into configuration coordinates by the element's placement.
//! The reasoning, since the two orders are not interchangeable in general:
//!
//! - PANGEN distributes nodes so that `(1 + C·κ·s_ref)·ds` is equal along the
//!   contour, with `s_ref` half the contour's own arc length and the curvature
//!   array normalised by its own maximum. Every factor is either dimensionless
//!   or a ratio of arc lengths on the same contour, so a rotation, a translation
//!   and a *uniform* scale all leave the arc-length distribution unchanged in
//!   exact arithmetic. Its Newton iteration, though, stops on an absolute test
//!   (`dmax < 1e-3` in arc-length units, as in XFOIL, which works in unit-chord
//!   coordinates), and that test is not scale-invariant: paneling a contour ten
//!   times unit chord stops relatively earlier than paneling the same shape at
//!   unit chord. So the two orders agree to roughly `1e-5` of the chord for an
//!   element larger than unit chord, and to machine precision for one smaller.
//!   `a_uniform_scale_changes_paneling_only_at_iteration_tolerance` measures it.
//! - What the order decides, then, is what the distribution is a *function* of.
//!   Paneling first makes it a function of the element's shape alone: a
//!   deflection or scale sweep re-places the same nodes instead of
//!   redistributing them, so a change in the answer across the sweep is not
//!   confounded with a change in paneling, and the iteration always runs in the
//!   element's own near-unit-chord coordinates where its absolute tolerance
//!   means what it does in XFOIL. Placing first would re-run the iteration on
//!   transformed coordinates and shift every node a little.
//! - The two orders would disagree outright for a non-uniform scale, which
//!   changes the curvature ratio between the leading edge and the rest of the
//!   contour and so would make the panel distribution depend on the placement.
//!   `Placement` has no non-uniform scale today; paneling first is the order
//!   that stays correct if one is ever added.
//! - With the identity placement the mapping is the identity, so a
//!   single-element configuration panels bit-for-bit identically to calling
//!   [`CubicSpline::resample_xfoil`] on the body directly. That is the
//!   backward-compatibility guarantee, pinned by
//!   `single_element_is_bit_identical_to_resample_xfoil`.
//!
//! Leading-edge indices are located in element coordinates for the same reason:
//! the minimum-x node of a flap deflected 30° down is not its leading edge.
//!
//! # Scope
//! Geometry only. Nothing here solves anything, and no existing solve path calls
//! it; the current solve path panels a single body. See the development phases in
//! `README.md` for status.

use crate::body::Body;
use crate::configuration::{Configuration, Element};
use crate::error::GeometryError;
use crate::layout::Layout;
use crate::point::Point;
use crate::spline::CubicSpline;

/// Fewest points an element may be paneled with.
///
/// A closed contour of `n` points has `n - 1` distinct nodes, so four points is
/// the smallest request that leaves three nodes — the fewest that enclose an
/// area. Below three, `resample_xfoil` also abandons the PANGEN distribution and
/// falls back to uniform spacing, which is not what a caller asking for
/// curvature-based paneling wants.
pub const MIN_POINTS_PER_ELEMENT: usize = 4;

/// XFOIL's default panel-node count (`PPAR`'s `N`), used per element by
/// [`PanelCounts::default`].
pub const DEFAULT_POINTS_PER_ELEMENT: usize = 160;

/// How many points each element is paneled with.
///
/// # Points, nodes and panels
/// A count here is the number of points handed to
/// [`CubicSpline::resample_xfoil`], which is XFOIL's `N`: the first and last
/// point both sit at the trailing edge. So an element asked for `n` points ends
/// up with
///
/// - `n - 1` nodes and `n - 1` panels if its contour is closed (sharp trailing
///   edge, where one node serves both surfaces);
/// - `n` nodes and `n - 1` panels if it is open (blunt trailing edge).
///
/// [`Layout`] counts nodes, so `layout.span(k).len` is the node count, not the
/// point count requested here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PanelCounts {
    /// The same number of points for every element.
    ///
    /// For a single-element configuration, `Each(160)` is XFOIL's default.
    Each(usize),

    /// One count per element, in configuration order.
    ///
    /// Must have exactly as many entries as the configuration has elements.
    PerElement(Vec<usize>),

    /// A total point budget shared out between the elements.
    ///
    /// Every element first gets `min_per_element` points (or
    /// [`MIN_POINTS_PER_ELEMENT`], whichever is larger); what is left over is
    /// then split in proportion to each element's **placed** arc length, that
    /// is, its contour arc length times its placement scale. Fractional shares
    /// are settled by largest remainder, ties going to the lowest element index,
    /// so the split is deterministic and the parts sum to `total` exactly.
    ///
    /// Sharing by arc length rather than by chord means a small element gets
    /// panels in proportion to how much surface it actually has, so every
    /// element ends up at a comparable number of panels per unit arc.
    Shared {
        /// Total points across all elements.
        total: usize,
        /// Floor for each element, before the remainder is shared out.
        min_per_element: usize,
    },
}

impl Default for PanelCounts {
    /// [`DEFAULT_POINTS_PER_ELEMENT`] points for every element — XFOIL's default
    /// applied to each element in turn.
    fn default() -> Self {
        Self::Each(DEFAULT_POINTS_PER_ELEMENT)
    }
}

impl PanelCounts {
    /// Work out the per-element point counts for one configuration.
    ///
    /// The result always has one entry per element, in configuration order, and
    /// every entry is at least [`MIN_POINTS_PER_ELEMENT`]. An empty
    /// configuration resolves to an empty vector whatever the variant.
    ///
    /// # Errors
    /// - `InvalidParameter { name: "element_point_count" }` if a requested count
    ///   is below [`MIN_POINTS_PER_ELEMENT`].
    /// - `InvalidParameter { name: "panel_counts_len" }` if a
    ///   [`PerElement`](Self::PerElement) list does not have one entry per
    ///   element. The reported value is the length that was supplied.
    /// - `InvalidParameter { name: "panel_counts_total" }` if a
    ///   [`Shared`](Self::Shared) budget cannot even cover the per-element floor.
    pub fn resolve(&self, config: &Configuration) -> Result<Vec<usize>, GeometryError> {
        let n_elements = config.len();
        match self {
            Self::Each(count) => {
                check_point_count(*count)?;
                Ok(vec![*count; n_elements])
            }
            Self::PerElement(counts) => {
                if counts.len() != n_elements {
                    return Err(GeometryError::InvalidParameter {
                        name: "panel_counts_len",
                        value: counts.len() as f64,
                    });
                }
                for &count in counts {
                    check_point_count(count)?;
                }
                Ok(counts.clone())
            }
            Self::Shared {
                total,
                min_per_element,
            } => share_points(config, *total, *min_per_element),
        }
    }
}

/// One point count, checked against [`MIN_POINTS_PER_ELEMENT`].
fn check_point_count(count: usize) -> Result<(), GeometryError> {
    if count < MIN_POINTS_PER_ELEMENT {
        return Err(GeometryError::InvalidParameter {
            name: "element_point_count",
            value: count as f64,
        });
    }
    Ok(())
}

/// The element's contour arc length in configuration coordinates.
///
/// A placement's rotation and translation do not change a length; its scale
/// does. `abs` because a negative scale is a reflection, which has the same
/// arc length as its mirror image.
fn placed_arc_length(element: &Element) -> f64 {
    element.body.arc_length() * element.placement.scale.abs()
}

/// Share a total point budget out between the elements — see
/// [`PanelCounts::Shared`].
fn share_points(
    config: &Configuration,
    total: usize,
    min_per_element: usize,
) -> Result<Vec<usize>, GeometryError> {
    let n_elements = config.len();
    if n_elements == 0 {
        return Ok(Vec::new());
    }

    let floor = min_per_element.max(MIN_POINTS_PER_ELEMENT);
    let required = floor.saturating_mul(n_elements);
    if total < required {
        return Err(GeometryError::InvalidParameter {
            name: "panel_counts_total",
            value: total as f64,
        });
    }

    let mut counts = vec![floor; n_elements];
    let spare = total - required;
    if spare == 0 {
        return Ok(counts);
    }

    let weights: Vec<f64> = config.iter().map(placed_arc_length).collect();
    let weight_sum: f64 = weights.iter().sum();

    // No usable arc lengths to weight by (an unmeasurable or non-finite
    // contour): fall back to an even split, remainder to the lowest indices.
    if !weight_sum.is_finite() || weight_sum <= 0.0 {
        let each = spare / n_elements;
        let remainder = spare % n_elements;
        for (i, count) in counts.iter_mut().enumerate() {
            *count += each + usize::from(i < remainder);
        }
        return Ok(counts);
    }

    // Whole shares first, then the largest remainders, so the parts sum to
    // `spare` exactly rather than to whatever the rounding happens to give.
    let mut fractions: Vec<(usize, f64)> = Vec::with_capacity(n_elements);
    let mut allocated = 0usize;
    for (i, &weight) in weights.iter().enumerate() {
        let exact = spare as f64 * weight / weight_sum;
        let whole = exact.floor();
        // Non-finite or negative would cast to 0, which the largest-remainder
        // pass below then tops up; it cannot over-allocate.
        let whole = whole as usize;
        counts[i] += whole;
        allocated += whole;
        fractions.push((i, exact - whole as f64));
    }

    let mut remainder = spare.saturating_sub(allocated);
    fractions.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    for &(i, _) in &fractions {
        if remainder == 0 {
            break;
        }
        counts[i] += 1;
        remainder -= 1;
    }

    Ok(counts)
}

/// A paneled configuration: one global node array plus the [`Layout`] that says
/// which nodes belong to which element.
///
/// `nodes` is in **configuration coordinates** — each element's placement has
/// been applied — and holds exactly [`Layout::total_nodes`] entries, so a global
/// node index from the layout indexes `nodes` directly.
///
/// # Nodes, not contour points
/// An element's slice holds its *distinct* nodes. A closed contour's duplicated
/// closing point is not among them: the element's last node is the one before
/// the trailing edge is reached again, and [`Layout::next_node`] closes the
/// contour by wrapping to the element's first node. So a slice is not a
/// ready-to-plot polyline for a sharp-trailing-edge element — repeat the first
/// node to close it. [`element_is_closed`](Self::element_is_closed) says which
/// case an element is in.
#[derive(Debug, Clone, PartialEq)]
pub struct PaneledConfiguration {
    /// Every element's nodes, concatenated in configuration order, in
    /// configuration coordinates.
    pub nodes: Vec<Point>,

    /// The node numbering of `nodes`, one span per element.
    pub layout: Layout,
}

impl PaneledConfiguration {
    /// Number of elements.
    #[inline]
    pub fn n_elements(&self) -> usize {
        self.layout.n_elements()
    }

    /// Total number of nodes — the length of `nodes`.
    #[inline]
    pub fn total_nodes(&self) -> usize {
        self.layout.total_nodes()
    }

    /// True if there are no elements, and therefore no nodes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.layout.n_elements() == 0
    }

    /// One element's nodes.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`, like indexing a slice.
    #[inline]
    pub fn element_nodes(&self, element: usize) -> &[Point] {
        let span = self.layout.span(element);
        &self.nodes[span.start..span.end()]
    }

    /// Each element's nodes in turn, in configuration order.
    pub fn contours(&self) -> impl Iterator<Item = &[Point]> + '_ {
        self.layout
            .spans()
            .iter()
            .map(move |span| &self.nodes[span.start..span.end()])
    }

    /// Whether one element's paneled contour closes on itself, i.e. whether its
    /// trailing edge is sharp.
    ///
    /// The same question as `!layout.span(element).has_blunt_te()`, which is
    /// meaningful here because [`Configuration::panel_all`] always fills the
    /// landmark indices in.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_is_closed(&self, element: usize) -> bool {
        !self.layout.span(element).has_blunt_te()
    }
}

impl Element {
    /// The element's contour as points, in its own coordinates.
    ///
    /// This is the point set the body was built from: one point per panel start,
    /// plus the last panel's end point.
    fn contour_points(&self) -> Vec<Point> {
        let panels = self.body.panels();
        let mut points = Vec::with_capacity(panels.len() + 1);
        points.extend(panels.iter().map(|panel| panel.p1));
        if let Some(last) = panels.last() {
            points.push(last.p2);
        }
        points
    }

    /// Re-panel the element's contour with `n_points` points, in the element's
    /// own coordinates.
    ///
    /// One [`CubicSpline::resample_xfoil`] call on this element alone, with this
    /// element's [`PanelingParams`](crate::spline::PanelingParams).
    fn repanel_unplaced(&self, n_points: usize) -> Result<Vec<Point>, GeometryError> {
        check_point_count(n_points)?;
        let spline = CubicSpline::from_points(&self.contour_points())?;
        Ok(spline.resample_xfoil(n_points, &self.paneling))
    }
}

impl Configuration {
    /// Panel every element independently and concatenate the results.
    ///
    /// Each element is paneled by one [`CubicSpline::resample_xfoil`] call on its
    /// own contour, with its own
    /// [`PanelingParams`](crate::spline::PanelingParams) and its own point count
    /// from `counts`, and the resulting nodes are then mapped into configuration
    /// coordinates by the element's placement. See the module documentation for
    /// why that is one call per element, and why the placement is applied
    /// afterwards.
    ///
    /// An empty configuration panels to an empty result rather than an error.
    ///
    /// # Errors
    /// - Whatever [`PanelCounts::resolve`] reports for an unusable point count.
    /// - `SplineInterpolationFailed` if an element's contour cannot be splined.
    /// - `DegeneratePanel` if an element's re-paneled contour has a zero-length
    ///   panel, which needs two coincident nodes.
    /// - `InvalidParameter { name: "element_span_le" }` if a re-paneled contour's
    ///   minimum-x point is its duplicated closing point, which a well-formed
    ///   airfoil does not have.
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::naca::naca4;
    /// use rustfoil_core::paneling::PanelCounts;
    /// use rustfoil_core::point::{point, Point};
    /// use rustfoil_core::{Body, Configuration, Element, Placement};
    ///
    /// // A NACA 2412 scaled to a main element and a 30% flap.
    /// let scaled = |chord: f64| -> Vec<Point> {
    ///     naca4(2412, Some(60))
    ///         .iter()
    ///         .map(|p| point(p.x * chord, p.y * chord))
    ///         .collect()
    /// };
    ///
    /// let main = Element::from_body(Body::from_points("main", &scaled(1.0)).unwrap());
    /// let mut flap = Element::from_body(Body::from_points("flap", &scaled(0.3)).unwrap());
    /// flap.placement = Placement::from_translation(0.95, -0.04);
    ///
    /// let config = Configuration::new(vec![main, flap]);
    /// let paneled = config.panel_all(&PanelCounts::Each(80)).unwrap();
    ///
    /// // The four-digit series has a finite trailing-edge thickness, so both
    /// // contours are open and all 80 points of each are distinct nodes.
    /// assert_eq!(paneled.total_nodes(), 160);
    /// assert_eq!(paneled.element_nodes(1).len(), 80);
    /// assert!(!paneled.element_is_closed(1));
    ///
    /// // Each element closes onto itself, never onto its neighbour.
    /// assert_eq!(paneled.layout.next_node(79), 0);
    /// assert_eq!(paneled.layout.next_node(159), 80);
    /// ```
    pub fn panel_all(&self, counts: &PanelCounts) -> Result<PaneledConfiguration, GeometryError> {
        let point_counts = counts.resolve(self)?;

        // Re-panel each element in its own coordinates, and keep the re-paneled
        // contour as a Body so the layout's node counts and landmark indices
        // come from `Layout::from_configuration` — one definition of those
        // conventions rather than a second copy here. The placements are left
        // off these bodies because the layout does not depend on them, and
        // because the leading edge has to be located before any rotation.
        let mut repaneled = Vec::with_capacity(self.len());
        let mut contours = Vec::with_capacity(self.len());
        for (element, &n_points) in self.iter().zip(&point_counts) {
            let points = element.repanel_unplaced(n_points)?;
            let body = Body::from_points(&element.body.name, &points)?;
            repaneled.push(Element::from_body(body));
            contours.push(points);
        }

        let layout = Layout::from_configuration(&Configuration::new(repaneled))?;

        // `Body::n_nodes` is the contour's point count, less one where the
        // contour closes, so a span is never longer than its contour.
        let mut nodes = Vec::with_capacity(layout.total_nodes());
        for (element, (contour, span)) in self.iter().zip(contours.iter().zip(layout.spans())) {
            let placement = &element.placement;
            if placement.is_identity() {
                // Exactly what `apply` would return, without the arithmetic —
                // this is the path a single-element configuration takes, and it
                // has to stay bit-for-bit equal to `resample_xfoil`.
                nodes.extend_from_slice(&contour[..span.len]);
            } else {
                nodes.extend(contour[..span.len].iter().map(|&p| placement.apply(p)));
            }
        }

        Ok(PaneledConfiguration { nodes, layout })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::naca::naca4;
    use crate::placement::Placement;
    use crate::point::{point, vec2};
    use crate::spline::PanelingParams;
    use approx::assert_relative_eq;

    /// A NACA 4412 contour scaled to `chord`, about the origin.
    ///
    /// The four-digit thickness distribution has a finite trailing-edge
    /// thickness, so this contour is **open**: `n` requested points give `n`
    /// nodes. [`sharp_foil`] is the closed counterpart.
    fn foil(chord: f64) -> Vec<Point> {
        naca4(4412, Some(80))
            .iter()
            .map(|p| point(p.x * chord, p.y * chord))
            .collect()
    }

    /// [`foil`] with its trailing edge closed, so the contour is closed and `n`
    /// requested points give `n - 1` nodes.
    fn sharp_foil(chord: f64) -> Vec<Point> {
        let mut points = foil(chord);
        let last = points.len() - 1;
        let te = point(chord, 0.0);
        points[0] = te;
        points[last] = te;
        points
    }

    fn element(name: &str, chord: f64, placement: Placement) -> Element {
        let body = Body::from_points(name, &foil(chord)).unwrap();
        Element::new(body, placement, PanelingParams::default(), name)
    }

    /// Slat, main and flap with chords 0.15, 1.0 and 0.3 — the ratio that makes
    /// a shared curvature scale a problem.
    fn slat_main_flap() -> Configuration {
        Configuration::new(vec![
            element("slat", 0.15, Placement::from_translation(-0.18, 0.02)),
            element("main", 1.0, Placement::identity()),
            element(
                "flap",
                0.3,
                Placement {
                    pivot: point(0.0, 0.0),
                    rotation_deg: -30.0,
                    translation: vec2(1.05, -0.06),
                    scale: 1.0,
                },
            ),
        ])
    }

    /// Panel one contour on its own, the way the single-element path does.
    fn alone(points: &[Point], n_points: usize, params: &PanelingParams) -> Vec<Point> {
        CubicSpline::from_points(points)
            .unwrap()
            .resample_xfoil(n_points, params)
    }

    // --- the backward-compatibility guarantee ----------------------------

    #[test]
    fn single_element_is_bit_identical_to_resample_xfoil() {
        // The whole of W1 rests on this: routing a single element through the
        // new path must not move a coordinate, or it moves a solver number.
        // Both trailing-edge cases, since they differ in how many of the
        // resampled points are distinct nodes.
        for points in [foil(1.0), sharp_foil(1.0)] {
            let closed = Body::from_points("main", &points).unwrap().is_closed();
            for n_points in [40usize, 61, 160, 201] {
                let direct = alone(&points, n_points, &PanelingParams::default());

                let config = Configuration::single(Body::from_points("main", &points).unwrap());
                let paneled = config.panel_all(&PanelCounts::Each(n_points)).unwrap();

                let expected = if closed { n_points - 1 } else { n_points };
                assert_eq!(paneled.total_nodes(), expected, "n = {n_points}");
                for (i, (&got, &want)) in paneled.nodes.iter().zip(direct.iter()).enumerate() {
                    assert_eq!(got.x, want.x, "n = {n_points}, node {i} x");
                    assert_eq!(got.y, want.y, "n = {n_points}, node {i} y");
                }
            }
        }
    }

    #[test]
    fn a_closed_contour_drops_only_its_closing_point() {
        // What `total_nodes == n_points - 1` leaves out: the closing point,
        // which is the first node again.
        let points = sharp_foil(1.0);
        let direct = alone(&points, 160, &PanelingParams::default());
        assert_relative_eq!(direct[159].x, direct[0].x, epsilon = 1e-12);
        assert_relative_eq!(direct[159].y, direct[0].y, epsilon = 1e-12);

        let config = Configuration::single(Body::from_points("main", &points).unwrap());
        let paneled = config.panel_all(&PanelCounts::Each(160)).unwrap();
        assert_eq!(paneled.total_nodes(), 159);
        assert!(paneled.element_is_closed(0));
        assert!(!paneled.layout.span(0).has_blunt_te());
    }

    #[test]
    fn a_blunt_trailing_edge_keeps_every_node() {
        // An open contour has no duplicated point to drop. Nodes run
        // TE(upper) → LE → TE(lower), so its lower trailing-edge node is the
        // last one and its upper trailing-edge node is node 0.
        let body = Body::from_points("blunt", &foil(1.0)).unwrap();
        assert!(!body.is_closed());

        let paneled = Configuration::single(body)
            .panel_all(&PanelCounts::Each(120))
            .unwrap();
        assert_eq!(paneled.total_nodes(), 120);
        assert!(!paneled.element_is_closed(0));
        assert_eq!(paneled.layout.span(0).te_upper, 0);
        assert_eq!(paneled.layout.span(0).te_lower, 119);
        assert!(paneled.layout.span(0).has_blunt_te());
    }

    #[test]
    fn default_counts_are_xfoils_default_per_element() {
        let config = Configuration::single(Body::from_points("main", &sharp_foil(1.0)).unwrap());
        let paneled = config.panel_all(&PanelCounts::default()).unwrap();
        assert_eq!(paneled.total_nodes(), DEFAULT_POINTS_PER_ELEMENT - 1);
    }

    // --- per-element independence ----------------------------------------

    #[test]
    fn every_element_is_paneled_as_if_it_were_alone() {
        // The point of the exercise: three elements of very different sizes,
        // each getting exactly the distribution it would get on its own. Equal
        // to the last bit, because the calls really are independent.
        let config = slat_main_flap();
        let counts = vec![70usize, 200, 100];
        let paneled = config
            .panel_all(&PanelCounts::PerElement(counts.clone()))
            .unwrap();

        for (k, element) in config.iter().enumerate() {
            let solo = alone(&element.contour_points(), counts[k], &element.paneling);
            let got = paneled.element_nodes(k);
            // Open contours, so every requested point is a distinct node.
            assert_eq!(got.len(), counts[k], "element {k} node count");
            for (&node, &want) in got.iter().zip(solo.iter()) {
                // The element's own nodes, mapped by its placement.
                let want = element.placement.apply(want);
                assert_eq!(node.x, want.x, "element {k}");
                assert_eq!(node.y, want.y, "element {k}");
            }
        }
    }

    #[test]
    fn a_small_element_is_not_starved_by_a_large_neighbour() {
        // What per-element paneling changes. Paneling the concatenated point set
        // with one PANGEN call shares one leading-edge search, one pair of
        // trailing-edge bunching points and one curvature maximum between all
        // three elements. Measured for a 0.15 / 1.0 / 0.3 chord configuration at
        // 480 nodes, that costs the small element node count and costs the large
        // element resolution:
        //
        //                       nodes        finest panel
        //   single contour   37 / 287 / 156   1.3e-3 / 4.9e-3 / 6.9e-4
        //   per element      77 / 288 / 115   5.8e-4 / 1.0e-3 / 7.6e-4
        //
        // Two separate effects, in opposite directions. The slat gets half the
        // nodes it should, because one arc-length distribution over the whole
        // concatenation does not know the slat is a separate body needing its own
        // leading-edge cluster. Meanwhile the *main* element is the one whose
        // leading edge is left five times coarser, because the shared curvature
        // maximum is set by the tightest leading edge anywhere in the
        // configuration — the slat's — and normalising by it flattens the main
        // element's own curvature peak. Both are the same root cause: one call
        // for three contours.
        let config = slat_main_flap();
        let total = 480usize;

        // Placed contours, in configuration order — what a single-contour
        // paneler would be handed.
        let placed: Vec<Vec<Point>> = config
            .iter()
            .map(|e| {
                e.contour_points()
                    .iter()
                    .map(|&p| e.placement.apply(p))
                    .collect()
            })
            .collect();
        let concatenated: Vec<Point> = placed.iter().flatten().copied().collect();
        let single_contour = alone(&concatenated, total, &PanelingParams::default());

        // Attribute each output node to the nearest input contour. The three
        // elements are well separated here, so this is unambiguous.
        let nearest = |p: Point| -> usize {
            placed
                .iter()
                .enumerate()
                .map(|(k, contour)| {
                    let d = contour
                        .iter()
                        .map(|&q| (q - p).norm())
                        .fold(f64::INFINITY, f64::min);
                    (k, d)
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|(k, _)| k)
                .unwrap()
        };
        let mut single_contour_by_element: Vec<Vec<Point>> = vec![Vec::new(); config.len()];
        for &p in &single_contour {
            single_contour_by_element[nearest(p)].push(p);
        }

        let paneled = config
            .panel_all(&PanelCounts::Shared {
                total,
                min_per_element: 40,
            })
            .unwrap();
        let per_element_slat = paneled.element_nodes(0);

        let min_spacing = |nodes: &[Point]| -> f64 {
            nodes
                .windows(2)
                .map(|w| (w[1] - w[0]).norm())
                .fold(f64::INFINITY, f64::min)
        };

        // The slat carries 0.15 / 1.45 of the configuration's chord and a
        // similar share of its arc length, so a proportional share of 480 nodes
        // is somewhere near 75. Per-element paneling gives it that; the
        // single-contour path gives it about half.
        let single_contour_slat = &single_contour_by_element[0];
        assert!(
            per_element_slat.len() >= 60,
            "per-element paneling gave the slat only {} of {total} nodes",
            per_element_slat.len()
        );
        assert!(
            per_element_slat.len() as f64 > 1.5 * single_contour_slat.len() as f64,
            "expected per-element paneling to give the slat more nodes than the \
             single-contour path: {} vs {}",
            per_element_slat.len(),
            single_contour_slat.len()
        );

        // And the nodes it gets are distributed for its own size, so its finest
        // panel is finer despite the slat being the smallest element.
        let per_element_min = min_spacing(per_element_slat);
        let single_contour_min = min_spacing(single_contour_slat);
        assert!(
            per_element_min < 0.75 * single_contour_min,
            "expected a finer minimum panel on the slat from per-element \
             paneling: {per_element_min:.3e} vs {single_contour_min:.3e}"
        );

        // The other half of the effect, and the larger one: the main element's
        // leading edge. Under a shared curvature maximum it is normalised
        // against the slat's tighter leading edge and comes out several times
        // coarser than the main element on its own would be.
        let main_per_element = min_spacing(paneled.element_nodes(1));
        let main_single_contour = min_spacing(&single_contour_by_element[1]);
        assert!(
            main_per_element < 0.5 * main_single_contour,
            "expected a finer minimum panel on the main element from \
             per-element paneling: {main_per_element:.3e} vs \
             {main_single_contour:.3e}"
        );
    }

    #[test]
    fn each_element_uses_its_own_paneling_params() {
        // Two elements of the same shape, one with curvature bunching and one
        // without, must come out differently — the params are per element, not
        // per configuration.
        let mut bunched = element("bunched", 1.0, Placement::identity());
        bunched.paneling = PanelingParams::default();
        let mut uniform = element("uniform", 1.0, Placement::identity());
        uniform.paneling = PanelingParams::uniform();

        let config = Configuration::new(vec![bunched, uniform]);
        let paneled = config.panel_all(&PanelCounts::Each(120)).unwrap();

        let spread = |nodes: &[Point]| -> f64 {
            let lengths: Vec<f64> = nodes.windows(2).map(|w| (w[1] - w[0]).norm()).collect();
            let max = lengths.iter().copied().fold(0.0_f64, f64::max);
            let min = lengths.iter().copied().fold(f64::INFINITY, f64::min);
            max / min
        };

        // Curvature bunching spreads panel lengths much more than uniform
        // spacing does.
        assert!(
            spread(paneled.element_nodes(0)) > 5.0 * spread(paneled.element_nodes(1)),
            "bunched spread {:.2} vs uniform spread {:.2}",
            spread(paneled.element_nodes(0)),
            spread(paneled.element_nodes(1))
        );
    }

    // --- placement -------------------------------------------------------

    #[test]
    fn placement_is_applied_to_the_paneled_nodes() {
        // Same shape, two placements: the nodes differ by the placement alone,
        // which is what "panel first, then place" means.
        let unplaced = Configuration::new(vec![element("flap", 0.3, Placement::identity())]);
        let placement = Placement {
            pivot: point(0.1, 0.0),
            rotation_deg: -25.0,
            translation: vec2(1.02, -0.05),
            scale: 1.0,
        };
        let placed = Configuration::new(vec![element("flap", 0.3, placement)]);

        let a = unplaced.panel_all(&PanelCounts::Each(90)).unwrap();
        let b = placed.panel_all(&PanelCounts::Each(90)).unwrap();
        assert_eq!(a.total_nodes(), b.total_nodes());

        for (&raw, &moved) in a.nodes.iter().zip(b.nodes.iter()) {
            let want = placement.apply(raw);
            assert_relative_eq!(moved.x, want.x, epsilon = 1e-15);
            assert_relative_eq!(moved.y, want.y, epsilon = 1e-15);
        }
    }

    #[test]
    fn paneling_is_independent_of_the_placement_scale() {
        // The property the "panel first, then place" order actually buys: the
        // paneling runs on the unplaced contour whatever the placement, so the
        // distribution is identical to the last bit across a scale sweep. Panel
        // lengths divided by the scale are the same numbers every time.
        let reference = Configuration::new(vec![element("flap", 1.0, Placement::identity())])
            .panel_all(&PanelCounts::Each(140))
            .unwrap();
        let reference_lengths: Vec<f64> = reference
            .nodes
            .windows(2)
            .map(|w| (w[1] - w[0]).norm())
            .collect();

        for scale in [0.15, 0.4, 2.5, 10.0] {
            let scaled = Configuration::new(vec![element(
                "flap",
                1.0,
                Placement {
                    pivot: point(0.3, 0.0),
                    scale,
                    ..Placement::identity()
                },
            )])
            .panel_all(&PanelCounts::Each(140))
            .unwrap();

            for (i, w) in scaled.nodes.windows(2).enumerate() {
                assert_relative_eq!(
                    (w[1] - w[0]).norm() / scale,
                    reference_lengths[i],
                    epsilon = 1e-15,
                    max_relative = 1e-14
                );
            }
        }
    }

    #[test]
    fn a_uniform_scale_changes_paneling_only_at_iteration_tolerance() {
        // The other half of the reasoning, measured rather than asserted.
        // PANGEN's spacing functional is scale-invariant in exact arithmetic, so
        // scaling a contour and then paneling it should give the same
        // distribution as paneling it and then scaling. It does — to the point
        // its Newton iteration is carried, which is an absolute test in
        // arc-length units and therefore stops relatively earlier on a larger
        // contour.
        for (scale, tolerance) in [(0.15, 1e-13), (0.4, 1e-13), (2.5, 1e-4), (10.0, 1e-4)] {
            let scale_first =
                Configuration::new(vec![element("scaled", scale, Placement::identity())])
                    .panel_all(&PanelCounts::Each(140))
                    .unwrap();
            let panel_first = Configuration::new(vec![element(
                "placed",
                1.0,
                Placement {
                    pivot: point(0.0, 0.0),
                    scale,
                    ..Placement::identity()
                },
            )])
            .panel_all(&PanelCounts::Each(140))
            .unwrap();

            assert_eq!(scale_first.total_nodes(), panel_first.total_nodes());
            let worst = scale_first
                .nodes
                .iter()
                .zip(panel_first.nodes.iter())
                .map(|(a, b)| (a - b).norm() / scale)
                .fold(0.0_f64, f64::max);
            assert!(
                worst < tolerance,
                "scale {scale}: worst node disagreement {worst:.3e} of the chord"
            );
        }
    }

    #[test]
    fn a_deflection_sweep_does_not_redistribute_nodes() {
        // Rotating an element re-places its nodes and nothing else, so a
        // deflection sweep varies one thing at a time.
        let reference = Configuration::new(vec![element("flap", 0.3, Placement::identity())])
            .panel_all(&PanelCounts::Each(100))
            .unwrap();

        for deflection in [-40.0, -20.0, 0.0, 15.0] {
            let deflected = Configuration::new(vec![element(
                "flap",
                0.3,
                Placement::rotation_about(point(0.05, 0.0), deflection),
            )])
            .panel_all(&PanelCounts::Each(100))
            .unwrap();

            // Panel lengths are unchanged: the distribution did not move.
            for (a, b) in reference
                .nodes
                .windows(2)
                .zip(deflected.nodes.windows(2))
            {
                assert_relative_eq!(
                    (b[1] - b[0]).norm(),
                    (a[1] - a[0]).norm(),
                    epsilon = 1e-12
                );
            }
        }
    }

    #[test]
    fn leading_edges_are_located_before_the_placement() {
        // A flap deflected 30° down: its minimum-x node in configuration
        // coordinates is not its leading edge, so the landmark has to come from
        // the element's own coordinates.
        let placement = Placement {
            pivot: point(0.0, 0.0),
            rotation_deg: -30.0,
            translation: vec2(1.05, -0.06),
            scale: 1.0,
        };
        let config = Configuration::new(vec![element("flap", 0.3, placement)]);
        let paneled = config.panel_all(&PanelCounts::Each(100)).unwrap();

        let span = paneled.layout.span(0);
        assert!(span.has_le());

        let nodes = paneled.element_nodes(0);
        let le = nodes[span.le];
        // The leading edge is the node furthest from the trailing edge, whatever
        // the deflection. Node 0 is the lower trailing edge.
        let te = nodes[0];
        let furthest = nodes
            .iter()
            .map(|&p| (p - te).norm())
            .fold(0.0_f64, f64::max);
        assert_relative_eq!((le - te).norm(), furthest, epsilon = 0.02 * 0.3);

        // Whereas the minimum-x node of the deflected contour is somewhere else
        // entirely, which is what makes locating it after the placement wrong.
        let min_x = nodes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert_ne!(min_x, span.le);
    }

    // --- layout agreement -------------------------------------------------

    #[test]
    fn the_layout_describes_the_returned_nodes() {
        let config = slat_main_flap();
        let paneled = config
            .panel_all(&PanelCounts::PerElement(vec![60, 180, 90]))
            .unwrap();

        assert_eq!(paneled.n_elements(), 3);
        assert_eq!(paneled.nodes.len(), paneled.total_nodes());
        // Open contours, so a node per requested point.
        assert_eq!(paneled.total_nodes(), 60 + 180 + 90);

        // The spans tile the node array, and each slice is its span's length.
        let mut expected_start = 0;
        for (k, contour) in paneled.contours().enumerate() {
            let span = paneled.layout.span(k);
            assert_eq!(span.start, expected_start);
            assert_eq!(contour.len(), span.len);
            assert_eq!(contour, paneled.element_nodes(k));
            expected_start = span.end();
        }
        assert_eq!(expected_start, paneled.total_nodes());

        // And no panel closes across the gap between two elements.
        for global in 0..paneled.total_nodes() {
            let owner = paneled.layout.element_of(global);
            assert_eq!(
                paneled.layout.element_of(paneled.layout.next_node(global)),
                owner
            );
        }
        assert_eq!(paneled.layout.next_node(59), 0);
        assert_eq!(paneled.layout.next_node(239), 60);
    }

    #[test]
    fn an_empty_configuration_panels_to_nothing() {
        let config = Configuration::new(vec![]);
        for counts in [
            PanelCounts::Each(160),
            PanelCounts::PerElement(vec![]),
            PanelCounts::Shared {
                total: 400,
                min_per_element: 40,
            },
        ] {
            let paneled = config.panel_all(&counts).unwrap();
            assert!(paneled.is_empty());
            assert_eq!(paneled.total_nodes(), 0);
            assert!(paneled.nodes.is_empty());
        }
    }

    // --- point counts -----------------------------------------------------

    #[test]
    fn shared_budget_follows_placed_arc_length() {
        let config = slat_main_flap();
        let counts = PanelCounts::Shared {
            total: 480,
            min_per_element: 40,
        }
        .resolve(&config)
        .unwrap();

        assert_eq!(counts.iter().sum::<usize>(), 480, "the parts must sum");
        assert!(counts.iter().all(|&c| c >= 40), "the floor must hold");

        // Ordered by arc length: main, then flap, then slat.
        assert!(counts[1] > counts[2] && counts[2] > counts[0], "{counts:?}");

        // The share above the floor tracks arc length. The slat is 0.15 chord
        // against the main's 1.0, so it takes roughly a seventh of the main's
        // spare share.
        let spare: Vec<f64> = counts.iter().map(|&c| (c - 40) as f64).collect();
        assert_relative_eq!(spare[0] / spare[1], 0.15, epsilon = 0.02);
        assert_relative_eq!(spare[2] / spare[1], 0.30, epsilon = 0.02);
    }

    #[test]
    fn shared_budget_scales_with_the_placement() {
        // Halving an element's placed size halves its share.
        let mut config = slat_main_flap();
        let before = PanelCounts::Shared {
            total: 480,
            min_per_element: 4,
        }
        .resolve(&config)
        .unwrap();

        config.elements[2].placement.scale = 0.5;
        let after = PanelCounts::Shared {
            total: 480,
            min_per_element: 4,
        }
        .resolve(&config)
        .unwrap();

        assert!(after[2] < before[2], "{before:?} -> {after:?}");
        assert_eq!(after.iter().sum::<usize>(), 480);
    }

    #[test]
    fn shared_budget_is_deterministic_for_equal_elements() {
        // Three identical elements and a budget that does not divide by three:
        // the remainder goes to the lowest indices, so the answer is stable
        // rather than dependent on float comparison order.
        let config = Configuration::new(vec![
            element("a", 1.0, Placement::identity()),
            element("b", 1.0, Placement::identity()),
            element("c", 1.0, Placement::identity()),
        ]);
        let counts = PanelCounts::Shared {
            total: 100,
            min_per_element: 4,
        }
        .resolve(&config)
        .unwrap();
        assert_eq!(counts, vec![34, 33, 33]);
    }

    #[test]
    fn shared_budget_rejects_a_total_below_the_floor() {
        let config = slat_main_flap();
        assert_eq!(
            PanelCounts::Shared {
                total: 100,
                min_per_element: 40,
            }
            .resolve(&config),
            Err(GeometryError::InvalidParameter {
                name: "panel_counts_total",
                value: 100.0,
            })
        );

        // The floor is never below the per-element minimum, even when asked.
        assert_eq!(
            PanelCounts::Shared {
                total: 11,
                min_per_element: 0,
            }
            .resolve(&config),
            Err(GeometryError::InvalidParameter {
                name: "panel_counts_total",
                value: 11.0,
            })
        );
        assert!(PanelCounts::Shared {
            total: 12,
            min_per_element: 0,
        }
        .resolve(&config)
        .is_ok());
    }

    #[test]
    fn per_element_counts_must_match_the_element_count() {
        let config = slat_main_flap();
        assert_eq!(
            config.panel_all(&PanelCounts::PerElement(vec![100, 100])),
            Err(GeometryError::InvalidParameter {
                name: "panel_counts_len",
                value: 2.0,
            })
        );
    }

    #[test]
    fn a_count_below_the_minimum_is_rejected() {
        let config = slat_main_flap();
        for counts in [
            PanelCounts::Each(MIN_POINTS_PER_ELEMENT - 1),
            PanelCounts::PerElement(vec![100, 3, 100]),
        ] {
            assert_eq!(
                config.panel_all(&counts),
                Err(GeometryError::InvalidParameter {
                    name: "element_point_count",
                    value: 3.0,
                }),
                "{counts:?}"
            );
        }
    }

    // --- the real three-element geometry ----------------------------------

    /// The McDonnell Douglas 30P-30N slat/main/flap fixture, as one element per
    /// `#`-headed block.
    fn mda_30p_30n() -> Vec<(String, Vec<Point>)> {
        let text = include_str!("../../../testdata/mda_30p_30n_trimmed.dat");
        let mut blocks: Vec<(String, Vec<Point>)> = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if let Some(header) = line.strip_prefix('#') {
                // Only a header immediately followed by coordinates starts a
                // block; the file's title lines are picked up and then replaced
                // by the next header before any point is read.
                match blocks.last() {
                    Some((_, points)) if points.is_empty() => {
                        blocks.pop();
                    }
                    _ => {}
                }
                blocks.push((header.trim().to_string(), Vec::new()));
                continue;
            }
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split_whitespace();
            let x: f64 = fields.next().unwrap().parse().unwrap();
            let y: f64 = fields.next().unwrap().parse().unwrap();
            blocks
                .last_mut()
                .expect("a coordinate line before any block header")
                .1
                .push(point(x, y));
        }
        blocks
    }

    #[test]
    fn the_real_three_element_fixture_panels_per_element() {
        let blocks = mda_30p_30n();
        let names: Vec<&str> = blocks.iter().map(|(name, _)| name.as_str()).collect();
        assert_eq!(names, ["Slat", "Main Element", "Flap"]);

        let config = Configuration::new(
            blocks
                .iter()
                .map(|(name, points)| {
                    Element::from_body(Body::from_points(name, points).unwrap())
                })
                .collect(),
        );

        // The elements are already positioned in the file, so the placements are
        // identities and the chords come straight from the coordinates.
        assert_eq!(config.main_element_index(), Some(1));

        let paneled = config
            .panel_all(&PanelCounts::Shared {
                total: 360,
                min_per_element: 60,
            })
            .unwrap();
        assert_eq!(paneled.n_elements(), 3);

        // Slat and main close on themselves in the file; the flap has a blunt
        // trailing edge, so it keeps every point as a node.
        assert!(paneled.element_is_closed(0));
        assert!(paneled.element_is_closed(1));
        assert!(!paneled.element_is_closed(2));

        // Each element's leading edge is its own minimum-x node, since nothing
        // here is deflected.
        for k in 0..3 {
            let nodes = paneled.element_nodes(k);
            let span = paneled.layout.span(k);
            let min_x = nodes
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap())
                .map(|(i, _)| i)
                .unwrap();
            assert_eq!(min_x, span.le, "element {k}");
        }

        // And the small slat still gets a share proportional to its surface
        // rather than being swamped by the main element.
        let slat_nodes = paneled.element_nodes(0).len();
        let main_nodes = paneled.element_nodes(1).len();
        assert!(slat_nodes >= 60, "slat got {slat_nodes} nodes");
        assert!(main_nodes > slat_nodes, "{main_nodes} vs {slat_nodes}");

        // Every element closes onto itself.
        for global in 0..paneled.total_nodes() {
            let next = paneled.layout.next_node(global);
            assert_eq!(
                paneled.layout.element_of(next),
                paneled.layout.element_of(global)
            );
        }
    }

    #[test]
    fn the_real_fixture_elements_are_paneled_independently() {
        // The same per-element equality as the synthetic case, on real
        // multi-element geometry.
        let blocks = mda_30p_30n();
        let config = Configuration::new(
            blocks
                .iter()
                .map(|(name, points)| {
                    Element::from_body(Body::from_points(name, points).unwrap())
                })
                .collect(),
        );

        let counts = vec![80usize, 200, 120];
        let paneled = config
            .panel_all(&PanelCounts::PerElement(counts.clone()))
            .unwrap();

        for (k, (_, points)) in blocks.iter().enumerate() {
            let solo = alone(points, counts[k], &PanelingParams::default());
            let got = paneled.element_nodes(k);
            for (&node, &want) in got.iter().zip(solo.iter()) {
                assert_eq!(node.x, want.x, "element {k}");
                assert_eq!(node.y, want.y, "element {k}");
            }
        }
    }
}
