//! Multi-element configurations: a set of placed elements plus the reference
//! quantities their coefficients are normalised by.
//!
//! An [`Element`] pairs one [`Body`] (a shape, in its own coordinates) with a
//! [`Placement`] (where it sits) and the [`PanelingParams`] it should be
//! paneled with. A [`Configuration`] is an ordered list of elements and the
//! reference chord and moment reference point used to non-dimensionalise the
//! configuration's total forces.
//!
//! # Scope
//! This module is geometry bookkeeping. It holds no flow state and performs no
//! aerodynamic calculation; the solve path is single-element today. See the
//! development phases in `README.md` for status.

use crate::body::Body;
use crate::placement::Placement;
use crate::point::{point, Point};
use crate::spline::PanelingParams;

/// One element of a configuration: a shape, where it sits, and how to panel it.
///
/// The element's `body` keeps the coordinates as imported. `placement` is
/// applied when the configuration is paneled and solved, so a gap, overlap or
/// deflection sweep varies the placement and leaves the contour alone.
#[derive(Debug, Clone)]
pub struct Element {
    /// The element's shape, in its own coordinates.
    pub body: Body,

    /// Where the element sits in configuration coordinates.
    pub placement: Placement,

    /// Panel distribution parameters for this element.
    ///
    /// Per-element, because a slat and the main element it sits ahead of do not
    /// want the same panel density.
    pub paneling: PanelingParams,

    /// Role of this element in the configuration — conventionally `"slat"`,
    /// `"main"` or `"flap"`.
    ///
    /// Distinct from `body.name`, which identifies the *source geometry* (an
    /// airfoil name, or a file's element label). Two configurations can reuse
    /// the same body under different roles.
    pub name: String,
}

impl Element {
    /// Build an element from its parts.
    pub fn new(body: Body, placement: Placement, paneling: PanelingParams, name: &str) -> Self {
        Self {
            body,
            placement,
            paneling,
            name: name.to_string(),
        }
    }

    /// Build an element that sits in its own coordinates, with default
    /// paneling, taking its role name from the body.
    pub fn from_body(body: Body) -> Self {
        let name = body.name.clone();
        Self {
            body,
            placement: Placement::identity(),
            paneling: PanelingParams::default(),
            name,
        }
    }

    /// The element's leading- and trailing-edge points, **in configuration
    /// coordinates** (that is, with the placement applied).
    ///
    /// `None` for a body with no panels.
    ///
    /// # The trailing edge is the trailing-edge midpoint
    /// The trailing-edge point is `0.5 * (first node + last node)` of the
    /// element's contour, which is the definition the solver uses:
    /// `xte = 0.5 * (x[0] + x[n-1])` in
    /// `rustfoil-inviscid`'s `AirfoilGeometry` (its `compute_te_geometry`,
    /// XFOIL's TECALC). A trailing-edge *corner* — either of the two nodes on
    /// its own — differs from the midpoint by half the trailing-edge gap
    /// whenever the trailing edge is blunt, and this length is the per-element
    /// coefficient normalisation and the input to
    /// [`main_element_index`](Configuration::main_element_index), so the two
    /// definitions have to be the same one. For a closed contour the first and
    /// last nodes coincide and the midpoint is that node exactly.
    ///
    /// # The leading edge is the minimum-x node
    /// Not the spline leading edge the solver's LEFIND locates, which is the
    /// authoritative one — see [`chord`](Self::chord) for how far apart the two
    /// are and why this one is still the definition here. The minimum-x node is
    /// taken in the body's **own** coordinates
    /// ([`Body::le_index`](crate::body::Body::le_index)) and then placed, so a
    /// flap whose deflection is carried in its [`Placement`] has its leading edge
    /// located before the rotation. A flap whose deflection is already baked into
    /// its coordinates — as it is in an imported multi-element coordinate file —
    /// does not get that, and its minimum-x node can be well away from its
    /// leading edge.
    pub fn chord_endpoints(&self) -> Option<(Point, Point)> {
        let panels = self.body.panels();
        if panels.is_empty() {
            return None;
        }
        let last = panels.len() - 1;
        let le = panels[self.body.le_index().min(last)].p1;
        // The two trailing-edge nodes are the ends of the contour: the first
        // panel's start and the last panel's end. Written in the solver's form
        // so the two agree bit for bit.
        let te_first = panels[0].p1;
        let te_last = panels[last].p2;
        let te = point(
            0.5 * (te_first.x + te_last.x),
            0.5 * (te_first.y + te_last.y),
        );
        Some((self.placement.apply(le), self.placement.apply(te)))
    }

    /// The element's chord length in configuration coordinates: a geometric
    /// length, not the chord reported coefficients are normalised by.
    ///
    /// The distance between the two points
    /// [`chord_endpoints`](Self::chord_endpoints) returns — the minimum-x node
    /// and the trailing-edge midpoint — so it is scaled by the placement's
    /// `scale` and unaffected by its rotation and translation. `0.0` for a body
    /// with no panels.
    ///
    /// This is *not* [`Body::chord`](crate::body::Body::chord), which measures
    /// to a trailing-edge corner instead of the midpoint. See
    /// [`chord_endpoints`](Self::chord_endpoints).
    ///
    /// # Which chord is authoritative
    /// `ElementGeometry::chord` in `rustfoil-inviscid` is: leading edge from
    /// LEFIND on the element's spline, trailing edge from TECALC, both of them
    /// the landmarks the solver uses for its own geometry. Per-element and
    /// configuration coefficients (decision D4) normalise by that one.
    ///
    /// This length differs from it, because a minimum-x node is not a spline
    /// leading edge. Measured on this repository's fixtures: 7.9e-5 relative on
    /// `naca2412.dat`, 2.6e-5 on `naca0012_xfoil_paneled.dat`, zero on the
    /// symmetric sections whose leading-edge node sits exactly at `x = 0`; and
    /// 3.3% on the slat and 0.8% on the flap of the full 30P-30N section, whose
    /// deflections are baked into their coordinates so that the minimum-x node is
    /// not near the leading edge at all.
    /// `the_two_chord_definitions_disagree_by_a_pinned_amount` in
    /// `rustfoil-inviscid`'s geometry module pins both.
    ///
    /// It is kept as it is, rather than moved onto the spline definition, for two
    /// reasons:
    ///
    /// - it has to answer for any [`Body`], including a coarse polygon or an
    ///   un-paneled import. A spline leading edge is not defined on a shape whose
    ///   own leading edge is a vertex the spline rounds off, and this crate holds
    ///   no inviscid geometry to ask;
    /// - it is what the D6 clearance thresholds are scaled by, through
    ///   [`Configuration::resolved_ref_chord`], so changing it moves published
    ///   clearance numbers. On the real 30P-30N the two candidate reference
    ///   chords are 0.831564097 here and 0.831574314 from the spline, a relative
    ///   difference of 1.2e-5 that moves the 0.5% floor from 4.157820e-3 to
    ///   4.157871e-3 against a smallest pair margin of 8.5e-3, so no verdict
    ///   changes — but that is a measurement on one section, not a licence.
    ///
    /// So: this for geometric bookkeeping (ranking elements, scaling clearance
    /// fractions), `ElementGeometry::chord` for anything aerodynamic.
    pub fn chord(&self) -> f64 {
        match self.chord_endpoints() {
            Some((le, te)) => (te - le).norm(),
            None => 0.0,
        }
    }

    /// The element's quarter-chord point in configuration coordinates.
    ///
    /// A quarter of the way from the leading edge to the trailing edge along
    /// the chord line. `None` for a body with no panels.
    pub fn quarter_chord(&self) -> Option<Point> {
        let (le, te) = self.chord_endpoints()?;
        Some(le + (te - le) * 0.25)
    }
}

/// A multi-element configuration and its reference quantities.
///
/// # Coefficient normalisation (MSES-style)
/// Two different reference lengths are in play, and which one applies depends on
/// what is being reported:
///
/// - **Per-element coefficients** normalise by that element's *own* chord. A
///   flap's `Cl` is the flap's load over the flap's chord, so it is comparable
///   with that flap run on its own.
/// - **Configuration totals** normalise by the configuration reference chord,
///   and take moments about the configuration reference point
///   ([`Configuration::resolved_ref_point`]). Totals summed from per-element
///   coefficients therefore have to be rescaled by each element's chord ratio
///   first.
///
/// This is the MSES convention. Mixing the two — summing per-element
/// coefficients directly into a total — gives a number that is not any
/// recognised coefficient.
///
/// An element's own chord there is `ElementGeometry::chord` in
/// `rustfoil-inviscid`, measured between the leading edge LEFIND locates and the
/// trailing-edge midpoint, and the configuration reference chord defaults to the
/// largest of them, `ConfigGeometry::default_ref_chord`. [`Element::chord`] and
/// [`Configuration::resolved_ref_chord`] are the geometry-time forms of the same
/// two quantities, from the contour's minimum-x node instead of the spline
/// leading edge; they rank elements and scale the clearance diagnostics, and
/// they are not what a reported coefficient divides by. [`Element::chord`]
/// records how far apart the two are.
///
/// # Default reference quantities (decision D4)
/// Both are `Option` and both have a documented default when unset:
///
/// - `ref_chord` defaults to the chord of the **largest-chord element**, which
///   in a slat/main/flap configuration is the main element. See
///   [`main_element_index`](Configuration::main_element_index).
/// - `ref_point` defaults to that same element's quarter-chord.
///
/// The defaults make a single-element configuration behave exactly like a plain
/// airfoil: its own chord and its own quarter-chord.
///
/// # Example
/// ```
/// use rustfoil_core::body::Body;
/// use rustfoil_core::configuration::{Configuration, Element};
/// use rustfoil_core::placement::Placement;
/// use rustfoil_core::point::point;
///
/// // Placeholder contours: a unit-chord main and a quarter-chord flap.
/// let contour = |chord: f64| vec![
///     point(chord, 0.0),
///     point(0.5 * chord, -0.05 * chord),
///     point(0.0, 0.0),
///     point(0.5 * chord, 0.05 * chord),
///     point(chord, 0.0),
/// ];
///
/// let main = Element::from_body(Body::from_points("main", &contour(1.0)).unwrap());
/// let mut flap = Element::from_body(Body::from_points("flap", &contour(0.25)).unwrap());
/// flap.placement = Placement::from_translation(0.95, -0.02);
///
/// let config = Configuration::new(vec![main, flap]);
///
/// assert_eq!(config.len(), 2);
/// // The main element is the largest-chord element, so it sets the reference.
/// assert_eq!(config.main_element_index(), Some(0));
/// assert!((config.resolved_ref_chord() - 1.0).abs() < 1e-12);
/// assert!((config.resolved_ref_point().x - 0.25).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct Configuration {
    /// The elements, in configuration order.
    ///
    /// Order is the caller's; nothing here requires slat-main-flap ordering.
    /// It is the order [`crate::layout::Layout`] concatenates nodes in, so
    /// changing it renumbers global node indices.
    pub elements: Vec<Element>,

    /// Reference chord for configuration totals. `None` uses the D4 default —
    /// see [`resolved_ref_chord`](Configuration::resolved_ref_chord).
    pub ref_chord: Option<f64>,

    /// Moment reference point for configuration totals. `None` uses the D4
    /// default — see [`resolved_ref_point`](Configuration::resolved_ref_point).
    pub ref_point: Option<Point>,
}

impl Configuration {
    /// Build a configuration from its elements, with default reference
    /// quantities.
    pub fn new(elements: Vec<Element>) -> Self {
        Self {
            elements,
            ref_chord: None,
            ref_point: None,
        }
    }

    /// Build a single-element configuration from one body, placed in its own
    /// coordinates.
    ///
    /// The reference quantities then resolve to the body's own chord and
    /// quarter-chord, so this is equivalent to treating the body as a plain
    /// airfoil.
    pub fn single(body: Body) -> Self {
        Self::new(vec![Element::from_body(body)])
    }

    /// Override the reference chord.
    #[must_use]
    pub fn with_ref_chord(mut self, ref_chord: f64) -> Self {
        self.ref_chord = Some(ref_chord);
        self
    }

    /// Override the moment reference point.
    #[must_use]
    pub fn with_ref_point(mut self, ref_point: Point) -> Self {
        self.ref_point = Some(ref_point);
        self
    }

    /// Number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// True if the configuration holds no elements.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Iterate over the elements in configuration order.
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, Element> {
        self.elements.iter()
    }

    /// Borrow one element by index, or `None` if out of range.
    #[inline]
    pub fn element(&self, index: usize) -> Option<&Element> {
        self.elements.get(index)
    }

    /// Index of the largest-chord element — the main element, by the D4 rule.
    ///
    /// Chords are compared in configuration coordinates, so an element scaled
    /// down by its placement is compared at its placed size. Ties go to the
    /// lowest index; `None` for an empty configuration.
    ///
    /// The rule is geometric rather than name-based on purpose: it gives the
    /// same answer for an imported multi-element file whose elements are
    /// unlabelled, or labelled in another language.
    pub fn main_element_index(&self) -> Option<usize> {
        self.elements
            .iter()
            .enumerate()
            .fold(None, |best, (i, element)| {
                let chord = element.chord();
                match best {
                    // Strictly greater, so the first of equal chords wins.
                    Some((_, best_chord)) if chord <= best_chord => best,
                    _ => Some((i, chord)),
                }
            })
            .map(|(i, _)| i)
    }

    /// The reference chord for configuration totals.
    ///
    /// `ref_chord` if set, otherwise the chord of the largest-chord element
    /// (D4). Falls back to `1.0` for an empty configuration, or for one whose
    /// main element has no measurable chord, so callers dividing by this never
    /// divide by zero.
    ///
    /// # What this length moves
    /// It carries [`Element::chord`]'s definition — the minimum-x node, not the
    /// spline leading edge — into everything that reads it, which today is:
    ///
    /// - [`main_element_index`](Self::main_element_index) and
    ///   [`resolved_ref_point`](Self::resolved_ref_point), which pick and measure
    ///   the same largest element;
    /// - [`clearance_with`](Self::clearance_with), where it sets the `floor` the
    ///   D6 hard minimum is taken as a fraction of, and divides every
    ///   `min_distance_fraction`, `gap_fraction` and `overlap_fraction` in the
    ///   report;
    /// - the summary written by [`crate::config_io`].
    ///
    /// So it is a threshold input for the clearance diagnostics, not only a
    /// reporting scale. [`Element::chord`] documents the size of the difference
    /// against the authoritative spline chord, and its effect on the 30P-30N
    /// verdicts.
    pub fn resolved_ref_chord(&self) -> f64 {
        if let Some(chord) = self.ref_chord {
            return chord;
        }
        let chord = self
            .main_element_index()
            .and_then(|i| self.elements.get(i))
            .map(Element::chord)
            .unwrap_or(0.0);

        if chord > 0.0 {
            chord
        } else {
            1.0
        }
    }

    /// The moment reference point for configuration totals.
    ///
    /// `ref_point` if set, otherwise the quarter-chord of the largest-chord
    /// element (D4). Falls back to the origin for an empty configuration.
    pub fn resolved_ref_point(&self) -> Point {
        if let Some(p) = self.ref_point {
            return p;
        }
        self.main_element_index()
            .and_then(|i| self.elements.get(i))
            .and_then(Element::quarter_chord)
            .unwrap_or_else(|| point(0.0, 0.0))
    }
}

impl<'a> IntoIterator for &'a Configuration {
    type Item = &'a Element;
    type IntoIter = core::slice::Iter<'a, Element>;

    fn into_iter(self) -> Self::IntoIter {
        self.elements.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::vec2;
    use approx::assert_relative_eq;

    /// A closed diamond of the given chord, leading edge at the origin.
    fn diamond(name: &str, chord: f64) -> Body {
        let pts = vec![
            point(chord, 0.0),
            point(0.5 * chord, -0.05 * chord),
            point(0.0, 0.0),
            point(0.5 * chord, 0.05 * chord),
            point(chord, 0.0),
        ];
        Body::from_points(name, &pts).unwrap()
    }

    fn slat_main_flap() -> Configuration {
        let mut slat = Element::from_body(diamond("slat", 0.15));
        slat.placement = Placement::from_translation(-0.12, 0.03);

        let main = Element::from_body(diamond("main", 1.0));

        let mut flap = Element::from_body(diamond("flap", 0.3));
        flap.placement = Placement {
            pivot: point(0.0, 0.0),
            rotation_deg: -35.0,
            translation: vec2(0.92, -0.02),
            scale: 1.0,
        };

        Configuration::new(vec![slat, main, flap])
    }

    #[test]
    fn element_from_body_is_unplaced_with_default_paneling() {
        let element = Element::from_body(diamond("main", 1.0));
        assert!(element.placement.is_identity());
        assert_eq!(element.name, "main");
        assert_relative_eq!(
            element.paneling.curv_param,
            PanelingParams::default().curv_param
        );
    }

    #[test]
    fn element_name_is_independent_of_body_name() {
        let element = Element::new(
            diamond("naca2412", 1.0),
            Placement::identity(),
            PanelingParams::default(),
            "main",
        );
        assert_eq!(element.name, "main");
        assert_eq!(element.body.name, "naca2412");
    }

    #[test]
    fn element_chord_is_measured_in_configuration_coordinates() {
        let mut element = Element::from_body(diamond("flap", 0.4));
        assert_relative_eq!(element.chord(), 0.4, epsilon = 1e-12);

        // Rotation and translation do not change a length.
        element.placement = Placement {
            pivot: point(0.1, 0.0),
            rotation_deg: -40.0,
            translation: vec2(0.9, -0.05),
            scale: 1.0,
        };
        assert_relative_eq!(element.chord(), 0.4, epsilon = 1e-12);

        // Scale does.
        element.placement.scale = 0.5;
        assert_relative_eq!(element.chord(), 0.2, epsilon = 1e-12);
    }

    /// An open contour of the given chord whose trailing edge is `gap` thick,
    /// split evenly about `y = 0`. Nodes run TE(upper) → LE → TE(lower).
    fn blunt(name: &str, chord: f64, gap: f64) -> Body {
        let pts = vec![
            point(chord, 0.5 * gap),
            point(0.5 * chord, 0.05 * chord),
            point(0.0, 0.0),
            point(0.5 * chord, -0.05 * chord),
            point(chord, -0.5 * gap),
        ];
        Body::from_points(name, &pts).unwrap()
    }

    #[test]
    fn element_chord_measures_to_the_trailing_edge_midpoint() {
        // The two trailing-edge nodes sit at y = ±0.05 on a unit chord, so a
        // corner is 0.05 off the chord line and the midpoint is on it. Measuring
        // to a corner would report sqrt(1 + 0.05^2), which is the error this
        // definition removes from every per-element coefficient.
        let element = Element::from_body(blunt("blunt", 1.0, 0.1));
        assert_relative_eq!(element.chord(), 1.0, epsilon = 1e-15);

        let (le, te) = element.chord_endpoints().unwrap();
        assert_relative_eq!(le.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(le.y, 0.0, epsilon = 1e-15);
        assert_relative_eq!(te.x, 1.0, epsilon = 1e-15);
        assert_relative_eq!(te.y, 0.0, epsilon = 1e-15);

        // Body::chord is the corner measurement and is deliberately left alone;
        // the two are different quantities.
        assert_relative_eq!(
            element.body.chord(),
            (1.0f64 + 0.05 * 0.05).sqrt(),
            epsilon = 1e-15
        );
    }

    #[test]
    fn element_chord_is_unchanged_for_a_closed_contour() {
        // A closed contour's two trailing-edge nodes are the same point, so the
        // midpoint is that point and the corrected definition agrees with the
        // old one exactly — which is why no single-element number moves.
        let element = Element::from_body(diamond("sharp", 1.0));
        let (_, te) = element.chord_endpoints().unwrap();
        let corner = element.body.panels()[0].p1;
        assert_eq!(te.x.to_bits(), corner.x.to_bits());
        assert_eq!(te.y.to_bits(), corner.y.to_bits());
        assert_eq!(element.chord().to_bits(), element.body.chord().to_bits());
    }

    #[test]
    fn trailing_edge_definition_can_decide_the_main_element() {
        // Two elements whose corner-measured chords rank one way and whose
        // midpoint-measured chords rank the other. The reference chord follows
        // the midpoint definition, so it is the flatter element that wins.
        let a = Element::from_body(blunt("a", 1.0, 0.6)); // corner: sqrt(1+0.09)
        let b = Element::from_body(blunt("b", 1.02, 0.0)); // corner: 1.02
        let config = Configuration::new(vec![a, b]);

        assert_eq!(config.main_element_index(), Some(1));
        assert_relative_eq!(config.resolved_ref_chord(), 1.02, epsilon = 1e-15);
    }

    #[test]
    fn element_quarter_chord_follows_the_placement() {
        let mut element = Element::from_body(diamond("main", 1.0));
        let unplaced = element.quarter_chord().unwrap();
        assert_relative_eq!(unplaced.x, 0.25, epsilon = 1e-12);
        assert_relative_eq!(unplaced.y, 0.0, epsilon = 1e-12);

        element.placement = Placement::from_translation(0.5, 0.1);
        let placed = element.quarter_chord().unwrap();
        assert_relative_eq!(placed.x, 0.75, epsilon = 1e-12);
        assert_relative_eq!(placed.y, 0.1, epsilon = 1e-12);
    }

    #[test]
    fn main_element_is_the_largest_chord() {
        let config = slat_main_flap();
        assert_eq!(config.main_element_index(), Some(1));
    }

    #[test]
    fn main_element_uses_placed_chords() {
        // Body chords say element 0 is largest; placement scale says element 1.
        let mut a = Element::from_body(diamond("a", 1.0));
        a.placement.scale = 0.4;
        let b = Element::from_body(diamond("b", 0.8));

        let config = Configuration::new(vec![a, b]);
        assert_eq!(config.main_element_index(), Some(1));
    }

    #[test]
    fn main_element_ties_go_to_the_lowest_index() {
        let config = Configuration::new(vec![
            Element::from_body(diamond("first", 1.0)),
            Element::from_body(diamond("second", 1.0)),
        ]);
        assert_eq!(config.main_element_index(), Some(0));
    }

    #[test]
    fn main_element_of_an_empty_configuration() {
        let config = Configuration::new(vec![]);
        assert_eq!(config.main_element_index(), None);
    }

    #[test]
    fn reference_quantities_default_to_the_main_element() {
        let config = slat_main_flap();
        assert_relative_eq!(config.resolved_ref_chord(), 1.0, epsilon = 1e-12);

        let ref_point = config.resolved_ref_point();
        assert_relative_eq!(ref_point.x, 0.25, epsilon = 1e-12);
        assert_relative_eq!(ref_point.y, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn reference_quantities_are_overridable() {
        let config = slat_main_flap()
            .with_ref_chord(1.37)
            .with_ref_point(point(0.4, -0.01));

        assert_relative_eq!(config.resolved_ref_chord(), 1.37);
        assert_relative_eq!(config.resolved_ref_point().x, 0.4);
        assert_relative_eq!(config.resolved_ref_point().y, -0.01);
    }

    #[test]
    fn single_element_configuration_references_its_own_geometry() {
        let config = Configuration::single(diamond("main", 2.0));
        assert_eq!(config.len(), 1);
        assert_relative_eq!(config.resolved_ref_chord(), 2.0, epsilon = 1e-12);
        assert_relative_eq!(config.resolved_ref_point().x, 0.5, epsilon = 1e-12);
    }

    #[test]
    fn empty_configuration_has_safe_reference_quantities() {
        let config = Configuration::new(vec![]);
        assert!(config.is_empty());
        assert_eq!(config.len(), 0);
        // Never zero: callers divide by this.
        assert_relative_eq!(config.resolved_ref_chord(), 1.0);
        assert_relative_eq!(config.resolved_ref_point().x, 0.0);
        assert_relative_eq!(config.resolved_ref_point().y, 0.0);
    }

    #[test]
    fn iteration_visits_elements_in_configuration_order() {
        let config = slat_main_flap();
        let names: Vec<&str> = config.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(names, ["slat", "main", "flap"]);

        let by_into_iter: Vec<&str> = (&config).into_iter().map(|e| e.name.as_str()).collect();
        assert_eq!(by_into_iter, names);
    }

    #[test]
    fn element_accessor_is_bounds_checked() {
        let config = slat_main_flap();
        assert_eq!(config.element(2).map(|e| e.name.as_str()), Some("flap"));
        assert!(config.element(3).is_none());
    }
}
