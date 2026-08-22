//! Clearance between the elements of a configuration.
//!
//! This module measures how close the elements of a
//! [`Configuration`] come to each other once their placements are applied, and
//! reports when a placement brings two elements closer than the geometry layer
//! is prepared to represent — including the case where two contours overlap.
//!
//! The measurement is geometric. Slot flow between closely spaced elements, and
//! the confluent boundary layers it produces, are out of scope for this engine;
//! that regime belongs to a Navier-Stokes solver. What this module does is keep
//! the boundary of that scope checkable: it is the test that says whether a
//! configuration is spaced widely enough for the engine's independent-wake
//! treatment to be a defensible approximation.
//!
//! # What is reported
//! [`Configuration::clearance`] returns a [`ClearanceReport`] holding, for every
//! pair of elements:
//!
//! - the **minimum surface-to-surface distance**, measured segment to segment
//!   over the placed contours (see [`PairClearance`]);
//! - whether the two contours **intersect** — cross, or one lies inside the
//!   other (see [`ContourRelation`]). An intersecting configuration is
//!   geometrically invalid and must not reach a solver;
//!
//! and, for each **adjacent** pair in chordwise order, the two rigging readouts
//! a high-lift engineer works in: **gap** and **overlap** (see
//! [`PairRigging`], which states the overlap sign convention precisely).
//!
//! A report rather than a bare pass/fail, so a caller can show *why* a
//! configuration was rejected. [`ClearanceReport::check`] reduces it to a
//! `Result` for callers that only need the gate.
//!
//! # Decision D6
//! The floor implemented here — reject below
//! [`DEFAULT_MIN_CLEARANCE_FRACTION`] of the reference chord — is the
//! *geometric* half of decision D6. The other half is a warning when a traced
//! wake passes close to a downstream element. That needs wake geometry, which
//! the geometry layer does not have; it belongs with the wake tracing work, and
//! is deliberately not attempted here.
//!
//! # Cost
//! Every element pair is swept segment against segment, so the work is
//! `sum over pairs of n_i * n_j` segment tests, each a handful of flops. For a
//! three-element configuration at a few hundred nodes per element that is of
//! order 10^5 tests — well under a millisecond, and not on any solver hot path.
//! No spatial index or bounding-box rejection is used: at these sizes the
//! straightforward sweep is fast enough, and being obviously correct matters
//! more here than being quick.

use core::cmp::Ordering;
use core::fmt;

use crate::configuration::{Configuration, Element};
use crate::error::GeometryError;
use crate::point::{cross_2d, point, Point};

/// Hard-fail floor for element-to-element clearance, as a fraction of the
/// configuration reference chord (decision D6, geometric half).
///
/// `0.005` is 0.5% of the reference chord.
///
/// # Why a floor at all
/// Its job is to keep degenerate geometry away from the panel method, not to
/// police aerodynamic validity. At the floor the two surfaces are already closer
/// together than the discretisation that represents them: on the 30P-30N section
/// below, the main element's 220 nodes give panels averaging 0.8% of the
/// retracted chord, so half a percent is under one panel length. Closer than
/// that the paneled representation of the space between the elements stops
/// meaning anything, and the step past it — contours touching or crossing — is a
/// configuration no panel method can be asked to solve.
///
/// # Why this value
/// It sits below every rigging gap in service use, so it rejects nothing a user
/// would legitimately ask for. Measured on the McDonnell Douglas 30P-30N
/// slat/main/flap section (`flexfoil-ui/public/airfoils/30p-30n.dat`), whose
/// published rigging table quotes gaps of 2.95% and 1.27% of the retracted
/// chord, the tightest element-to-element clearance is 1.27% — about 2.5 times
/// this floor.
///
/// It is a provisional value, and a starting point rather than a settled one:
/// the number that ends up shipping should be chosen against a wider set of
/// measured geometries. Override it per call with
/// [`ClearanceCriteria::with_min_clearance_fraction`] rather than editing this
/// constant.
pub const DEFAULT_MIN_CLEARANCE_FRACTION: f64 = 0.005;

/// What a pair of contours does to each other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContourRelation {
    /// The contours are disjoint: neither crosses the other and neither is
    /// inside the other. The only valid relation.
    Separate,

    /// The contours cross — some segment of one meets some segment of the
    /// other. Geometrically invalid.
    Crossing,

    /// Neither contour crosses the other, but one lies wholly inside it.
    ///
    /// Geometrically invalid, and not detectable from distances alone: a small
    /// element sitting well inside a large one can be a long way from its
    /// walls, so this case can pass a distance floor comfortably.
    Contained {
        /// Index of the element that lies inside the other.
        inside: usize,
    },
}

impl ContourRelation {
    /// True unless the relation is [`Separate`](ContourRelation::Separate).
    #[inline]
    pub fn intersects(&self) -> bool {
        !matches!(self, ContourRelation::Separate)
    }
}

/// Outcome of the clearance test, for one pair or for a whole configuration.
///
/// The variants are ordered by severity — `Pass < BelowFloor < Intersecting` —
/// so a configuration's verdict is the maximum over its pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ClearanceVerdict {
    /// Clear of the floor, and no intersection.
    Pass,

    /// Separate contours, but closer than the floor allows.
    BelowFloor,

    /// The contours intersect: they cross, or one is inside the other.
    Intersecting,
}

impl ClearanceVerdict {
    /// True for [`Pass`](ClearanceVerdict::Pass) only.
    #[inline]
    pub fn is_pass(&self) -> bool {
        matches!(self, ClearanceVerdict::Pass)
    }
}

impl fmt::Display for ClearanceVerdict {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            ClearanceVerdict::Pass => "pass",
            ClearanceVerdict::BelowFloor => "below floor",
            ClearanceVerdict::Intersecting => "intersecting",
        };
        f.write_str(text)
    }
}

/// What the clearance test accepts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceCriteria {
    /// Smallest element-to-element clearance that passes, as a fraction of the
    /// configuration reference chord.
    ///
    /// Defaults to [`DEFAULT_MIN_CLEARANCE_FRACTION`]. A non-positive value
    /// disables the floor, leaving only intersection detection.
    pub min_clearance_fraction: f64,
}

impl Default for ClearanceCriteria {
    fn default() -> Self {
        Self {
            min_clearance_fraction: DEFAULT_MIN_CLEARANCE_FRACTION,
        }
    }
}

impl ClearanceCriteria {
    /// Criteria with an explicit floor, as a fraction of the reference chord.
    #[inline]
    pub fn with_min_clearance_fraction(min_clearance_fraction: f64) -> Self {
        Self {
            min_clearance_fraction,
        }
    }

    /// The floor as an absolute distance, given the reference chord it is a
    /// fraction of.
    #[inline]
    pub fn floor(&self, ref_chord: f64) -> f64 {
        self.min_clearance_fraction * ref_chord
    }
}

/// How close two elements come, and whether their contours intersect.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairClearance {
    /// The two element indices, lower first.
    pub elements: (usize, usize),

    /// Minimum surface-to-surface distance between the two placed contours, in
    /// configuration units.
    ///
    /// Measured segment to segment: the closest approach of any surface segment
    /// of one element to any surface segment of the other, not merely the
    /// closest approach of their nodes. Node-to-node can only ever *over*state
    /// the clearance, which is the direction that hides a violation, so the
    /// segment measure is the one reported.
    ///
    /// Zero for crossing contours.
    pub min_distance: f64,

    /// [`min_distance`](PairClearance::min_distance) as a fraction of the
    /// configuration reference chord.
    pub min_distance_fraction: f64,

    /// The two points where the minimum is attained — the first on element
    /// `elements.0`, the second on element `elements.1`, both in configuration
    /// coordinates.
    ///
    /// A point may lie part-way along a segment rather than on a node. For
    /// crossing contours both points are a crossing point.
    pub witness: (Point, Point),

    /// Whether the contours are disjoint, cross, or one contains the other.
    pub relation: ContourRelation,

    /// Outcome for this pair.
    pub verdict: ClearanceVerdict,
}

impl PairClearance {
    /// True if this pair alone would not reject the configuration.
    #[inline]
    pub fn passes(&self) -> bool {
        self.verdict.is_pass()
    }

    /// True if the two contours intersect — see [`ContourRelation`].
    #[inline]
    pub fn intersects(&self) -> bool {
        self.relation.intersects()
    }

    /// True if this pair holds `element`.
    #[inline]
    pub fn involves(&self, element: usize) -> bool {
        self.elements.0 == element || self.elements.1 == element
    }
}

/// The rigging of one adjacent pair of elements: gap and overlap.
///
/// These are the two numbers a high-lift rigging table is written in, and both
/// are stated here in the vocabulary that table uses.
///
/// # Adjacency and ordering
/// Elements are ordered by the **x of their leading edge**
/// ([`Element::chord_endpoints`], placement applied), and adjacent pairs are the
/// consecutive ones in that order. For a slat/main/flap section that gives
/// slat → main → flap whatever order the elements were listed in. `upstream` is
/// the forward element of the pair, `downstream` the aft one.
///
/// # GAP
/// The minimum surface-to-surface distance for the pair — the same number as
/// [`PairClearance::min_distance`], restated here because it is the readout an
/// engineer rigging the configuration reads, and the one the floor is tested
/// against.
///
/// Classical rigging tables define the gap slightly more narrowly, as the
/// distance from the *upstream element's trailing edge* to the downstream
/// element's surface. The minimum over both whole surfaces can only be smaller,
/// so it is the conservative form of the same measurement. On the 30P-30N the
/// minimum is in fact attained at the upstream trailing edge and the two
/// definitions agree to the precision the published table is quoted at.
///
/// # OVERLAP — sign convention
/// Chordwise overlap, measured along the configuration x axis (the chordwise
/// direction; see the crate's coordinate convention):
///
/// ```text
/// overlap = x_TE(upstream) - x_LE(downstream)
/// ```
///
/// where `x_TE` is the x of the upstream element's trailing edge and `x_LE` the
/// x of the downstream element's leading edge, both from
/// [`Element::chord_endpoints`] with the placement applied.
///
/// - **Positive** overlap: the downstream element's leading edge lies *ahead
///   of* (upstream of) the upstream element's trailing edge, so the two elements
///   overlap chordwise. This is the normal state of a single-slotted flap tucked
///   under the main element's trailing edge.
/// - **Negative** overlap: the downstream leading edge lies *behind* the
///   upstream trailing edge — clear air between them in the chordwise sense.
///   This is the normal state of a deployed slat, and it is why published slat
///   rigging is quoted as a *negative* overhang.
/// - **Zero**: trailing edge and leading edge at the same chordwise station.
///
/// Getting this backwards is a classic error, so it is pinned by test against a
/// published table: the 30P-30N file quotes overhangs of -2.50% and +0.25% of
/// the retracted chord for its slat and flap, and this definition reproduces
/// both, sign included.
///
/// Measuring along x assumes the configuration is held in body axes, with the
/// reference chord line along x — the convention the rest of the crate uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairRigging {
    /// Index of the forward element of the pair.
    pub upstream: usize,

    /// Index of the aft element of the pair.
    pub downstream: usize,

    /// Minimum surface-to-surface distance, in configuration units. See the
    /// type documentation.
    pub gap: f64,

    /// [`gap`](PairRigging::gap) as a fraction of the configuration reference
    /// chord.
    pub gap_fraction: f64,

    /// Signed chordwise overlap, in configuration units. Positive when the
    /// downstream leading edge lies ahead of the upstream trailing edge — see
    /// the type documentation for the full convention.
    pub overlap: f64,

    /// [`overlap`](PairRigging::overlap) as a fraction of the configuration
    /// reference chord, carrying the same sign.
    pub overlap_fraction: f64,
}

/// Clearance diagnostics for a whole configuration.
///
/// # Which reference chord the fractions are against
/// Every `_fraction` field divides by [`ref_chord`](ClearanceReport::ref_chord),
/// which is [`Configuration::resolved_ref_chord`] — by default (decision D4) the
/// chord of the largest *placed* element.
///
/// That is not always the chord a published rigging table is quoted against. A
/// high-lift table normally normalises by the **retracted** (cruise) chord of the
/// whole section, which is longer than the main element alone: for the 30P-30N,
/// 1.0 against a main element chord of 0.832 in the same units, a factor of 1.20
/// between the two sets of percentages. Set
/// [`Configuration::with_ref_chord`] to the retracted chord to read fractions
/// that compare directly with such a table.
#[derive(Debug, Clone, PartialEq)]
pub struct ClearanceReport {
    /// The reference chord every `_fraction` in this report is against.
    pub ref_chord: f64,

    /// The criteria the report was measured against.
    pub criteria: ClearanceCriteria,

    /// The floor as an absolute distance:
    /// `criteria.min_clearance_fraction * ref_chord`.
    pub floor: f64,

    /// Every unordered pair of elements, in `(0,1), (0,2), … (1,2), …` order.
    /// Empty for a configuration of fewer than two elements.
    pub pairs: Vec<PairClearance>,

    /// Rigging readouts for the adjacent pairs, in chordwise order — one fewer
    /// than the number of elements that have a chord line to measure between,
    /// and empty for a configuration of fewer than two.
    pub rigging: Vec<PairRigging>,
}

impl ClearanceReport {
    /// The worst outcome over all pairs — the configuration's verdict.
    ///
    /// [`Pass`](ClearanceVerdict::Pass) for a configuration with no pairs at
    /// all: a single element has nothing to clear.
    pub fn verdict(&self) -> ClearanceVerdict {
        self.pairs
            .iter()
            .map(|pair| pair.verdict)
            .max()
            .unwrap_or(ClearanceVerdict::Pass)
    }

    /// True if the configuration is fit to panel and solve as far as element
    /// clearance is concerned.
    #[inline]
    pub fn passes(&self) -> bool {
        self.verdict().is_pass()
    }

    /// Reduce the report to a `Result`, for a caller that only needs the gate.
    ///
    /// Errors:
    /// - contours intersecting →
    ///   `InvalidParameter { name: "element_intersection", value: <lower index
    ///   of the first intersecting pair> }`
    /// - clearance below the floor →
    ///   `InvalidParameter { name: "element_clearance", value: <the measured
    ///   clearance as a fraction of the reference chord> }`
    ///
    /// Intersection is reported first when both apply, being the stronger
    /// statement. Keep the report itself when the message matters: it says which
    /// pair, how close, and where.
    pub fn check(&self) -> Result<(), GeometryError> {
        if let Some(pair) = self
            .pairs
            .iter()
            .find(|pair| pair.verdict == ClearanceVerdict::Intersecting)
        {
            return Err(GeometryError::InvalidParameter {
                name: "element_intersection",
                value: pair.elements.0 as f64,
            });
        }
        if let Some(pair) = self
            .pairs
            .iter()
            .find(|pair| pair.verdict == ClearanceVerdict::BelowFloor)
        {
            return Err(GeometryError::InvalidParameter {
                name: "element_clearance",
                value: pair.min_distance_fraction,
            });
        }
        Ok(())
    }

    /// The smallest element-to-element clearance in the configuration, or `None`
    /// if there are no pairs.
    ///
    /// Note that this is not on its own a validity test: a contained element can
    /// be far from the contour that holds it. Use [`verdict`](Self::verdict).
    pub fn min_clearance(&self) -> Option<f64> {
        self.tightest_pair().map(|pair| pair.min_distance)
    }

    /// [`min_clearance`](Self::min_clearance) as a fraction of the reference
    /// chord.
    pub fn min_clearance_fraction(&self) -> Option<f64> {
        self.tightest_pair().map(|pair| pair.min_distance_fraction)
    }

    /// The pair that comes closest to touching. Ties go to the earlier pair.
    pub fn tightest_pair(&self) -> Option<&PairClearance> {
        self.pairs.iter().fold(None, |best, pair| match best {
            Some(b) if b.min_distance <= pair.min_distance => Some(b),
            _ => Some(pair),
        })
    }

    /// The pairs that do not pass, in report order.
    pub fn failures(&self) -> impl Iterator<Item = &PairClearance> {
        self.pairs.iter().filter(|pair| !pair.passes())
    }

    /// The pairs whose contours intersect, in report order.
    pub fn intersections(&self) -> impl Iterator<Item = &PairClearance> {
        self.pairs.iter().filter(|pair| pair.intersects())
    }

    /// The entry for one pair of elements, in either order.
    pub fn pair(&self, a: usize, b: usize) -> Option<&PairClearance> {
        self.pairs
            .iter()
            .find(|pair| pair.elements == (a.min(b), a.max(b)))
    }

    /// The rigging entry for one adjacent pair, in either order. `None` if the
    /// two elements are not adjacent in chordwise order.
    pub fn rigging_between(&self, a: usize, b: usize) -> Option<&PairRigging> {
        self.rigging
            .iter()
            .find(|r| (r.upstream, r.downstream) == (a, b) || (r.upstream, r.downstream) == (b, a))
    }
}

impl fmt::Display for ClearanceReport {
    /// A short human-readable summary, for a CLI to print when a configuration
    /// is rejected.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "clearance {}: reference chord {:.6}, floor {:.3}% chord",
            self.verdict(),
            self.ref_chord,
            100.0 * self.criteria.min_clearance_fraction
        )?;

        for pair in &self.pairs {
            let (a, b) = pair.elements;
            let note = match pair.relation {
                ContourRelation::Crossing => "contours cross".to_string(),
                ContourRelation::Contained { inside } => {
                    let other = if inside == a { b } else { a };
                    format!("element {inside} lies inside element {other}")
                }
                ContourRelation::Separate => match pair.verdict {
                    ClearanceVerdict::Pass => "pass".to_string(),
                    _ => format!(
                        "below the {:.3}% chord floor",
                        100.0 * self.criteria.min_clearance_fraction
                    ),
                },
            };
            writeln!(
                f,
                "  elements {a}-{b}: min distance {:.6} ({:.3}% chord) - {note}",
                pair.min_distance,
                100.0 * pair.min_distance_fraction
            )?;
        }

        for r in &self.rigging {
            writeln!(
                f,
                "  rigging {} -> {}: gap {:.3}% chord, overlap {:+.3}% chord",
                r.upstream,
                r.downstream,
                100.0 * r.gap_fraction,
                100.0 * r.overlap_fraction
            )?;
        }

        Ok(())
    }
}

impl Configuration {
    /// Clearance diagnostics with the default criteria.
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::body::Body;
    /// use rustfoil_core::configuration::{Configuration, Element};
    /// use rustfoil_core::placement::Placement;
    /// use rustfoil_core::point::point;
    ///
    /// let contour = |chord: f64| vec![
    ///     point(chord, 0.0),
    ///     point(0.5 * chord, -0.05 * chord),
    ///     point(0.0, 0.0),
    ///     point(0.5 * chord, 0.05 * chord),
    ///     point(chord, 0.0),
    /// ];
    ///
    /// let main = Element::from_body(Body::from_points("main", &contour(1.0)).unwrap());
    /// let mut flap = Element::from_body(Body::from_points("flap", &contour(0.3)).unwrap());
    /// // Tucked under the main element's trailing edge and below it.
    /// flap.placement = Placement::from_translation(0.98, -0.06);
    ///
    /// let report = Configuration::new(vec![main, flap]).clearance();
    ///
    /// assert!(report.passes());
    /// // One pair, one adjacent-pair rigging readout.
    /// assert_eq!(report.pairs.len(), 1);
    /// assert_eq!(report.rigging.len(), 1);
    ///
    /// // The main element is the forward one of the pair.
    /// let station = report.rigging[0];
    /// assert_eq!((station.upstream, station.downstream), (0, 1));
    /// // The flap's leading edge is at x = 0.98, ahead of the main element's
    /// // trailing edge at x = 1.0, so the two overlap chordwise by 0.02 and
    /// // the overlap is positive.
    /// assert!((station.overlap - 0.02).abs() < 1e-12);
    /// ```
    pub fn clearance(&self) -> ClearanceReport {
        self.clearance_with(ClearanceCriteria::default())
    }

    /// Clearance diagnostics against explicit criteria.
    ///
    /// A clearance exactly equal to the floor passes; only a strictly smaller
    /// one fails.
    pub fn clearance_with(&self, criteria: ClearanceCriteria) -> ClearanceReport {
        let ref_chord = self.resolved_ref_chord();
        let floor = criteria.floor(ref_chord);
        let contours: Vec<Vec<Point>> = self.iter().map(element_contour).collect();

        let mut pairs = Vec::new();
        for i in 0..contours.len() {
            for j in (i + 1)..contours.len() {
                let approach = contour_approach(&contours[i], &contours[j]);
                let relation =
                    contour_relation(i, j, &contours[i], &contours[j], approach.crossing);

                let verdict = if relation.intersects() {
                    ClearanceVerdict::Intersecting
                } else if approach.distance < floor {
                    ClearanceVerdict::BelowFloor
                } else {
                    ClearanceVerdict::Pass
                };

                pairs.push(PairClearance {
                    elements: (i, j),
                    min_distance: approach.distance,
                    min_distance_fraction: approach.distance / ref_chord,
                    witness: (approach.on_a, approach.on_b),
                    relation,
                    verdict,
                });
            }
        }

        let rigging = rigging_stations(self, &pairs, ref_chord);

        ClearanceReport {
            ref_chord,
            criteria,
            floor,
            pairs,
            rigging,
        }
    }
}

/// Rigging readouts for the adjacent pairs of `config`, given the pair
/// clearances the gaps are taken from.
///
/// Adjacency is chordwise, by leading-edge x — see [`PairRigging`].
fn rigging_stations(
    config: &Configuration,
    pairs: &[PairClearance],
    ref_chord: f64,
) -> Vec<PairRigging> {
    // Elements with no panels have no chord line to measure against, so they
    // take no part in a rigging readout.
    let mut order: Vec<(usize, Point, Point)> = config
        .iter()
        .enumerate()
        .filter_map(|(i, element)| {
            let (le, te) = element.chord_endpoints()?;
            Some((i, le, te))
        })
        .collect();

    // Stable, so equal leading-edge stations keep configuration order.
    order.sort_by(|(_, le_a, _), (_, le_b, _)| {
        le_a.x.partial_cmp(&le_b.x).unwrap_or(Ordering::Equal)
    });

    order
        .windows(2)
        .map(|w| {
            let (upstream, _, te_upstream) = w[0];
            let (downstream, le_downstream, _) = w[1];

            // The gap is the pair's measured minimum distance, looked up rather
            // than measured again so the two readouts cannot disagree. Every
            // pair of existing elements is in `pairs`, so the fallback is
            // unreachable.
            let gap = pairs
                .iter()
                .find(|pair| pair.elements == (upstream.min(downstream), upstream.max(downstream)))
                .map_or(f64::NAN, |pair| pair.min_distance);

            let overlap = te_upstream.x - le_downstream.x;

            PairRigging {
                upstream,
                downstream,
                gap,
                gap_fraction: gap / ref_chord,
                overlap,
                overlap_fraction: overlap / ref_chord,
            }
        })
        .collect()
}

/// The element's contour as a cyclic sequence of nodes, in configuration
/// coordinates.
///
/// One node per distinct point of the contour — [`crate::Body::n_nodes`] of
/// them — with the surface understood to run from each node to the next and
/// from the last back to the first. That last, wrapping segment is what closes
/// the contour: for a body with a blunt trailing edge it is the trailing-edge
/// base, a real surface that has to be measured like any other, and which the
/// body's panel list does not contain.
///
/// The node ordering is the one [`crate::Layout`] numbers in, so node `k` here
/// is the element's local node `k` there.
pub fn element_contour(element: &Element) -> Vec<Point> {
    let panels = element.body.panels();
    let mut nodes = Vec::with_capacity(element.body.n_nodes());

    for panel in panels {
        nodes.push(element.placement.apply(panel.p1));
    }
    if !element.body.is_closed() {
        if let Some(last) = panels.last() {
            nodes.push(element.placement.apply(last.p2));
        }
    }

    nodes
}

/// Closest approach of two contours.
#[derive(Debug, Clone, Copy)]
struct Approach {
    /// Minimum distance; `0.0` when the two cross.
    distance: f64,
    /// Where the minimum is attained on the first contour.
    on_a: Point,
    /// Where the minimum is attained on the second contour.
    on_b: Point,
    /// Whether any segment of one met any segment of the other.
    crossing: bool,
}

impl Approach {
    /// A starting value for a minimisation: nothing measured yet.
    fn none_yet() -> Self {
        Self {
            distance: f64::INFINITY,
            on_a: point(0.0, 0.0),
            on_b: point(0.0, 0.0),
            crossing: false,
        }
    }
}

/// The segments of a cyclic contour, including the one that wraps from the last
/// node back to the first. Empty for a contour of fewer than two nodes.
fn contour_segments(contour: &[Point]) -> impl Iterator<Item = (Point, Point)> + '_ {
    let n = contour.len();
    let count = if n < 2 { 0 } else { n };
    (0..count).map(move |i| (contour[i], contour[(i + 1) % n]))
}

/// Minimum surface-to-surface distance between two contours, and whether they
/// cross.
///
/// Every segment of `a` against every segment of `b`: `n_a * n_b` tests. See the
/// module documentation on cost.
fn contour_approach(a: &[Point], b: &[Point]) -> Approach {
    let mut best = Approach::none_yet();

    for (a0, a1) in contour_segments(a) {
        for (b0, b1) in contour_segments(b) {
            let approach = segment_approach(a0, a1, b0, b1);

            // The first crossing found settles the pair at zero distance;
            // nothing later can beat it.
            if approach.crossing {
                if !best.crossing {
                    best = approach;
                }
            } else if !best.crossing && approach.distance < best.distance {
                best = approach;
            }
        }
    }

    best
}

/// How two contours are placed relative to each other, given whether they were
/// found to cross.
///
/// Without a crossing, either the two are disjoint or one is wholly inside the
/// other, so testing a single node of each against the other contour settles it.
fn contour_relation(
    i: usize,
    j: usize,
    a: &[Point],
    b: &[Point],
    crossing: bool,
) -> ContourRelation {
    if crossing {
        return ContourRelation::Crossing;
    }
    if a.first().is_some_and(|p| point_inside_contour(*p, b)) {
        return ContourRelation::Contained { inside: i };
    }
    if b.first().is_some_and(|p| point_inside_contour(*p, a)) {
        return ContourRelation::Contained { inside: j };
    }
    ContourRelation::Separate
}

/// Closest approach of two line segments.
fn segment_approach(a0: Point, a1: Point, b0: Point, b1: Point) -> Approach {
    if let Some(meeting) = segment_crossing(a0, a1, b0, b1) {
        return Approach {
            distance: 0.0,
            on_a: meeting,
            on_b: meeting,
            crossing: true,
        };
    }

    // Two segments that do not meet attain their minimum separation at an
    // endpoint of one of them, so these four tests are exact rather than a
    // sampling.
    let mut best = Approach::none_yet();
    for p in [a0, a1] {
        let (distance, q) = closest_point_on_segment(p, b0, b1);
        if distance < best.distance {
            best = Approach {
                distance,
                on_a: p,
                on_b: q,
                crossing: false,
            };
        }
    }
    for p in [b0, b1] {
        let (distance, q) = closest_point_on_segment(p, a0, a1);
        if distance < best.distance {
            best = Approach {
                distance,
                on_a: q,
                on_b: p,
                crossing: false,
            };
        }
    }
    best
}

/// Distance from `p` to the segment `a`-`b`, and the point on the segment where
/// it is attained.
fn closest_point_on_segment(p: Point, a: Point, b: Point) -> (f64, Point) {
    let d = b - a;
    let length_squared = d.norm_squared();
    let q = if length_squared > 0.0 {
        let t = ((p - a).dot(&d) / length_squared).clamp(0.0, 1.0);
        a + d * t
    } else {
        a
    };
    ((p - q).norm(), q)
}

/// Where two segments meet, if they meet at all.
///
/// The generic case solves for the two segment parameters. Exactly parallel
/// segments are handled separately: they meet only if they are collinear and
/// their spans overlap, in which case the start of the shared span is returned.
///
/// Near-parallel segments go through the generic case, where the parameters lose
/// precision; a crossing so shallow that it is missed there has its surfaces
/// within round-off of each other, and is caught by the distance floor instead.
fn segment_crossing(a0: Point, a1: Point, b0: Point, b1: Point) -> Option<Point> {
    let da = a1 - a0;
    let db = b1 - b0;
    let r = b0 - a0;

    let denominator = cross_2d(&da, &db);
    if denominator != 0.0 {
        let t = cross_2d(&r, &db) / denominator;
        let u = cross_2d(&r, &da) / denominator;
        if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
            return Some(a0 + da * t);
        }
        return None;
    }

    // Parallel. Only a collinear pair can meet.
    if cross_2d(&r, &da) != 0.0 {
        return None;
    }
    let length_squared = da.norm_squared();
    if length_squared <= 0.0 {
        return None;
    }
    let t0 = r.dot(&da) / length_squared;
    let t1 = (b1 - a0).dot(&da) / length_squared;
    let low = t0.min(t1).max(0.0);
    let high = t0.max(t1).min(1.0);
    if low <= high {
        Some(a0 + da * low)
    } else {
        None
    }
}

/// Whether `p` lies inside a closed contour.
///
/// Crossing number against a ray in +x: count the contour segments that straddle
/// the point's y and pass to its right. An odd count means inside. Points
/// exactly on the contour may report either way, which does not matter here
/// because a contour that passes through another contour's node has already been
/// found to cross.
fn point_inside_contour(p: Point, contour: &[Point]) -> bool {
    if contour.len() < 3 {
        return false;
    }

    let mut inside = false;
    for (a, b) in contour_segments(contour) {
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if p.x < x {
                inside = !inside;
            }
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Body;
    use crate::placement::Placement;
    use crate::point::point;
    use crate::spline::PanelingParams;
    use approx::assert_relative_eq;
    use std::path::PathBuf;

    // --- Synthetic geometry ----------------------------------------------

    /// A closed diamond of the given chord, leading edge at the origin.
    fn diamond(name: &str, chord: f64) -> Body {
        Body::from_points(
            name,
            &[
                point(chord, 0.0),
                point(0.5 * chord, -0.05 * chord),
                point(0.0, 0.0),
                point(0.5 * chord, 0.05 * chord),
                point(chord, 0.0),
            ],
        )
        .unwrap()
    }

    /// A closed axis-aligned rectangle, nodes at the corners only, so that any
    /// close approach along an edge is mid-segment rather than at a node.
    fn rectangle(name: &str, x0: f64, y0: f64, width: f64, height: f64) -> Body {
        Body::from_points(
            name,
            &[
                point(x0 + width, y0),
                point(x0, y0),
                point(x0, y0 + height),
                point(x0 + width, y0 + height),
                point(x0 + width, y0),
            ],
        )
        .unwrap()
    }

    fn element(body: Body) -> Element {
        Element::from_body(body)
    }

    /// Smallest node-to-node distance between two elements — the measure this
    /// module deliberately does *not* use, kept so tests can show the
    /// difference.
    fn node_to_node(a: &Element, b: &Element) -> f64 {
        let (ca, cb) = (element_contour(a), element_contour(b));
        ca.iter()
            .flat_map(|p| cb.iter().map(move |q| (p - q).norm()))
            .fold(f64::INFINITY, f64::min)
    }

    // --- Trivial cases ---------------------------------------------------

    #[test]
    fn a_single_element_has_nothing_to_clear() {
        let report = Configuration::single(diamond("main", 1.0)).clearance();

        assert!(report.pairs.is_empty());
        assert!(report.rigging.is_empty());
        assert_eq!(report.verdict(), ClearanceVerdict::Pass);
        assert!(report.passes());
        assert_eq!(report.min_clearance(), None);
        assert!(report.check().is_ok());
    }

    #[test]
    fn an_empty_configuration_passes() {
        let report = Configuration::new(vec![]).clearance();

        assert!(report.pairs.is_empty());
        assert!(report.rigging.is_empty());
        assert!(report.passes());
        // resolved_ref_chord never returns zero, so the fractions stay finite.
        assert_relative_eq!(report.ref_chord, 1.0);
    }

    #[test]
    fn elements_far_apart_pass_with_a_large_clearance() {
        let mut aft = element(diamond("aft", 1.0));
        aft.placement = Placement::from_translation(5.0, 0.0);

        let report = Configuration::new(vec![element(diamond("fore", 1.0)), aft]).clearance();

        assert!(report.passes());
        assert_eq!(report.pairs.len(), 1);
        assert_eq!(report.pairs[0].relation, ContourRelation::Separate);
        // Fore trailing edge at x = 1, aft leading edge at x = 5.
        assert_relative_eq!(report.pairs[0].min_distance, 4.0, epsilon = 1e-12);
        assert_relative_eq!(report.pairs[0].min_distance_fraction, 4.0, epsilon = 1e-12);
    }

    // --- Intersection ----------------------------------------------------

    #[test]
    fn crossing_contours_are_rejected() {
        let mut aft = element(diamond("aft", 1.0));
        // Half a chord downstream: the two diamonds interpenetrate.
        aft.placement = Placement::from_translation(0.5, 0.0);

        let report = Configuration::new(vec![element(diamond("fore", 1.0)), aft]).clearance();

        assert_eq!(report.pairs[0].relation, ContourRelation::Crossing);
        assert_eq!(report.pairs[0].verdict, ClearanceVerdict::Intersecting);
        assert_relative_eq!(report.pairs[0].min_distance, 0.0);
        assert!(!report.passes());
        assert_eq!(report.intersections().count(), 1);
        assert!(matches!(
            report.check(),
            Err(GeometryError::InvalidParameter {
                name: "element_intersection",
                ..
            })
        ));
    }

    #[test]
    fn an_element_inside_another_is_rejected() {
        // A small diamond well inside a large rectangle: no contour crosses the
        // other, and the surfaces are nowhere near each other, so a distance
        // floor alone would pass this.
        let mut inner = element(diamond("inner", 0.2));
        inner.placement = Placement::from_translation(0.4, 0.5);

        let report =
            Configuration::new(vec![element(rectangle("outer", 0.0, 0.0, 2.0, 1.0)), inner])
                .clearance();

        assert_eq!(
            report.pairs[0].relation,
            ContourRelation::Contained { inside: 1 }
        );
        assert_eq!(report.pairs[0].verdict, ClearanceVerdict::Intersecting);
        assert!(report.pairs[0].min_distance > report.floor);
        assert!(!report.passes());
    }

    #[test]
    fn containment_is_reported_whichever_element_is_the_smaller() {
        // Same geometry, elements listed the other way round.
        let mut inner = element(diamond("inner", 0.2));
        inner.placement = Placement::from_translation(0.4, 0.5);

        let report =
            Configuration::new(vec![inner, element(rectangle("outer", 0.0, 0.0, 2.0, 1.0))])
                .clearance();

        assert_eq!(
            report.pairs[0].relation,
            ContourRelation::Contained { inside: 0 }
        );
    }

    #[test]
    fn touching_contours_count_as_crossing() {
        // The aft rectangle's leading edge sits exactly on the fore
        // rectangle's trailing edge.
        let report = Configuration::new(vec![
            element(rectangle("fore", 0.0, 0.0, 1.0, 0.1)),
            element(rectangle("aft", 1.0, 0.0, 1.0, 0.1)),
        ])
        .clearance();

        assert_eq!(report.pairs[0].relation, ContourRelation::Crossing);
        assert!(!report.passes());
    }

    // --- The floor -------------------------------------------------------

    /// Two rectangles separated by a controlled vertical distance, with the
    /// reference chord pinned to 1.0 so the fraction is the distance itself.
    fn stacked_at(separation: f64) -> Configuration {
        Configuration::new(vec![
            element(rectangle("lower", 0.0, 0.0, 1.0, 0.1)),
            element(rectangle("upper", 0.0, 0.1 + separation, 1.0, 0.1)),
        ])
        .with_ref_chord(1.0)
    }

    #[test]
    fn just_above_the_floor_passes() {
        let separation = DEFAULT_MIN_CLEARANCE_FRACTION * 1.2;
        let report = stacked_at(separation).clearance();

        assert_relative_eq!(report.floor, DEFAULT_MIN_CLEARANCE_FRACTION);
        assert_relative_eq!(report.pairs[0].min_distance, separation, epsilon = 1e-12);
        assert_eq!(report.pairs[0].verdict, ClearanceVerdict::Pass);
        assert!(report.passes());
    }

    #[test]
    fn just_below_the_floor_is_rejected() {
        let separation = DEFAULT_MIN_CLEARANCE_FRACTION * 0.8;
        let report = stacked_at(separation).clearance();

        assert_relative_eq!(report.pairs[0].min_distance, separation, epsilon = 1e-12);
        assert_eq!(report.pairs[0].verdict, ClearanceVerdict::BelowFloor);
        assert_eq!(report.pairs[0].relation, ContourRelation::Separate);
        assert!(!report.passes());
        assert_eq!(report.failures().count(), 1);
        assert!(matches!(
            report.check(),
            Err(GeometryError::InvalidParameter {
                name: "element_clearance",
                ..
            })
        ));
    }

    #[test]
    fn the_floor_is_a_fraction_of_the_reference_chord() {
        // The same absolute separation, read against two reference chords.
        let separation = 0.01;
        let tight = stacked_at(separation).with_ref_chord(4.0).clearance();
        let loose = stacked_at(separation).with_ref_chord(1.0).clearance();

        assert_relative_eq!(tight.floor, 0.02);
        assert_eq!(tight.pairs[0].verdict, ClearanceVerdict::BelowFloor);

        assert_relative_eq!(loose.floor, 0.005);
        assert_eq!(loose.pairs[0].verdict, ClearanceVerdict::Pass);
    }

    #[test]
    fn the_floor_is_overridable() {
        let separation = DEFAULT_MIN_CLEARANCE_FRACTION * 0.8;
        let config = stacked_at(separation);

        assert!(!config.clearance().passes());
        assert!(config
            .clearance_with(ClearanceCriteria::with_min_clearance_fraction(0.001))
            .passes());
        // A non-positive floor leaves only intersection detection.
        assert!(config
            .clearance_with(ClearanceCriteria::with_min_clearance_fraction(0.0))
            .passes());
    }

    // --- Segment-to-segment, not node-to-node ----------------------------

    #[test]
    fn the_minimum_is_measured_between_segments_not_between_nodes() {
        // A wide rectangle whose upper surface is a single segment from (0, 0)
        // to (1, 0), and a tall narrow tab whose lowest vertex hangs just above
        // the middle of it. The closest nodes are half a chord apart; the
        // surfaces are 0.4% of chord apart.
        let separation = DEFAULT_MIN_CLEARANCE_FRACTION * 0.8;
        let plate = element(rectangle("plate", 0.0, -0.1, 1.0, 0.1));
        let tab = element(
            Body::from_points(
                "tab",
                &[
                    point(0.55, 1.0),
                    point(0.5, separation),
                    point(0.45, 1.0),
                    point(0.55, 1.0),
                ],
            )
            .unwrap(),
        );

        let node_gap = node_to_node(&plate, &tab);
        let report = Configuration::new(vec![plate, tab])
            .with_ref_chord(1.0)
            .clearance();

        // Node to node would have called this configuration clear by a wide
        // margin; the surfaces are in fact below the floor.
        assert!(node_gap > 0.49, "node-to-node gap = {node_gap}");
        assert!(node_gap > report.floor);
        assert_relative_eq!(report.pairs[0].min_distance, separation, epsilon = 1e-12);
        assert_eq!(report.pairs[0].verdict, ClearanceVerdict::BelowFloor);
    }

    #[test]
    fn the_witness_points_lie_on_the_two_contours() {
        let separation = 0.02;
        let report = stacked_at(separation).clearance();

        let (on_lower, on_upper) = report.pairs[0].witness;
        assert_relative_eq!(on_lower.y, 0.1, epsilon = 1e-12);
        assert_relative_eq!(on_upper.y, 0.1 + separation, epsilon = 1e-12);
        assert_relative_eq!(
            (on_upper - on_lower).norm(),
            report.pairs[0].min_distance,
            epsilon = 1e-12
        );
    }

    #[test]
    fn the_trailing_edge_base_of_a_blunt_element_is_measured() {
        // An open contour: its last node and its first do not meet, so the
        // segment between them is a real surface the body's panel list does not
        // contain. Placed so that base is the closest thing to the other
        // element, it has to be measured or the clearance comes out too large.
        let blunt = Body::from_points(
            "blunt",
            &[
                point(1.0, -0.02),
                point(0.5, -0.05),
                point(0.0, 0.0),
                point(0.5, 0.05),
                point(1.0, 0.02),
            ],
        )
        .unwrap();
        assert!(!blunt.is_closed());

        let element_count = element_contour(&element(blunt.clone())).len();
        assert_eq!(element_count, blunt.n_nodes());

        // A tab facing the middle of the trailing-edge base, level with y = 0,
        // where no node of the blunt body sits.
        let mut tab = element(rectangle("tab", 0.0, -0.01, 0.2, 0.02));
        tab.placement = Placement::from_translation(1.03, 0.0);

        let report = Configuration::new(vec![element(blunt), tab])
            .with_ref_chord(1.0)
            .clearance();

        // The base spans x = 1.0; the tab's leading face is at x = 1.03.
        assert_relative_eq!(report.pairs[0].min_distance, 0.03, epsilon = 1e-12);
    }

    // --- Rigging: order and sign ----------------------------------------

    #[test]
    fn rigging_pairs_are_chordwise_adjacent_whatever_the_element_order() {
        // Listed flap, main, slat; rigged slat -> main -> flap.
        let mut slat = element(diamond("slat", 0.15));
        slat.placement = Placement::from_translation(-0.12, 0.03);
        let main = element(diamond("main", 1.0));
        let mut flap = element(diamond("flap", 0.3));
        flap.placement = Placement::from_translation(0.95, -0.08);

        let report = Configuration::new(vec![flap, main, slat]).clearance();

        // Three pairs, but only two adjacent stations.
        assert_eq!(report.pairs.len(), 3);
        let stations: Vec<(usize, usize)> = report
            .rigging
            .iter()
            .map(|r| (r.upstream, r.downstream))
            .collect();
        assert_eq!(stations, [(2, 1), (1, 0)]);

        // The gap at a station is the pair's minimum distance.
        let station = report.rigging_between(1, 0).unwrap();
        assert_relative_eq!(station.gap, report.pair(0, 1).unwrap().min_distance);
        assert!(report.rigging_between(2, 0).is_none(), "not adjacent");
    }

    #[test]
    fn overlap_is_positive_when_the_downstream_leading_edge_is_ahead() {
        // Upstream trailing edge at x = 1.0; downstream leading edge at
        // x = 0.9, tucked under it. Overlap = 1.0 - 0.9 = +0.1.
        let mut aft = element(rectangle("aft", 0.0, -0.2, 1.0, 0.1));
        aft.placement = Placement::from_translation(0.9, 0.0);

        let report = Configuration::new(vec![element(rectangle("fore", 0.0, 0.0, 1.0, 0.1)), aft])
            .with_ref_chord(1.0)
            .clearance();

        let station = report.rigging[0];
        assert_eq!((station.upstream, station.downstream), (0, 1));
        assert_relative_eq!(station.overlap, 0.1, epsilon = 1e-12);
        assert_relative_eq!(station.overlap_fraction, 0.1, epsilon = 1e-12);
    }

    #[test]
    fn overlap_is_negative_when_the_downstream_leading_edge_is_behind() {
        // Downstream leading edge at x = 1.1, aft of the upstream trailing edge
        // at x = 1.0. Overlap = 1.0 - 1.1 = -0.1.
        let mut aft = element(rectangle("aft", 0.0, -0.2, 1.0, 0.1));
        aft.placement = Placement::from_translation(1.1, 0.0);

        let report = Configuration::new(vec![element(rectangle("fore", 0.0, 0.0, 1.0, 0.1)), aft])
            .with_ref_chord(1.0)
            .clearance();

        assert_relative_eq!(report.rigging[0].overlap, -0.1, epsilon = 1e-12);
    }

    #[test]
    fn overlap_follows_the_placement_not_the_stored_coordinates() {
        let fore = element(rectangle("fore", 0.0, 0.0, 1.0, 0.1));
        let aft_body = rectangle("aft", 0.0, -0.2, 1.0, 0.1);

        let overlap_of = |dx: f64| {
            let mut aft = element(aft_body.clone());
            aft.placement = Placement::from_translation(dx, 0.0);
            Configuration::new(vec![fore.clone(), aft])
                .with_ref_chord(1.0)
                .clearance()
                .rigging[0]
                .overlap
        };

        assert_relative_eq!(overlap_of(0.8), 0.2, epsilon = 1e-12);
        assert_relative_eq!(overlap_of(1.0), 0.0, epsilon = 1e-12);
        assert_relative_eq!(overlap_of(1.2), -0.2, epsilon = 1e-12);
    }

    #[test]
    fn a_placement_scale_is_carried_into_the_measurement() {
        // Halving the aft element's size halves its chord, so its leading edge
        // and the clearance both move.
        let mut aft = element(rectangle("aft", 0.0, 0.4, 1.0, 0.1));
        aft.placement = Placement {
            pivot: point(0.0, 0.4),
            rotation_deg: 0.0,
            translation: crate::point::vec2(0.0, 0.0),
            scale: 0.5,
        };

        let report = Configuration::new(vec![element(rectangle("fore", 0.0, 0.0, 1.0, 0.1)), aft])
            .with_ref_chord(1.0)
            .clearance();

        // Scaled about (0, 0.4): the aft rectangle now spans y = 0.4 to 0.45.
        assert_relative_eq!(report.pairs[0].min_distance, 0.3, epsilon = 1e-12);
    }

    // --- Verdict plumbing ------------------------------------------------

    #[test]
    fn verdicts_are_ordered_by_severity() {
        assert!(ClearanceVerdict::Pass < ClearanceVerdict::BelowFloor);
        assert!(ClearanceVerdict::BelowFloor < ClearanceVerdict::Intersecting);
    }

    #[test]
    fn the_configuration_verdict_is_the_worst_of_its_pairs() {
        // Element 1 lies inside element 0; element 2 is clear of both.
        let mut inner = element(diamond("inner", 0.5));
        inner.placement = Placement::from_translation(0.25, 0.0);
        let mut far = element(diamond("far", 0.5));
        far.placement = Placement::from_translation(8.0, 0.0);

        let report =
            Configuration::new(vec![element(diamond("outer", 1.0)), inner, far]).clearance();

        assert_eq!(report.pairs.len(), 3);
        assert_eq!(report.verdict(), ClearanceVerdict::Intersecting);
        assert_eq!(report.pair(0, 2).unwrap().verdict, ClearanceVerdict::Pass);
        assert_eq!(report.pair(1, 2).unwrap().verdict, ClearanceVerdict::Pass);
    }

    #[test]
    fn the_tightest_pair_is_the_closest_one() {
        let mut near = element(diamond("near", 0.3));
        near.placement = Placement::from_translation(1.1, 0.0);
        let mut far = element(diamond("far", 0.3));
        far.placement = Placement::from_translation(4.0, 0.0);

        let report = Configuration::new(vec![element(diamond("main", 1.0)), near, far]).clearance();

        assert_eq!(report.tightest_pair().unwrap().elements, (0, 1));
        assert_relative_eq!(
            report.min_clearance().unwrap(),
            report.pair(0, 1).unwrap().min_distance
        );
    }

    #[test]
    fn the_report_prints_the_reason_a_configuration_was_rejected() {
        let separation = DEFAULT_MIN_CLEARANCE_FRACTION * 0.8;
        let text = stacked_at(separation).clearance().to_string();

        assert!(text.contains("clearance below floor"), "{text}");
        assert!(text.contains("elements 0-1"), "{text}");
        assert!(text.contains("floor"), "{text}");
        assert!(text.contains("rigging"), "{text}");
    }

    #[test]
    fn pair_lookup_accepts_either_order() {
        let mut aft = element(diamond("aft", 0.3));
        aft.placement = Placement::from_translation(2.0, 0.0);
        let report = Configuration::new(vec![element(diamond("fore", 1.0)), aft]).clearance();

        assert_eq!(report.pair(0, 1), report.pair(1, 0));
        assert!(report.pair(0, 2).is_none());
        assert!(report.pairs[0].involves(1));
        assert!(!report.pairs[0].involves(2));
    }

    #[test]
    fn paneling_parameters_do_not_affect_the_measurement() {
        // Clearance is measured on the contour as it stands. Re-paneling is a
        // separate step, and until it has run the paneling parameters cannot
        // change a distance.
        let base = Configuration::new(vec![
            element(rectangle("lower", 0.0, 0.0, 1.0, 0.1)),
            element(rectangle("upper", 0.0, 0.15, 1.0, 0.1)),
        ]);
        let mut respecified = base.clone();
        respecified.elements[0].paneling = PanelingParams::uniform();

        assert_eq!(
            base.clearance().pairs[0].min_distance,
            respecified.clearance().pairs[0].min_distance
        );
    }

    // --- Real geometry: McDonnell Douglas 30P-30N -------------------------

    fn repo_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    /// Read a multi-block Selig/XFOIL `.dat` file into one element per block.
    ///
    /// Blank lines, comment lines and `999.0` sentinels separate blocks. A
    /// deliberately minimal reader: the production import lives in the CLI and
    /// the UI, and this only has to open the two fixtures below.
    fn load_configuration(relative: &str) -> Configuration {
        let path = repo_root().join(relative);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));

        let mut blocks: Vec<Vec<Point>> = Vec::new();
        let mut current: Vec<Point> = Vec::new();
        for line in text.lines() {
            let trimmed = line.trim();
            let coordinate = if trimmed.is_empty() || trimmed.starts_with('#') {
                None
            } else {
                let mut parts = trimmed.split_whitespace();
                match (
                    parts.next().and_then(|s| s.parse::<f64>().ok()),
                    parts.next().and_then(|s| s.parse::<f64>().ok()),
                ) {
                    (Some(x), Some(y)) if x < 999.0 && y < 999.0 => Some(point(x, y)),
                    _ => None,
                }
            };

            match coordinate {
                Some(p) => current.push(p),
                None => {
                    if current.len() >= 3 {
                        blocks.push(std::mem::take(&mut current));
                    } else {
                        current.clear();
                    }
                }
            }
        }
        if current.len() >= 3 {
            blocks.push(current);
        }

        let elements = blocks
            .iter()
            .enumerate()
            .map(|(i, pts)| {
                let name = ["slat", "main", "flap"]
                    .get(i)
                    .copied()
                    .unwrap_or("element");
                Element::new(
                    Body::from_points(name, pts).unwrap(),
                    Placement::identity(),
                    PanelingParams::default(),
                    name,
                )
            })
            .collect();

        // The published 30P-30N rigging table is quoted against the retracted
        // chord, which these coordinates are normalised to — not against the
        // main element's own chord, which is the D4 default. Setting it
        // explicitly is what makes the fractions below comparable with the
        // table.
        Configuration::new(elements).with_ref_chord(1.0)
    }

    /// The full 201/221/242-point section, coordinates as shipped in the
    /// airfoil library.
    #[test]
    fn real_30p_30n_reproduces_its_published_rigging_table() {
        let config = load_configuration("flexfoil-ui/public/airfoils/30p-30n.dat");
        assert_eq!(config.len(), 3, "slat, main, flap");

        let report = config.clearance();

        // Geometrically valid: nothing crosses, nothing is contained.
        assert_eq!(report.intersections().count(), 0);
        for pair in &report.pairs {
            assert_eq!(pair.relation, ContourRelation::Separate, "{pair:?}");
        }

        // And inside the regime: the tightest clearance is 1.27% of chord,
        // about 2.5 times the 0.5% floor.
        assert!(report.passes(), "{report}");
        assert_relative_eq!(
            report.min_clearance_fraction().unwrap(),
            0.0127,
            epsilon = 1e-4
        );

        // The file's own header quotes, for the slat and the flap:
        //   Gap        2.95       1.27      (% retracted chord)
        //   Overhang  -2.50       0.25
        // Both are reproduced, the overhang sign included. Note that the
        // published overhang is this module's overlap: negative for the
        // deployed slat, positive for the flap tucked under the main element.
        let slat = report.rigging_between(0, 1).unwrap();
        assert_eq!((slat.upstream, slat.downstream), (0, 1));
        assert_relative_eq!(slat.gap_fraction, 0.0295, epsilon = 1e-4);
        assert_relative_eq!(slat.overlap_fraction, -0.0250, epsilon = 1e-4);

        let flap = report.rigging_between(1, 2).unwrap();
        assert_eq!((flap.upstream, flap.downstream), (1, 2));
        assert_relative_eq!(flap.gap_fraction, 0.0127, epsilon = 1e-4);
        assert_relative_eq!(flap.overlap_fraction, 0.0025, epsilon = 1e-4);

        // The slat and the flap are a chord apart and not an adjacent pair.
        assert!(report.rigging_between(0, 2).is_none());
        assert!(report.pair(0, 2).unwrap().min_distance > 0.8);
    }

    /// The trimmed fixture is a coarse subsample of the same section — about a
    /// tenth of the points. Its rigging is the same to within the resolution it
    /// has left, and it is where the difference between measuring segments and
    /// measuring nodes shows up on real coordinates.
    #[test]
    fn real_30p_30n_trimmed_is_measured_between_segments() {
        let config = load_configuration("testdata/mda_30p_30n_trimmed.dat");
        assert_eq!(config.len(), 3);

        let report = config.clearance();
        assert!(report.passes(), "{report}");
        assert_eq!(report.intersections().count(), 0);

        // Same rigging as the full section, to the accuracy a tenth of the
        // points supports.
        let slat = report.rigging_between(0, 1).unwrap();
        assert_relative_eq!(slat.gap_fraction, 0.0296, epsilon = 5e-4);
        assert_relative_eq!(slat.overlap_fraction, -0.0251, epsilon = 5e-4);

        let flap = report.rigging_between(1, 2).unwrap();
        assert_relative_eq!(flap.gap_fraction, 0.0138, epsilon = 5e-4);
        assert_relative_eq!(flap.overlap_fraction, 0.0024, epsilon = 5e-4);

        // At this resolution the main-to-flap clearance measured node to node
        // is 6% larger than the true surface-to-surface distance — an error in
        // the direction that reports a configuration as clearer than it is.
        let main_to_flap = report.pair(1, 2).unwrap().min_distance;
        let nodes = node_to_node(&config.elements[1], &config.elements[2]);
        assert!(
            nodes > main_to_flap * 1.05,
            "node-to-node {nodes}, segment-to-segment {main_to_flap}"
        );
    }
}
