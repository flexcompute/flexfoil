//! Node numbering for a multi-element configuration.
//!
//! A configuration's elements are solved in one linear system, so their nodes
//! are concatenated into one array: element 0's nodes, then element 1's, and so
//! on in configuration order. A [`Layout`] is the table that says which stretch
//! of that array belongs to which element, and — the part that matters — how a
//! panel closes.
//!
//! # Why the closure rule needs a table
//! A single-element panel method closes each panel with `jp = (jo + 1) % n`:
//! node `n - 1` connects back to node `0`, which is exactly the assumption that
//! the node array is one closed contour. Concatenate three elements into that
//! array and the same expression connects element 0's last node to element 1's
//! first node, adding a panel that spans the gap between two separate bodies.
//! The resulting influence matrix is well-formed and solvable, so the failure
//! shows up as a plausible-looking pressure distribution rather than as an
//! error.
//!
//! [`Layout::next_node`] is the replacement: it wraps within the owning
//! element's span, so element `k`'s last node closes onto element `k`'s first
//! node and never onto element `k + 1`. For a single-element layout it reduces
//! to `(global + 1) % n` exactly, which is the backward-compatibility guarantee
//! the existing kernels rely on.
//!
//! The kernels themselves still use the modulo form
//! (`rustfoil-inviscid/src/influence.rs`, `rustfoil-inviscid/src/velocity.rs`);
//! converting them is solver work, not geometry work.
//!
//! # Index vocabulary
//! Two kinds of index appear here and they are not interchangeable:
//!
//! - a **global** index numbers a node in the concatenated array,
//!   `0 .. total_nodes`;
//! - a **local** index numbers a node within one element, `0 .. span.len`.
//!
//! [`ElementSpan`]'s `le`, `te_upper` and `te_lower` fields are **local**. See
//! [`ElementSpan`] for why.

use crate::configuration::Configuration;
use crate::error::GeometryError;

/// The stretch of the global node array belonging to one element, and that
/// element's landmark nodes.
///
/// # `le`, `te_upper` and `te_lower` are element-local
/// They are indices into this element's own nodes, `0 .. len`, not into the
/// global array. That choice keeps them meaningful on their own: an element's
/// leading edge is a property of the element, and reordering a configuration or
/// repaneling a neighbour must not silently invalidate it. Convert with
/// [`ElementSpan::global`] or the named helpers
/// ([`le_global`](ElementSpan::le_global) and friends) — never by adding
/// `start` by hand at the call site, which is where the two conventions get
/// mixed up.
///
/// # Node ordering within a span
/// Nodes run in Selig/XFOIL order: upper-surface trailing edge → leading edge →
/// lower-surface trailing edge. So `te_upper` is 0, `le` is somewhere in the
/// middle, and `te_lower` is either the last node (blunt trailing edge) or 0
/// again (sharp trailing edge, where the two trailing-edge nodes are one node).
///
/// That direction is set by the geometry the crate actually produces and
/// consumes — [`crate::naca::naca4`] emits "upper surface (TE→LE), then lower
/// surface (LE→TE)", and the surface extraction in rustfoil-solver
/// (`viscous::setup`) documents the same order for Selig-format input. Note
/// that [`crate::body`]'s module documentation describes the opposite walk; it
/// does not match the emitted data, and this convention follows the data.
///
/// # These are node indices, not panel indices
/// [`crate::body::Body`]'s `te_upper_index` and `te_lower_index` are **panel**
/// indices — the two panels that meet at the trailing edge, which is what the
/// Kutta condition is written in terms of. The fields here are **node**
/// indices. The two do not coincide: for a sharp trailing edge, `Body`'s upper
/// trailing-edge panel is the last panel of the contour, whose *end* node is
/// node 0. [`Layout::from_configuration`] does that conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ElementSpan {
    /// Global index of this element's first node.
    pub start: usize,

    /// Number of nodes belonging to this element.
    pub len: usize,

    /// Local index of the leading-edge node.
    ///
    /// `0` means *not located*: node 0 is a trailing-edge node by the ordering
    /// convention, so a real leading edge is never at local index 0. Spans
    /// built by [`Layout::from_node_counts`] carry `le == 0` because a node
    /// count says nothing about where the leading edge is;
    /// [`Layout::from_configuration`] fills it in from the element's geometry.
    /// [`ElementSpan::has_le`] reports which case a span is in.
    pub le: usize,

    /// Local index of the upper-surface trailing-edge node.
    ///
    /// `0` under the ordering convention: the contour walk starts at the upper
    /// trailing edge.
    pub te_upper: usize,

    /// Local index of the lower-surface trailing-edge node.
    ///
    /// `len - 1` for a blunt trailing edge, whose upper and lower trailing-edge
    /// nodes are distinct. `0`, equal to `te_upper`, for a sharp trailing edge,
    /// where one node serves both surfaces — so there is no base panel to
    /// close. [`ElementSpan::has_blunt_te`] distinguishes the two.
    pub te_lower: usize,
}

impl ElementSpan {
    /// A connectivity-only span of `len` nodes starting at `start`.
    ///
    /// `start` and `len` are all the index arithmetic and
    /// [`Layout::next_node`] need. The landmark fields are set to their
    /// "nothing known" values, because a node count carries no geometry:
    /// `le` is `0` (not located) and `te_upper` is `0` (no separate upper
    /// trailing-edge node). Use [`with_indices`](Self::with_indices), or build
    /// the layout with [`Layout::from_configuration`], when the landmarks
    /// matter.
    ///
    /// # Errors
    /// `InvalidParameter { name: "element_span_len" }` if `len` is zero: an
    /// element with no nodes has no panels and no trailing edge, and would make
    /// [`Layout::next_node`] meaningless for it.
    pub fn new(start: usize, len: usize) -> Result<Self, GeometryError> {
        if len == 0 {
            return Err(GeometryError::InvalidParameter {
                name: "element_span_len",
                value: 0.0,
            });
        }
        Ok(Self {
            start,
            len,
            le: 0,
            te_upper: 0,
            te_lower: 0,
        })
    }

    /// A span with all landmark indices stated explicitly, in **local**
    /// coordinates.
    ///
    /// # Errors
    /// - `InvalidParameter { name: "element_span_len" }` if `len` is zero.
    /// - `InvalidParameter { name: "element_span_le" | "element_span_te_upper"
    ///   | "element_span_te_lower" }` if that local index is not less than
    ///   `len`. Passing a global index by mistake is the usual cause.
    pub fn with_indices(
        start: usize,
        len: usize,
        le: usize,
        te_upper: usize,
        te_lower: usize,
    ) -> Result<Self, GeometryError> {
        if len == 0 {
            return Err(GeometryError::InvalidParameter {
                name: "element_span_len",
                value: 0.0,
            });
        }
        for (name, local) in [
            ("element_span_le", le),
            ("element_span_te_upper", te_upper),
            ("element_span_te_lower", te_lower),
        ] {
            if local >= len {
                return Err(GeometryError::InvalidParameter {
                    name,
                    value: local as f64,
                });
            }
        }
        Ok(Self {
            start,
            len,
            le,
            te_upper,
            te_lower,
        })
    }

    /// Global index one past this element's last node.
    #[inline]
    pub fn end(&self) -> usize {
        self.start + self.len
    }

    /// Whether `global` is one of this element's nodes.
    #[inline]
    pub fn contains(&self, global: usize) -> bool {
        global >= self.start && global < self.end()
    }

    /// Convert a local index in this span to a global index.
    ///
    /// # Panics
    /// If `local >= self.len`.
    #[inline]
    pub fn global(&self, local: usize) -> usize {
        assert!(
            local < self.len,
            "local node index {local} is out of range for a span of {} nodes",
            self.len
        );
        self.start + local
    }

    /// Convert a global index in this span to a local index.
    ///
    /// # Panics
    /// If `global` does not belong to this span.
    #[inline]
    pub fn local(&self, global: usize) -> usize {
        assert!(
            self.contains(global),
            "global node index {global} is outside the span [{}, {})",
            self.start,
            self.end()
        );
        global - self.start
    }

    /// Whether a leading edge has been located for this element.
    ///
    /// `false` for a span built from a node count alone; see the `le` field.
    #[inline]
    pub fn has_le(&self) -> bool {
        self.le != 0
    }

    /// Whether this element's trailing edge is blunt, i.e. whether its upper
    /// and lower trailing-edge nodes are distinct.
    ///
    /// `false` for a sharp trailing edge, and also for a connectivity-only span
    /// built from a node count alone.
    #[inline]
    pub fn has_blunt_te(&self) -> bool {
        self.te_upper != self.te_lower
    }

    /// Global index of the leading-edge node.
    #[inline]
    pub fn le_global(&self) -> usize {
        self.start + self.le
    }

    /// Global index of the upper-surface trailing-edge node.
    #[inline]
    pub fn te_upper_global(&self) -> usize {
        self.start + self.te_upper
    }

    /// Global index of the lower-surface trailing-edge node.
    #[inline]
    pub fn te_lower_global(&self) -> usize {
        self.start + self.te_lower
    }
}

/// The node numbering of a whole configuration: one [`ElementSpan`] per
/// element, in configuration order.
///
/// # Invariants, checked at construction
/// - There is at least as much structure as the spans claim: the first span
///   starts at 0, each subsequent span starts where the previous one ended, and
///   `total_nodes` is the sum of the lengths. The spans therefore tile
///   `0 .. total_nodes` exactly, with no gap and no overlap.
/// - Every span has at least one node.
///
/// Those invariants are what let [`element_of`](Layout::element_of) locate a
/// node by binary search, and what make [`next_node`](Layout::next_node) total
/// on `0 .. total_nodes`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Layout {
    /// One span per element, in configuration order, tiling `0 .. total_nodes`.
    pub spans: Vec<ElementSpan>,

    /// Total number of nodes across all elements.
    pub total_nodes: usize,
}

impl Layout {
    /// Build a layout from spans that already carry their `start` offsets.
    ///
    /// # Errors
    /// `InvalidParameter { name: "layout_spans" }` if the spans do not tile
    /// `0 .. total`: the first must start at 0 and each subsequent one must
    /// start where the previous ended. The reported value is the index of the
    /// first span that breaks the chain.
    pub fn from_spans(spans: Vec<ElementSpan>) -> Result<Self, GeometryError> {
        let mut expected_start = 0usize;
        for (i, span) in spans.iter().enumerate() {
            if span.len == 0 {
                return Err(GeometryError::InvalidParameter {
                    name: "element_span_len",
                    value: 0.0,
                });
            }
            if span.start != expected_start {
                return Err(GeometryError::InvalidParameter {
                    name: "layout_spans",
                    value: i as f64,
                });
            }
            if span.le >= span.len || span.te_upper >= span.len || span.te_lower >= span.len {
                return Err(GeometryError::InvalidParameter {
                    name: "layout_spans",
                    value: i as f64,
                });
            }
            expected_start = span.end();
        }
        Ok(Self {
            spans,
            total_nodes: expected_start,
        })
    }

    /// Build a layout from per-element node counts.
    ///
    /// Connectivity only: the spans get their trailing-edge indices from the
    /// ordering convention and no leading edge, because a node count does not
    /// say where a leading edge is. Use [`from_configuration`](Self::from_configuration)
    /// when the leading edges matter.
    ///
    /// # Errors
    /// `InvalidParameter { name: "element_span_len" }` if any count is zero.
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::layout::Layout;
    ///
    /// let layout = Layout::from_node_counts(&[4, 6, 5]).unwrap();
    /// assert_eq!(layout.total_nodes, 15);
    ///
    /// // Element 0 closes onto itself, not onto element 1's first node.
    /// assert_eq!(layout.next_node(3), 0);
    /// // Element 1 spans 4..10 and closes the same way.
    /// assert_eq!(layout.next_node(9), 4);
    /// ```
    pub fn from_node_counts(counts: &[usize]) -> Result<Self, GeometryError> {
        let mut spans = Vec::with_capacity(counts.len());
        let mut start = 0usize;
        for &len in counts {
            let span = ElementSpan::new(start, len)?;
            start = span.end();
            spans.push(span);
        }
        Ok(Self {
            spans,
            total_nodes: start,
        })
    }

    /// Build a layout from a configuration, taking each element's node count
    /// and landmark indices from its geometry.
    ///
    /// Node counts come from [`crate::body::Body::n_nodes`], so a blunt
    /// trailing edge is counted as the extra distinct node it is.
    ///
    /// # Converting the body's indices
    /// A `Body` numbers its trailing edge in *panels* and its leading edge in
    /// *points*, while a span numbers everything in nodes, so the conversion is
    /// explicit here rather than left to each caller:
    ///
    /// - `te_lower` is `0`. The body's lower trailing-edge panel is panel 0,
    ///   whose first node is node 0.
    /// - `te_upper` is the *end* node of the body's upper trailing-edge panel
    ///   (its last panel): node 0 for a closed contour, where that panel runs
    ///   back to the start, and `len - 1` for an open one.
    /// - `le` is [`crate::body::Body::le_index`] unchanged. It is a point index,
    ///   and point `i` is node `i` for every node of the contour.
    ///
    /// Placement is irrelevant here: moving an element does not change how many
    /// nodes it has, or which of them is its leading edge.
    ///
    /// # Errors
    /// - `InvalidParameter { name: "element_span_len" }` for an element with no
    ///   nodes.
    /// - `InvalidParameter { name: "element_span_le" }` if a body's leading-edge
    ///   point index falls outside its own node range. That needs a contour
    ///   whose minimum-x point is the duplicated closing point, which a
    ///   well-formed airfoil does not have.
    pub fn from_configuration(config: &Configuration) -> Result<Self, GeometryError> {
        let mut spans = Vec::with_capacity(config.len());
        let mut start = 0usize;
        for element in config.iter() {
            let body = &element.body;
            let len = body.n_nodes();
            // Node 0 is the UPPER trailing edge: the contour runs
            // TE(upper) → upper surface → LE → lower surface → TE(lower), which
            // is Selig/XFOIL order. See `crate::naca::naca4` and the surface
            // extraction in rustfoil-solver (`viscous::setup`), which both
            // depend on that direction.
            let te_upper = 0;
            let te_lower = if body.is_closed() {
                0
            } else {
                len.saturating_sub(1)
            };
            let span = ElementSpan::with_indices(start, len, body.le_index(), te_upper, te_lower)?;
            start = span.end();
            spans.push(span);
        }
        Ok(Self {
            spans,
            total_nodes: start,
        })
    }

    /// Number of elements.
    #[inline]
    pub fn n_elements(&self) -> usize {
        self.spans.len()
    }

    /// Total number of nodes across all elements.
    #[inline]
    pub fn total_nodes(&self) -> usize {
        self.total_nodes
    }

    /// The spans, in configuration order.
    #[inline]
    pub fn spans(&self) -> &[ElementSpan] {
        &self.spans
    }

    /// The span of one element.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn span(&self, element: usize) -> &ElementSpan {
        &self.spans[element]
    }

    /// The span of one element, or `None` if out of range.
    #[inline]
    pub fn try_span(&self, element: usize) -> Option<&ElementSpan> {
        self.spans.get(element)
    }

    /// Convert an element-local node index to a global one.
    ///
    /// # Panics
    /// If `element` is out of range, or `local` is outside that element's span.
    #[inline]
    pub fn to_global(&self, element: usize, local: usize) -> usize {
        self.spans[element].global(local)
    }

    /// Convert a global node index to `(element, local)`.
    ///
    /// # Panics
    /// If `global >= self.total_nodes`.
    #[inline]
    pub fn from_global(&self, global: usize) -> (usize, usize) {
        let element = self.element_of(global);
        (element, global - self.spans[element].start)
    }

    /// Which element owns a global node index.
    ///
    /// # Panics
    /// If `global >= self.total_nodes`.
    #[inline]
    pub fn element_of(&self, global: usize) -> usize {
        assert!(
            global < self.total_nodes,
            "global node index {global} is out of range for {} nodes",
            self.total_nodes
        );
        // The spans tile 0..total_nodes in increasing order (checked at
        // construction), so the owner is the last span that starts at or before
        // `global`.
        self.spans.partition_point(|span| span.start <= global) - 1
    }

    /// The node that closes the panel starting at `global` — the next node
    /// **within the same element**.
    ///
    /// The last node of element `k` wraps to the *first* node of element `k`,
    /// never to element `k + 1`. For a single-element layout this is exactly
    /// `(global + 1) % total_nodes`.
    ///
    /// A one-node element is its own successor.
    ///
    /// # Panics
    /// If `global >= self.total_nodes`.
    #[inline]
    pub fn next_node(&self, global: usize) -> usize {
        let span = &self.spans[self.element_of(global)];
        let local = global - span.start;
        if local + 1 == span.len {
            span.start
        } else {
            global + 1
        }
    }

    /// The node before `global` within the same element — the inverse of
    /// [`next_node`](Self::next_node).
    ///
    /// The first node of element `k` wraps to the *last* node of element `k`.
    /// For a single-element layout this is exactly
    /// `(global + total_nodes - 1) % total_nodes`.
    ///
    /// # Panics
    /// If `global >= self.total_nodes`.
    #[inline]
    pub fn prev_node(&self, global: usize) -> usize {
        let span = &self.spans[self.element_of(global)];
        if global == span.start {
            span.end() - 1
        } else {
            global - 1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Body;
    use crate::configuration::{Configuration, Element};
    use crate::point::point;

    /// Three elements of deliberately different sizes.
    fn three_elements() -> Layout {
        Layout::from_node_counts(&[4, 7, 3]).unwrap()
    }

    #[test]
    fn node_counts_build_contiguous_spans() {
        let layout = three_elements();
        assert_eq!(layout.n_elements(), 3);
        assert_eq!(layout.total_nodes(), 14);

        assert_eq!(layout.span(0).start, 0);
        assert_eq!(layout.span(0).len, 4);
        assert_eq!(layout.span(1).start, 4);
        assert_eq!(layout.span(1).len, 7);
        assert_eq!(layout.span(2).start, 11);
        assert_eq!(layout.span(2).len, 3);
        assert_eq!(layout.span(2).end(), 14);
    }

    #[test]
    fn node_counts_reject_an_empty_element() {
        assert_eq!(
            Layout::from_node_counts(&[4, 0, 3]),
            Err(GeometryError::InvalidParameter {
                name: "element_span_len",
                value: 0.0,
            })
        );
    }

    #[test]
    fn an_empty_layout_is_allowed_and_has_no_nodes() {
        let layout = Layout::from_node_counts(&[]).unwrap();
        assert_eq!(layout.n_elements(), 0);
        assert_eq!(layout.total_nodes(), 0);
    }

    // --- next_node: the closure rule -------------------------------------

    #[test]
    fn next_node_wraps_within_every_element() {
        let layout = three_elements();

        // Element 0: nodes 0..4.
        assert_eq!(layout.next_node(0), 1);
        assert_eq!(layout.next_node(1), 2);
        assert_eq!(layout.next_node(2), 3);
        assert_eq!(layout.next_node(3), 0, "element 0 must close onto itself");

        // Element 1: nodes 4..11.
        assert_eq!(layout.next_node(4), 5);
        assert_eq!(layout.next_node(9), 10);
        assert_eq!(layout.next_node(10), 4, "element 1 must close onto itself");

        // Element 2: nodes 11..14.
        assert_eq!(layout.next_node(11), 12);
        assert_eq!(layout.next_node(12), 13);
        assert_eq!(layout.next_node(13), 11, "element 2 must close onto itself");
    }

    #[test]
    fn next_node_never_crosses_an_element_boundary() {
        // The property the modulo-n form gets wrong: no successor may leave the
        // element it started in.
        let layout = three_elements();
        for global in 0..layout.total_nodes() {
            let owner = layout.element_of(global);
            let successor = layout.next_node(global);
            assert_eq!(
                layout.element_of(successor),
                owner,
                "next_node({global}) left element {owner}"
            );
        }
    }

    #[test]
    fn next_node_differs_from_modulo_n_exactly_at_the_seams() {
        let layout = three_elements();
        let n = layout.total_nodes();

        // The last node of every element is where plain modulo-n would connect
        // two separate bodies. Everywhere else the two agree.
        let seams = [3usize, 10, 13];
        for global in 0..n {
            let modulo = (global + 1) % n;
            let ours = layout.next_node(global);
            if seams.contains(&global) {
                assert_ne!(ours, modulo, "expected a difference at seam {global}");
            } else {
                assert_eq!(ours, modulo, "unexpected difference at {global}");
            }
        }
    }

    #[test]
    fn single_element_layout_is_plain_modulo_n() {
        // The backward-compatibility guarantee: with one element, next_node is
        // the modulo-n closure the existing kernels already use, so converting
        // them cannot change a single-element result.
        for n in [1usize, 2, 3, 8, 61, 160] {
            let layout = Layout::from_node_counts(&[n]).unwrap();
            for global in 0..n {
                assert_eq!(
                    layout.next_node(global),
                    (global + 1) % n,
                    "n = {n}, global = {global}"
                );
                assert_eq!(
                    layout.prev_node(global),
                    (global + n - 1) % n,
                    "n = {n}, global = {global}"
                );
            }
        }
    }

    #[test]
    fn a_one_node_element_is_its_own_successor() {
        // Degenerate but well-defined: a single node closes onto itself rather
        // than onto a neighbouring element.
        let layout = Layout::from_node_counts(&[3, 1, 2]).unwrap();
        assert_eq!(layout.next_node(3), 3);
        assert_eq!(layout.prev_node(3), 3);
        assert_eq!(layout.element_of(3), 1);
    }

    #[test]
    fn prev_node_inverts_next_node() {
        let layout = three_elements();
        for global in 0..layout.total_nodes() {
            assert_eq!(layout.prev_node(layout.next_node(global)), global);
            assert_eq!(layout.next_node(layout.prev_node(global)), global);
        }
    }

    #[test]
    fn walking_next_node_visits_one_element_and_returns() {
        let layout = three_elements();
        for element in 0..layout.n_elements() {
            let span = *layout.span(element);
            let mut seen = Vec::new();
            let mut node = span.start;
            for _ in 0..span.len {
                seen.push(node);
                node = layout.next_node(node);
            }
            // Back where we started, having visited every node of this element
            // exactly once.
            assert_eq!(node, span.start);
            let mut expected: Vec<usize> = (span.start..span.end()).collect();
            seen.sort_unstable();
            expected.sort_unstable();
            assert_eq!(seen, expected);
        }
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn next_node_rejects_an_out_of_range_index() {
        let layout = three_elements();
        layout.next_node(layout.total_nodes());
    }

    // --- index conversion ------------------------------------------------

    #[test]
    fn global_and_local_round_trip_over_every_node() {
        let layout = three_elements();
        for global in 0..layout.total_nodes() {
            let (element, local) = layout.from_global(global);
            assert!(local < layout.span(element).len);
            assert_eq!(layout.to_global(element, local), global);
        }
    }

    #[test]
    fn local_and_global_round_trip_over_every_element() {
        let layout = three_elements();
        for element in 0..layout.n_elements() {
            for local in 0..layout.span(element).len {
                let global = layout.to_global(element, local);
                assert_eq!(layout.from_global(global), (element, local));
            }
        }
    }

    #[test]
    fn element_of_matches_the_span_boundaries() {
        let layout = three_elements();
        let expected = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2];
        for (global, &element) in expected.iter().enumerate() {
            assert_eq!(layout.element_of(global), element, "global = {global}");
            assert!(layout.span(element).contains(global));
        }
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn element_of_rejects_an_out_of_range_index() {
        three_elements().element_of(14);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn to_global_rejects_an_out_of_range_local_index() {
        // Element 0 has 4 nodes, so local index 4 does not exist — even though
        // global index 4 does (it is element 1's first node).
        three_elements().to_global(0, 4);
    }

    // --- landmark indices ------------------------------------------------

    #[test]
    fn node_count_spans_carry_no_landmarks() {
        // A node count says nothing about geometry, so neither landmark is
        // claimed rather than being guessed at.
        let layout = three_elements();
        let span = layout.span(1);
        assert_eq!(span.te_lower, 0);
        assert_eq!(span.te_upper, 0);
        assert!(!span.has_le());
        assert!(!span.has_blunt_te());
        // The last node of the element is still available from the span itself.
        assert_eq!(span.end() - 1, 10);
    }

    #[test]
    fn landmark_helpers_convert_local_to_global() {
        let span = ElementSpan::with_indices(11, 5, 2, 4, 0).unwrap();
        assert!(span.has_le());
        assert!(span.has_blunt_te());
        assert_eq!(span.le_global(), 13);
        assert_eq!(span.te_upper_global(), 15);
        assert_eq!(span.te_lower_global(), 11);
        assert_eq!(span.global(3), 14);
        assert_eq!(span.local(14), 3);
    }

    #[test]
    fn with_indices_rejects_a_local_index_outside_the_span() {
        // The usual mistake: passing a global index where a local one is wanted.
        assert_eq!(
            ElementSpan::with_indices(11, 5, 13, 4, 0),
            Err(GeometryError::InvalidParameter {
                name: "element_span_le",
                value: 13.0,
            })
        );
        assert!(ElementSpan::with_indices(0, 5, 2, 5, 0).is_err());
        assert!(ElementSpan::with_indices(0, 5, 2, 4, 9).is_err());
    }

    #[test]
    fn spans_must_tile_the_node_array() {
        let good = vec![
            ElementSpan::new(0, 4).unwrap(),
            ElementSpan::new(4, 3).unwrap(),
        ];
        assert_eq!(Layout::from_spans(good).unwrap().total_nodes(), 7);

        // A gap between the spans.
        let gapped = vec![
            ElementSpan::new(0, 4).unwrap(),
            ElementSpan::new(5, 3).unwrap(),
        ];
        assert_eq!(
            Layout::from_spans(gapped),
            Err(GeometryError::InvalidParameter {
                name: "layout_spans",
                value: 1.0,
            })
        );

        // An overlap.
        let overlapping = vec![
            ElementSpan::new(0, 4).unwrap(),
            ElementSpan::new(3, 3).unwrap(),
        ];
        assert!(Layout::from_spans(overlapping).is_err());

        // Not starting at zero.
        let offset = vec![ElementSpan::new(1, 4).unwrap()];
        assert!(Layout::from_spans(offset).is_err());
    }

    // --- from a configuration --------------------------------------------

    /// A closed diamond of the given chord: 5 points, 4 panels, 4 nodes.
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

    #[test]
    fn configuration_layout_takes_counts_and_landmarks_from_geometry() {
        let config = Configuration::new(vec![
            Element::from_body(diamond("slat", 0.2)),
            Element::from_body(diamond("main", 1.0)),
        ]);
        let layout = Layout::from_configuration(&config).unwrap();

        assert_eq!(layout.n_elements(), 2);
        assert_eq!(layout.total_nodes(), 8);

        // The diamond's leading edge is its third point, so local index 2. Its
        // trailing edge is sharp, so one node serves both surfaces.
        for element in 0..2 {
            let span = layout.span(element);
            assert_eq!(span.len, 4);
            assert_eq!(span.le, 2);
            assert_eq!(span.te_lower, 0);
            assert_eq!(span.te_upper, 0);
            assert!(span.has_le());
            assert!(!span.has_blunt_te());
        }

        // And the closure rule still holds per element.
        assert_eq!(layout.next_node(3), 0);
        assert_eq!(layout.next_node(7), 4);
    }

    #[test]
    fn configuration_layout_counts_a_blunt_trailing_edge_node() {
        // An open contour has one more distinct node than it has panels.
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

        let config = Configuration::new(vec![
            Element::from_body(diamond("sharp", 1.0)),
            Element::from_body(blunt),
        ]);
        let layout = Layout::from_configuration(&config).unwrap();

        assert_eq!(layout.span(0).len, 4);
        assert_eq!(layout.span(1).len, 5);
        assert_eq!(layout.total_nodes(), 9);
        assert_eq!(layout.next_node(8), 4);

        // The sharp element has one trailing-edge node; the blunt one has two.
        // Nodes run TE(upper) → LE → TE(lower), so the upper trailing edge is
        // the element's first node and the lower one is its last.
        assert!(!layout.span(0).has_blunt_te());
        assert!(layout.span(1).has_blunt_te());
        assert_eq!(layout.span(1).te_upper, 0);
        assert_eq!(layout.span(1).te_lower, 4);
        assert_eq!(layout.span(1).te_upper_global(), 4);
        assert_eq!(layout.span(1).te_lower_global(), 8);
    }

    #[test]
    fn span_trailing_edge_is_a_node_index_not_a_panel_index() {
        // The conversion that is easy to get wrong. For the sharp diamond,
        // Body reports its upper trailing-edge *panel* as panel 3, but that
        // panel ends on node 0 — so the span's upper trailing-edge *node* is 0,
        // not 3. Copying the panel index across would name an ordinary
        // upper-surface node as the trailing edge.
        let body = diamond("main", 1.0);
        assert_eq!(body.te_upper_index(), 3, "panel index");

        let layout = Layout::from_configuration(&Configuration::single(body)).unwrap();
        assert_eq!(layout.span(0).te_upper, 0, "node index");
    }

    #[test]
    fn configuration_layout_is_unaffected_by_placement() {
        let unplaced = Configuration::new(vec![Element::from_body(diamond("flap", 0.3))]);
        let mut placed = unplaced.clone();
        placed.elements[0].placement =
            crate::placement::Placement::rotation_about(point(0.0, 0.0), -30.0);

        assert_eq!(
            Layout::from_configuration(&unplaced).unwrap(),
            Layout::from_configuration(&placed).unwrap()
        );
    }

    #[test]
    fn an_empty_configuration_gives_an_empty_layout() {
        let layout = Layout::from_configuration(&Configuration::new(vec![])).unwrap();
        assert_eq!(layout.n_elements(), 0);
        assert_eq!(layout.total_nodes(), 0);
    }

    #[test]
    fn closure_holds_at_the_seams_of_a_real_three_element_geometry() {
        // Node counts of the McDonnell Douglas 30P-30N slat/main/flap in
        // flexfoil-ui/public/airfoils/30p-30n.dat: blocks of 201, 221 and 242
        // points, of which the first two are closed contours (one point fewer
        // than that in distinct nodes) and the flap is a blunt trailing edge
        // (all 242 points are distinct nodes).
        let layout = Layout::from_node_counts(&[200, 220, 242]).unwrap();
        assert_eq!(layout.total_nodes(), 662);

        // Each element's last node closes onto its own first node. Plain
        // modulo-n would instead join the slat to the main, and the main to the
        // flap, with panels spanning the slots between them.
        assert_eq!(layout.next_node(199), 0);
        assert_eq!(layout.next_node(419), 200);
        assert_eq!(layout.next_node(661), 420);

        // And no successor anywhere crosses into another element.
        for global in 0..layout.total_nodes() {
            assert_eq!(
                layout.element_of(layout.next_node(global)),
                layout.element_of(global)
            );
        }
    }

    #[test]
    fn single_element_configuration_layout_is_plain_modulo_n() {
        let config = Configuration::single(diamond("main", 1.0));
        let layout = Layout::from_configuration(&config).unwrap();
        let n = layout.total_nodes();
        for global in 0..n {
            assert_eq!(layout.next_node(global), (global + 1) % n);
        }
    }
}
