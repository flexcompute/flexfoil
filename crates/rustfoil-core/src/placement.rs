//! Rigid placement of one element inside a multi-element configuration.
//!
//! A [`Placement`] says *where an element sits*, separately from *what shape it
//! is*. The element's coordinates stay as imported; the placement is applied
//! when the configuration is paneled and solved, and is never baked back into
//! the stored coordinates. That split is what lets gap, overlap and deflection
//! be swept without re-importing geometry: the sweep varies the placement and
//! leaves the contour untouched.
//!
//! # Relationship to `flap`
//! [`crate::flap`] is XFOIL's plain flap: it rotates a *region of one contour*
//! about a hinge and the result is still one body. A `Placement` moves a whole
//! element and leaves its shape alone. The two are independent.
//!
//! # Cross-language agreement
//! The same four fields, in the same order, are the placement half of the
//! geometry cache key in `packages/flexfoil-python/src/flexfoil/airfoil.py`
//! (`Placement`) and `flexfoil-ui/src/lib/airfoilHash.ts` (`Placement`).
//! [`Placement::is_identity`] deliberately uses the same exact-comparison rule
//! as those two, so a placement that hashes as absent also reports as identity
//! here.

use crate::point::{point, vec2, Point, Vec2};

/// Rigid-body placement of one element: uniform scale and rotation about a
/// pivot, then a translation.
///
/// # Apply order
/// [`apply`](Placement::apply) evaluates, in this order:
///
/// 1. **scale** about `pivot` — `pivot + scale * (p - pivot)`
/// 2. **rotate** about `pivot` by `rotation_deg`, counter-clockwise positive
/// 3. **translate** by `translation`
///
/// so that
///
/// ```text
/// apply(p) = pivot + R(rotation_deg) * (scale * (p - pivot)) + translation
/// ```
///
/// The order matters and is fixed: scaling and rotating both happen about
/// `pivot` (and therefore commute with each other), but the translation is
/// applied last and is *not* scaled or rotated. A placement built to put a flap
/// at a given gap and overlap describes that gap in configuration coordinates,
/// not in the flap's own pre-scale coordinates.
///
/// `pivot` is expressed in the element's own (unplaced) coordinates.
/// `translation` and the result of `apply` are in configuration coordinates.
///
/// # Sign convention
/// `rotation_deg` is counter-clockwise positive, consistent with the crate's
/// coordinate convention (x downstream, y up). A trailing-edge flap deflected
/// *down* — the usual high-lift sense — therefore has a **negative**
/// `rotation_deg`.
///
/// # Example
/// ```
/// use rustfoil_core::placement::Placement;
/// use rustfoil_core::point::{point, vec2};
///
/// // A flap hinged at (0.7, 0.0), deflected 30° down and shifted aft.
/// let p = Placement {
///     pivot: point(0.7, 0.0),
///     rotation_deg: -30.0,
///     translation: vec2(0.02, -0.01),
///     scale: 1.0,
/// };
///
/// // The pivot itself only sees the translation.
/// let moved = p.apply(point(0.7, 0.0));
/// assert!((moved.x - 0.72).abs() < 1e-12);
/// assert!((moved.y + 0.01).abs() < 1e-12);
///
/// // `inverse` undoes `apply`.
/// let there_and_back = p.inverse(p.apply(point(1.0, 0.05)));
/// assert!((there_and_back.x - 1.0).abs() < 1e-12);
/// assert!((there_and_back.y - 0.05).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Placement {
    /// Centre of the scale and rotation, in the element's own coordinates.
    ///
    /// A pivot on its own moves nothing: with zero rotation, unit scale and no
    /// translation the placement is the identity whatever the pivot is.
    pub pivot: Point,

    /// Rotation about `pivot` in **degrees**, counter-clockwise positive.
    pub rotation_deg: f64,

    /// Translation applied after the scale and rotation, in configuration
    /// coordinates.
    pub translation: Vec2,

    /// Uniform scale about `pivot`. `1.0` leaves the element's size alone.
    pub scale: f64,
}

impl Default for Placement {
    /// The identity placement — see [`Placement::identity`].
    fn default() -> Self {
        Self::identity()
    }
}

impl Placement {
    /// The identity placement: no scale, no rotation, no translation.
    ///
    /// [`apply`](Placement::apply) returns its argument unchanged (bit for bit,
    /// since the arithmetic reduces to `p + 0.0` on each component).
    #[inline]
    pub fn identity() -> Self {
        Self {
            pivot: point(0.0, 0.0),
            rotation_deg: 0.0,
            translation: vec2(0.0, 0.0),
            scale: 1.0,
        }
    }

    /// A placement that only translates.
    #[inline]
    pub fn from_translation(dx: f64, dy: f64) -> Self {
        Self {
            translation: vec2(dx, dy),
            ..Self::identity()
        }
    }

    /// A placement that only rotates, by `rotation_deg` about `pivot`.
    #[inline]
    pub fn rotation_about(pivot: Point, rotation_deg: f64) -> Self {
        Self {
            pivot,
            rotation_deg,
            ..Self::identity()
        }
    }

    /// True if this placement moves nothing.
    ///
    /// The test is `rotation_deg == 0`, `translation == (0, 0)` and
    /// `scale == 1`; `pivot` is ignored because a pivot on its own is a no-op.
    ///
    /// # Exact comparison, on purpose
    /// The comparisons are exact rather than tolerant, matching
    /// `Placement.is_identity` in `packages/flexfoil-python` and
    /// `isIdentityPlacement` in `flexfoil-ui`. Those two decide whether a
    /// placement contributes to the geometry cache key at all, so a tolerance
    /// here and an exact test there would disagree about which configurations
    /// share a digest.
    ///
    /// One consequence: a rotation of `360.0` degrees is geometrically the
    /// identity but does not report as one. That is the same on all three
    /// sides.
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.rotation_deg == 0.0
            && self.translation.x == 0.0
            && self.translation.y == 0.0
            && self.scale == 1.0
    }

    /// Map a point from the element's own coordinates into configuration
    /// coordinates.
    ///
    /// See the type-level documentation for the apply order.
    #[inline]
    pub fn apply(&self, p: Point) -> Point {
        let (sin_t, cos_t) = self.rotation_deg.to_radians().sin_cos();
        let d = (p - self.pivot) * self.scale;
        point(
            self.pivot.x + cos_t * d.x - sin_t * d.y + self.translation.x,
            self.pivot.y + sin_t * d.x + cos_t * d.y + self.translation.y,
        )
    }

    /// Map a point from configuration coordinates back into the element's own
    /// coordinates.
    ///
    /// The inverse of [`apply`](Placement::apply): undo the translation, then
    /// the rotation, then the scale. `inverse(apply(p)) == p` to floating-point
    /// round-off (exactly, for the identity placement).
    ///
    /// # Degenerate scale
    /// A `scale` of zero collapses the element to a point and has no inverse;
    /// this returns `pivot` in that case rather than dividing by zero. Callers
    /// that care should reject a zero scale when the placement is built.
    #[inline]
    pub fn inverse(&self, p: Point) -> Point {
        let (sin_t, cos_t) = self.rotation_deg.to_radians().sin_cos();

        // Undo the translation, then measure from the pivot.
        let d = vec2(
            p.x - self.translation.x - self.pivot.x,
            p.y - self.translation.y - self.pivot.y,
        );

        // Undo the rotation (transpose of R).
        let rx = cos_t * d.x + sin_t * d.y;
        let ry = -sin_t * d.x + cos_t * d.y;

        if self.scale == 0.0 {
            return self.pivot;
        }

        point(
            self.pivot.x + rx / self.scale,
            self.pivot.y + ry / self.scale,
        )
    }

    /// Map a *direction* (tangent, normal, velocity) into configuration
    /// coordinates.
    ///
    /// Directions transform by the linear part of the placement only — scale
    /// and rotation, no translation and no pivot. Use this rather than
    /// differencing two [`apply`](Placement::apply) results when transforming
    /// panel tangents and normals.
    ///
    /// Note that the result is scaled by `scale`; re-normalise if a unit vector
    /// is required.
    #[inline]
    pub fn apply_direction(&self, v: Vec2) -> Vec2 {
        let (sin_t, cos_t) = self.rotation_deg.to_radians().sin_cos();
        let s = self.scale;
        vec2(
            s * (cos_t * v.x - sin_t * v.y),
            s * (sin_t * v.x + cos_t * v.y),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// A placement exercising all four fields at once.
    fn nontrivial() -> Placement {
        Placement {
            pivot: point(0.7, -0.03),
            rotation_deg: -27.5,
            translation: vec2(0.031, -0.017),
            scale: 0.85,
        }
    }

    #[test]
    fn identity_is_identity() {
        let id = Placement::identity();
        assert!(id.is_identity());
        assert_eq!(Placement::default(), id);
    }

    #[test]
    fn identity_leaves_points_untouched() {
        let id = Placement::identity();
        for p in [point(0.0, 0.0), point(1.0, 0.0), point(0.3, -0.11)] {
            // Exact, not approximate: the identity must not perturb coordinates.
            assert_eq!(id.apply(p), p);
            assert_eq!(id.inverse(p), p);
        }
    }

    #[test]
    fn pivot_alone_is_still_the_identity() {
        // W0's hashing rule: a pivot with no rotation, scale or translation
        // moves nothing and must not read as a real placement.
        let p = Placement {
            pivot: point(0.42, 0.13),
            ..Placement::identity()
        };
        assert!(p.is_identity());
        assert_eq!(p.apply(point(1.0, 0.05)), point(1.0, 0.05));
    }

    #[test]
    fn real_placements_are_not_identity() {
        let id = Placement::identity();

        assert!(!Placement {
            rotation_deg: 1e-9,
            ..id
        }
        .is_identity());
        assert!(!Placement {
            translation: vec2(1e-9, 0.0),
            ..id
        }
        .is_identity());
        assert!(!Placement {
            translation: vec2(0.0, -1e-9),
            ..id
        }
        .is_identity());
        assert!(!Placement {
            scale: 1.0 + 1e-9,
            ..id
        }
        .is_identity());

        // Documented consequence of the exact comparison.
        assert!(!Placement {
            rotation_deg: 360.0,
            ..id
        }
        .is_identity());
    }

    #[test]
    fn translation_only() {
        let p = Placement::from_translation(0.25, -0.1);
        let moved = p.apply(point(1.0, 0.0));
        assert_relative_eq!(moved.x, 1.25);
        assert_relative_eq!(moved.y, -0.1);
    }

    #[test]
    fn rotation_is_counter_clockwise_positive() {
        let p = Placement::rotation_about(point(0.0, 0.0), 90.0);
        let rotated = p.apply(point(1.0, 0.0));
        assert_relative_eq!(rotated.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(rotated.y, 1.0, epsilon = 1e-15);
    }

    #[test]
    fn rotation_leaves_the_pivot_where_it_is() {
        let pivot = point(0.7, -0.03);
        let p = Placement::rotation_about(pivot, -35.0);
        let moved = p.apply(pivot);
        assert_relative_eq!(moved.x, pivot.x, epsilon = 1e-15);
        assert_relative_eq!(moved.y, pivot.y, epsilon = 1e-15);
    }

    #[test]
    fn scale_is_about_the_pivot() {
        let pivot = point(0.5, 0.0);
        let p = Placement {
            pivot,
            scale: 2.0,
            ..Placement::identity()
        };
        // A point one unit downstream of the pivot ends up two units away.
        let moved = p.apply(point(1.5, 0.0));
        assert_relative_eq!(moved.x, 2.5);
        assert_relative_eq!(moved.y, 0.0);
        // The pivot itself does not move.
        assert_relative_eq!(p.apply(pivot).x, pivot.x);
    }

    #[test]
    fn apply_order_is_scale_then_rotate_then_translate() {
        // Hand-computed reference: pivot at the origin, scale 2, rotate 90°,
        // then translate. Under any other order the answer differs, so this
        // pins the documented order rather than just self-consistency.
        let p = Placement {
            pivot: point(0.0, 0.0),
            rotation_deg: 90.0,
            translation: vec2(1.0, 0.0),
            scale: 2.0,
        };
        // (1,0) -> scale -> (2,0) -> rotate 90° -> (0,2) -> translate -> (1,2)
        let moved = p.apply(point(1.0, 0.0));
        assert_relative_eq!(moved.x, 1.0, epsilon = 1e-15);
        assert_relative_eq!(moved.y, 2.0, epsilon = 1e-15);
    }

    #[test]
    fn translation_is_not_scaled_or_rotated() {
        // If the translation were applied first (and so picked up the scale and
        // rotation), the result would be (0, 2) rather than (1, 0).
        let p = Placement {
            pivot: point(0.0, 0.0),
            rotation_deg: 90.0,
            translation: vec2(1.0, 0.0),
            scale: 2.0,
        };
        let moved = p.apply(point(0.0, 0.0));
        assert_relative_eq!(moved.x, 1.0, epsilon = 1e-15);
        assert_relative_eq!(moved.y, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn inverse_round_trips() {
        let p = nontrivial();
        for q in [
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(0.5, 0.06),
            point(0.5, -0.06),
            point(-2.3, 11.7),
        ] {
            let back = p.inverse(p.apply(q));
            assert_relative_eq!(back.x, q.x, epsilon = 1e-12);
            assert_relative_eq!(back.y, q.y, epsilon = 1e-12);
        }
    }

    #[test]
    fn apply_round_trips_the_other_way() {
        let p = nontrivial();
        let q = point(0.31, -0.04);
        let there = p.apply(p.inverse(q));
        assert_relative_eq!(there.x, q.x, epsilon = 1e-12);
        assert_relative_eq!(there.y, q.y, epsilon = 1e-12);
    }

    #[test]
    fn placement_preserves_distance_up_to_scale() {
        let p = nontrivial();
        let a = point(0.0, 0.0);
        let b = point(1.0, 0.0);
        let placed = (p.apply(b) - p.apply(a)).norm();
        assert_relative_eq!(placed, (b - a).norm() * p.scale, epsilon = 1e-12);
    }

    #[test]
    fn directions_ignore_the_translation() {
        let p = nontrivial();
        let a = point(0.2, 0.01);
        let b = point(0.9, -0.02);
        let by_difference = p.apply(b) - p.apply(a);
        let by_direction = p.apply_direction(b - a);
        assert_relative_eq!(by_direction.x, by_difference.x, epsilon = 1e-12);
        assert_relative_eq!(by_direction.y, by_difference.y, epsilon = 1e-12);
    }

    #[test]
    fn zero_scale_inverse_returns_the_pivot() {
        let p = Placement {
            pivot: point(0.3, 0.02),
            scale: 0.0,
            ..Placement::identity()
        };
        assert_eq!(p.inverse(point(9.0, 9.0)), p.pivot);
    }
}
