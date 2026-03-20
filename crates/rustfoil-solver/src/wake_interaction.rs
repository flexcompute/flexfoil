//! Wake-body interaction detection for multi-element configurations (Phase 2).
//!
//! Detects where the wake from an upstream element impinges on a downstream
//! element's surface. At impingement points, the BL transitions from single-layer
//! to confluent (two-layer) mode.
//!
//! # Algorithm
//!
//! 1. For each body, compute a straight-line wake trajectory from the TE
//! 2. For each downstream body, check if the wake intersects its surface
//! 3. If intersection found, record the impingement location and upstream wake state

use rustfoil_core::{Body, Point};

/// Describes where an upstream body's wake hits a downstream body.
#[derive(Debug, Clone)]
pub struct WakeImpingement {
    /// Index of the body that sheds the wake
    pub source_body: usize,
    /// Index of the body that receives the wake
    pub target_body: usize,
    /// Station index on the target body where impingement occurs (arc-length based)
    pub target_station: usize,
    /// x-coordinate of the impingement point
    pub x_impingement: f64,
    /// y-coordinate of the impingement point
    pub y_impingement: f64,
    /// Wake theta at impingement (from upstream body's TE)
    pub wake_theta: f64,
    /// Wake delta_star at impingement
    pub wake_delta_star: f64,
    /// Wake ctau at impingement
    pub wake_ctau: f64,
}

/// Straight-line wake trajectory from a body's TE.
#[derive(Debug, Clone)]
pub struct WakeTrajectory {
    /// Body index
    pub body_idx: usize,
    /// TE midpoint (start of wake)
    pub start: Point,
    /// Wake direction (unit vector)
    pub direction: Point,
    /// Wake length (in chord units from source body)
    pub length: f64,
}

/// Compute wake trajectories for all bodies.
///
/// Each body sheds a straight-line wake from its TE along the bisector
/// of the upper and lower TE tangent directions.
pub fn compute_wake_trajectories(bodies: &[Body]) -> Vec<WakeTrajectory> {
    bodies.iter().enumerate().map(|(idx, body)| {
        let panels = body.panels();
        let n = panels.len();
        let te_upper = panels[0].p1;
        let te_lower = panels[n - 1].p2;

        let te_mid = Point::new(
            0.5 * (te_upper.x + te_lower.x),
            0.5 * (te_upper.y + te_lower.y),
        );

        let dxte = te_lower.x - te_upper.x;
        let dyte = te_lower.y - te_upper.y;

        // Wake direction: perpendicular to TE gap, pointing downstream
        let mut wake_dx = -dyte;
        let mut wake_dy = dxte;
        if wake_dx < 0.0 {
            wake_dx = -wake_dx;
            wake_dy = -wake_dy;
        }

        let len = (wake_dx * wake_dx + wake_dy * wake_dy).sqrt();
        if len > 1e-30 {
            wake_dx /= len;
            wake_dy /= len;
        } else {
            wake_dx = 1.0;
            wake_dy = 0.0;
        }

        WakeTrajectory {
            body_idx: idx,
            start: te_mid,
            direction: Point::new(wake_dx, wake_dy),
            length: 10.0 * body.chord(),
        }
    }).collect()
}

/// Detect wake impingements: where does each body's wake hit other bodies.
///
/// Uses ray-segment intersection to find where the wake ray from body_i
/// crosses the surface panels of body_j.
pub fn detect_wake_impingements(
    bodies: &[Body],
    wakes: &[WakeTrajectory],
) -> Vec<WakeImpingement> {
    let mut impingements = Vec::new();

    for wake in wakes {
        let src = wake.body_idx;

        for (tgt, body) in bodies.iter().enumerate() {
            if tgt == src { continue; }

            let panels = body.panels();
            for (panel_idx, panel) in panels.iter().enumerate() {
                let p1 = panel.p1;
                let p2 = panel.p2;

                // Ray-segment intersection
                if let Some(t) = ray_segment_intersect(
                    wake.start, wake.direction,
                    p1, p2,
                ) {
                    if t > 0.0 && t < wake.length {
                        let x_hit = wake.start.x + t * wake.direction.x;
                        let y_hit = wake.start.y + t * wake.direction.y;

                        impingements.push(WakeImpingement {
                            source_body: src,
                            target_body: tgt,
                            target_station: panel_idx,
                            x_impingement: x_hit,
                            y_impingement: y_hit,
                            wake_theta: 0.0, // will be set from upstream BL solution
                            wake_delta_star: 0.0,
                            wake_ctau: 0.0,
                        });
                        break; // first intersection only
                    }
                }
            }
        }
    }

    impingements
}

/// Ray-segment intersection: find parameter t where ray(t) = origin + t*dir
/// intersects the segment p1-p2. Returns None if no intersection.
fn ray_segment_intersect(
    origin: Point,
    dir: Point,
    p1: Point,
    p2: Point,
) -> Option<f64> {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;

    let denom = dir.x * dy - dir.y * dx;
    if denom.abs() < 1e-15 {
        return None; // parallel
    }

    let t = ((p1.x - origin.x) * dy - (p1.y - origin.y) * dx) / denom;
    let u = ((p1.x - origin.x) * dir.y - (p1.y - origin.y) * dir.x) / denom;

    if t >= 0.0 && u >= 0.0 && u <= 1.0 {
        Some(t)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    #[test]
    fn test_ray_segment_intersect_hit() {
        let origin = point(0.0, 0.0);
        let dir = point(1.0, 0.0);
        let p1 = point(5.0, -1.0);
        let p2 = point(5.0, 1.0);

        let t = ray_segment_intersect(origin, dir, p1, p2);
        assert!(t.is_some());
        let t = t.unwrap();
        assert!((t - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_segment_intersect_miss() {
        let origin = point(0.0, 0.0);
        let dir = point(1.0, 0.0);
        let p1 = point(5.0, 2.0);
        let p2 = point(5.0, 3.0);

        let t = ray_segment_intersect(origin, dir, p1, p2);
        assert!(t.is_none());
    }

    #[test]
    fn test_ray_segment_behind() {
        let origin = point(10.0, 0.0);
        let dir = point(1.0, 0.0);
        let p1 = point(5.0, -1.0);
        let p2 = point(5.0, 1.0);

        let t = ray_segment_intersect(origin, dir, p1, p2);
        // t would be -5, which is < 0, so no intersection
        assert!(t.is_none());
    }
}
