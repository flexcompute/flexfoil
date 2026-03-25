//! Structured O-grid mesh generation for CFD.
//!
//! Generates a body-fitted O-topology mesh around an airfoil by:
//! 1. Distributing points along the airfoil surface (circumferential direction)
//! 2. Algebraically growing layers outward (radial direction)
//! 3. Optionally smoothing with Poisson iteration

/// Result of mesh generation: flat arrays ready for GPU upload.
#[derive(Debug, Clone)]
pub struct MeshResult {
    /// Node x-coordinates, length ni*nj (row-major: j*ni + i)
    pub x: Vec<f32>,
    /// Node y-coordinates, length ni*nj
    pub y: Vec<f32>,
    /// Grid dimensions
    pub ni: u32,
    pub nj: u32,
}

/// Generate an O-type structured mesh around an airfoil.
///
/// # Arguments
/// * `airfoil_x` - Airfoil x-coordinates (CCW from TE lower, around LE, to TE upper)
/// * `airfoil_y` - Airfoil y-coordinates
/// * `ni` - Number of circumferential points (should match airfoil point count or resample)
/// * `nj` - Number of radial layers
/// * `far_field` - Far-field radius in chord lengths
/// * `ds0` - First cell wall spacing (normalized by chord)
///
/// # Returns
/// `MeshResult` with flat coordinate arrays
pub fn generate_o_mesh(
    airfoil_x: &[f64],
    airfoil_y: &[f64],
    ni: u32,
    nj: u32,
    far_field: f64,
    ds0: f64,
) -> MeshResult {
    let ni = ni as usize;
    let nj = nj as usize;
    // Resample airfoil to ni points using linear interpolation along arc length
    let (body_x, body_y) = resample_airfoil(airfoil_x, airfoil_y, ni);

    // Compute airfoil centroid for radial direction
    let cx: f64 = body_x.iter().sum::<f64>() / ni as f64;
    let cy: f64 = body_y.iter().sum::<f64>() / ni as f64;

    // Compute outward normal directions at each surface point
    let mut nx = vec![0.0f64; ni];
    let mut ny = vec![0.0f64; ni];
    for i in 0..ni {
        // Use centroid-to-point direction as approximate outward normal
        let dx = body_x[i] - cx;
        let dy = body_y[i] - cy;
        let mag = (dx * dx + dy * dy).sqrt().max(1e-15);
        nx[i] = dx / mag;
        ny[i] = dy / mag;
    }

    // Generate radial distribution with geometric stretching
    let radial_dist = generate_radial_distribution(nj, ds0, far_field);

    // Build mesh: algebraic growth along normals
    let mut x = vec![0.0f32; ni * nj];
    let mut y = vec![0.0f32; ni * nj];

    for j in 0..nj {
        let r = radial_dist[j];
        for i in 0..ni {
            let idx = j * ni + i;
            x[idx] = (body_x[i] + r * nx[i]) as f32;
            y[idx] = (body_y[i] + r * ny[i]) as f32;
        }
    }

    // Laplacian smoothing passes for interior points
    laplacian_smooth(&mut x, &mut y, ni, nj, 50);

    MeshResult {
        x,
        y,
        ni: ni as u32,
        nj: nj as u32,
    }
}

/// Resample airfoil coordinates to exactly `n` points using arc-length parameterization.
fn resample_airfoil(ax: &[f64], ay: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let m = ax.len();
    if m == n {
        return (ax.to_vec(), ay.to_vec());
    }

    // Compute cumulative arc length
    let mut s = vec![0.0f64; m];
    for i in 1..m {
        let dx = ax[i] - ax[i - 1];
        let dy = ay[i] - ay[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    let total_s = s[m - 1];

    // Interpolate at uniform arc-length intervals
    let mut rx = vec![0.0f64; n];
    let mut ry = vec![0.0f64; n];
    let mut j = 0usize;
    for i in 0..n {
        let target = total_s * (i as f64) / ((n - 1) as f64);
        // Advance to the segment containing target
        while j < m - 2 && s[j + 1] < target {
            j += 1;
        }
        let seg_len = s[j + 1] - s[j];
        let t = if seg_len > 1e-15 {
            (target - s[j]) / seg_len
        } else {
            0.0
        };
        rx[i] = ax[j] + t * (ax[j + 1] - ax[j]);
        ry[i] = ay[j] + t * (ay[j + 1] - ay[j]);
    }

    (rx, ry)
}

/// Generate radial distance distribution with geometric stretching.
///
/// The first cell has height `ds0`, and cells grow geometrically until
/// the total radial extent reaches `far_field`.
fn generate_radial_distribution(nj: usize, ds0: f64, far_field: f64) -> Vec<f64> {
    let mut dist = vec![0.0f64; nj];
    // j=0 is the body surface
    dist[0] = 0.0;

    if nj < 2 {
        return dist;
    }

    // Find geometric stretch ratio using bisection
    // Total distance: ds0 * (ratio^(nj-1) - 1) / (ratio - 1) = far_field
    let ratio = find_stretch_ratio(nj - 1, ds0, far_field);

    let mut ds = ds0;
    for j in 1..nj {
        dist[j] = dist[j - 1] + ds;
        ds *= ratio;
    }

    dist
}

/// Find geometric stretch ratio via bisection.
fn find_stretch_ratio(n: usize, ds0: f64, total: f64) -> f64 {
    let mut lo = 1.0001f64;
    let mut hi = 2.0f64;

    // Ensure hi is large enough
    while geometric_sum(n, ds0, hi) < total {
        hi *= 2.0;
    }

    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let sum = geometric_sum(n, ds0, mid);
        if sum < total {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    0.5 * (lo + hi)
}

fn geometric_sum(n: usize, ds0: f64, ratio: f64) -> f64 {
    if (ratio - 1.0).abs() < 1e-10 {
        return ds0 * n as f64;
    }
    ds0 * (ratio.powi(n as i32) - 1.0) / (ratio - 1.0)
}

/// Simple Laplacian smoothing for interior mesh points.
fn laplacian_smooth(x: &mut [f32], y: &mut [f32], ni: usize, nj: usize, iters: usize) {
    let omega = 0.3f32; // Under-relaxation factor

    for _ in 0..iters {
        for j in 1..nj - 1 {
            for i in 0..ni {
                let idx = j * ni + i;
                // Periodic in i-direction (O-grid wraps around)
                let ip = if i + 1 < ni { i + 1 } else { 0 };
                let im = if i > 0 { i - 1 } else { ni - 1 };

                let avg_x = 0.25
                    * (x[j * ni + ip] + x[j * ni + im] + x[(j + 1) * ni + i] + x[(j - 1) * ni + i]);
                let avg_y = 0.25
                    * (y[j * ni + ip] + y[j * ni + im] + y[(j + 1) * ni + i] + y[(j - 1) * ni + i]);

                x[idx] += omega * (avg_x - x[idx]);
                y[idx] += omega * (avg_y - y[idx]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_radial_distribution() {
        let dist = generate_radial_distribution(64, 0.001, 20.0);
        assert_eq!(dist.len(), 64);
        assert!((dist[0]).abs() < 1e-10);
        assert!((dist[1] - 0.001).abs() < 1e-6);
        assert!((dist[63] - 20.0).abs() / 20.0 < 0.01);
    }

    #[test]
    fn test_resample_identity() {
        let x = vec![0.0, 0.5, 1.0];
        let y = vec![0.0, 0.1, 0.0];
        let (rx, ry) = resample_airfoil(&x, &y, 3);
        for i in 0..3 {
            assert!((rx[i] - x[i]).abs() < 1e-10);
            assert!((ry[i] - y[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_generate_o_mesh() {
        // Simple circle as "airfoil" for testing
        let n = 64;
        let mut ax = vec![0.0; n];
        let mut ay = vec![0.0; n];
        for i in 0..n {
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            ax[i] = 0.5 + 0.5 * theta.cos();
            ay[i] = 0.5 * theta.sin();
        }

        let mesh = generate_o_mesh(&ax, &ay, 64, 32, 10.0, 0.01);
        assert_eq!(mesh.x.len(), 64 * 32);
        assert_eq!(mesh.y.len(), 64 * 32);
    }
}
