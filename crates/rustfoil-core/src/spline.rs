//! Cubic spline interpolation for airfoil geometry.
//!
//! Provides spline interpolation for smoothing and repaneling airfoil coordinates.
//!
//! Uses natural cubic splines with arc-length parameterization:
//! - Parameter `s` = cumulative arc length along the airfoil
//! - `x(s)` and `y(s)` are each represented by cubic splines
//!
//! Includes XFOIL-exact repaneling via `resample_xfoil()` (PANGEN algorithm).

use crate::error::GeometryError;
use crate::point::Point;

/// A parametric cubic spline representing an airfoil contour.
///
/// The spline is parameterized by arc length `s ∈ [0, s_max]`, where
/// `s_max` is the total arc length of the original point set.
///
/// # Construction
/// Use `CubicSpline::from_points()` to create a spline from discrete points.
/// The spline will pass exactly through all input points.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    /// Arc-length parameter values at each input point
    s_values: Vec<f64>,
    /// Original input points (stored for XFOIL-exact paneling)
    points: Vec<Point>,
    /// Spline coefficients for x(s)
    x_coeffs: SplineCoeffs,
    /// Spline coefficients for y(s)
    y_coeffs: SplineCoeffs,
}

/// Coefficients for a 1D cubic spline.
///
/// For segment `i`, the spline is:
/// ```text
/// f(t) = a[i] + b[i]*t + c[i]*t² + d[i]*t³
/// ```
/// where `t = s - s[i]` is the local parameter within the segment.
#[derive(Debug, Clone)]
struct SplineCoeffs {
    a: Vec<f64>, // Function values at knots
    b: Vec<f64>, // First derivative coefficients
    c: Vec<f64>, // Second derivative coefficients  
    d: Vec<f64>, // Third derivative coefficients
}

impl CubicSpline {
    /// Construct a cubic spline from a sequence of points.
    ///
    /// # Arguments
    /// * `points` - Ordered points along the airfoil contour
    ///
    /// # Errors
    /// - `InsufficientPoints` if fewer than 2 points provided
    /// - `SplineInterpolationFailed` if points are degenerate
    ///
    /// # Boundary Conditions
    /// Uses "natural" boundary conditions (zero second derivative at endpoints).
    /// This is appropriate for closed airfoil contours where the actual
    /// boundary condition is periodic.
    pub fn from_points(points: &[Point]) -> Result<Self, GeometryError> {
        if points.len() < 2 {
            return Err(GeometryError::InsufficientPoints {
                required: 2,
                provided: points.len(),
            });
        }

        // Compute arc-length parameterization
        let s_values = compute_arc_lengths(points);

        // Extract x and y coordinates
        let x_vals: Vec<f64> = points.iter().map(|p| p.x).collect();
        let y_vals: Vec<f64> = points.iter().map(|p| p.y).collect();

        // Build spline coefficients for x(s) and y(s)
        let x_coeffs = build_natural_spline(&s_values, &x_vals)?;
        let y_coeffs = build_natural_spline(&s_values, &y_vals)?;

        Ok(Self {
            s_values,
            points: points.to_vec(),
            x_coeffs,
            y_coeffs,
        })
    }

    /// Total arc length of the spline.
    pub fn total_arc_length(&self) -> f64 {
        *self.s_values.last().unwrap_or(&0.0)
    }

    /// Evaluate the spline at arc-length parameter `s`.
    ///
    /// # Arguments
    /// * `s` - Arc-length parameter in `[0, total_arc_length()]`
    ///
    /// # Returns
    /// The interpolated point `(x(s), y(s))`.
    ///
    /// # Panics
    /// May panic if `s` is significantly outside the valid range.
    pub fn evaluate(&self, s: f64) -> Point {
        let x = evaluate_spline(&self.s_values, &self.x_coeffs, s);
        let y = evaluate_spline(&self.s_values, &self.y_coeffs, s);
        Point::new(x, y)
    }

    /// Evaluate the first derivative (tangent) at arc-length `s`.
    ///
    /// Returns `(dx/ds, dy/ds)`, which is the unit tangent vector
    /// (since `s` is arc-length parameterized).
    pub fn derivative(&self, s: f64) -> (f64, f64) {
        let dx = evaluate_spline_derivative(&self.s_values, &self.x_coeffs, s);
        let dy = evaluate_spline_derivative(&self.s_values, &self.y_coeffs, s);
        (dx, dy)
    }

    /// Evaluate the second derivative (curvature-related) at arc-length `s`.
    pub fn second_derivative(&self, s: f64) -> (f64, f64) {
        let d2x = evaluate_spline_second_derivative(&self.s_values, &self.x_coeffs, s);
        let d2y = evaluate_spline_second_derivative(&self.s_values, &self.y_coeffs, s);
        (d2x, d2y)
    }

    /// Compute curvature at arc-length `s`.
    ///
    /// # Mathematical Background
    /// For a parametric curve `(x(s), y(s))`:
    /// ```text
    /// κ = (x'*y'' - y'*x'') / (x'² + y'²)^(3/2)
    /// ```
    /// Since `s` is arc-length, `x'² + y'² = 1`, so:
    /// ```text
    /// κ = x'*y'' - y'*x''
    /// ```
    pub fn curvature(&self, s: f64) -> f64 {
        let (dx, dy) = self.derivative(s);
        let (d2x, d2y) = self.second_derivative(s);

        // For arc-length parameterization, denominator ≈ 1
        let denom = (dx * dx + dy * dy).powf(1.5);
        if denom < 1e-12 {
            return 0.0;
        }

        (dx * d2y - dy * d2x) / denom
    }

    /// Resample the spline at `n` uniformly-spaced arc-length positions.
    ///
    /// This is useful for creating evenly-spaced panels (though cosine
    /// spacing is usually preferred for airfoils).
    pub fn resample_uniform(&self, n: usize) -> Vec<Point> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.evaluate(0.0)];
        }

        let s_max = self.total_arc_length();
        let ds = s_max / (n - 1) as f64;

        (0..n).map(|i| self.evaluate(i as f64 * ds)).collect()
    }

    /// Resample using cosine spacing (clusters points at LE and TE).
    ///
    /// # Cosine Spacing
    /// Instead of uniform `s` values, we use:
    /// ```text
    /// s_i = s_max * (1 - cos(π * i / (n-1))) / 2
    /// ```
    /// This produces dense spacing near `s = 0` and `s = s_max` (the trailing
    /// edge) and at the leading edge (if properly set up).
    ///
    /// # Arguments
    /// * `n` - Number of output points
    pub fn resample_cosine(&self, n: usize) -> Vec<Point> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.evaluate(0.0)];
        }

        let s_max = self.total_arc_length();
        let pi = std::f64::consts::PI;

        (0..n)
            .map(|i| {
                let theta = pi * i as f64 / (n - 1) as f64;
                let s = s_max * (1.0 - theta.cos()) / 2.0;
                self.evaluate(s)
            })
            .collect()
    }

    /// Resample using XFOIL-style curvature-based paneling.
    ///
    /// This implements Mark Drela's PANGEN algorithm from XFOIL, which
    /// distributes panel nodes based on local curvature. Panels are bunched
    /// in regions of high curvature (leading edge) and at the trailing edge.
    ///
    /// # Algorithm
    /// 1. Compute curvature at dense set of points
    /// 2. Smooth the curvature distribution
    /// 3. Add artificial curvature at TE for bunching
    /// 4. Newton iteration to find node positions where `(1 + C*κ)*ds` is equal
    ///
    /// # Arguments
    /// * `n` - Number of output points (XFOIL default: 160)
    /// * `params` - Paneling parameters (use `PanelingParams::default()` for XFOIL defaults)
    ///
    /// # Returns
    /// Resampled points with curvature-based distribution
    pub fn resample_xfoil(&self, n: usize, params: &PanelingParams) -> Vec<Point> {
        if n < 3 {
            return self.resample_uniform(n.max(2));
        }

        let s_max = self.total_arc_length();
        let sbref = s_max / 2.0;  // Normalizing length (~ chord)

        // Step 1: Compute curvature at EXACT buffer points (like XFOIL)
        // XFOIL computes CURV at each of the NB buffer points
        // This is critical for exact matching
        use crate::xfoil_spline::XfoilSpline;
        let xfoil_spline = XfoilSpline::from_points(&self.points)
            .expect("Failed to create XFOIL spline");
        
        let n_buffer = self.points.len();
        let s_buffer = self.s_values.clone();  // Use exact arc-length values
        let mut curv_buffer: Vec<f64> = s_buffer.iter()
            .map(|&s| xfoil_spline.curvature(s).abs() * sbref)
            .collect();

        // Force symmetry in curvature buffer (average symmetric pairs)
        // This ensures symmetric airfoils produce symmetric paneling
        for i in 0..n_buffer / 2 {
            let j = n_buffer - 1 - i;
            let avg = (curv_buffer[i] + curv_buffer[j]) / 2.0;
            curv_buffer[i] = avg;
            curv_buffer[j] = avg;
        }

        // Find LE using XFOIL's LEFIND algorithm (Newton iteration)
        // This finds the exact arc-length where tangent is normal to chord
        let s_le = xfoil_spline.lefind();
        
        // Find the buffer index closest to LE
        let le_idx = s_buffer.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - s_le).abs()).partial_cmp(&((**b - s_le).abs())).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(n_buffer / 2);

        // Compute LE curvature at exact LE position (like XFOIL)
        let cv_le = xfoil_spline.curvature(s_le).abs() * sbref;

        // Average curvature near LE (XFOIL lines 1697-1705)
        // XFOIL samples 2*NK+1 = 7 points around the LE using the spline
        let nk = 3;
        let mut cv_sum = 0.0;
        for k in -nk..=nk {
            let frac = k as f64 / nk as f64;
            let sbk = s_le + frac * sbref / cv_le.max(20.0);
            let cvk = xfoil_spline.curvature(sbk).abs() * sbref;
            cv_sum += cvk;
        }
        let cv_avg = (cv_sum / (2 * nk + 1) as f64).max(10.0);

        // Step 2: Smooth curvature and set TE bunching
        let cv_te = cv_avg * params.te_le_ratio;
        curv_buffer[0] = cv_te;
        curv_buffer[n_buffer - 1] = cv_te;

        // Smooth the curvature array
        // Pass exact LE position for XFOIL-style smoothing
        let smool = (1.0 / cv_avg.max(20.0)).max(0.25 / (n as f64 / 2.0));
        let mut smoothed = smooth_curvature_with_le_fix(&s_buffer, &curv_buffer, smool * sbref, s_le, cv_le);

        // Force symmetry again after smoothing
        for i in 0..n_buffer / 2 {
            let j = n_buffer - 1 - i;
            let avg = (smoothed[i] + smoothed[j]) / 2.0;
            smoothed[i] = avg;
            smoothed[j] = avg;
        }

        // Normalize curvature
        let cv_max = smoothed.iter().fold(0.0_f64, |acc, &x| acc.max(x.abs()));
        let normalized: Vec<f64> = smoothed.iter().map(|&x| x / cv_max.max(1e-10)).collect();

        // Create a spline of curvature values (XFOIL's W5/W6 arrays)
        // This is the key to exact XFOIL matching - XFOIL uses SEGSPL on curvature
        // then SEVAL/DEVAL to interpolate during Newton iteration
        use crate::xfoil_spline::Spline1D;
        let curv_spline = Spline1D::new(&s_buffer, &normalized)
            .expect("Failed to create curvature spline");

        // Step 3: Newton iteration for panel node positions
        // XFOIL uses a denser temporary grid (IPFAC times more nodes) for iteration
        let ipfac = 5;
        let nn = ipfac * (n - 1) + 1;  // Number of temporary nodes
        
        // cc controls how much curvature affects panel density (XFOIL: CC = 6.0 * CVPAR)
        let cc = 6.0 * params.curv_param;

        // XFOIL: RTF = (RDSTE-1.0)*2.0 + 1.0
        // With RDSTE=0.667: RTF = (0.667-1.0)*2.0 + 1.0 = 0.334
        let rtf = (params.te_spacing_ratio - 1.0) * 2.0 + 1.0;
        
        // XFOIL initial guess: uniform spacing with TE ratio adjustment
        // DSAVG = (SB(NB)-SB(1)) / (FLOAT(NN-3) + 2.0*RTF)
        let ds_avg = s_max / ((nn - 3) as f64 + 2.0 * rtf);
        
        let mut s_new = Vec::with_capacity(nn);
        s_new.push(0.0);  // SNEW(1) = SB(1) = 0
        for i in 1..nn - 1 {
            // XFOIL (1-indexed): SNEW(I) = SB(1) + DSAVG * (FLOAT(I-2) + RTF)
            s_new.push(ds_avg * ((i - 1) as f64 + rtf));
        }
        s_new.push(s_max);  // SNEW(NN) = SB(NB)

        // Newton iteration (20 iterations max, like XFOIL)
        for _iter in 0..20 {
            // Build tridiagonal system for nn temporary nodes
            let mut w1 = vec![0.0; nn];  // Lower diagonal
            let mut w2 = vec![0.0; nn];  // Main diagonal
            let mut w3 = vec![0.0; nn];  // Upper diagonal
            let mut w4 = vec![0.0; nn];  // RHS / solution

            // Interpolate curvature and its derivative using XFOIL's spline
            // XFOIL: CV = SEVAL(SNEW,W5,W6,SB,NB), CVS = DEVAL(SNEW,W5,W6,SB,NB)
            let get_curv_and_deriv = |s: f64| -> (f64, f64) {
                (curv_spline.eval(s), curv_spline.deriv(s))
            };

            let (mut cv1, mut cvs1) = get_curv_and_deriv(s_new[0]);
            let (mut cv2, mut cvs2) = get_curv_and_deriv(s_new[1]);
            
            // XFOIL: CAVM = SQRT(CV1**2 + CV2**2)
            let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
            let (mut cavm_s1, mut cavm_s2) = if cavm > 1e-12 {
                (cvs1 * cv1 / cavm, cvs2 * cv2 / cavm)
            } else {
                (0.0, 0.0)
            };

            for i in 1..nn - 1 {
                let (cv3, cvs3) = get_curv_and_deriv(s_new[(i + 1).min(nn - 1)]);

                // Match XFOIL's sign convention:
                // DSM = S(I) - S(I-1) > 0  (panel to the left)
                // DSP = S(I) - S(I+1) < 0  (panel to the right, negative!)
                let dsm = s_new[i] - s_new[i - 1];  // positive
                let dsp = s_new[i] - s_new[i + 1];  // negative (XFOIL convention)

                // XFOIL: CAVP = SQRT(CV3**2 + CV2**2)
                let cavp = (cv3 * cv3 + cv2 * cv2).sqrt();
                let (cavp_s2, cavp_s3) = if cavp > 1e-12 {
                    (cvs2 * cv2 / cavp, cvs3 * cv3 / cavp)
                } else {
                    (0.0, 0.0)
                };

                let fm = cc * cavm + 1.0;
                let fp = cc * cavp + 1.0;

                // XFOIL residual: REZ = DSP*FP + DSM*FM
                let rez = dsp * fp + dsm * fm;

                // XFOIL Jacobian WITH curvature derivatives (lines 1912-1914)
                // W1(I) = -FM + CC*DSM*CAVM_S1
                // W2(I) = FP + FM + CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
                // W3(I) = -FP + CC*DSP*CAVP_S3
                w1[i] = -fm + cc * dsm * cavm_s1;
                w2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
                w3[i] = -fp + cc * dsp * cavp_s3;
                w4[i] = -rez;

                // Shift variables for next iteration (XFOIL lines 1920-1926)
                cv1 = cv2;
                cv2 = cv3;
                cvs1 = cvs2;
                cvs2 = cvs3;
                cavm = cavp;
                cavm_s1 = cavp_s2;
                cavm_s2 = cavp_s3;
            }

            // Fix endpoints
            w2[0] = 1.0;
            w3[0] = 0.0;
            w4[0] = 0.0;
            w1[nn - 1] = 0.0;
            w2[nn - 1] = 1.0;
            w4[nn - 1] = 0.0;

            // Apply TE spacing ratio constraint
            if (rtf - 1.0).abs() > 1e-10 {
                // Near upper TE (i=1 in XFOIL, which is index 1 here)
                w4[1] = -((s_new[1] - s_new[0]) + rtf * (s_new[1] - s_new[2]));
                w1[1] = -1.0;
                w2[1] = 1.0 + rtf;
                w3[1] = -rtf;

                // Near lower TE (i=NN-1 in XFOIL, which is nn-2 here)
                let i = nn - 2;
                w4[i] = -((s_new[i] - s_new[i + 1]) + rtf * (s_new[i] - s_new[i - 1]));
                w3[i] = -1.0;
                w2[i] = 1.0 + rtf;
                w1[i] = -rtf;
            }

            // Solve tridiagonal system (inline to avoid slice issues)
            // Forward sweep
            let mut c_prime = vec![0.0; nn];
            let mut d_prime = vec![0.0; nn];

            c_prime[0] = w3[0] / w2[0];
            d_prime[0] = w4[0] / w2[0];

            for ii in 1..nn {
                let denom = w2[ii] - w1[ii] * c_prime[ii - 1];
                if denom.abs() < 1e-14 {
                    break;
                }
                c_prime[ii] = w3[ii] / denom;
                d_prime[ii] = (w4[ii] - w1[ii] * d_prime[ii - 1]) / denom;
            }

            // Back substitution
            let mut delta = vec![0.0; nn];
            delta[nn - 1] = d_prime[nn - 1];

            for ii in (0..nn - 1).rev() {
                delta[ii] = d_prime[ii] - c_prime[ii] * delta[ii + 1];
            }

            // Under-relaxation to prevent node crossing
            let mut rlx = 1.0;
            let mut dmax = 0.0_f64;
            for i in 0..nn - 1 {
                let ds = s_new[i + 1] - s_new[i];
                let dds = delta[i + 1] - delta[i];
                if ds.abs() > 1e-12 {
                    let dsrat = 1.0 + rlx * dds / ds;
                    if dsrat > 4.0 {
                        rlx = (4.0 - 1.0) * ds / dds;
                    }
                    if dsrat < 0.2 {
                        rlx = (0.2 - 1.0) * ds / dds;
                    }
                }
                dmax = dmax.max(delta[i].abs());
            }

            // Update node positions
            for i in 1..nn - 1 {
                s_new[i] += rlx * delta[i];
            }

            // Check convergence
            if dmax < 1e-3 {
                break;
            }
        }

        // Subsample: take every IPFAC-th point to get final n points
        // XFOIL: IND = IPFAC*(I-1) + 1, S(I) = SNEW(IND)
        let mut s_final = Vec::with_capacity(n);
        for i in 0..n {
            let ind = ipfac * i;
            s_final.push(s_new[ind]);
        }

        // NOTE: XFOIL does NOT enforce symmetry in the paneling - it relies purely on 
        // the curvature distribution. We previously had symmetry enforcement here which
        // caused issues for cambered airfoils. Now removed to match XFOIL behavior.

        // Evaluate spline at final node positions
        s_final.iter().map(|&s| self.evaluate(s)).collect()
    }

    /// Get the arc-length parameter values at original knots.
    pub fn arc_length_params(&self) -> &[f64] {
        &self.s_values
    }
}

/// Parameters for XFOIL-style panel distribution.
///
/// These control how panel nodes are distributed based on local curvature.
#[derive(Debug, Clone, Copy)]
pub struct PanelingParams {
    /// Curvature attraction parameter (XFOIL's CVPAR).
    /// - 0 = uniform panel spacing
    /// - 1 = nodes bunched in high-curvature regions (default)
    pub curv_param: f64,

    /// TE/LE panel density ratio (XFOIL's CTERAT).
    /// Ratio of panel density at TE to panel density at LE.
    /// Default: 0.15 (TE panels are ~15% as dense as LE panels)
    pub te_le_ratio: f64,

    /// TE panel spacing ratio (XFOIL's RDSTE).
    /// Ratio of TE panel length to adjacent panel length.
    /// Default: 0.667
    pub te_spacing_ratio: f64,
}

impl Default for PanelingParams {
    fn default() -> Self {
        Self {
            curv_param: 1.0,       // XFOIL default (CVPAR)
            te_le_ratio: 0.15,     // XFOIL default (CTERAT)
            te_spacing_ratio: 0.667, // XFOIL default (RDSTE)
        }
    }
}

impl PanelingParams {
    /// Create paneling parameters with custom curvature bunching.
    pub fn with_curvature(curv_param: f64) -> Self {
        Self {
            curv_param,
            ..Default::default()
        }
    }

    /// Create uniform paneling (no curvature bunching).
    pub fn uniform() -> Self {
        Self {
            curv_param: 0.0,
            te_le_ratio: 1.0,
            te_spacing_ratio: 1.0,
        }
    }
}

/// Smooth a curvature array using tridiagonal diffusion with LE curvature preservation.
/// This matches XFOIL's PANGEN smoothing exactly (xfoil.f lines 1730-1784).
///
/// # Arguments
/// * `s` - Arc-length values at buffer points
/// * `curv` - Curvature values at buffer points  
/// * `smoo_len` - Smoothing length scale
/// * `s_le` - Exact LE arc-length position (from LEFIND)
/// * `cv_le` - Curvature value at LE
fn smooth_curvature_with_le_fix(s: &[f64], curv: &[f64], smoo_len: f64, s_le: f64, cv_le: f64) -> Vec<f64> {
    let n = curv.len();
    if n < 3 {
        return curv.to_vec();
    }

    let smoosq = smoo_len * smoo_len;

    // Build tridiagonal system (XFOIL: W1=lower, W2=diag, W3=upper, W5=rhs)
    let mut lower = vec![0.0; n];
    let mut diag = vec![0.0; n];
    let mut upper = vec![0.0; n];
    let mut rhs = curv.to_vec();

    // First point: W2(1) = 1.0, W3(1) = 0.0 (identity)
    diag[0] = 1.0;
    upper[0] = 0.0;

    // Interior points: diffusion smoothing (XFOIL lines 1733-1747)
    for i in 1..n - 1 {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];
        let dso = 0.5 * (s[i + 1] - s[i - 1]);

        if dsm.abs() < 1e-12 || dsp.abs() < 1e-12 {
            // Corner point - leave unchanged
            lower[i] = 0.0;
            diag[i] = 1.0;
            upper[i] = 0.0;
        } else {
            lower[i] = smoosq * (-1.0 / dsm) / dso;
            diag[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            upper[i] = smoosq * (-1.0 / dsp) / dso;
        }
    }

    // Last point: W1(NB) = 0.0, W2(NB) = 1.0 (identity)
    lower[n - 1] = 0.0;
    diag[n - 1] = 1.0;

    // XFOIL lines 1754-1783: Fix curvature at/near LE point
    // Handle three cases:
    // 1. Node falls exactly on LE (s[i] == s_le)
    // 2. LE is between nodes s[i-1] and s[i]
    // 3. Sharp LE (doubled point) - handled by corner point logic above
    
    for i in 1..n - 1 {
        if (s[i] - s_le).abs() < 1e-10 {
            // Case 1: Node falls right on LE - fix curvature there (XFOIL lines 1756-1760)
            lower[i] = 0.0;
            diag[i] = 1.0;
            upper[i] = 0.0;
            rhs[i] = cv_le;
        } else if s[i - 1] < s_le && s[i] > s_le {
            // Case 2: LE is between s[i-1] and s[i]
            // Modify equation at node just BEFORE LE (i-1) - XFOIL lines 1763-1770
            if i >= 2 {
                let dsm = s[i - 1] - s[i - 2];
                let dsp = s_le - s[i - 1];  // Distance to exact LE
                let dso = 0.5 * (s_le - s[i - 2]);

                if dsm.abs() > 1e-12 && dsp.abs() > 1e-12 {
                    lower[i - 1] = smoosq * (-1.0 / dsm) / dso;
                    diag[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                    upper[i - 1] = 0.0;  // No coupling past LE
                    rhs[i - 1] += smoosq * cv_le / (dsp * dso);
                }
            }

            // Modify equation at node just AFTER LE (i) - XFOIL lines 1772-1779
            if i + 1 < n {
                let dsm = s[i] - s_le;  // Distance from exact LE
                let dsp = s[i + 1] - s[i];
                let dso = 0.5 * (s[i + 1] - s_le);

                if dsm.abs() > 1e-12 && dsp.abs() > 1e-12 {
                    lower[i] = 0.0;  // No coupling past LE
                    diag[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                    upper[i] = smoosq * (-1.0 / dsp) / dso;
                    rhs[i] += smoosq * cv_le / (dsm * dso);
                }
            }
            break;  // Only process LE once
        }
    }

    // Solve tridiagonal system using Thomas algorithm (XFOIL: TRISOL)
    let mut c_prime = vec![0.0; n];
    let mut d_prime = vec![0.0; n];

    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];

    for i in 1..n {
        let denom = diag[i] - lower[i] * c_prime[i - 1];
        if denom.abs() < 1e-14 {
            return curv.to_vec();
        }
        c_prime[i] = upper[i] / denom;
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    let mut result = vec![0.0; n];
    result[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        result[i] = d_prime[i] - c_prime[i] * result[i + 1];
    }

    result
}

/// Compute cumulative arc lengths for a sequence of points.
fn compute_arc_lengths(points: &[Point]) -> Vec<f64> {
    let mut s = Vec::with_capacity(points.len());
    s.push(0.0);

    for i in 1..points.len() {
        let ds = (points[i] - points[i - 1]).norm();
        s.push(s[i - 1] + ds);
    }

    s
}

/// Build natural cubic spline coefficients using the Thomas algorithm.
///
/// Natural spline: second derivative is zero at both endpoints.
fn build_natural_spline(s: &[f64], f: &[f64]) -> Result<SplineCoeffs, GeometryError> {
    let n = s.len();
    if n < 2 {
        return Err(GeometryError::SplineInterpolationFailed {
            reason: "need at least 2 points",
        });
    }

    if n == 2 {
        // Linear interpolation
        let h = s[1] - s[0];
        if h.abs() < 1e-14 {
            return Err(GeometryError::SplineInterpolationFailed {
                reason: "coincident points",
            });
        }
        return Ok(SplineCoeffs {
            a: vec![f[0]],
            b: vec![(f[1] - f[0]) / h],
            c: vec![0.0],
            d: vec![0.0],
        });
    }

    // Number of segments
    let n_seg = n - 1;

    // Interval sizes
    let h: Vec<f64> = (0..n_seg).map(|i| s[i + 1] - s[i]).collect();

    // Check for zero intervals
    for &hi in h.iter() {
        if hi.abs() < 1e-14 {
            return Err(GeometryError::SplineInterpolationFailed {
                reason: "coincident knots in spline",
            });
        }
    }

    // Build tridiagonal system for second derivatives (c values)
    // Natural BCs: c[0] = 0, c[n-1] = 0
    let mut c = vec![0.0; n];

    if n > 2 {
        // Interior equations
        let n_interior = n - 2;
        let mut lower = vec![0.0; n_interior];
        let mut diag = vec![0.0; n_interior];
        let mut upper = vec![0.0; n_interior];
        let mut rhs = vec![0.0; n_interior];

        for i in 0..n_interior {
            let j = i + 1; // Index in original arrays
            lower[i] = h[j - 1];
            diag[i] = 2.0 * (h[j - 1] + h[j]);
            upper[i] = h[j];
            rhs[i] = 3.0 * ((f[j + 1] - f[j]) / h[j] - (f[j] - f[j - 1]) / h[j - 1]);
        }

        // Solve tridiagonal system (Thomas algorithm)
        let c_interior = solve_tridiagonal(&lower, &diag, &upper, &rhs);

        for (i, &ci) in c_interior.iter().enumerate() {
            c[i + 1] = ci;
        }
    }

    // Compute remaining coefficients
    let a: Vec<f64> = f[..n_seg].to_vec();

    let b: Vec<f64> = (0..n_seg)
        .map(|i| (f[i + 1] - f[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0)
        .collect();

    let d: Vec<f64> = (0..n_seg)
        .map(|i| (c[i + 1] - c[i]) / (3.0 * h[i]))
        .collect();

    let c_seg: Vec<f64> = c[..n_seg].to_vec();

    Ok(SplineCoeffs {
        a,
        b,
        c: c_seg,
        d,
    })
}

/// Solve a tridiagonal system using the Thomas algorithm.
fn solve_tridiagonal(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    if n == 0 {
        return vec![];
    }

    // Forward elimination
    let mut c_prime = vec![0.0; n];
    let mut d_prime = vec![0.0; n];

    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];

    for i in 1..n {
        let denom = diag[i] - lower[i] * c_prime[i - 1];
        c_prime[i] = if i < n - 1 {
            upper[i] / denom
        } else {
            0.0
        };
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    let mut x = vec![0.0; n];
    x[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    x
}

/// Find the segment index for parameter value `s`.
fn find_segment(s_values: &[f64], s: f64) -> usize {
    // Binary search for the segment containing s
    match s_values.binary_search_by(|&si| si.partial_cmp(&s).unwrap()) {
        Ok(i) => i.min(s_values.len() - 2),
        Err(i) => (i.saturating_sub(1)).min(s_values.len() - 2),
    }
}

/// Evaluate spline at parameter `s`.
fn evaluate_spline(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    coeffs.a[i] + coeffs.b[i] * t + coeffs.c[i] * t * t + coeffs.d[i] * t * t * t
}

/// Evaluate first derivative at parameter `s`.
fn evaluate_spline_derivative(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    coeffs.b[i] + 2.0 * coeffs.c[i] * t + 3.0 * coeffs.d[i] * t * t
}

/// Evaluate second derivative at parameter `s`.
fn evaluate_spline_second_derivative(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    2.0 * coeffs.c[i] + 6.0 * coeffs.d[i] * t
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;
    use approx::assert_relative_eq;

    #[test]
    fn test_linear_interpolation() {
        let points = vec![point(0.0, 0.0), point(1.0, 1.0)];
        let spline = CubicSpline::from_points(&points).unwrap();

        let mid = spline.evaluate(spline.total_arc_length() / 2.0);
        assert_relative_eq!(mid.x, 0.5, epsilon = 0.01);
        assert_relative_eq!(mid.y, 0.5, epsilon = 0.01);
    }

    #[test]
    fn test_spline_passes_through_points() {
        let points = vec![
            point(0.0, 0.0),
            point(1.0, 1.0),
            point(2.0, 0.0),
            point(3.0, 1.0),
        ];
        let spline = CubicSpline::from_points(&points).unwrap();

        // Should pass through first and last points
        let p0 = spline.evaluate(0.0);
        assert_relative_eq!(p0.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(p0.y, 0.0, epsilon = 1e-10);

        let p_end = spline.evaluate(spline.total_arc_length());
        assert_relative_eq!(p_end.x, 3.0, epsilon = 1e-10);
        assert_relative_eq!(p_end.y, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cosine_spacing() {
        let points = vec![point(0.0, 0.0), point(1.0, 0.0), point(2.0, 0.0)];
        let spline = CubicSpline::from_points(&points).unwrap();

        let resampled = spline.resample_cosine(5);
        assert_eq!(resampled.len(), 5);

        // First and last should match original endpoints
        assert_relative_eq!(resampled[0].x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(resampled[4].x, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_insufficient_points() {
        let points = vec![point(0.0, 0.0)];
        let result = CubicSpline::from_points(&points);
        assert!(matches!(
            result,
            Err(GeometryError::InsufficientPoints { .. })
        ));
    }

    #[test]
    fn test_xfoil_paneling_preserves_endpoints() {
        // Create a simple airfoil-like shape
        let n = 50;
        let pi = std::f64::consts::PI;
        let points: Vec<_> = (0..n)
            .map(|i| {
                let theta = 2.0 * pi * i as f64 / (n - 1) as f64;
                // Ellipse with 4:1 aspect ratio
                point(0.5 * (1.0 + theta.cos()), 0.125 * theta.sin())
            })
            .collect();

        let spline = CubicSpline::from_points(&points).unwrap();
        let params = PanelingParams::default();
        let resampled = spline.resample_xfoil(100, &params);

        // Should have 100 points
        assert_eq!(resampled.len(), 100);

        // First and last should be near original endpoints
        assert_relative_eq!(resampled[0].x, points[0].x, epsilon = 0.01);
        assert_relative_eq!(resampled[0].y, points[0].y, epsilon = 0.01);
        assert_relative_eq!(resampled[99].x, points[n - 1].x, epsilon = 0.01);
        assert_relative_eq!(resampled[99].y, points[n - 1].y, epsilon = 0.01);
    }

    #[test]
    fn test_xfoil_paneling_symmetric() {
        // Create a symmetric airfoil (NACA 0012-like ellipse)
        let n = 100;
        let pi = std::f64::consts::PI;
        let points: Vec<_> = (0..n)
            .map(|i| {
                let theta = 2.0 * pi * i as f64 / (n - 1) as f64;
                // Symmetric ellipse - starts and ends at TE (x=1)
                point(0.5 * (1.0 + theta.cos()), 0.06 * theta.sin())
            })
            .collect();

        let spline = CubicSpline::from_points(&points).unwrap();
        let params = PanelingParams::default();
        let resampled = spline.resample_xfoil(80, &params);

        // Check symmetry: point i should mirror point (n-1-i) about y=0
        let n_out = resampled.len();
        for i in 0..n_out / 2 {
            let j = n_out - 1 - i;
            let p_i = &resampled[i];
            let p_j = &resampled[j];
            
            // x coordinates should be equal
            let x_diff = (p_i.x - p_j.x).abs();
            assert!(x_diff < 1e-6, 
                "Point {} x={:.6} should equal point {} x={:.6}, diff={:.2e}", 
                i, p_i.x, j, p_j.x, x_diff);
            
            // y coordinates should be opposite
            let y_diff = (p_i.y + p_j.y).abs();
            assert!(y_diff < 1e-6,
                "Point {} y={:.6} should equal -point {} y={:.6}, diff={:.2e}", 
                i, p_i.y, j, p_j.y, y_diff);
        }
    }

    #[test]
    fn test_xfoil_paneling_non_uniform() {
        // Create ellipse shape (airfoil-like)
        let n = 100;
        let pi = std::f64::consts::PI;
        let points: Vec<_> = (0..n)
            .map(|i| {
                let theta = 2.0 * pi * i as f64 / (n - 1) as f64;
                point(0.5 * (1.0 + theta.cos()), 0.1 * theta.sin())
            })
            .collect();

        let spline = CubicSpline::from_points(&points).unwrap();
        let params = PanelingParams::default();
        let resampled = spline.resample_xfoil(50, &params);

        // Compute all panel lengths
        let panel_lengths: Vec<f64> = (0..resampled.len() - 1)
            .map(|i| {
                let p1 = &resampled[i];
                let p2 = &resampled[i + 1];
                ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
            })
            .collect();

        // With curvature-based paneling, panels should NOT be uniform
        let min_len = panel_lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_len = panel_lengths.iter().fold(0.0_f64, |a, &b| a.max(b));

        // The ratio of max to min panel length should be > 2 for non-uniform paneling
        let ratio = max_len / min_len;
        assert!(
            ratio > 1.5,
            "Expected non-uniform paneling (ratio {}), got nearly uniform",
            ratio
        );

        // Compare with uniform paneling
        let uniform = spline.resample_uniform(50);
        let uniform_lengths: Vec<f64> = (0..uniform.len() - 1)
            .map(|i| {
                let p1 = &uniform[i];
                let p2 = &uniform[i + 1];
                ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
            })
            .collect();

        let uniform_min = uniform_lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let uniform_max = uniform_lengths.iter().fold(0.0_f64, |a, &b| a.max(b));
        let uniform_ratio = uniform_max / uniform_min;

        // XFOIL paneling should have more variation than uniform
        assert!(
            ratio > uniform_ratio,
            "XFOIL paneling ratio ({}) should exceed uniform ratio ({})",
            ratio,
            uniform_ratio
        );
    }
}

#[cfg(test)]
mod comparison_tests {
    use super::*;
    use crate::naca::naca4;
    use crate::point::point;
    use crate::xfoil_spline::XfoilSpline;
    
    /// Load coordinates from XFOIL's PSAV output file.
    fn load_xfoil_dat(content: &str) -> Vec<Point> {
        content
            .lines()
            .filter_map(|line| {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let x = parts[0].parse::<f64>().ok()?;
                    let y = parts[1].parse::<f64>().ok()?;
                    Some(point(x, y))
                } else {
                    None
                }
            })
            .collect()
    }
    
    
    #[test]
    fn test_naca4_xfoil_symmetry() {
        // Test that XFOIL's NACA generator produces symmetric coordinates
        let coords = naca4(12, Some(100));  // NACA 0012 with 100 points per side
        
        let n = coords.len();
        let nside = (n + 1) / 2;  // number of points per side
        
        println!("\nNACA 0012 from XFOIL generator (nside={}, total={}):", nside, n);
        println!("First 5 points (upper TE):");
        for i in 0..5 {
            println!("  [{:3}] x={:.8}, y={:+.8}", i, coords[i].x, coords[i].y);
        }
        println!("Around LE (nside-1={}):", nside - 1);
        for i in (nside - 3)..=(nside + 1).min(n - 1) {
            println!("  [{:3}] x={:.8}, y={:+.8}", i, coords[i].x, coords[i].y);
        }
        println!("Last 5 points (lower TE):");
        for i in (n - 5)..n {
            println!("  [{:3}] x={:.8}, y={:+.8}", i, coords[i].x, coords[i].y);
        }
        
        // Verify symmetry: point i on upper should match point n-1-i on lower
        // (with same x, opposite y)
        for i in 0..nside {
            let j = n - 1 - i;
            if i >= j { break; }
            
            let x_diff = (coords[i].x - coords[j].x).abs();
            let y_sum = (coords[i].y + coords[j].y).abs();
            
            assert!(x_diff < 1e-12, 
                "Point {} x={:.10} != point {} x={:.10}", i, coords[i].x, j, coords[j].x);
            assert!(y_sum < 1e-12,
                "Point {} y={:.10} != -point {} y={:.10}", i, coords[i].y, j, coords[j].y);
        }
        
        println!("\nSymmetry verified!");
    }
    
    #[test]
    fn compare_spline_curvatures() {
        // Test that XfoilSpline curvature matches our CubicSpline
        let coords = naca4(12, Some(100));  // Use XFOIL generator
        
        let cubic_spline = CubicSpline::from_points(&coords).unwrap();
        let xfoil_spline = XfoilSpline::from_points(&coords).unwrap();
        
        let s_max_cubic = cubic_spline.total_arc_length();
        let s_max_xfoil = xfoil_spline.total_arc_length();
        
        println!("\nComparing spline curvatures (CubicSpline vs XfoilSpline):");
        println!("Arc lengths: CubicSpline={:.6}, XfoilSpline={:.6}", s_max_cubic, s_max_xfoil);
        
        // Compare curvature at LE (s=0.5)
        let s_le_cubic = s_max_cubic * 0.5;
        let s_le_xfoil = s_max_xfoil * 0.5;
        
        let k_cubic = cubic_spline.curvature(s_le_cubic);
        let k_xfoil = xfoil_spline.curvature(s_le_xfoil);
        
        println!("LE curvature: CubicSpline={:.2}, XfoilSpline={:.2}", k_cubic, k_xfoil);
        
        // Compare at several points
        println!("\nCurvature comparison:");
        for i in 0..=10 {
            let s_norm = i as f64 / 10.0;
            let k_c = cubic_spline.curvature(s_norm * s_max_cubic);
            let k_x = xfoil_spline.curvature(s_norm * s_max_xfoil);
            let diff = if k_c.abs() > 1e-6 { ((k_x - k_c) / k_c * 100.0).abs() } else { 0.0 };
            println!("  s={:.1}: cubic={:8.2}, xfoil={:8.2}, diff={:5.1}%", 
                s_norm, k_c, k_x, diff);
        }
    }
    
    #[test]
    fn compare_with_xfoil_paneling() {
        // Use XFOIL's exact NACA generator with same NSIDE as XFOIL
        // XFOIL: NSIDE = IQX/3 = 370/3 = 123 points per side = 245 total
        let coords = naca4(12, Some(123));  // 245 total points, matches XFOIL exactly
        
        let spline = CubicSpline::from_points(&coords).unwrap();
        
        // Debug curvature
        let s_max = spline.total_arc_length();
        println!("\nCurvature along arc-length (using XFOIL NACA4):");
        for i in 0..=10 {
            let s = s_max * i as f64 / 10.0;
            let k = spline.curvature(s);
            let pt = spline.evaluate(s);
            println!("  s={:.3} (x={:.4}, y={:+.4}): κ={:.2}", 
                s/s_max, pt.x, pt.y, k);
        }
        
        // XFOIL-style paneling
        let params = PanelingParams::default();
        let resampled = spline.resample_xfoil(160, &params);
        
        // Verify output is symmetric
        let n = resampled.len();
        println!("\nVerifying output symmetry:");
        for i in 0..5 {
            let j = n - 1 - i;
            println!("  [{:3}] x={:.8}, y={:+.8}  <=>  [{:3}] x={:.8}, y={:+.8}",
                i, resampled[i].x, resampled[i].y,
                j, resampled[j].x, resampled[j].y);
        }
        
        let panel_lengths: Vec<f64> = (0..resampled.len() - 1)
            .map(|i| {
                let dx = resampled[i+1].x - resampled[i].x;
                let dy = resampled[i+1].y - resampled[i].y;
                (dx*dx + dy*dy).sqrt()
            })
            .collect();
        
        let min_len = panel_lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_len = panel_lengths.iter().fold(0.0_f64, |a, &b| a.max(b));
        
        let min_idx = panel_lengths.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        
        println!("\nXFOIL-style panel lengths:");
        println!("  Min: {:.6}, Max: {:.6}, Ratio: {:.2} (target: 9.28)", min_len, max_len, max_len/min_len);
        println!("  Smallest at idx {}: {:.6} (target: idx ~79, 0.001811)", min_idx, panel_lengths[min_idx]);
        
        println!("\nPanel lengths at key positions:");
        println!("  [0-1] TE upper: {:.6} (target: 0.0084)", panel_lengths[0]);
        println!("  [79-80] LE: {:.6} (target: 0.0018)", panel_lengths[79]);
        println!("  [158-159] TE lower: {:.6} (target: 0.0084)", panel_lengths[158]);
        
        // Verify we're within acceptable tolerance
        let le_diff = (panel_lengths[79] - 0.001811).abs() / 0.001811 * 100.0;
        let te_diff = (panel_lengths[0] - 0.0084).abs() / 0.0084 * 100.0;
        
        println!("\nDifference from XFOIL:");
        println!("  LE panel: {:.1}%", le_diff);
        println!("  TE panel: {:.1}%", te_diff);
        
        // Verify symmetry of output panels
        for i in 0..80 {
            let j = 158 - i;
            let len_diff = (panel_lengths[i] - panel_lengths[j]).abs();
            assert!(len_diff < 1e-10, 
                "Panel {} length {} != panel {} length {}", 
                i, panel_lengths[i], j, panel_lengths[j]);
        }
        println!("Output panel symmetry verified!");
    }
    
    #[test]
    fn test_exact_xfoil_match() {
        // Load actual XFOIL output and compare
        let xfoil_dat = include_str!("../../../naca0012_xfoil_paneled.dat");
        let xfoil_coords = load_xfoil_dat(xfoil_dat);
        
        println!("\nLoaded {} points from XFOIL .dat file", xfoil_coords.len());
        println!("First point: ({:.8}, {:.8})", xfoil_coords[0].x, xfoil_coords[0].y);
        println!("Last point: ({:.8}, {:.8})", 
            xfoil_coords.last().unwrap().x, xfoil_coords.last().unwrap().y);
        
        // Verify XFOIL's output is symmetric
        let n = xfoil_coords.len();
        println!("\nXFOIL output symmetry check:");
        for i in 0..5 {
            let j = n - 1 - i;
            let x_diff = (xfoil_coords[i].x - xfoil_coords[j].x).abs();
            let y_sum = (xfoil_coords[i].y + xfoil_coords[j].y).abs();
            println!("  [{:3}] ({:.6}, {:+.6}) <=> [{:3}] ({:.6}, {:+.6})  x_diff={:.2e}, y_sum={:.2e}",
                i, xfoil_coords[i].x, xfoil_coords[i].y,
                j, xfoil_coords[j].x, xfoil_coords[j].y,
                x_diff, y_sum);
        }
        
        // Calculate XFOIL's panel lengths
        let xfoil_lengths: Vec<f64> = (0..xfoil_coords.len() - 1)
            .map(|i| {
                let dx = xfoil_coords[i+1].x - xfoil_coords[i].x;
                let dy = xfoil_coords[i+1].y - xfoil_coords[i].y;
                (dx*dx + dy*dy).sqrt()
            })
            .collect();
        
        let xfoil_min = xfoil_lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let xfoil_max = xfoil_lengths.iter().fold(0.0_f64, |a, &b| a.max(b));
        
        println!("\nXFOIL panel lengths:");
        println!("  Min: {:.6}, Max: {:.6}, Ratio: {:.2}", xfoil_min, xfoil_max, xfoil_max/xfoil_min);
        println!("  TE: {:.6}, LE: {:.6}", xfoil_lengths[0], xfoil_lengths[79]);
    }
    
    #[test]
    fn test_point_by_point_comparison() {
        // Load actual XFOIL output
        let xfoil_dat = include_str!("../../../naca0012_xfoil_paneled.dat");
        let xfoil_coords = load_xfoil_dat(xfoil_dat);
        
        // IMPORTANT: XFOIL loads the NACA coordinates, then repanels them.
        // We need to load the SAME input that XFOIL used.
        // XFOIL's NACA command generates 245 buffer points (NSIDE=123), 
        // then PANGEN repanels to NPAN=160.
        
        // Try using XFOIL's output directly as input to see if spline is the issue
        // First, let's use our NACA generator
        let input_coords = naca4(12, Some(123));  // 245 buffer points
        
        // Create spline and repanel to 160 points (same as XFOIL output)
        let spline = CubicSpline::from_points(&input_coords).unwrap();
        let params = PanelingParams::default();
        let our_coords = spline.resample_xfoil(160, &params);
        
        assert_eq!(our_coords.len(), xfoil_coords.len(), 
            "Point count mismatch: ours={}, XFOIL={}", our_coords.len(), xfoil_coords.len());
        
        // Calculate point-by-point errors
        let mut x_errors: Vec<f64> = Vec::new();
        let mut y_errors: Vec<f64> = Vec::new();
        let mut total_errors: Vec<f64> = Vec::new();
        
        println!("\nPoint-by-point comparison (first 10 and last 10):");
        println!("{:>4} {:>12} {:>12} {:>12} {:>12} {:>10}", 
            "idx", "our_x", "xfoil_x", "our_y", "xfoil_y", "error");
        
        for i in 0..our_coords.len() {
            let dx = our_coords[i].x - xfoil_coords[i].x;
            let dy = our_coords[i].y - xfoil_coords[i].y;
            let err = (dx*dx + dy*dy).sqrt();
            
            x_errors.push(dx.abs());
            y_errors.push(dy.abs());
            total_errors.push(err);
            
            // Print first 10 and last 10
            if i < 10 || i >= our_coords.len() - 10 {
                println!("{:4} {:12.8} {:12.8} {:12.8} {:12.8} {:10.2e}",
                    i, our_coords[i].x, xfoil_coords[i].x, 
                    our_coords[i].y, xfoil_coords[i].y, err);
            } else if i == 10 {
                println!("  ... (140 more points) ...");
            }
        }
        
        // Calculate RMS errors
        let n = total_errors.len() as f64;
        let rms_x = (x_errors.iter().map(|e| e*e).sum::<f64>() / n).sqrt();
        let rms_y = (y_errors.iter().map(|e| e*e).sum::<f64>() / n).sqrt();
        let rms_total = (total_errors.iter().map(|e| e*e).sum::<f64>() / n).sqrt();
        
        let max_error = total_errors.iter().fold(0.0_f64, |a, &b| a.max(b));
        let max_idx = total_errors.iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        
        println!("\n=== ERROR SUMMARY ===");
        println!("RMS X error:     {:.6e}", rms_x);
        println!("RMS Y error:     {:.6e}", rms_y);
        println!("RMS total error: {:.6e}", rms_total);
        println!("Max error:       {:.6e} at point {}", max_error, max_idx);
        println!("Max error point: ours=({:.8}, {:.8}), XFOIL=({:.8}, {:.8})",
            our_coords[max_idx].x, our_coords[max_idx].y,
            xfoil_coords[max_idx].x, xfoil_coords[max_idx].y);
        
        // Also show LE comparison (point 79-80)
        println!("\nLE region comparison (points 78-81):");
        for i in 78..=81 {
            let dx = our_coords[i].x - xfoil_coords[i].x;
            let dy = our_coords[i].y - xfoil_coords[i].y;
            let err = (dx*dx + dy*dy).sqrt();
            println!("  [{:3}] ours=({:.8}, {:+.8}), XFOIL=({:.8}, {:+.8}), err={:.2e}",
                i, our_coords[i].x, our_coords[i].y,
                xfoil_coords[i].x, xfoil_coords[i].y, err);
        }
    }
    
    #[test]
    fn test_with_xfoil_buffer_input() {
        // Load XFOIL's buffer coordinates (the input to PANGEN)
        let xfoil_buffer = include_str!("../../../naca0012_buffer_real.dat");
        let buffer_coords = load_xfoil_dat(xfoil_buffer);
        
        // Load XFOIL's paneled output
        let xfoil_dat = include_str!("../../../naca0012_xfoil_paneled.dat");
        let xfoil_coords = load_xfoil_dat(xfoil_dat);
        
        println!("\nUsing XFOIL's actual buffer as input:");
        println!("Buffer points: {}", buffer_coords.len());
        println!("Target points: {}", xfoil_coords.len());
        
        // Compare our NACA generator with XFOIL's buffer
        let our_buffer = naca4(12, Some(123));
        println!("\nComparing buffer coordinates:");
        println!("  Our first 3: ({:.8}, {:.8}), ({:.8}, {:.8}), ({:.8}, {:.8})",
            our_buffer[0].x, our_buffer[0].y,
            our_buffer[1].x, our_buffer[1].y,
            our_buffer[2].x, our_buffer[2].y);
        println!("  XFOIL first 3: ({:.8}, {:.8}), ({:.8}, {:.8}), ({:.8}, {:.8})",
            buffer_coords[0].x, buffer_coords[0].y,
            buffer_coords[1].x, buffer_coords[1].y,
            buffer_coords[2].x, buffer_coords[2].y);
        
        // Use XFOIL's buffer to create spline and repanel
        let spline = CubicSpline::from_points(&buffer_coords).unwrap();
        let params = PanelingParams::default();
        let our_coords = spline.resample_xfoil(160, &params);
        
        // Calculate errors
        let mut total_errors: Vec<f64> = Vec::new();
        for i in 0..our_coords.len().min(xfoil_coords.len()) {
            let dx = our_coords[i].x - xfoil_coords[i].x;
            let dy = our_coords[i].y - xfoil_coords[i].y;
            total_errors.push((dx*dx + dy*dy).sqrt());
        }
        
        let n = total_errors.len() as f64;
        let rms_total = (total_errors.iter().map(|e| e*e).sum::<f64>() / n).sqrt();
        let max_error = total_errors.iter().fold(0.0_f64, |a, &b| a.max(b));
        
        println!("\n=== WITH XFOIL BUFFER INPUT (CubicSpline) ===");
        println!("RMS error: {:.6e}", rms_total);
        println!("Max error: {:.6e}", max_error);
        
        // Also try with XfoilSpline
        let xfoil_spline = XfoilSpline::from_points(&buffer_coords).unwrap();
        println!("\n=== XFOIL Spline comparison ===");
        println!("Arc length: CubicSpline={:.6}, XfoilSpline={:.6}", 
            spline.total_arc_length(), xfoil_spline.total_arc_length());
        
        // Check curvature at LE
        let s_le = spline.total_arc_length() / 2.0;
        let k_cubic = spline.curvature(s_le);
        let k_xfoil = xfoil_spline.curvature(s_le);
        println!("LE curvature: CubicSpline={:.2}, XfoilSpline={:.2}", k_cubic, k_xfoil);
    }
}
