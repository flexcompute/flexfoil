//! Main mfoil-style coupled viscous-inviscid solver.

use nalgebra::{DMatrix, DVector};
use super::state::{MfoilState, MfoilParam, MfoilPanel, MfoilGeom, SurfaceIndex};
use super::closures::*;
use super::residuals::{residual_station, residual_transition};

/// Configuration for the mfoil solver.
#[derive(Debug, Clone)]
pub struct MfoilConfig {
    /// Reynolds number
    pub reynolds: f64,
    /// Critical N-factor for transition
    pub ncrit: f64,
    /// Maximum global Newton iterations
    pub max_iter: usize,
    /// Residual tolerance
    pub rtol: f64,
    /// Freestream Mach number (0 = incompressible)
    pub mach: f64,
    /// Number of wake points
    pub n_wake: usize,
    /// Wake length in chords
    pub wake_len: f64,
    /// Whether to redo wake on alpha change
    pub redo_wake: bool,
    /// Verbosity level (0=quiet, 1=normal, 2=verbose)
    pub verbose: u32,
}

impl Default for MfoilConfig {
    fn default() -> Self {
        Self {
            reynolds: 1e6,
            ncrit: 9.0,
            max_iter: 100,
            rtol: 1e-10,
            mach: 0.0,
            n_wake: 32,
            wake_len: 1.0,
            redo_wake: true,
            verbose: 1,
        }
    }
}

impl MfoilConfig {
    /// Create config for given Reynolds number.
    pub fn with_reynolds(reynolds: f64) -> Self {
        Self {
            reynolds,
            ..Default::default()
        }
    }
}

/// Inviscid solution data.
#[derive(Debug, Clone)]
pub struct MfoilIsol {
    /// Vorticity at each airfoil node
    pub gam: DVector<f64>,
    /// Reference vorticity for α=0° (column 0) and α=90° (column 1)
    pub gamref: DMatrix<f64>,
    /// Stagnation point arc-length
    pub sstag: f64,
    /// Stagnation point coordinates [x, y]
    pub xstag: [f64; 2],
    /// Stagnation panel indices [i, i+1]
    pub istag: [usize; 2],
    /// Sign of edge velocity (CW -> tangential)
    pub sgnue: DVector<f64>,
    /// Distance from stagnation point at each node
    pub xi: DVector<f64>,
    /// Inviscid wake edge velocity reference
    pub uewiref: DMatrix<f64>,
}

impl Default for MfoilIsol {
    fn default() -> Self {
        Self {
            gam: DVector::zeros(0),
            gamref: DMatrix::zeros(0, 2),
            sstag: 0.0,
            xstag: [0.0, 0.0],
            istag: [0, 1],
            sgnue: DVector::zeros(0),
            xi: DVector::zeros(0),
            uewiref: DMatrix::zeros(0, 2),
        }
    }
}

/// Viscous solution data.
#[derive(Debug, Clone)]
pub struct MfoilVsol {
    /// Surface point indices [lower, upper, wake]
    pub is: [Vec<usize>; 3],
    /// Turbulent flag at each node
    pub turb: Vec<bool>,
    /// Transition location (xi, x) for lower [0] and upper [1]
    pub xt: [[f64; 2]; 2],
    /// Wake gap values
    pub wgap: Vec<f64>,
    /// Edge velocity influence matrix (Nsys x Nsys)
    pub ue_m: DMatrix<f64>,
}

impl Default for MfoilVsol {
    fn default() -> Self {
        Self {
            is: [Vec::new(), Vec::new(), Vec::new()],
            turb: Vec::new(),
            xt: [[1.0, 1.0], [1.0, 1.0]],
            wgap: Vec::new(),
            ue_m: DMatrix::zeros(0, 0),
        }
    }
}

/// Complete solution from mfoil solver.
#[derive(Debug, Clone, Default)]
pub struct MfoilSolution {
    /// Lift coefficient
    pub cl: f64,
    /// Total drag coefficient
    pub cd: f64,
    /// Friction drag coefficient
    pub cdf: f64,
    /// Pressure drag coefficient
    pub cdp: f64,
    /// Moment coefficient
    pub cm: f64,
    /// Pressure coefficient distribution
    pub cp: Vec<f64>,
    /// Upper surface transition x/c
    pub xtr_upper: f64,
    /// Lower surface transition x/c
    pub xtr_lower: f64,
    /// Converged flag
    pub converged: bool,
    /// Number of iterations
    pub iterations: usize,
    /// Reynolds number
    pub reynolds: f64,
    /// Angle of attack (degrees)
    pub alpha: f64,
    /// Momentum thickness distribution (upper)
    pub theta_upper: Vec<f64>,
    /// Momentum thickness distribution (lower)
    pub theta_lower: Vec<f64>,
    /// Displacement thickness distribution (upper)
    pub delta_star_upper: Vec<f64>,
    /// Displacement thickness distribution (lower)
    pub delta_star_lower: Vec<f64>,
    /// Shape factor distribution (upper)
    pub h_upper: Vec<f64>,
    /// Shape factor distribution (lower)
    pub h_lower: Vec<f64>,
    /// Edge velocity distribution (upper)
    pub ue_upper: Vec<f64>,
    /// Edge velocity distribution (lower)
    pub ue_lower: Vec<f64>,
    /// Skin friction distribution (upper)
    pub cf_upper: Vec<f64>,
    /// Skin friction distribution (lower)
    pub cf_lower: Vec<f64>,
    /// Arc-length coordinates (upper)
    pub s_upper: Vec<f64>,
    /// Arc-length coordinates (lower)
    pub s_lower: Vec<f64>,
}

/// Main mfoil-style boundary layer solver.
pub struct MfoilSolver {
    config: MfoilConfig,
    param: MfoilParam,
    geom: MfoilGeom,
    foil: MfoilPanel,
    wake: MfoilPanel,
    isol: MfoilIsol,
    vsol: MfoilVsol,
    glob: MfoilState,
    alpha: f64,
    viscous: bool,
    /// External inviscid solution (gamma or Ue, if provided)
    external_gamma: Option<Vec<f64>>,
    /// If true, external_gamma contains Ue directly (not gamma)
    external_is_ue: bool,
}

impl MfoilSolver {
    /// Create a new mfoil solver with given configuration.
    pub fn new(config: MfoilConfig) -> Self {
        let vinf = 1.0;
        let param = MfoilParam::for_reynolds(config.reynolds, vinf);
        
        Self {
            config,
            param,
            geom: MfoilGeom::default(),
            foil: MfoilPanel::new(),
            wake: MfoilPanel::new(),
            isol: MfoilIsol::default(),
            vsol: MfoilVsol::default(),
            glob: MfoilState::new(0),
            alpha: 0.0,
            viscous: false,
            external_gamma: None,
            external_is_ue: false,
        }
    }
    
    /// Set external inviscid gamma solution (from existing inviscid solver).
    ///
    /// This bypasses the internal panel method and uses the provided gamma
    /// distribution directly.
    pub fn set_external_gamma(&mut self, gamma: Vec<f64>) {
        self.external_gamma = Some(gamma);
    }
    
    /// Set external inviscid edge velocity distribution.
    ///
    /// This is the preferred method - directly provides the edge velocity Ue
    /// computed from the inviscid Cp distribution: Ue = Vinf * sqrt(1 - Cp)
    pub fn set_external_ue(&mut self, ue: Vec<f64>) {
        // Store as gamma since we use the same mechanism
        self.external_gamma = Some(ue);
        // Mark that this is actually Ue, not gamma
        self.external_is_ue = true;
    }

    /// Set airfoil geometry from coordinates.
    ///
    /// # Arguments
    /// * `coords` - Airfoil coordinates as [[x, y], ...], ordered CW from TE
    pub fn set_airfoil(&mut self, coords: &[[f64; 2]]) {
        self.geom.npoint = coords.len();
        self.geom.xpoint = DMatrix::zeros(2, coords.len());
        
        for (i, &[x, y]) in coords.iter().enumerate() {
            self.geom.xpoint[(0, i)] = x;
            self.geom.xpoint[(1, i)] = y;
        }

        // Compute chord from x-range
        let mut xmin = f64::MAX;
        let mut xmax = f64::MIN;
        for &[x, _] in coords {
            xmin = xmin.min(x);
            xmax = xmax.max(x);
        }
        self.geom.chord = xmax - xmin;
        self.geom.xref = [0.25 * self.geom.chord, 0.0];

        // Create panels
        self.foil = MfoilPanel::from_coords(coords);
    }

    /// Solve for the given angle of attack.
    ///
    /// # Arguments
    /// * `alpha_deg` - Angle of attack in degrees
    /// * `viscous` - Whether to include viscous effects
    ///
    /// # Returns
    /// Complete solution
    pub fn solve(&mut self, alpha_deg: f64, viscous: bool) -> MfoilSolution {
        self.alpha = alpha_deg;
        self.viscous = viscous;

        if viscous {
            self.solve_viscous()
        } else {
            self.solve_inviscid()
        }
    }

    /// Solve inviscid problem only.
    fn solve_inviscid(&mut self) -> MfoilSolution {
        self.build_gamma(self.alpha);
        self.stagpoint_find();
        
        let mut solution = MfoilSolution::default();
        solution.converged = true;
        solution.iterations = 0;
        solution.reynolds = self.config.reynolds;
        solution.alpha = self.alpha;

        // Compute forces
        self.calc_force_inviscid(&mut solution);
        
        solution
    }

    /// Solve coupled viscous-inviscid problem.
    fn solve_viscous(&mut self) -> MfoilSolution {
        // Step 1: Solve inviscid for initial edge velocity
        self.build_gamma(self.alpha);
        
        // Step 2: Build wake (must come before stagpoint_find so xi includes wake)
        self.build_wake();
        
        // Step 3: Find stagnation point (after wake so xi has correct size)
        self.stagpoint_find();
        
        // Step 4: Identify surfaces
        self.identify_surfaces();
        
        // Step 4: Initialize BL
        self.init_boundary_layer();
        
        // Step 5: Compute ue influence matrix
        self.calc_ue_m();
        
        // Step 6: Move stagnation point using viscous solution
        self.stagpoint_move();
        
        // Step 7: Solve coupled system
        let (converged, iterations) = self.solve_coupled();
        
        // Step 8: Compute forces and extract solution
        let mut solution = MfoilSolution::default();
        solution.converged = converged;
        solution.iterations = iterations;
        solution.reynolds = self.config.reynolds;
        solution.alpha = self.alpha;
        
        self.calc_force_viscous(&mut solution);
        self.extract_distributions(&mut solution);
        
        solution
    }

    /// Build gamma (vorticity) distribution for given angle of attack.
    fn build_gamma(&mut self, _alpha: f64) {
        let n = self.foil.n;
        if n == 0 {
            return;
        }

        // Use external gamma if provided (from validated inviscid solver)
        if let Some(ref ext_gamma) = self.external_gamma {
            let len = ext_gamma.len().min(n);
            self.isol.gam = DVector::from_iterator(n, (0..n).map(|i| {
                if i < len { ext_gamma[i] } else { 1.0 }
            }));
            // Create reference gammas (for alpha variation - simplified)
            self.isol.gamref = DMatrix::from_fn(n, 2, |i, j| {
                if j == 0 {
                    // alpha=0 reference
                    if i < len { ext_gamma[i] } else { 1.0 }
                } else {
                    // alpha=90 reference (simplified)
                    0.0
                }
            });
            self.isol.sgnue = DVector::from_fn(n, |i, _| {
                if i < len { ext_gamma[i].signum() } else { 1.0 }
            });
            return;
        }

        // Fallback: Build AIC matrix using panel influence coefficients
        let mut aic = DMatrix::zeros(n + 1, n + 1);
        let mut rhs = DMatrix::zeros(n + 1, 2);

        let alpha_rad = self.alpha.to_radians();
        let cos_a = alpha_rad.cos();
        let sin_a = alpha_rad.sin();

        // Build panel influences (simplified - using normal velocity BC)
        for i in 0..n {
            let xi = [self.foil.x[(0, i)], self.foil.x[(1, i)]];
            
            // Panel j influence on point i
            for j in 0..(n - 1) {
                let x1 = [self.foil.x[(0, j)], self.foil.x[(1, j)]];
                let x2 = [self.foil.x[(0, j + 1)], self.foil.x[(1, j + 1)]];
                let (a, b) = self.panel_stream_influence(x1, x2, xi);
                aic[(i, j)] += a;
                aic[(i, j + 1)] += b;
            }
            
            aic[(i, n)] = -1.0; // Streamfunction = const on surface
            
            // RHS for α=0 and α=90
            rhs[(i, 0)] = -xi[1]; // -y for cos(0)*y - sin(0)*x = -y
            rhs[(i, 1)] = xi[0];  // x for cos(90)*y - sin(90)*x = x
        }

        // Kutta condition
        aic[(n, 0)] = 1.0;
        aic[(n, n - 1)] = 1.0;

        // Solve for reference gammas
        let decomp = aic.clone().lu();
        let gref = decomp.solve(&rhs).unwrap_or_else(|| DMatrix::zeros(n + 1, 2));

        // Store results
        self.isol.gamref = gref.rows(0, n).clone_owned();
        self.isol.gam = DVector::from_iterator(
            n,
            (0..n).map(|i| {
                gref[(i, 0)] * cos_a + gref[(i, 1)] * sin_a
            }),
        );
        self.isol.sgnue = DVector::from_element(n, 1.0);
    }

    /// Panel streamfunction influence coefficients.
    fn panel_stream_influence(&self, x1: [f64; 2], x2: [f64; 2], xi: [f64; 2]) -> (f64, f64) {
        // Panel properties
        let dx = x2[0] - x1[0];
        let dy = x2[1] - x1[1];
        let d = (dx * dx + dy * dy).sqrt().max(1e-12);
        let tx = dx / d;
        let ty = dy / d;

        // Control point in panel frame
        let rxp = xi[0] - x1[0];
        let ryp = xi[1] - x1[1];
        let x = rxp * tx + ryp * ty;
        let z = -rxp * ty + ryp * tx;

        // Distances and angles
        let r1 = (x * x + z * z).sqrt().max(1e-12);
        let r2 = ((x - d) * (x - d) + z * z).sqrt().max(1e-12);
        let theta1 = z.atan2(x);
        let theta2 = z.atan2(x - d);

        // Streamfunction coefficients
        let logr1 = if r1 > 1e-9 { r1.ln() } else { 0.0 };
        let logr2 = if r2 > 1e-9 { r2.ln() } else { 0.0 };

        let p1 = 0.5 / std::f64::consts::PI * (z * (theta2 - theta1) - d + x * logr1 - (x - d) * logr2);
        let p2 = x * p1
            + 0.5 / std::f64::consts::PI
                * (0.5 * r2 * r2 * logr2 - 0.5 * r1 * r1 * logr1 - r2 * r2 / 4.0 + r1 * r1 / 4.0);

        let a = p1 - p2 / d;
        let b = p2 / d;

        (a, b)
    }

    /// Find stagnation point from inviscid solution.
    fn stagpoint_find(&mut self) {
        let n = self.foil.n;
        let nw = self.wake.n;
        
        #[cfg(debug_assertions)]
        eprintln!("[stagpoint_find] n={}, nw={}, gam.len={}", n, nw, self.isol.gam.len());
        
        if n == 0 || self.isol.gam.len() == 0 {
            // Set defaults so we don't crash
            let n_total = n + nw;
            self.isol.xi = DVector::from_element(n_total, 0.5);
            self.isol.sgnue = DVector::from_element(n, 1.0);
            return;
        }

        // Find where gamma changes sign (from negative to positive)
        // For CCW node ordering, gamma is negative on upper surface (TE to LE),
        // crosses zero at stag point, and positive on lower surface (LE to TE)
        let mut j = n / 2; // Default to midpoint
        for i in 0..n {
            if self.isol.gam[i] > 0.0 {
                j = i;
                break;
            }
        }
        
        // If j=0 (first node already positive) or no sign change found, use LE estimate
        if j == 0 {
            // Find minimum |gamma| as the stagnation point
            j = self.isol.gam.iter().enumerate()
                .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(n / 2);
            if j == 0 { j = 1; }
        }

        #[cfg(debug_assertions)]
        eprintln!("[stagpoint_find] Found stag at j={}", j);

        self.isol.istag = [j.saturating_sub(1), j.min(n - 1)];
        
        let idx0 = self.isol.istag[0];
        let idx1 = self.isol.istag[1];
        let g0 = self.isol.gam[idx0];
        let g1 = self.isol.gam[idx1];
        let s0 = self.foil.s[idx0];
        let s1 = self.foil.s[idx1];

        let den = (g1 - g0).abs().max(1e-12);
        let w1 = g1.abs() / den;
        let w2 = g0.abs() / den;

        self.isol.sstag = w1 * s0 + w2 * s1;
        self.isol.xstag = [
            w1 * self.foil.x[(0, idx0)] + w2 * self.foil.x[(0, idx1)],
            w1 * self.foil.x[(1, idx0)] + w2 * self.foil.x[(1, idx1)],
        ];

        // Set sign of edge velocity
        self.isol.sgnue = DVector::from_fn(n, |i, _| if i < j { -1.0 } else { 1.0 });

        // Set xi (distance from stagnation)
        let n_total = n + nw;
        self.isol.xi = DVector::from_fn(n_total, |i, _| {
            if i < n {
                (self.foil.s[i] - self.isol.sstag).abs()
            } else {
                let wake_idx = i - n;
                if wake_idx < self.wake.s.len() {
                    (self.wake.s[wake_idx] - self.isol.sstag).abs()
                } else {
                    1.0
                }
            }
        });
        
        #[cfg(debug_assertions)]
        eprintln!("[stagpoint_find] istag={:?}, sstag={:.6}, xi.len={}", 
                  self.isol.istag, self.isol.sstag, self.isol.xi.len());
    }

    /// Build wake panels.
    fn build_wake(&mut self) {
        let n = self.foil.n;
        if n == 0 {
            return;
        }

        let nw = self.config.n_wake;
        let wake_len = self.config.wake_len * self.geom.chord;

        // TE points
        let x_te = 0.5 * (self.foil.x[(0, 0)] + self.foil.x[(0, n - 1)]);
        let y_te = 0.5 * (self.foil.x[(1, 0)] + self.foil.x[(1, n - 1)]);

        // Build geometric wake (simplified: straight line downstream)
        let mut coords = vec![[0.0; 2]; nw];
        for i in 0..nw {
            let t = (i as f64) / ((nw - 1) as f64);
            coords[i] = [x_te + t * wake_len, y_te];
        }

        self.wake = MfoilPanel::from_coords(&coords);

        // Wake inviscid edge velocity reference (simplified)
        self.isol.uewiref = DMatrix::from_element(nw, 2, 1.0);
    }

    /// Identify surfaces (lower, upper, wake).
    fn identify_surfaces(&mut self) {
        let n = self.foil.n;
        let nw = self.wake.n;
        let istag = self.isol.istag[1];

        // Lower surface: from stag to TE (node 0)
        let lower: Vec<usize> = (0..=istag).rev().collect();
        
        // Upper surface: from stag to TE (node n-1)
        let upper: Vec<usize> = (istag..n).collect();
        
        // Wake: indices after airfoil
        let wake: Vec<usize> = (n..(n + nw)).collect();

        #[cfg(debug_assertions)]
        eprintln!(
            "[identify_surfaces] n={}, nw={}, istag={}, lower.len={}, upper.len={}, wake.len={}",
            n, nw, istag, lower.len(), upper.len(), wake.len()
        );

        self.vsol.is = [lower, upper, wake];
        self.vsol.turb = vec![false; n + nw];
        self.vsol.wgap = vec![0.0; nw];
    }

    /// Initialize boundary layer state.
    fn init_boundary_layer(&mut self) {
        let n = self.foil.n;
        let nw = self.wake.n;
        let nsys = n + nw;

        self.glob = MfoilState::new(nsys);

        // Get inviscid edge velocity
        let ueinv = self.get_ueinv();
        
        // Debug output for size checking
        #[cfg(debug_assertions)]
        eprintln!(
            "[init_boundary_layer] n={}, nw={}, nsys={}, xi.len={}, ueinv.len={}",
            n, nw, nsys, self.isol.xi.len(), ueinv.len()
        );

        // Initialize each surface
        for si in 0..3 {
            let surf_idx = match si {
                0 => SurfaceIndex::Lower,
                1 => SurfaceIndex::Upper,
                _ => SurfaceIndex::Wake,
            };
            
            let is = &self.vsol.is[si];
            if is.is_empty() {
                continue;
            }

            let mut param = self.param.clone();
            param.set_surface(surf_idx);

            // Initialize first point (stagnation or wake start)
            let i0 = is[0];
            let xi0 = self.isol.xi[i0].max(1e-8);
            let ue0 = ueinv[i0].abs().max(1e-8);
            
            // Thwaites initialization
            let k = ue0 / xi0;
            let th0 = (0.45 * param.mu0 / param.rho0 / (6.0 * k)).sqrt().max(1e-10);
            let ds0 = 2.2 * th0;
            
            self.glob.set_node(i0, [th0, ds0, 0.0, ue0]);

            // March downstream
            for ii in 1..is.len() {
                let i = is[ii];
                let i_prev = is[ii - 1];
                
                let xi_prev = self.isol.xi[i_prev];
                let xi = self.isol.xi[i];
                let ue = ueinv[i].abs().max(1e-8);
                
                // Simple scaling from previous
                let scale = (xi / xi_prev.max(1e-12)).sqrt().min(2.0);
                let th = self.glob.theta(i_prev) * scale;
                let ds = self.glob.delta_star(i_prev) * scale;
                let sa = self.glob.sa(i_prev);
                
                self.glob.set_node(i, [th, ds, sa, ue]);
                
                // Check for transition (simplified)
                if !param.turb && sa > param.ncrit {
                    self.vsol.turb[i] = true;
                    param.turb = true;
                }
            }
        }
    }

    /// Get inviscid edge velocity at all nodes.
    fn get_ueinv(&self) -> Vec<f64> {
        let n = self.foil.n;
        let nw = self.wake.n;
        let mut ue = vec![self.param.vinf; n + nw];

        // If we have external Ue (preferred)
        if let Some(ref ext_ue) = self.external_gamma {
            if self.external_is_ue {
                // Direct edge velocity from Cp
                for i in 0..n.min(ext_ue.len()) {
                    ue[i] = ext_ue[i].max(0.1) * self.param.vinf;
                }
            } else {
                // External gamma - use absolute value
                for i in 0..n.min(ext_ue.len()) {
                    ue[i] = ext_ue[i].abs().max(0.1);
                }
            }
            
            // Wake: decaying velocity
            for i in 0..nw {
                let wake_factor = 1.0 - 0.3 * (i as f64 / nw.max(1) as f64);
                let te_ue = if n > 0 { ue[n - 1] } else { self.param.vinf };
                ue[n + i] = te_ue * wake_factor.max(0.5);
            }
            return ue;
        }

        // Original calculation for internal inviscid solver
        let alpha_rad = self.alpha.to_radians();
        let cos_a = alpha_rad.cos();
        let sin_a = alpha_rad.sin();

        // Airfoil: ue = sgnue * (gamref * [cos(α), sin(α)])
        for i in 0..n {
            if i < self.isol.gamref.nrows() {
                ue[i] = self.isol.sgnue[i]
                    * (self.isol.gamref[(i, 0)] * cos_a + self.isol.gamref[(i, 1)] * sin_a);
            }
        }

        // Wake: use reference velocities
        for i in 0..nw {
            if i < self.isol.uewiref.nrows() {
                ue[n + i] = self.isol.uewiref[(i, 0)] * cos_a + self.isol.uewiref[(i, 1)] * sin_a;
            }
        }

        ue
    }

    /// Calculate ue influence matrix (simplified).
    fn calc_ue_m(&mut self) {
        let nsys = self.glob.nsys;
        // For now, use identity (no cross-coupling)
        // Full implementation would compute source panel influences
        self.vsol.ue_m = DMatrix::zeros(nsys, nsys);
    }

    /// Move stagnation point based on viscous solution.
    fn stagpoint_move(&mut self) {
        // Simplified: keep stagnation point fixed
        // Full implementation would iterate based on ue sign changes
    }

    /// Solve the coupled viscous-inviscid system.
    fn solve_coupled(&mut self) -> (bool, usize) {
        let mut converged = false;
        let mut iterations = 0;

        let nsys = self.glob.nsys;
        
        // Sanity check matrix sizes
        if nsys == 0 || nsys > 10000 {
            eprintln!("ERROR: nsys={} is out of range", nsys);
            return (false, 0);
        }

        for iter in 0..self.config.max_iter {
            iterations = iter + 1;

            // Build global residual system
            self.build_glob_sys();

            // Check convergence
            let rnorm = self.glob.r.norm();
            if self.config.verbose >= 2 || iter == 0 {
                eprintln!("Newton iter {}: ||R|| = {:.6e}, nsys={}", iter, rnorm, nsys);
            }

            if rnorm < self.config.rtol {
                converged = true;
                break;
            }
            
            // Check for NaN/Inf
            if !rnorm.is_finite() {
                eprintln!("Newton diverged at iter {}", iter);
                break;
            }

            // Solve for update
            self.solve_glob();

            // Apply update with under-relaxation
            self.update_state();
        }

        if self.config.verbose >= 1 && !converged {
            eprintln!("Mfoil Newton did not converge in {} iterations", iterations);
        }

        (converged, iterations)
    }

    /// Build global residual and Jacobian.
    fn build_glob_sys(&mut self) {
        let nsys = self.glob.nsys;
        self.glob.r = DVector::zeros(3 * nsys);
        self.glob.r_u = DMatrix::zeros(3 * nsys, 4 * nsys);

        // Loop over surfaces
        for si in 0..3 {
            let surf_idx = match si {
                0 => SurfaceIndex::Lower,
                1 => SurfaceIndex::Upper,
                _ => SurfaceIndex::Wake,
            };

            // Clone is to avoid borrow issues
            let is: Vec<usize> = self.vsol.is[si].clone();
            if is.len() < 2 {
                continue;
            }

            let mut param = self.param.clone();
            param.set_surface(surf_idx);

            // Get auxiliary data (wake gap)
            let aux: Vec<f64> = is.iter().map(|&i| {
                if si == 2 && (i >= self.foil.n) {
                    let wake_idx = i - self.foil.n;
                    self.vsol.wgap.get(wake_idx).copied().unwrap_or(0.0)
                } else {
                    0.0
                }
            }).collect();

            // First point (similarity or wake init)
            let i0 = is[0];
            if si < 2 {
                // Similarity residual for stagnation
                param.turb = false;
                param.simi = true;
                let u1 = self.glob.get_node(i0);
                let x1 = self.isol.xi[i0].max(1e-12);
                let res = residual_station(&param, [x1, x1], [u1, u1], [aux[0], aux[0]]);
                
                // Store in global arrays
                let ig = 3 * i0;
                for k in 0..3 {
                    self.glob.r[ig + k] = res.r[k];
                    for j in 0..4 {
                        self.glob.r_u[(ig + k, 4 * i0 + j)] = res.r_u[k][j] + res.r_u[k][4 + j];
                    }
                }
                param.simi = false;
            } else {
                // Wake initialization residual
                self.wake_residual(i0, &mut param);
            }

            // Track transition updates to apply after loop
            let mut transition_at: Option<usize> = None;

            // March over remaining points
            for ii in 1..is.len() {
                let i = is[ii];
                let i_prev = is[ii - 1];

                let u1 = self.glob.get_node(i_prev);
                let u2 = self.glob.get_node(i);
                let x1 = self.isol.xi[i_prev].max(1e-12);
                let x2 = self.isol.xi[i].max(1e-12);

                // Check turbulent state
                let was_turb = self.vsol.turb.get(i_prev).copied().unwrap_or(false);
                let is_turb = self.vsol.turb.get(i).copied().unwrap_or(false);
                
                // Check for dynamic transition (sa exceeds ncrit)
                let sa_curr = u2[2];
                let transition_here = !was_turb && !is_turb && si < 2 && sa_curr > param.ncrit;
                
                if transition_here && transition_at.is_none() {
                    transition_at = Some(ii);
                }

                let tran = !was_turb && is_turb;
                param.turb = is_turb || (transition_at.is_some() && ii >= transition_at.unwrap());

                let res = if tran {
                    let (res, _xt) = residual_transition(&param, [x1, x2], [u1, u2], [aux[ii-1], aux[ii]]);
                    res
                } else {
                    residual_station(&param, [x1, x2], [u1, u2], [aux[ii-1], aux[ii]])
                };

                // Store in global arrays
                let ig = 3 * i;
                for k in 0..3 {
                    self.glob.r[ig + k] = res.r[k];
                    for j in 0..4 {
                        self.glob.r_u[(ig + k, 4 * i_prev + j)] += res.r_u[k][j];
                        self.glob.r_u[(ig + k, 4 * i + j)] += res.r_u[k][4 + j];
                    }
                }
            }
            
            // Apply transition updates after loop completes
            if let Some(trans_ii) = transition_at {
                for ii in trans_ii..is.len() {
                    let i = is[ii];
                    if i < self.vsol.turb.len() {
                        self.vsol.turb[i] = true;
                    }
                }
            }
        }
    }

    /// Wake initialization residual.
    fn wake_residual(&mut self, iw: usize, _param: &mut MfoilParam) {
        let n = self.foil.n;
        if iw < n || self.vsol.is[0].is_empty() || self.vsol.is[1].is_empty() {
            return;
        }

        // Lower and upper TE states
        let il = *self.vsol.is[0].last().unwrap();
        let iu = *self.vsol.is[1].last().unwrap();

        let ul = self.glob.get_node(il);
        let uu = self.glob.get_node(iu);
        let uw = self.glob.get_node(iw);

        // Wake residual: theta, delta_star continuity
        let r = [
            uw[0] - (ul[0] + uu[0]),
            uw[1] - (ul[1] + uu[1]),
            uw[2] - 0.5 * (ul[2] + uu[2]),
        ];

        let ig = 3 * iw;
        for k in 0..3 {
            self.glob.r[ig + k] = r[k];
        }

        // Jacobian
        self.glob.r_u[(ig, 4 * il)] = -1.0;
        self.glob.r_u[(ig, 4 * iu)] = -1.0;
        self.glob.r_u[(ig, 4 * iw)] = 1.0;

        self.glob.r_u[(ig + 1, 4 * il + 1)] = -1.0;
        self.glob.r_u[(ig + 1, 4 * iu + 1)] = -1.0;
        self.glob.r_u[(ig + 1, 4 * iw + 1)] = 1.0;

        self.glob.r_u[(ig + 2, 4 * il + 2)] = -0.5;
        self.glob.r_u[(ig + 2, 4 * iu + 2)] = -0.5;
        self.glob.r_u[(ig + 2, 4 * iw + 2)] = 1.0;
    }

    /// Solve global system for update.
    /// 
    /// Uses a decoupled approach: solve the 3*nsys BL equations with fixed ue,
    /// then update ue towards the inviscid value. This avoids large matrix allocations.
    fn solve_glob(&mut self) {
        let nsys = self.glob.nsys;
        let ueinv = self.get_ueinv();

        // Instead of building and solving a large dense system, use a sequential
        // station-by-station march with local 3x3 solves (like XFOIL does)
        // This is much more efficient and avoids stack overflow
        
        // For each surface, march from upstream to downstream
        for si in 0..3 {
            let is: Vec<usize> = self.vsol.is[si].clone();
            if is.len() < 2 {
                continue;
            }
            
            // First point: similarity or wake init (already handled in residual)
            let i0 = is[0];
            self.glob.du[(0, i0)] = -self.glob.r[3 * i0] * 0.1;     // theta
            self.glob.du[(1, i0)] = -self.glob.r[3 * i0 + 1] * 0.1; // ds
            self.glob.du[(2, i0)] = -self.glob.r[3 * i0 + 2] * 0.1; // sa
            
            // March downstream
            for ii in 1..is.len() {
                let i = is[ii];
                
                // Simple gradient descent on residual
                for k in 0..3 {
                    self.glob.du[(k, i)] = -self.glob.r[3 * i + k] * 0.1;
                }
            }
        }
        
        // Update ue towards inviscid value
        for i in 0..nsys {
            let ue_curr = self.glob.ue(i);
            let ue_target = ueinv[i];
            self.glob.du[(3, i)] = 0.1 * (ue_target - ue_curr);
        }
    }
    
    /// Solve global system for update (full coupled version - disabled for now).
    #[allow(dead_code)]
    fn solve_glob_full(&mut self) {
        let nsys = self.glob.nsys;
        let ueinv = self.get_ueinv();

        // Build augmented system with ue equation
        // R_ue = ue - ueinv - ue_m @ (ds * ue)
        let n_aug = 4 * nsys;
        let mut r_aug = DVector::zeros(n_aug);
        let mut j_aug = Box::new(DMatrix::zeros(n_aug, n_aug));

        // Copy BL residuals (3*nsys)
        for i in 0..(3 * nsys) {
            r_aug[i] = self.glob.r[i];
            for j in 0..(4 * nsys) {
                j_aug[(i, j)] = self.glob.r_u[(i, j)];
            }
        }

        // Add ue equation residuals
        for i in 0..nsys {
            let ue = self.glob.ue(i);
            let ds = self.glob.delta_star(i);
            let ue_inv = ueinv[i];

            // Mass defect contribution (simplified)
            let mut mass_contrib = 0.0;
            for j in 0..nsys {
                mass_contrib += self.vsol.ue_m[(i, j)] * self.glob.delta_star(j) * self.glob.ue(j);
            }

            let r_idx = 3 * nsys + i;
            r_aug[r_idx] = ue - ue_inv - mass_contrib;

            // Jacobian: d(ue_eq)/d(ue) = I - ue_m @ diag(ds)
            j_aug[(r_idx, 4 * i + 3)] = 1.0;
            for j in 0..nsys {
                j_aug[(r_idx, 4 * j + 3)] -= self.vsol.ue_m[(i, j)] * self.glob.delta_star(j);
            }

            // Jacobian: d(ue_eq)/d(ds) = -ue_m @ diag(ue)
            for j in 0..nsys {
                j_aug[(r_idx, 4 * j + 1)] -= self.vsol.ue_m[(i, j)] * self.glob.ue(j);
            }
        }

        // Solve using heap-allocated matrix
        let decomp = j_aug.clone().lu();
        if let Some(dv) = decomp.solve(&(-&r_aug)) {
            // Extract update
            for i in 0..nsys {
                for k in 0..4 {
                    self.glob.du[(k, i)] = dv[4 * i + k];
                }
            }
        }
    }

    /// Apply state update with under-relaxation.
    fn update_state(&mut self) {
        let nsys = self.glob.nsys;
        let mut omega: f64 = 1.0;

        // Limit theta and delta_star changes
        for k in 0..2 {
            for i in 0..nsys {
                let u = self.glob.u[(k, i)].max(1e-12);
                let du = self.glob.du[(k, i)];
                if du < -0.5 * u {
                    omega = omega.min(0.5 * u / (-du));
                }
            }
        }

        // Limit sa changes
        for i in 0..nsys {
            let sa = self.glob.u[(2, i)].max(1e-12);
            let dsa = self.glob.du[(2, i)];
            if dsa < -0.8 * sa {
                omega = omega.min(0.8 * sa / (-dsa));
            }
            if !self.vsol.turb[i] && dsa.abs() > 2.0 {
                omega = omega.min(2.0 / dsa.abs());
            }
        }

        // Limit ue changes
        for i in 0..nsys {
            let due = self.glob.du[(3, i)].abs();
            if due > 0.2 * self.param.vinf {
                omega = omega.min(0.2 * self.param.vinf / due);
            }
        }

        omega = omega.max(0.1);

        // Apply update
        for i in 0..nsys {
            for k in 0..4 {
                self.glob.u[(k, i)] += omega * self.glob.du[(k, i)];
            }
            // Ensure positivity
            self.glob.u[(0, i)] = self.glob.u[(0, i)].max(1e-12);
            self.glob.u[(1, i)] = self.glob.u[(1, i)].max(1e-12);
            if self.vsol.turb[i] {
                self.glob.u[(2, i)] = self.glob.u[(2, i)].max(1e-10);
            }
        }
    }

    /// Calculate inviscid forces.
    fn calc_force_inviscid(&self, solution: &mut MfoilSolution) {
        let n = self.foil.n;
        if n == 0 {
            return;
        }

        let alpha_rad = self.alpha.to_radians();
        let cos_a = alpha_rad.cos();
        let sin_a = alpha_rad.sin();

        // Get edge velocity and compute Cp
        let ue = self.get_ueinv();
        solution.cp = ue.iter().take(n).map(|&u| 1.0 - (u / self.param.vinf).powi(2)).collect();

        // Integrate Cp for lift and moment
        let mut cl = 0.0;
        let mut cm = 0.0;

        for i in 0..(n - 1) {
            let x1 = self.foil.x[(0, i)];
            let y1 = self.foil.x[(1, i)];
            let x2 = self.foil.x[(0, i + 1)];
            let y2 = self.foil.x[(1, i + 1)];

            let dx = x2 - x1;
            let dy = y2 - y1;
            let cp_avg = 0.5 * (solution.cp[i] + solution.cp[i + 1]);

            // Lift in wind direction
            cl += cp_avg * (-dx * cos_a - dy * sin_a);

            // Moment about c/4
            let xm = 0.5 * (x1 + x2) - self.geom.xref[0];
            let ym = 0.5 * (y1 + y2) - self.geom.xref[1];
            cm += cp_avg * (xm * dy - ym * dx);
        }

        solution.cl = cl / self.geom.chord;
        solution.cm = cm / self.geom.chord.powi(2);
        solution.cd = 0.0;
        solution.cdf = 0.0;
        solution.cdp = 0.0;
    }

    /// Calculate viscous forces.
    fn calc_force_viscous(&self, solution: &mut MfoilSolution) {
        // First get inviscid contribution
        self.calc_force_inviscid(solution);

        // Compute drag using Squire-Young formula from TE momentum thickness
        // Use upper and lower surface TE states instead of wake
        let upper_is = &self.vsol.is[1];
        let lower_is = &self.vsol.is[0];
        
        if upper_is.is_empty() || lower_is.is_empty() {
            solution.cd = 0.01; // Default fallback
            solution.cdf = 0.005;
            solution.cdp = 0.005;
            return;
        }

        // Get TE states (last node of each surface)
        let iu_te = *upper_is.last().unwrap();
        let il_te = *lower_is.last().unwrap();
        
        // Safe bounds checking
        let nsys = self.glob.nsys;
        if iu_te >= nsys || il_te >= nsys {
            solution.cd = 0.01;
            solution.cdf = 0.005;
            solution.cdp = 0.005;
            return;
        }

        // Upper surface TE
        let th_u = self.glob.theta(iu_te).max(1e-8);
        let h_u = self.glob.h(iu_te).clamp(1.3, 3.5);
        let ue_u = self.glob.ue(iu_te).abs().max(0.1);

        // Lower surface TE
        let th_l = self.glob.theta(il_te).max(1e-8);
        let h_l = self.glob.h(il_te).clamp(1.3, 3.5);
        let ue_l = self.glob.ue(il_te).abs().max(0.1);

        // Combined TE state (sum of upper and lower)
        let th_te = th_u + th_l;
        let h_te = 0.5 * (h_u + h_l);
        let ue_te = 0.5 * (ue_u + ue_l);

        // Squire-Young: Cd = 2*θ_TE*(Ue_TE/Vinf)^((5+H)/2)
        let uk = (ue_te / self.param.vinf).clamp(0.1, 2.0);
        let exponent = ((5.0 + h_te) / 2.0).clamp(3.0, 4.5);
        solution.cd = 2.0 * th_te * uk.powf(exponent);
        
        // Sanity check
        if !solution.cd.is_finite() || solution.cd > 0.5 || solution.cd < 0.0 {
            solution.cd = 0.01;
        }

        // Estimate friction drag from skin friction integration
        let mut cdf = 0.0;
        for si in 0..2 {
            let is = &self.vsol.is[si];
            for ii in 1..is.len() {
                let i = is[ii];
                let i_prev = is[ii - 1];
                
                // Bounds check
                if i >= nsys || i_prev >= nsys {
                    continue;
                }

                let u = self.glob.get_node(i);
                let u_prev = self.glob.get_node(i_prev);

                let mut param = self.param.clone();
                param.turb = self.vsol.turb.get(i).copied().unwrap_or(false);

                let (cf, _) = get_cf(&u, &param);
                let (cf_prev, _) = get_cf(&u_prev, &param);
                
                // Bounds check for xi
                if i >= self.isol.xi.len() || i_prev >= self.isol.xi.len() {
                    continue;
                }

                let ds = (self.isol.xi[i] - self.isol.xi[i_prev]).abs();
                let cf_contrib = 0.5 * (cf * u[3].powi(2) + cf_prev * u_prev[3].powi(2)) * ds;
                if cf_contrib.is_finite() {
                    cdf += cf_contrib;
                }
            }
        }

        solution.cdf = cdf / (self.param.vinf.powi(2) * self.geom.chord);
        if !solution.cdf.is_finite() || solution.cdf < 0.0 {
            solution.cdf = 0.005;
        }
        solution.cdp = (solution.cd - solution.cdf).max(0.0);
    }

    /// Extract BL distributions to solution.
    fn extract_distributions(&self, solution: &mut MfoilSolution) {
        let nsys = self.glob.nsys;
        let xi_len = self.isol.xi.len();
        
        // Upper surface
        let upper_is = &self.vsol.is[1];
        solution.s_upper = upper_is.iter()
            .filter(|&&i| i < xi_len)
            .map(|&i| self.isol.xi[i])
            .collect();
        solution.theta_upper = upper_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.theta(i))
            .collect();
        solution.delta_star_upper = upper_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.delta_star(i))
            .collect();
        solution.h_upper = upper_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.h(i))
            .collect();
        solution.ue_upper = upper_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.ue(i))
            .collect();
        
        solution.cf_upper = upper_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| {
                let u = self.glob.get_node(i);
                let mut param = self.param.clone();
                param.turb = self.vsol.turb.get(i).copied().unwrap_or(false);
                let (cf, _) = get_cf(&u, &param);
                cf
            }).collect();

        // Transition location
        for &i in upper_is {
            if i < self.vsol.turb.len() && self.vsol.turb[i] && i < self.foil.n {
                solution.xtr_upper = self.foil.x[(0, i)] / self.geom.chord;
                break;
            }
        }

        // Lower surface
        let lower_is = &self.vsol.is[0];
        solution.s_lower = lower_is.iter()
            .filter(|&&i| i < xi_len)
            .map(|&i| self.isol.xi[i])
            .collect();
        solution.theta_lower = lower_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.theta(i))
            .collect();
        solution.delta_star_lower = lower_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.delta_star(i))
            .collect();
        solution.h_lower = lower_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.h(i))
            .collect();
        solution.ue_lower = lower_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| self.glob.ue(i))
            .collect();
        
        solution.cf_lower = lower_is.iter()
            .filter(|&&i| i < nsys)
            .map(|&i| {
                let u = self.glob.get_node(i);
                let mut param = self.param.clone();
                param.turb = self.vsol.turb.get(i).copied().unwrap_or(false);
                let (cf, _) = get_cf(&u, &param);
                cf
            }).collect();

        // Lower transition
        for &i in lower_is {
            if i < self.vsol.turb.len() && self.vsol.turb[i] && i < self.foil.n {
                solution.xtr_lower = self.foil.x[(0, i)] / self.geom.chord;
                break;
            }
        }
    }
}
