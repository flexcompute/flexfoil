//! State structures for mfoil-style boundary layer solver.

use nalgebra::{DMatrix, DVector};

/// Surface index for mfoil surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceIndex {
    /// Lower surface (from stagnation to lower TE)
    Lower = 0,
    /// Upper surface (from stagnation to upper TE)
    Upper = 1,
    /// Wake (from TE downstream)
    Wake = 2,
}

/// Parameters for mfoil closure relations and BL calculations.
#[derive(Debug, Clone)]
pub struct MfoilParam {
    // Operating conditions
    /// Freestream velocity
    pub vinf: f64,
    /// Freestream Mach number
    pub minf: f64,
    /// Viscosity at freestream
    pub mu0: f64,
    /// Density at stagnation
    pub rho0: f64,
    /// Stagnation enthalpy
    pub h0: f64,
    /// Sutherland temperature ratio Ts/Tref
    pub tsrat: f64,
    /// Ratio of specific heats
    pub gam: f64,
    /// Karman-Tsien beta = sqrt(1 - Minf^2)
    pub ktb: f64,
    /// Karman-Tsien lambda = Minf^2 / (1 + KTb)^2
    pub ktl: f64,

    // Viscous parameters
    /// Critical amplification factor
    pub ncrit: f64,
    /// Scales the uq term in shear lag equation
    pub cuq: f64,
    /// Wall/wake dissipation length ratio
    pub dlr: f64,
    /// Shear lag constant
    pub slag_k: f64,

    // Initial Ctau after transition
    /// Ctau constant
    pub ctau_c: f64,
    /// Ctau exponent
    pub ctau_e: f64,

    // G-beta locus constants
    /// G-beta A constant
    pub ga: f64,
    /// G-beta B constant
    pub gb: f64,
    /// G-beta C constant
    pub gc: f64,

    // Station flags
    /// True if at a similarity station
    pub simi: bool,
    /// True if station is turbulent
    pub turb: bool,
    /// True if in wake
    pub wake: bool,
}

impl Default for MfoilParam {
    fn default() -> Self {
        Self {
            vinf: 1.0,
            minf: 0.0,
            mu0: 1.0,
            rho0: 1.0,
            h0: 0.0,
            tsrat: 0.35,
            gam: 1.4,
            ktb: 1.0,
            ktl: 0.0,

            ncrit: 9.0,
            cuq: 1.0,
            dlr: 0.9,
            slag_k: 5.6,

            ctau_c: 1.8,
            ctau_e: 3.3,

            ga: 6.7,
            gb: 0.75,
            gc: 18.0,

            simi: false,
            turb: false,
            wake: false,
        }
    }
}

impl MfoilParam {
    /// Create parameters for a given Reynolds number.
    pub fn for_reynolds(re: f64, vinf: f64) -> Self {
        // mu0 = rho0 * vinf * chord / Re (assuming chord = 1, rho0 = 1)
        let mu0 = vinf / re;
        Self {
            vinf,
            mu0,
            ..Default::default()
        }
    }

    /// Create parameters with compressibility.
    pub fn with_mach(re: f64, vinf: f64, mach: f64) -> Self {
        let mut p = Self::for_reynolds(re, vinf);
        p.minf = mach;
        if mach > 0.0 {
            p.ktb = (1.0 - mach * mach).sqrt();
            p.ktl = mach * mach / (1.0 + p.ktb).powi(2);
            p.h0 = vinf * vinf / (p.gam - 1.0); // Simplified stagnation enthalpy
        }
        p
    }

    /// Set station flags for a given surface.
    pub fn set_surface(&mut self, si: SurfaceIndex) {
        self.wake = si == SurfaceIndex::Wake;
        self.turb = self.wake; // Wake is always turbulent
        self.simi = false;
    }
}

/// Mfoil global state vector.
///
/// The state at all nodes is stored as a 4 x Nsys matrix where:
/// - Row 0: θ (momentum thickness)
/// - Row 1: δ* (displacement thickness)
/// - Row 2: sa (amplification or sqrt(ctau))
/// - Row 3: Ue (edge velocity)
#[derive(Debug, Clone)]
pub struct MfoilState {
    /// Primary states: [θ, δ*, sa, Ue] at each node (4 x Nsys)
    pub u: DMatrix<f64>,
    /// Proposed state update (4 x Nsys)
    pub du: DMatrix<f64>,
    /// Angle of attack update (for cl-constrained mode)
    pub dalpha: f64,
    /// Total number of nodes (airfoil + wake)
    pub nsys: usize,
    /// Converged flag
    pub conv: bool,
    /// Global residual vector (3 x Nsys, flattened)
    pub r: DVector<f64>,
    /// Residual Jacobian w.r.t. U (sparse, 3*Nsys x 4*Nsys)
    pub r_u: DMatrix<f64>,
    /// Residual Jacobian w.r.t. x (sparse, 3*Nsys x Nsys)
    pub r_x: DMatrix<f64>,
    /// Global Jacobian including Ue equation (4*Nsys x 4*Nsys)
    pub r_v: DMatrix<f64>,
}

impl MfoilState {
    /// Create a new state with the given number of nodes.
    pub fn new(nsys: usize) -> Self {
        Self {
            u: DMatrix::zeros(4, nsys),
            du: DMatrix::zeros(4, nsys),
            dalpha: 0.0,
            nsys,
            conv: false,
            r: DVector::zeros(3 * nsys),
            r_u: DMatrix::zeros(3 * nsys, 4 * nsys),
            r_x: DMatrix::zeros(3 * nsys, nsys),
            r_v: DMatrix::zeros(4 * nsys, 4 * nsys),
        }
    }

    /// Get state at node i as [θ, δ*, sa, Ue].
    pub fn get_node(&self, i: usize) -> [f64; 4] {
        [
            self.u[(0, i)],
            self.u[(1, i)],
            self.u[(2, i)],
            self.u[(3, i)],
        ]
    }

    /// Set state at node i from [θ, δ*, sa, Ue].
    pub fn set_node(&mut self, i: usize, state: [f64; 4]) {
        self.u[(0, i)] = state[0];
        self.u[(1, i)] = state[1];
        self.u[(2, i)] = state[2];
        self.u[(3, i)] = state[3];
    }

    /// Get theta at node i.
    #[inline]
    pub fn theta(&self, i: usize) -> f64 {
        self.u[(0, i)]
    }

    /// Get delta_star at node i.
    #[inline]
    pub fn delta_star(&self, i: usize) -> f64 {
        self.u[(1, i)]
    }

    /// Get sa (amplification or sqrt(ctau)) at node i.
    #[inline]
    pub fn sa(&self, i: usize) -> f64 {
        self.u[(2, i)]
    }

    /// Get edge velocity at node i.
    #[inline]
    pub fn ue(&self, i: usize) -> f64 {
        self.u[(3, i)]
    }

    /// Get H = delta_star / theta at node i.
    pub fn h(&self, i: usize) -> f64 {
        let th = self.theta(i).max(1e-12);
        self.delta_star(i) / th
    }
}

/// Geometry representation for mfoil.
#[derive(Debug, Clone)]
pub struct MfoilGeom {
    /// Chord length
    pub chord: f64,
    /// Wake extent length in chords
    pub wakelen: f64,
    /// Number of geometry representation points
    pub npoint: usize,
    /// Airfoil name
    pub name: String,
    /// Point coordinates [2 x npoint] (x in row 0, y in row 1)
    pub xpoint: DMatrix<f64>,
    /// Moment reference point [x, y]
    pub xref: [f64; 2],
}

impl Default for MfoilGeom {
    fn default() -> Self {
        Self {
            chord: 1.0,
            wakelen: 1.0,
            npoint: 0,
            name: String::from("noname"),
            xpoint: DMatrix::zeros(2, 0),
            xref: [0.25, 0.0],
        }
    }
}

/// Panel representation for mfoil.
#[derive(Debug, Clone)]
pub struct MfoilPanel {
    /// Number of nodes
    pub n: usize,
    /// Node coordinates [2 x N]
    pub x: DMatrix<f64>,
    /// Arc-length values at nodes
    pub s: DVector<f64>,
    /// Tangent vectors [2 x N] (dx/ds, dy/ds)
    pub t: DMatrix<f64>,
}

impl MfoilPanel {
    /// Create empty panel.
    pub fn new() -> Self {
        Self {
            n: 0,
            x: DMatrix::zeros(2, 0),
            s: DVector::zeros(0),
            t: DMatrix::zeros(2, 0),
        }
    }

    /// Create panel from node coordinates.
    pub fn from_coords(coords: &[[f64; 2]]) -> Self {
        let n = coords.len();
        let mut x = DMatrix::zeros(2, n);
        let mut s = DVector::zeros(n);
        let mut t = DMatrix::zeros(2, n);

        // Set coordinates
        for (i, &[xi, yi]) in coords.iter().enumerate() {
            x[(0, i)] = xi;
            x[(1, i)] = yi;
        }

        // Compute arc-lengths
        let mut total_s = 0.0;
        s[0] = 0.0;
        for i in 1..n {
            let dx = x[(0, i)] - x[(0, i - 1)];
            let dy = x[(1, i)] - x[(1, i - 1)];
            total_s += (dx * dx + dy * dy).sqrt();
            s[i] = total_s;
        }

        // Compute tangents (central differences for interior, one-sided at ends)
        for i in 0..n {
            let (ti, ti1) = if i == 0 {
                (0, 1)
            } else if i == n - 1 {
                (n - 2, n - 1)
            } else {
                (i - 1, i + 1)
            };

            let ds = (s[ti1] - s[ti]).max(1e-12);
            t[(0, i)] = (x[(0, ti1)] - x[(0, ti)]) / ds;
            t[(1, i)] = (x[(1, ti1)] - x[(1, ti)]) / ds;

            // Normalize
            let tnorm = (t[(0, i)].powi(2) + t[(1, i)].powi(2)).sqrt().max(1e-12);
            t[(0, i)] /= tnorm;
            t[(1, i)] /= tnorm;
        }

        Self { n, x, s, t }
    }
}

impl Default for MfoilPanel {
    fn default() -> Self {
        Self::new()
    }
}
