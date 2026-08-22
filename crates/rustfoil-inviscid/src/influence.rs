//! Influence coefficient calculations (XFOIL's PSILIN).
//!
//! This module computes the stream function ψ at a field point due to all vortex panels,
//! and the sensitivity ∂ψ/∂γ for each node. This forms the core of the panel method.
//!
//! # Linear Vorticity Panels
//!
//! XFOIL uses panels with linearly-varying vorticity (node-based unknowns).
//! For a panel from node JO to JP with γ varying linearly:
//!
//! ```text
//! γ(s) = γ_JO + (γ_JP - γ_JO) * (s - s_JO) / (s_JP - s_JO)
//! ```
//!
//! # Sum/Difference Formulation
//!
//! The key to XFOIL's efficiency is the PSIS/PSID formulation:
//! - PSIS: coefficient of (γ_JO + γ_JP), the "sum" term
//! - PSID: coefficient of (γ_JP - γ_JO), the "difference" term
//!
//! This allows influence coefficients to be computed incrementally.
//!
//! # One contour per element
//!
//! Panel JP is not `JO + 1` modulo the node count. It is the next node **within
//! one element**, so the panel starting at an element's last node ends at that
//! element's own first node. With several elements concatenated into one node
//! array, modulo-n closure would instead join one element's last node to the
//! next element's first, putting a panel across the physical gap between them.
//! Every trailing-edge quantity here — SCS, SDS, SEPS, the sharpness test — is
//! likewise the owning element's own.
//!
//! # XFOIL Reference
//!
//! - `xpanel.f`: PSILIN subroutine (lines 99-800)

use crate::geometry::{AirfoilGeometry, ConfigGeometry};
use crate::{QOPI, HOPI};
use core::ops::Range;
use nalgebra::DMatrix;
use std::f64::consts::PI;

/// Result from computing influence coefficients at a single field point.
#[derive(Debug, Clone)]
pub struct PsilinResult {
    /// Stream function value at the field point
    pub psi: f64,
    /// Tangential velocity at α=0° (contribution to QTAN1)
    pub qtan1: f64,
    /// Tangential velocity at α=90° (contribution to QTAN2)
    pub qtan2: f64,
    /// ∂ψ/∂γⱼ for each node j (influence coefficient array)
    pub dzdg: Vec<f64>,
    /// ∂ψ/∂σⱼ for each node j (source influence, if computed)
    pub dzdm: Vec<f64>,
    /// ∂Qtan/∂σⱼ for each node j (tangential velocity derivative w.r.t. source)
    pub dqdm: Vec<f64>,
    /// ∂Qtan/∂γⱼ for each node j (tangential velocity derivative w.r.t. vorticity)
    /// Used for CIJ in QDCALC and for computing QTAN1/QTAN2.
    pub dqdg: Vec<f64>,
}

/// Intermediate values from a single panel's contribution.
/// Used for testing against XFOIL's internal values.
#[derive(Debug, Clone, Default)]
pub struct PanelContribution {
    /// Sum term (coefficient of γ_JO + γ_JP)
    pub psis: f64,
    /// Difference term (coefficient of γ_JP - γ_JO)
    pub psid: f64,
    /// Log term for endpoint 1
    pub g1: f64,
    /// Log term for endpoint 2
    pub g2: f64,
    /// Angle term for endpoint 1
    pub t1: f64,
    /// Angle term for endpoint 2
    pub t2: f64,
    /// Local x-coordinate of point relative to panel start
    pub x1: f64,
    /// Local x-coordinate of point relative to panel end
    pub x2: f64,
    /// Local y-coordinate (perpendicular distance to panel line)
    pub yy: f64,
}

// ===========================================================================
// Panel connectivity
// ===========================================================================

/// The node arrays the influence kernels read.
///
/// Borrowed rather than owned, so one kernel serves a single
/// [`AirfoilGeometry`] and a multi-element [`ConfigGeometry`] without copying
/// either.
#[derive(Debug, Clone, Copy)]
struct Nodes<'a> {
    x: &'a [f64],
    y: &'a [f64],
    apanel: &'a [f64],
}

impl Nodes<'_> {
    /// Total number of nodes across every element.
    #[inline]
    fn len(&self) -> usize {
        self.x.len()
    }
}

/// One element's panel connectivity and trailing-edge data.
///
/// # Closure is per element
/// Panels close through [`next_node`](Self::next_node), which wraps **within
/// this element**: the panel starting at element k's last node ends at element
/// k's own first node, never at element k+1's. That is the rule
/// [`Layout::next_node`](rustfoil_core::layout::Layout::next_node) implements,
/// and `element_panels_match_the_layout_closure_rule` holds the two against each
/// other node by node on a real three-element geometry. For a single element it
/// reduces to `(jo + 1) % n`.
///
/// # Trailing edge
/// By the node ordering convention an element's first node is its
/// upper-surface trailing edge and its last node its lower-surface trailing
/// edge, so `first` and `last` are that element's own TE node pair and the panel
/// between them is that element's own TE panel. `sharp`, `scs`, `sds` and `seps`
/// all describe that trailing edge and no other.
#[derive(Debug, Clone, Copy, PartialEq)]
struct ElementPanels {
    /// Global index of this element's first node (its upper-surface TE).
    first: usize,
    /// Global index of this element's last node (its lower-surface TE).
    last: usize,
    /// Whether this element's own trailing edge is sharp.
    sharp: bool,
    /// This element's TE panel source coefficient (SCS in XFOIL).
    scs: f64,
    /// This element's TE panel vortex coefficient (SDS in XFOIL).
    sds: f64,
    /// TE gap below which this element's TE panel is skipped as closed
    /// (SEPS in XFOIL), from this element's own arc length.
    seps: f64,
}

impl ElementPanels {
    /// The one element of a single-body geometry.
    fn single(geom: &AirfoilGeometry) -> Self {
        let (scs, sds) = geom.te_coefficients();
        Self {
            first: 0,
            last: geom.n - 1,
            sharp: geom.sharp,
            scs,
            sds,
            seps: geom.total_arc_length() * 1e-5,
        }
    }

    /// One element of a configuration, with its own trailing edge and its own
    /// arc length.
    fn of_configuration(config: &ConfigGeometry, element: usize) -> Self {
        let geom = config.element(element);
        let (scs, sds) = geom.te_coefficients();
        Self {
            first: geom.start,
            last: geom.end() - 1,
            sharp: geom.sharp,
            scs,
            sds,
            seps: config.element_arc_length(element) * 1e-5,
        }
    }

    /// Every element of a configuration, in configuration order.
    fn all_of(config: &ConfigGeometry) -> Vec<Self> {
        (0..config.n_elements())
            .map(|element| Self::of_configuration(config, element))
            .collect()
    }

    /// The node that closes the panel starting at `global`.
    ///
    /// The replacement for `(jo + 1) % n`.
    #[inline]
    fn next_node(&self, global: usize) -> usize {
        if global == self.last {
            self.first
        } else {
            global + 1
        }
    }

    /// The nodes starting an ordinary surface panel: every node of this element
    /// except its last, whose panel is the TE panel and is treated separately.
    #[inline]
    fn surface_panels(&self) -> Range<usize> {
        self.first..self.last
    }

    /// The node before `global`, held at this element's first node so the
    /// source-gradient stencil never reaches into another element.
    #[inline]
    fn stencil_back(&self, global: usize) -> usize {
        if global == self.first {
            global
        } else {
            global - 1
        }
    }

    /// The node after `global`, held at this element's last node so the
    /// source-gradient stencil never reaches into another element.
    #[inline]
    fn stencil_forward(&self, global: usize) -> usize {
        if global == self.last {
            global
        } else {
            global + 1
        }
    }
}

/// Every element paired with each of its surface-panel start nodes, in
/// configuration order.
///
/// The panel loops iterate over this rather than `0..n` so that the element a
/// panel belongs to — and therefore which node closes it — travels with the
/// panel index.
fn surface_panels(
    elements: &[ElementPanels],
) -> impl Iterator<Item = (&ElementPanels, usize)> + '_ {
    elements
        .iter()
        .flat_map(|element| element.surface_panels().map(move |jo| (element, jo)))
}

/// Compute influence coefficients at a field point (i, xi, yi).
///
/// This is XFOIL's PSILIN subroutine, computing:
/// - The stream function ψ at (xi, yi) due to all panels
/// - The derivative ∂ψ/∂γⱼ for each node j
///
/// # Arguments
///
/// * `geom` - Airfoil geometry
/// * `i` - Index of the field point (for singularity handling)
/// * `xi`, `yi` - Field point coordinates
///
/// # Returns
///
/// A `PsilinResult` containing the stream function and influence coefficients.
pub fn psilin(geom: &AirfoilGeometry, i: usize, xi: f64, yi: f64) -> PsilinResult {
    psilin_internal(geom, i, xi, yi, false)
}

/// Compute influence coefficients including source sensitivity (SIGLIN).
pub fn psilin_with_sources(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinResult {
    psilin_internal(geom, i, xi, yi, true)
}

/// Compute influence coefficients with DQDM (tangential velocity derivatives).
///
/// This is XFOIL's PSILIN with SIGLIN=.TRUE., computing:
/// - dzdm: ∂ψ/∂σⱼ (streamfunction derivatives w.r.t. source strength)
/// - dqdm: ∂Qtan/∂σⱼ (tangential velocity derivatives w.r.t. source strength)
///
/// The DQDM computation uses half-panel decomposition with PSNI/PDNI normal
/// derivatives, matching XFOIL's QDCALC implementation.
///
/// # Arguments
///
/// * `geom` - Airfoil geometry
/// * `i` - Index of the field point (for singularity handling)
/// * `xi`, `yi` - Field point coordinates
/// * `nxi`, `nyi` - Field point normal vector (for tangential velocity computation)
///
/// # Returns
///
/// A `PsilinResult` containing dzdg, dzdm, and dqdm arrays.
pub fn psilin_with_dqdm(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
) -> PsilinResult {
    psilin_dqdm_kernel(
        Nodes { x: &geom.x, y: &geom.y, apanel: &geom.apanel },
        &[ElementPanels::single(geom)],
        i,
        xi,
        yi,
        nxi,
        nyi,
    )
}

/// [`psilin_with_dqdm`] over a whole configuration.
///
/// `i` and the returned arrays are indexed globally, spanning every element in
/// configuration order. Each element contributes its own panels closed onto its
/// own nodes and its own trailing edge, so no panel crosses the gap between two
/// elements.
pub fn psilin_config_with_dqdm(
    config: &ConfigGeometry,
    i: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
) -> PsilinResult {
    psilin_dqdm_kernel(
        Nodes { x: config.x(), y: config.y(), apanel: config.apanel() },
        &ElementPanels::all_of(config),
        i,
        xi,
        yi,
        nxi,
        nyi,
    )
}

fn psilin_dqdm_kernel(
    nodes: Nodes<'_>,
    elements: &[ElementPanels],
    i: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
) -> PsilinResult {
    let n = nodes.len();

    // Initialize influence coefficient arrays
    let mut dzdg = vec![0.0; n];
    let mut dzdm = vec![0.0; n];
    let mut dqdm = vec![0.0; n];
    let mut dqdg = vec![0.0; n];

    // Initialize accumulated values
    let psi = 0.0;
    let qtan1 = 0.0;
    let qtan2 = 0.0;

    // Determine if field point is on a surface (affects SGN)
    let sgn_surface = if i < n { 1.0 } else { 0.0 };

    // Loop over every element's surface panels. An element's last node starts
    // that element's TE panel, handled separately after this loop.
    for (element, jo) in surface_panels(elements) {
        let jp = element.next_node(jo); // closes within this element

        // Panel endpoints
        let x_jo = nodes.x[jo];
        let y_jo = nodes.y[jo];
        let x_jp = nodes.x[jp];
        let y_jp = nodes.y[jp];

        // Panel vector and length
        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        // Skip zero-length panels
        if ds_sq < 1e-24 {
            continue;
        }

        let dso = ds_sq.sqrt();
        let dsio = 1.0 / dso;

        // Unit tangent vector along panel
        let sx = dx * dsio;
        let sy = dy * dsio;

        // Vector from panel endpoints to field point
        let rx1 = xi - x_jo;
        let ry1 = yi - y_jo;
        let rx2 = xi - x_jp;
        let ry2 = yi - y_jp;

        // Transform to panel-local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        // Squared distances to endpoints
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // SGN reflection for branch cuts
        let sgn = if sgn_surface != 0.0 {
            1.0
        } else if yy >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let pi_offset = (0.5 - 0.5 * sgn) * PI;

        let (g1, t1) = if i != jo && rs1 > 1e-20 {
            (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if i != jp && rs2 > 1e-20 {
            (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };

        // Panel angle
        let apan = nodes.apanel[jo];

        // Midpoint quantities for half-panel decomposition
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = if rs0 > 0.0 { rs0.ln() } else { 0.0 };
        let t0 = (sgn * x0).atan2(sgn * yy) + pi_offset;

        // Transform normal to panel coordinates (for DQDM)
        let x1i = sx * nxi + sy * nyi;
        let x2i = x1i; // Same for both endpoints
        let yyi = sx * nyi - sy * nxi;

        // Neighbouring nodes for the source gradient stencil, held inside this
        // element at its first and last node
        let jm = element.stencil_back(jo);
        let jq = element.stencil_forward(jp);

        // ============ First half-panel (1-0) ============
        {
            let dxinv = safe_inv(x1 - x0);

            // Panel integrals for streamfunction (PSUM, PDIF)
            let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif = ((x1 + x0) * psum
                + rs1 * (t1 - apan)
                - rs0 * (t0 - apan)
                + (x0 - x1) * yy)
                * dxinv;

            // Source strength derivatives
            let dsm =
                ((nodes.x[jp] - nodes.x[jm]).powi(2) + (nodes.y[jp] - nodes.y[jm]).powi(2)).sqrt();
            let dsim = safe_inv(dsm);

            // Accumulate dPsi/dm (DZDM)
            dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
            dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
            dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            // Normal derivatives for tangential velocity (PSNI, PDNI)
            let psx1 = -(t1 - apan);
            let psx0 = t0 - apan;
            let psyy = 0.5 * (g1 - g0);

            let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
            let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
            let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

            let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

            // Accumulate dQtan/dm (DQDM)
            dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
            dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
            dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));
        }

        // ============ Second half-panel (0-2) ============
        {
            let dxinv = safe_inv(x0 - x2);

            let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif = ((x0 + x2) * psum
                + rs0 * (t0 - apan)
                - rs2 * (t2 - apan)
                + (x2 - x0) * yy)
                * dxinv;

            let dsp =
                ((nodes.x[jq] - nodes.x[jo]).powi(2) + (nodes.y[jq] - nodes.y[jo]).powi(2)).sqrt();
            let dsip = safe_inv(dsp);

            // Accumulate dPsi/dm (DZDM)
            dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
            dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);

            // Normal derivatives for tangential velocity
            let psx0 = -(t0 - apan);
            let psx2 = t2 - apan;
            let psyy = 0.5 * (g0 - g2);

            let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
            let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
            let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

            let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

            // Accumulate dQtan/dm (DQDM)
            dqdm[jo] += QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            dqdm[jp] += QOPI * (psni * dsio - pdni * dsio);
            dqdm[jq] += QOPI * (psni * dsip + pdni * dsip);
        }

        // Compute PSIS and PSID for vortex contribution (same as psilin_internal)
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);
        dzdg[jo] += QOPI * (contrib.psis - contrib.psid);
        dzdg[jp] += QOPI * (contrib.psis + contrib.psid);

        // ============ Vortex PSNI/PDNI for DQDG (Fortran PSILIN lines 348-384) ============
        // These are the normal derivatives of the vortex panel contribution,
        // distinct from the source normal derivatives used for DQDM above.
        {
            let dxinv_v = safe_inv(x1 - x2);

            // Vortex panel x/y derivatives (Fortran lines 353-359)
            let psx1_v = 0.5 * g1;
            let psx2_v = -0.5 * g2;
            let psyy_v = t1 - t2;

            let pdx1_v = ((x1 + x2) * psx1_v + contrib.psis - x1 * g1 - contrib.psid) * dxinv_v;
            let pdx2_v = ((x1 + x2) * psx2_v + contrib.psis + x2 * g2 + contrib.psid) * dxinv_v;
            let pdyy_v = ((x1 + x2) * psyy_v - yy * (g1 - g2)) * dxinv_v;

            // Normal derivatives (Fortran lines 376-377)
            let psni_v = psx1_v * x1i + psx2_v * x2i + psyy_v * yyi;
            let pdni_v = pdx1_v * x1i + pdx2_v * x2i + pdyy_v * yyi;

            // Accumulate dQtan/dGam (Fortran lines 383-384)
            dqdg[jo] += QOPI * (psni_v - pdni_v);
            dqdg[jp] += QOPI * (psni_v + pdni_v);
        }
    }

    // TE panel special treatment (XFOIL lines 407-438), once per element on
    // that element's own TE node pair: its last node closing onto its first.
    for element in elements {
        if element.sharp {
            continue;
        }

        let jo = element.last;
        let jp = element.first;
        let (scs, sds) = (element.scs, element.sds);
        let seps = element.seps;

        let x_jo = nodes.x[jo];
        let y_jo = nodes.y[jo];
        let x_jp = nodes.x[jp];
        let y_jp = nodes.y[jp];

        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        if ds_sq > seps * seps {
            let dso = ds_sq.sqrt();
            let dsio = 1.0 / dso;

            let sx = dx * dsio;
            let sy = dy * dsio;

            let rx1 = xi - x_jo;
            let ry1 = yi - y_jo;
            let rx2 = xi - x_jp;
            let ry2 = yi - y_jp;

            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;

            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;

            let apan = nodes.apanel[jo];

            let sgn = if sgn_surface != 0.0 {
                1.0
            } else if yy >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let pi_offset = (0.5 - 0.5 * sgn) * PI;

            let (g1, t1_raw) = if i != jo && rs1 > 1e-20 {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
            } else {
                (0.0, 0.0)
            };

            let (g2, t2_raw) = if i != jp && rs2 > 1e-20 {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
            } else {
                (0.0, 0.0)
            };

            // PSIG/PGAM formulas for TE panel
            let psig = 0.5 * yy * (g1 - g2) + x2 * (t2_raw - apan) - x1 * (t1_raw - apan);
            let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1_raw - t2_raw);

            // TE panel dPsi/dGam (Fortran lines 434-438)
            dzdg[jo] += HOPI * (-psig * scs + pgam * sds) * 0.5;
            dzdg[jp] += HOPI * (psig * scs - pgam * sds) * 0.5;

            // TE panel DQDG contribution (Fortran lines 411-447)
            // Transform normal to panel coordinates
            let x1i = sx * nxi + sy * nyi;
            let x2i = x1i;
            let yyi = sx * nyi - sy * nxi;

            // Normal derivatives of PSIG and PGAM for TE panel (Fortran lines 411-416)
            let psig_x1 = -(t1_raw - apan);
            let psig_x2 = t2_raw - apan;
            let psig_yy = 0.5 * (g1 - g2);
            let pgam_x1 = 0.5 * g1;
            let pgam_x2 = -0.5 * g2;
            let pgam_yy = t1_raw - t2_raw;

            // Fortran lines 418-419
            let psig_ni = psig_x1 * x1i + psig_x2 * x2i + psig_yy * yyi;
            let pgam_ni = pgam_x1 * x1i + pgam_x2 * x2i + pgam_yy * yyi;

            // DQDG(JO) -= HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
            // DQDG(JP) += HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
            let te_dqdg = HOPI * (psig_ni * 0.5 * scs - pgam_ni * 0.5 * sds);
            dqdg[jo] -= te_dqdg;
            dqdg[jp] += te_dqdg;
        }
    }

    PsilinResult {
        psi,
        qtan1,
        qtan2,
        dzdg,
        dzdm,
        dqdm,
        dqdg,
    }
}

fn psilin_internal(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
    compute_sources: bool,
) -> PsilinResult {
    psilin_kernel(
        Nodes { x: &geom.x, y: &geom.y, apanel: &geom.apanel },
        &[ElementPanels::single(geom)],
        i,
        xi,
        yi,
        compute_sources,
    )
}

/// [`psilin`] over a whole configuration.
///
/// `i` and the returned arrays are indexed globally, spanning every element in
/// configuration order. Each element contributes its own panels closed onto its
/// own nodes and its own trailing edge, so no panel crosses the gap between two
/// elements.
pub fn psilin_config(config: &ConfigGeometry, i: usize, xi: f64, yi: f64) -> PsilinResult {
    psilin_config_internal(config, i, xi, yi, false)
}

/// [`psilin_with_sources`] over a whole configuration.
///
/// Indexed as [`psilin_config`].
pub fn psilin_config_with_sources(
    config: &ConfigGeometry,
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinResult {
    psilin_config_internal(config, i, xi, yi, true)
}

fn psilin_config_internal(
    config: &ConfigGeometry,
    i: usize,
    xi: f64,
    yi: f64,
    compute_sources: bool,
) -> PsilinResult {
    psilin_kernel(
        Nodes { x: config.x(), y: config.y(), apanel: config.apanel() },
        &ElementPanels::all_of(config),
        i,
        xi,
        yi,
        compute_sources,
    )
}

fn psilin_kernel(
    nodes: Nodes<'_>,
    elements: &[ElementPanels],
    i: usize,
    xi: f64,
    yi: f64,
    compute_sources: bool,
) -> PsilinResult {
    let n = nodes.len();

    // Initialize influence coefficient arrays
    let mut dzdg = vec![0.0; n];
    let mut dzdm = vec![0.0; n];

    // Initialize accumulated values (TODO: compute these in full PSILIN)
    let psi = 0.0;
    let qtan1 = 0.0;
    let qtan2 = 0.0;

    let sgn_surface = if i < n { 1.0 } else { 0.0 };

    // Loop over every element's surface panels. An element's last node starts
    // that element's TE panel, handled separately after this loop.
    // XFOIL line 245: IF(JO.EQ.N) GO TO 11 - skips regular vortex calculation for TE
    for (element, jo) in surface_panels(elements) {
        let jp = element.next_node(jo); // closes within this element

        // Panel endpoints
        let x_jo = nodes.x[jo];
        let y_jo = nodes.y[jo];
        let x_jp = nodes.x[jp];
        let y_jp = nodes.y[jp];

        // Panel vector and length
        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        // Skip zero-length panels
        if ds_sq < 1e-24 {
            continue;
        }

        let dso = ds_sq.sqrt();
        let dsio = 1.0 / dso;

        // Unit tangent vector along panel
        let sx = dx * dsio;
        let sy = dy * dsio;

        // Vector from panel endpoints to field point
        let rx1 = xi - x_jo;
        let ry1 = yi - y_jo;
        let rx2 = xi - x_jp;
        let ry2 = yi - y_jp;

        // Transform to panel-local coordinates
        // x1, x2: along panel direction
        // yy: perpendicular to panel
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        // Squared distances to endpoints
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // Logarithm and arctangent terms
        // XFOIL uses SGN=1 on the airfoil, and sign(yy) in the wake.
        let sgn = if sgn_surface != 0.0 {
            1.0
        } else if yy >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let pi_offset = (0.5 - 0.5 * sgn) * PI; // 0 when sgn=1, PI when sgn=-1

        let (g1, t1) = if i != jo && rs1 > 1e-20 {
            (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if i != jp && rs2 > 1e-20 {
            (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };

        if compute_sources {
            let apan = nodes.apanel[jo];

            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = if rs0 > 0.0 { rs0.ln() } else { 0.0 };
            let t0 = (sgn * x0).atan2(sgn * yy) + pi_offset;

            let dxinv = safe_inv(x1 - x0);
            let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif = ((x1 + x0) * psum
                + rs1 * (t1 - apan)
                - rs0 * (t0 - apan)
                + (x0 - x1) * yy)
                * dxinv;

            // Held inside this element at its first and last node
            let jm = element.stencil_back(jo);
            let jq = element.stencil_forward(jp);

            let dsm =
                ((nodes.x[jp] - nodes.x[jm]).powi(2) + (nodes.y[jp] - nodes.y[jm]).powi(2)).sqrt();
            let dsim = safe_inv(dsm);

            dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
            dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
            dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            let dxinv = safe_inv(x0 - x2);
            let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif = ((x0 + x2) * psum
                + rs0 * (t0 - apan)
                - rs2 * (t2 - apan)
                + (x2 - x0) * yy)
                * dxinv;

            let dsp =
                ((nodes.x[jq] - nodes.x[jo]).powi(2) + (nodes.y[jq] - nodes.y[jo]).powi(2)).sqrt();
            let dsip = safe_inv(dsp);

            dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
            dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);
        }

        // Compute PSIS and PSID (sum/difference formulation)
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);

        // Accumulate influence coefficients
        // ψ += QOPI * (PSIS*(γ_JO+γ_JP) + PSID*(γ_JP-γ_JO))
        // Rearranging: ψ += QOPI * ((PSIS-PSID)*γ_JO + (PSIS+PSID)*γ_JP)
        dzdg[jo] += QOPI * (contrib.psis - contrib.psid);
        dzdg[jp] += QOPI * (contrib.psis + contrib.psid);
    }

    // TE panel special treatment (XFOIL lines 407-438), once per element.
    // This handles each element's trailing edge "wake" panel, running from that
    // element's last node to its own first node, with the HOPI/SCS/SDS formula.
    for element in elements {
        if element.sharp {
            continue;
        }

        let jo = element.last;
        let jp = element.first;
        let (scs, sds) = (element.scs, element.sds);
        let seps = element.seps;

        let x_jo = nodes.x[jo];
        let y_jo = nodes.y[jo];
        let x_jp = nodes.x[jp];
        let y_jp = nodes.y[jp];

        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        if ds_sq > seps * seps {
            let dso = ds_sq.sqrt();
            let dsio = 1.0 / dso;

            let sx = dx * dsio;
            let sy = dy * dsio;

            let rx1 = xi - x_jo;
            let ry1 = yi - y_jo;
            let rx2 = xi - x_jp;
            let ry2 = yi - y_jp;

            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;

            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;

            // TE panel angle (APAN in XFOIL)
            // XFOIL xpanel.f line 184: APAN = APANEL(JO)
            let apan = nodes.apanel[jo];

            // SGN reflection for TE panel
            let sgn = if sgn_surface != 0.0 {
                1.0
            } else if yy >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let pi_offset = (0.5 - 0.5 * sgn) * PI;

            let (g1, t1_raw) = if i != jo && rs1 > 1e-20 {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
            } else {
                (0.0, 0.0)
            };

            let (g2, t2_raw) = if i != jp && rs2 > 1e-20 {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
            } else {
                (0.0, 0.0)
            };

            // XFOIL formulas (lines 408-409):
            // PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
            // PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
            let psig = 0.5 * yy * (g1 - g2) + x2 * (t2_raw - apan) - x1 * (t1_raw - apan);
            let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1_raw - t2_raw);

            // TE panel dPsi/dGam (XFOIL lines 434-438):
            // DZDG(JO) += -HOPI*PSIG*SCS*0.5 + HOPI*PGAM*SDS*0.5
            // DZDG(JP) += +HOPI*PSIG*SCS*0.5 - HOPI*PGAM*SDS*0.5
            dzdg[jo] += HOPI * (-psig * scs + pgam * sds) * 0.5;
            dzdg[jp] += HOPI * (psig * scs - pgam * sds) * 0.5;
        }
    }

    PsilinResult {
        psi,
        qtan1,
        qtan2,
        dzdg,
        dzdm,
        dqdm: vec![0.0; n], // Empty for legacy functions
        dqdg: vec![0.0; n], // Empty for legacy functions (needs NXI/NYI)
    }
}

fn safe_inv(value: f64) -> f64 {
    if value.abs() > 1e-20 {
        1.0 / value
    } else {
        0.0
    }
}

/// Build the source influence matrix BIJ (dPsi/dSig) for airfoil nodes.
///
/// One body: the row layout is `n` surface rows plus a single Kutta row, and the
/// sharp-TE substitution replaces row `n - 1`. A configuration needs one Kutta
/// row per element, which is a system-assembly change rather than an influence
/// kernel one, so this stays single-element for now. The kernels it calls are
/// already element-aware.
pub fn build_source_influence_matrix(geom: &AirfoilGeometry) -> DMatrix<f64> {
    let n = geom.n;
    let mut bij = DMatrix::zeros(n + 1, n);

    for i in 0..n {
        let xi = geom.x[i];
        let yi = geom.y[i];
        let result = psilin_with_sources(geom, i, xi, yi);

        for j in 0..n {
            bij[(i, j)] = -result.dzdm[j];
        }
    }

    // XFOIL replaces the sharp-TE row with a tangential-velocity equation
    // evaluated at a bisector control point just inside the corner.
    if let Some((xbis, ybis, nxbis, nybis)) = geom.sharp_te_bisector_control() {
        let result = psilin_with_dqdm(geom, n, xbis, ybis, nxbis, nybis);
        for j in 0..n {
            bij[(n - 1, j)] = -result.dqdm[j];
        }
    }

    // Kutta row has no direct source influence
    for j in 0..n {
        bij[(n, j)] = 0.0;
    }

    bij
}

/// Compute PSIS and PSID for a single panel contribution.
///
/// This implements the core formulas from XFOIL (lines 194-220).
///
/// # Arguments
///
/// * `x1`, `x2` - Local x-coordinates of field point relative to panel endpoints
/// * `yy` - Local y-coordinate (perpendicular distance to panel line)
/// * `rs1`, `rs2` - Squared distances to panel endpoints
/// * `g1`, `g2` - Log terms: ln(r1²), ln(r2²)
/// * `t1`, `t2` - Angle terms: atan2(x1, yy), atan2(x2, yy)
pub fn compute_psis_psid(
    x1: f64, x2: f64, yy: f64,
    rs1: f64, rs2: f64,
    g1: f64, g2: f64, t1: f64, t2: f64,
) -> PanelContribution {
    // PSIS: sum term (XFOIL line ~215)
    // PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
    let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
    
    // PSID: difference term (XFOIL line ~216)
    // PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2 - RS1*G1 + X1*X1 - X2*X2)) / (X1-X2)
    let dxinv = if (x1 - x2).abs() > 1e-20 {
        1.0 / (x1 - x2)
    } else {
        0.0
    };
    
    let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
    
    PanelContribution {
        psis,
        psid,
        g1, g2, t1, t2, x1, x2, yy,
    }
}

/// Compute influence coefficients for a single panel.
///
/// This is a helper for testing that returns intermediate values.
///
/// The panel starting at `jo` closes onto the next node in `jo`'s own element,
/// which for a single body is `(jo + 1) % n`.
pub fn psilin_single_panel(
    geom: &AirfoilGeometry,
    i: usize,
    jo: usize,
) -> PanelContribution {
    let n = geom.n;
    let jp = ElementPanels::single(geom).next_node(jo);

    // Field point
    let xi = geom.x[i];
    let yi = geom.y[i];

    // Panel endpoints
    let x_jo = geom.x[jo];
    let y_jo = geom.y[jo];
    let x_jp = geom.x[jp];
    let y_jp = geom.y[jp];

    // Panel vector
    let dx = x_jp - x_jo;
    let dy = y_jp - y_jo;
    let ds_sq = dx * dx + dy * dy;
    
    if ds_sq < 1e-24 {
        return PanelContribution::default();
    }
    
    let dso = ds_sq.sqrt();
    let dsio = 1.0 / dso;
    
    // Unit tangent
    let sx = dx * dsio;
    let sy = dy * dsio;
    
    // Vectors to field point
    let rx1 = xi - x_jo;
    let ry1 = yi - y_jo;
    let rx2 = xi - x_jp;
    let ry2 = yi - y_jp;
    
    // Local coordinates
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;
    
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;
    
    // SGN reflection - match main psilin function
    // XFOIL: IF(IO.GE.1 .AND. IO.LE.N) THEN SGN = 1.0 (on surface)
    let sgn = if i < n {
        1.0  // On airfoil surface - no reflection needed
    } else if yy >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let pi_offset = (0.5 - 0.5 * sgn) * PI;
    
    let (g1, t1) = if i != jo && rs1 > 1e-20 {
        (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
    } else {
        (0.0, 0.0)
    };
    
    let (g2, t2) = if i != jp && rs2 > 1e-20 {
        (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
    } else {
        (0.0, 0.0)
    };
    
    compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2)
}

/// Debug structure containing all intermediate values for a single panel influence.
#[derive(Debug, Clone)]
pub struct PsilinDebugPanel {
    /// Panel index (jo)
    pub jo: usize,
    /// Closing node of the panel: the next node in jo's own element
    pub jp: usize,
    /// Panel start point
    pub x_jo: f64,
    pub y_jo: f64,
    /// Panel end point
    pub x_jp: f64,
    pub y_jp: f64,
    /// Panel length
    pub dso: f64,
    /// Unit tangent
    pub sx: f64,
    pub sy: f64,
    /// All local coordinate values
    pub contrib: PanelContribution,
    /// DZDG contribution to jo
    pub dzdg_jo: f64,
    /// DZDG contribution to jp
    pub dzdg_jp: f64,
}

/// Complete debug output for a psilin call.
#[derive(Debug, Clone)]
pub struct PsilinDebug {
    /// Field point index
    pub i: usize,
    /// Field point coordinates
    pub xi: f64,
    pub yi: f64,
    /// Per-panel debug info
    pub panels: Vec<PsilinDebugPanel>,
    /// Final DZDG array
    pub dzdg: Vec<f64>,
}

/// Compute influence coefficients with full debug output.
///
/// This is the same as `psilin` but returns all intermediate values
/// for comparison with XFOIL.
pub fn psilin_debug(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinDebug {
    psilin_debug_kernel(
        Nodes { x: &geom.x, y: &geom.y, apanel: &geom.apanel },
        &[ElementPanels::single(geom)],
        i,
        xi,
        yi,
    )
}

/// [`psilin_debug`] over a whole configuration.
///
/// The reported `jo`/`jp` pairs are global node indices, and each pair lies
/// inside one element.
pub fn psilin_config_debug(
    config: &ConfigGeometry,
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinDebug {
    psilin_debug_kernel(
        Nodes { x: config.x(), y: config.y(), apanel: config.apanel() },
        &ElementPanels::all_of(config),
        i,
        xi,
        yi,
    )
}

fn psilin_debug_kernel(
    nodes: Nodes<'_>,
    elements: &[ElementPanels],
    i: usize,
    xi: f64,
    yi: f64,
) -> PsilinDebug {
    let n = nodes.len();

    // Initialize
    let mut dzdg = vec![0.0; n];
    let mut panels = Vec::new();

    // Loop over every element's surface panels.
    // CRITICAL: an element's TE panel is excluded from the vortex calculation
    // XFOIL line 245: IF(JO.EQ.N) GO TO 11
    for (element, jo) in surface_panels(elements) {
        let jp = element.next_node(jo); // closes within this element

        // Panel endpoints
        let x_jo = nodes.x[jo];
        let y_jo = nodes.y[jo];
        let x_jp = nodes.x[jp];
        let y_jp = nodes.y[jp];

        // Panel vector and length
        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        // Skip zero-length panels
        if ds_sq < 1e-24 {
            continue;
        }

        let dso = ds_sq.sqrt();
        let dsio = 1.0 / dso;
        
        // Unit tangent vector
        let sx = dx * dsio;
        let sy = dy * dsio;
        
        // Vectors from panel endpoints to field point
        let rx1 = xi - x_jo;
        let ry1 = yi - y_jo;
        let rx2 = xi - x_jp;
        let ry2 = yi - y_jp;
        
        // Local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;
        
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;
        
        // Log and angle terms with singularity handling
        // SGN reflection - match main psilin function
        // XFOIL: IF(IO.GE.1 .AND. IO.LE.N) THEN SGN = 1.0 (on surface)
        let sgn = if i < n {
            1.0  // On airfoil surface - no reflection needed
        } else if yy >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let pi_offset = (0.5 - 0.5 * sgn) * PI;
        
        let (g1, t1) = if i != jo && rs1 > 1e-20 {
            (rs1.ln(), (sgn * x1).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };
        
        let (g2, t2) = if i != jp && rs2 > 1e-20 {
            (rs2.ln(), (sgn * x2).atan2(sgn * yy) + pi_offset)
        } else {
            (0.0, 0.0)
        };
        
        // Compute PSIS and PSID
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);
        
        // Accumulate influence coefficients
        let dzdg_jo_contrib = QOPI * (contrib.psis - contrib.psid);
        let dzdg_jp_contrib = QOPI * (contrib.psis + contrib.psid);
        
        dzdg[jo] += dzdg_jo_contrib;
        dzdg[jp] += dzdg_jp_contrib;
        
        panels.push(PsilinDebugPanel {
            jo,
            jp,
            x_jo,
            y_jo,
            x_jp,
            y_jp,
            dso,
            sx,
            sy,
            contrib,
            dzdg_jo: dzdg_jo_contrib,
            dzdg_jp: dzdg_jp_contrib,
        });
    }
    
    PsilinDebug {
        i,
        xi,
        yi,
        panels,
        dzdg,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_geometry() -> AirfoilGeometry {
        // Simple lens shape for testing
        let n = 20;
        let mut pts = Vec::with_capacity(n);
        for i in 0..n/2 {
            let t = i as f64 / (n/2 - 1) as f64;
            let x = 1.0 - t;
            let y = 0.1 * (1.0 - (2.0 * t - 1.0).powi(2));
            pts.push((x, y));
        }
        for i in 1..n/2 {
            let t = i as f64 / (n/2 - 1) as f64;
            let x = t;
            let y = -0.1 * (1.0 - (2.0 * t - 1.0).powi(2));
            pts.push((x, y));
        }
        AirfoilGeometry::from_points(&pts).unwrap()
    }

    // --- multi-element fixtures ------------------------------------------

    /// A NACA 0012 contour with a sharp trailing edge, `n_panels` nodes, scaled
    /// by `scale` and placed with its leading edge at `(dx, dy)`.
    fn naca0012_at(n_panels: usize, scale: f64, dx: f64, dy: f64) -> Vec<(f64, f64)> {
        let n_half = n_panels / 2;
        let thickness = |x: f64| -> f64 {
            0.6 * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4))
        };
        let station = |i: usize| -> f64 {
            let beta = PI * (i as f64) / (n_half as f64);
            0.5 * (1.0 - beta.cos())
        };
        let place = |x: f64, y: f64| (dx + scale * x, dy + scale * y);

        let mut points = Vec::with_capacity(2 * n_half);
        for i in (0..=n_half).rev() {
            let x = station(i);
            points.push(place(x, thickness(x)));
        }
        for i in 1..=n_half {
            let x = station(i);
            points.push(place(x, -thickness(x)));
        }
        points
    }

    /// The same contour opened out to a blunt trailing edge of `gap` chords.
    fn naca0012_blunt_at(
        n_panels: usize,
        gap: f64,
        scale: f64,
        dx: f64,
        dy: f64,
    ) -> Vec<(f64, f64)> {
        let mut points = naca0012_at(n_panels, scale, dx, dy);
        let last = points.len() - 1;
        points[0].1 += 0.5 * gap * scale;
        points[last].1 -= 0.5 * gap * scale;
        points
    }

    /// Every element of the real McDonnell Douglas 30P-30N slat/main/flap
    /// fixture, in configuration coordinates. Blank and comment lines separate
    /// the blocks.
    fn mda_elements() -> Vec<Vec<(f64, f64)>> {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata/mda_30p_30n_trimmed.dat");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));

        let mut blocks: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut current: Vec<(f64, f64)> = Vec::new();
        for line in text.lines() {
            let mut parts = line.split_whitespace();
            match (
                parts.next().and_then(|s| s.parse::<f64>().ok()),
                parts.next().and_then(|s| s.parse::<f64>().ok()),
            ) {
                (Some(x), Some(y)) => current.push((x, y)),
                _ if !current.is_empty() => blocks.push(std::mem::take(&mut current)),
                _ => {}
            }
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        assert_eq!(blocks.len(), 3, "expected slat, main and flap");
        blocks
    }

    /// Compare two f64 arrays bit for bit.
    fn assert_bits_eq(actual: &[f64], expected: &[f64], what: &str) {
        assert_eq!(actual.len(), expected.len(), "{what}: length");
        for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{what}[{i}]: {a:?} ({:#018x}) vs {b:?} ({:#018x})",
                a.to_bits(),
                b.to_bits()
            );
        }
    }

    // --- connectivity ----------------------------------------------------

    /// The reduction the whole refactor rests on: with one element `next_node`
    /// is exactly the modulo-n closure it replaced, so no single-element
    /// influence coefficient moves.
    #[test]
    fn single_element_next_node_reduces_to_modulo_n() {
        for points in [
            naca0012_at(80, 1.0, 0.0, 0.0),
            naca0012_blunt_at(80, 0.004, 1.0, 0.0, 0.0),
            naca0012_at(160, 0.3, 1.2, -0.1),
        ] {
            let geom = AirfoilGeometry::from_points(&points).unwrap();
            let element = ElementPanels::single(&geom);

            assert_eq!(element.first, 0);
            assert_eq!(element.last, geom.n - 1);
            for jo in 0..geom.n {
                assert_eq!(element.next_node(jo), (jo + 1) % geom.n, "node {jo}");
            }

            // And the panel loop still walks every panel but the TE panel.
            let walked: Vec<usize> = element.surface_panels().collect();
            assert_eq!(walked, (0..geom.n - 1).collect::<Vec<_>>());
        }
    }

    /// `ElementPanels` closes contours by the same rule as
    /// `Layout::next_node`, which `ConfigGeometry::next_node` forwards to.
    #[test]
    fn element_panels_match_the_layout_closure_rule() {
        let fixtures: Vec<Vec<Vec<(f64, f64)>>> = vec![
            mda_elements(),
            vec![
                naca0012_at(80, 1.0, 0.0, 0.0),
                naca0012_blunt_at(60, 0.004, 0.3, 1.05, -0.08),
            ],
            vec![naca0012_at(120, 1.0, 0.0, 0.0)],
        ];

        for elements in fixtures {
            let config = ConfigGeometry::from_element_points(&elements).unwrap();
            let panels = ElementPanels::all_of(&config);
            assert_eq!(panels.len(), config.n_elements());

            for global in 0..config.total_nodes() {
                let element = &panels[config.element_of(global)];
                assert!(
                    element.first <= global && global <= element.last,
                    "node {global} is outside its own element"
                );
                assert_eq!(
                    element.next_node(global),
                    config.next_node(global),
                    "node {global} of a {}-element configuration",
                    config.n_elements()
                );
            }
        }
    }

    /// The failure mode this refactor exists to remove. With several elements
    /// concatenated into one node array, `(jo + 1) % n` closes element 0's last
    /// panel onto element 1's *first* node, laying a panel across the physical
    /// gap between them. Nothing crashes when that happens — the Cp
    /// distribution is simply wrong — so it has to be asserted.
    #[test]
    fn no_panel_endpoint_pair_straddles_an_element_boundary() {
        let elements = vec![
            naca0012_at(80, 1.0, 0.0, 0.0),
            naca0012_blunt_at(60, 0.004, 0.3, 1.05, -0.08),
        ];
        let config = ConfigGeometry::from_element_points(&elements).unwrap();
        assert_eq!(config.n_elements(), 2);

        let total = config.total_nodes();
        let first = *config.element(0);
        let second = *config.element(1);
        let boundary = first.end() - 1; // element 0's last node

        // Every surface panel the kernel walks stays inside one element.
        let debug = psilin_config_debug(&config, total, 2.0, 0.4);
        for panel in &debug.panels {
            assert_eq!(
                config.element_of(panel.jo),
                config.element_of(panel.jp),
                "panel {} -> {} crosses an element boundary",
                panel.jo,
                panel.jp
            );
        }
        // One TE panel is held back per element, not one for the whole array.
        assert_eq!(debug.panels.len(), total - config.n_elements());

        // So does each element's TE panel, which is where modulo-n closure
        // reached into the neighbouring element.
        let panels = ElementPanels::all_of(&config);
        assert_eq!(panels[0].next_node(boundary), first.start);
        assert_eq!(panels[1].next_node(total - 1), second.start);

        // Stated against the arithmetic it replaces: at a boundary the two
        // answers differ, and the modulo-n one names the wrong element.
        let modulo_n = |global: usize| (global + 1) % total;
        assert_eq!(modulo_n(boundary), second.start);
        assert_ne!(panels[0].next_node(boundary), modulo_n(boundary));
        assert_ne!(panels[1].next_node(total - 1), modulo_n(total - 1));

        // The source-gradient stencil is held inside the element too.
        assert_eq!(panels[0].stencil_forward(boundary), boundary);
        assert_eq!(panels[1].stencil_back(second.start), second.start);
    }

    /// Each element's trailing edge is treated with its own coefficients,
    /// tolerance and sharpness verdict rather than the first element's.
    #[test]
    fn each_element_carries_its_own_trailing_edge_data() {
        let elements = vec![
            naca0012_blunt_at(80, 0.01, 1.0, 0.0, 0.0),
            naca0012_at(60, 0.3, 1.05, -0.08),
        ];
        let config = ConfigGeometry::from_element_points(&elements).unwrap();
        let panels = ElementPanels::all_of(&config);

        assert!(!panels[0].sharp, "the opened trailing edge is blunt");
        assert!(panels[1].sharp, "the closed trailing edge is sharp");

        assert_eq!((panels[0].scs, panels[0].sds), config.te_coefficients(0));
        assert_eq!((panels[1].scs, panels[1].sds), (1.0, 0.0));

        // SEPS comes from each element's own arc length, so the smaller element
        // gets the smaller tolerance.
        assert_eq!(panels[0].seps, config.element_arc_length(0) * 1e-5);
        assert_eq!(panels[1].seps, config.element_arc_length(1) * 1e-5);
        assert!(panels[1].seps < panels[0].seps);
    }

    // --- influence coefficients ------------------------------------------

    /// No element influences another through the panel loop: a configuration's
    /// coefficients are its elements' own, concatenated. Under modulo-n closure
    /// the two spurious gap-spanning panels show up here as differences at the
    /// nodes either side of each boundary.
    #[test]
    fn configuration_influence_is_the_sum_of_its_independent_elements() {
        let fixtures: Vec<Vec<Vec<(f64, f64)>>> = vec![
            vec![
                naca0012_blunt_at(80, 0.004, 1.0, 0.0, 0.0),
                naca0012_at(60, 0.3, 1.05, -0.08),
            ],
            mda_elements(),
        ];

        for elements in fixtures {
            let config = ConfigGeometry::from_element_points(&elements).unwrap();
            let total = config.total_nodes();

            // A field point off every surface, indexed past the last node, so
            // each element meets the same conditions alone as it does in the
            // configuration.
            let (xi, yi) = (0.42, 0.63);
            let joint = psilin_config_with_sources(&config, total, xi, yi);
            let joint_dqdm = psilin_config_with_dqdm(&config, total, xi, yi, 0.6, 0.8);

            for (k, points) in elements.iter().enumerate() {
                let alone = ConfigGeometry::from_points(points).unwrap();
                let n = alone.total_nodes();
                let range = config.element_range(k);
                let what = format!("element {k} of {}", config.n_elements());

                let solo = psilin_config_with_sources(&alone, n, xi, yi);
                assert_bits_eq(&joint.dzdg[range.clone()], &solo.dzdg, &format!("{what} dzdg"));
                assert_bits_eq(&joint.dzdm[range.clone()], &solo.dzdm, &format!("{what} dzdm"));

                let solo = psilin_config_with_dqdm(&alone, n, xi, yi, 0.6, 0.8);
                assert_bits_eq(
                    &joint_dqdm.dzdg[range.clone()],
                    &solo.dzdg,
                    &format!("{what} dqdm dzdg"),
                );
                assert_bits_eq(
                    &joint_dqdm.dzdm[range.clone()],
                    &solo.dzdm,
                    &format!("{what} dqdm dzdm"),
                );
                assert_bits_eq(
                    &joint_dqdm.dqdm[range.clone()],
                    &solo.dqdm,
                    &format!("{what} dqdm"),
                );
                assert_bits_eq(&joint_dqdm.dqdg[range], &solo.dqdg, &format!("{what} dqdg"));
            }
        }
    }

    /// The single-element path through the element-aware kernel is the
    /// single-body path: same nodes, same closure, same bits. This is the
    /// in-repo statement of the parity requirement — a single element must
    /// reproduce the existing numbers exactly, not closely.
    #[test]
    fn single_element_configuration_matches_the_single_body_kernel_bit_for_bit() {
        for points in [
            naca0012_blunt_at(80, 0.004, 1.0, 0.0, 0.0),
            naca0012_at(80, 1.0, 0.0, 0.0),
        ] {
            let geom = AirfoilGeometry::from_points(&points).unwrap();
            let config = ConfigGeometry::from_points(&points).unwrap();
            let n = geom.n;

            let probes: Vec<(usize, f64, f64)> = (0..n)
                .map(|i| (i, geom.x[i], geom.y[i]))
                .chain([(n, 1.4, 0.02), (n + 1, 0.5, 0.25)])
                .collect();

            for (i, xi, yi) in probes {
                let expected = psilin(&geom, i, xi, yi);
                let actual = psilin_config(&config, i, xi, yi);
                assert_bits_eq(&actual.dzdg, &expected.dzdg, "psilin dzdg");

                let expected = psilin_with_sources(&geom, i, xi, yi);
                let actual = psilin_config_with_sources(&config, i, xi, yi);
                assert_bits_eq(&actual.dzdg, &expected.dzdg, "sources dzdg");
                assert_bits_eq(&actual.dzdm, &expected.dzdm, "sources dzdm");

                let expected = psilin_with_dqdm(&geom, i, xi, yi, 0.6, 0.8);
                let actual = psilin_config_with_dqdm(&config, i, xi, yi, 0.6, 0.8);
                assert_bits_eq(&actual.dzdg, &expected.dzdg, "dqdm dzdg");
                assert_bits_eq(&actual.dzdm, &expected.dzdm, "dqdm dzdm");
                assert_bits_eq(&actual.dqdm, &expected.dqdm, "dqdm");
                assert_bits_eq(&actual.dqdg, &expected.dqdg, "dqdg");

                let expected = psilin_debug(&geom, i, xi, yi);
                let actual = psilin_config_debug(&config, i, xi, yi);
                assert_bits_eq(&actual.dzdg, &expected.dzdg, "debug dzdg");
                assert_eq!(actual.panels.len(), expected.panels.len());
                for (a, b) in actual.panels.iter().zip(&expected.panels) {
                    assert_eq!((a.jo, a.jp), (b.jo, b.jp));
                    assert_eq!(a.dzdg_jo.to_bits(), b.dzdg_jo.to_bits());
                    assert_eq!(a.dzdg_jp.to_bits(), b.dzdg_jp.to_bits());
                }
            }
        }
    }

    #[test]
    fn test_psilin_produces_finite_values() {
        let geom = make_simple_geometry();
        
        // Test at a few interior points
        for i in 1..geom.n - 1 {
            let result = psilin(&geom, i, geom.x[i], geom.y[i]);
            
            assert!(result.psi.is_finite(), "psi should be finite at node {}", i);
            for (j, &dz) in result.dzdg.iter().enumerate() {
                assert!(dz.is_finite(), "dzdg[{}] should be finite at node {}", j, i);
            }
        }
    }

    #[test]
    fn test_psis_psid_symmetry() {
        // For a point directly above the panel midpoint, certain symmetries hold
        let x1: f64 = 0.5;
        let x2: f64 = -0.5;
        let yy: f64 = 1.0;
        let rs1: f64 = x1 * x1 + yy * yy;
        let rs2: f64 = x2 * x2 + yy * yy;
        let g1 = rs1.ln();
        let g2 = rs2.ln();
        let t1 = x1.atan2(yy);
        let t2 = x2.atan2(yy);
        
        let contrib = compute_psis_psid(x1, x2, yy, rs1, rs2, g1, g2, t1, t2);
        
        // rs1 == rs2 and g1 == g2 due to symmetry
        assert!((rs1 - rs2).abs() < 1e-10);
        assert!((g1 - g2).abs() < 1e-10);
        
        // PSIS and PSID should be finite
        assert!(contrib.psis.is_finite());
        assert!(contrib.psid.is_finite());
    }
}
