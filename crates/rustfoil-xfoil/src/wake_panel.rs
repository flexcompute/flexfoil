use rustfoil_bl::{add_event, is_debug_active, DebugEvent};
use rustfoil_inviscid::geometry::AirfoilGeometry;
use rustfoil_inviscid::influence::psilin_with_dqdm;
use rustfoil_inviscid::system::{compute_wake_count, setexp};
use rustfoil_inviscid::FactorizedSystem;

use crate::{state::XfoilState, Result};

/// Faithful port of XFOIL's XYWAKE (xpanel.f lines 1326-1409).
///
/// Sets wake coordinate array by following streamlines from the trailing edge.
/// At each wake point, PSILIN is called to get the streamfunction gradient ∇ψ,
/// and the wake normal is set to -∇ψ/|∇ψ|. The next point is placed DS
/// downstream along the local flow direction (perpendicular to the normal).
///
/// This replaces the previous straight-line wake with a streamline-following wake
/// that matches XFOIL's behavior exactly.
pub fn xywake(state: &mut XfoilState, factorized: &FactorizedSystem, wake_length_chords: f64) {
    let geom = factorized.geometry();
    let n = geom.n;

    if n < 2 || wake_length_chords <= 0.0 {
        return;
    }

    // Number of wake points (Fortran line 1336: NW = N/12 + 10*INT(WAKLEN))
    let nw = compute_wake_count(n, wake_length_chords);
    if nw == 0 {
        return;
    }

    // Exponential spacing (Fortran lines 1343-1344)
    let ds1 = 0.5 * ((geom.s[1] - geom.s[0]) + (geom.s[n - 1] - geom.s[n - 2]));
    let snew = setexp(ds1, wake_length_chords * geom.chord, nw);

    // Trailing edge midpoint (Fortran lines 1353-1354)
    let xte = 0.5 * (geom.x[0] + geom.x[n - 1]);
    let yte = 0.5 * (geom.y[0] + geom.y[n - 1]);

    // Compute GAM at current alpha for PSILIN streamfunction gradient.
    // In Fortran XYWAKE, GAM is already set from SPECAL before XYWAKE is called.
    // Here we compute it from the decomposed gamu_0/gamu_90 to avoid depending
    // on prior specal call order.
    let cosa = state.alpha_rad.cos();
    let sina = state.alpha_rad.sin();
    let gam: Vec<f64> = (0..n)
        .map(|j| cosa * factorized.gamu_0[j] + sina * factorized.gamu_90[j])
        .collect();

    // First wake point: initial normal from TE tangents (Fortran lines 1358-1364)
    let sx = 0.5 * (geom.yp[n - 1] - geom.yp[0]);
    let sy = 0.5 * (geom.xp[0] - geom.xp[n - 1]);
    let smod = (sx * sx + sy * sy).sqrt().max(1e-12);
    let nx_init = sx / smod;
    let ny_init = sy / smod;

    let mut wake_x = Vec::with_capacity(nw);
    let mut wake_y = Vec::with_capacity(nw);

    // Place first point slightly behind TE (Fortran lines 1363-1364)
    wake_x.push(xte - 0.0001 * ny_init);
    wake_y.push(yte + 0.0001 * nx_init);

    // Calculate streamfunction gradient at first point (Fortran lines 1368-1369)
    // CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
    // CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
    let (psi_x, psi_y) = compute_psi_gradient(geom, n, wake_x[0], wake_y[0], &gam, cosa, sina);
    let grad_mag = (psi_x * psi_x + psi_y * psi_y).sqrt().max(1e-12);

    // Set unit vector normal to wake at next point (Fortran lines 1372-1373)
    let mut nx_next = -psi_x / grad_mag;
    let mut ny_next = -psi_y / grad_mag;

    // Set rest of wake points (Fortran lines 1379-1399)
    for i in 1..nw {
        // Use the normal computed from previous point's PSILIN call
        let nx_curr = nx_next;
        let ny_curr = ny_next;

        // Arc length step (Fortran line 1380: DS = SNEW(I) - SNEW(I-1))
        let ds = snew[i] - snew[i - 1];

        // Set new point DS downstream of last point (Fortran lines 1383-1385)
        // X(I) = X(I-1) - DS*NY(I),  Y(I) = Y(I-1) + DS*NX(I)
        let xi = wake_x[i - 1] - ds * ny_curr;
        let yi = wake_y[i - 1] + ds * nx_curr;
        wake_x.push(xi);
        wake_y.push(yi);

        // Skip PSILIN for last point (Fortran line 1387: IF(I.EQ.N+NW) GO TO 10)
        if i < nw - 1 {
            // Calculate normal vector for next point (Fortran lines 1390-1394)
            let (psi_x, psi_y) =
                compute_psi_gradient(geom, n + i, xi, yi, &gam, cosa, sina);
            let grad_mag = (psi_x * psi_x + psi_y * psi_y).sqrt().max(1e-12);
            nx_next = -psi_x / grad_mag;
            ny_next = -psi_y / grad_mag;
        }
    }

    // Set state wake arrays
    state.wake_x = wake_x;
    state.wake_y = wake_y;
    state.wake_s = compute_wake_arc_lengths(xte, yte, &state.wake_x, &state.wake_y);
}

/// Faithful port of XFOIL's QWCALC (xpanel.f lines 1142-1161).
///
/// Sets inviscid tangential velocity for alpha = 0, 90 on wake due to
/// freestream and airfoil surface vorticity.
///
/// - First wake point: copies from trailing edge (QINVU(N+1,:) = QINVU(N,:))
/// - Remaining wake points: calls PSILIN to get QTAN1, QTAN2 per-point
///
/// The per-point velocities are then combined at the current alpha:
///   wake_qinv  =  cos(α)*qinvu_0 + sin(α)*qinvu_90
///   wake_qinv_a = -sin(α)*qinvu_0 + cos(α)*qinvu_90
pub fn qwcalc(state: &mut XfoilState, factorized: &FactorizedSystem) {
    let n = state.n_panel_nodes();
    let nw = state.wake_x.len();
    if nw == 0 || n < 2 {
        return;
    }

    // Compute wake normals from geometry (finite-difference approximation).
    // In the faithful XYWAKE port, these will come from PSILIN streamfunction gradients.
    let (wake_nx, wake_ny) = compute_wake_normals(&state.wake_x, &state.wake_y);

    let geom = factorized.geometry();

    // Allocate per-wake-point QINVU arrays (alpha=0 and alpha=90 decomposition).
    let mut wake_qinvu_0 = vec![0.0; nw];
    let mut wake_qinvu_90 = vec![0.0; nw];

    // First wake point: same as trailing edge (Fortran QWCALC lines 1150-1151).
    // QINVU(N+1,1) = QINVU(N,1), QINVU(N+1,2) = QINVU(N,2)
    // In Rust: gamu_0[n-1] and gamu_90[n-1] are the last airfoil node values.
    wake_qinvu_0[0] = factorized.gamu_0[n - 1];
    wake_qinvu_90[0] = factorized.gamu_90[n - 1];

    // Rest of wake: call PSILIN at each point (Fortran QWCALC lines 1154-1158).
    // PSILIN returns QTAN1/QTAN2 = tangential velocity at alpha=0/90.
    // QTAN1 = sum_j DQDG(j)*GAMU(j,1) + QINF*NYI
    // QTAN2 = sum_j DQDG(j)*GAMU(j,2) - QINF*NXI
    for i in 1..nw {
        let xi = state.wake_x[i];
        let yi = state.wake_y[i];
        let nxi = wake_nx[i];
        let nyi = wake_ny[i];

        // Field point index is n + i (global index in combined airfoil+wake space).
        let psilin_result = psilin_with_dqdm(geom, n + i, xi, yi, nxi, nyi);

        // Compute QTAN1/QTAN2 from DQDG and gamu arrays.
        // Freestream contribution: QINF = 1.0
        let mut qtan1: f64 = nyi;   // QINF * NYI
        let mut qtan2: f64 = -nxi;  // -QINF * NXI
        for j in 0..n {
            qtan1 += psilin_result.dqdg[j] * factorized.gamu_0[j];
            qtan2 += psilin_result.dqdg[j] * factorized.gamu_90[j];
        }

        wake_qinvu_0[i] = qtan1;
        wake_qinvu_90[i] = qtan2;
    }

    // Combine into wake_qinv/wake_qinv_a at current alpha (XFOIL's QISET for wake).
    let cosa = state.alpha_rad.cos();
    let sina = state.alpha_rad.sin();
    state.wake_qinv = wake_qinvu_0
        .iter()
        .zip(wake_qinvu_90.iter())
        .map(|(&q0, &q90)| cosa * q0 + sina * q90)
        .collect();
    state.wake_qinv_a = wake_qinvu_0
        .iter()
        .zip(wake_qinvu_90.iter())
        .map(|(&q0, &q90)| -sina * q0 + cosa * q90)
        .collect();
}

pub fn qdcalc(state: &mut XfoilState, factorized: &FactorizedSystem) -> Result<()> {
    state.dij = factorized.build_dij_with_wake(&state.wake_x, &state.wake_y)?;
    if is_debug_active() {
        let nsys = state.dij.nrows();
        let diag_sample: Vec<f64> = (0..nsys.min(20)).map(|i| state.dij[(i, i)]).collect();
        let row1_sample: Vec<f64> = (0..nsys.min(20)).map(|j| state.dij[(0, j)]).collect();
        add_event(DebugEvent::qdcalc(
            state.n_panel_nodes(),
            state.wake_x.len(),
            nsys,
            diag_sample,
            row1_sample,
        ));
    }
    Ok(())
}

fn compute_wake_arc_lengths(te_x: f64, te_y: f64, wake_x: &[f64], wake_y: &[f64]) -> Vec<f64> {
    let mut s = Vec::with_capacity(wake_x.len());
    let mut prev_x = te_x;
    let mut prev_y = te_y;
    let mut accum = 0.0;
    for (&x, &y) in wake_x.iter().zip(wake_y.iter()) {
        let dx = x - prev_x;
        let dy = y - prev_y;
        accum += (dx * dx + dy * dy).sqrt();
        s.push(accum);
        prev_x = x;
        prev_y = y;
    }
    s
}

/// Compute the streamfunction gradient (∂ψ/∂x, ∂ψ/∂y) at a field point.
///
/// This replicates XFOIL's two PSILIN calls in XYWAKE (lines 1368-1369, 1390-1391):
/// - CALL PSILIN(I,XI,YI,1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.) → ∂ψ/∂x
/// - CALL PSILIN(I,XI,YI,0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.) → ∂ψ/∂y
///
/// The gradient is computed as:
///   PSI_NI = sum_j DQDG(j) * GAM(j) + QINF*(cos(α)*NYI - sin(α)*NXI)
///
/// where DQDG comes from psilin_with_dqdm and the freestream term matches
/// Fortran PSILIN line 477.
fn compute_psi_gradient(
    geom: &AirfoilGeometry,
    i: usize,
    xi: f64,
    yi: f64,
    gam: &[f64],
    cosa: f64,
    sina: f64,
) -> (f64, f64) {
    let n = geom.n;

    // PSILIN with NXI=1.0, NYI=0.0 → PSI_X (= ∂ψ/∂x)
    let result_x = psilin_with_dqdm(geom, i, xi, yi, 1.0, 0.0);
    // Freestream: QINF*(cos(α)*NYI - sin(α)*NXI) = QINF*(0 - sin(α)) = -sin(α)
    let mut psi_x: f64 = -sina;
    for j in 0..n {
        psi_x += result_x.dqdg[j] * gam[j];
    }

    // PSILIN with NXI=0.0, NYI=1.0 → PSI_Y (= ∂ψ/∂y)
    let result_y = psilin_with_dqdm(geom, i, xi, yi, 0.0, 1.0);
    // Freestream: QINF*(cos(α)*NYI - sin(α)*NXI) = QINF*(cos(α) - 0) = cos(α)
    let mut psi_y: f64 = cosa;
    for j in 0..n {
        psi_y += result_y.dqdg[j] * gam[j];
    }

    (psi_x, psi_y)
}

/// Compute wake panel normals from wake coordinates.
///
/// Uses finite-difference approximation of the wake tangent direction.
/// For a streamline-following wake from faithful XYWAKE, these normals
/// closely match the PSILIN-derived normals since the wake follows the flow.
fn compute_wake_normals(wake_x: &[f64], wake_y: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let nw = wake_x.len();
    if nw == 0 {
        return (Vec::new(), Vec::new());
    }

    let mut nx = vec![0.0; nw];
    let mut ny = vec![0.0; nw];

    for i in 0..nw {
        let (tx, ty) = if i < nw - 1 {
            // Forward difference for interior and first points
            let dx = wake_x[i + 1] - wake_x[i];
            let dy = wake_y[i + 1] - wake_y[i];
            (dx, dy)
        } else if nw >= 2 {
            // Backward difference for last point
            let dx = wake_x[i] - wake_x[i - 1];
            let dy = wake_y[i] - wake_y[i - 1];
            (dx, dy)
        } else {
            (1.0, 0.0)
        };

        let ds = (tx * tx + ty * ty).sqrt().max(1e-12);
        // Normal is perpendicular to tangent: n = (ty/ds, -tx/ds)
        nx[i] = ty / ds;
        ny[i] = -tx / ds;
    }

    (nx, ny)
}
