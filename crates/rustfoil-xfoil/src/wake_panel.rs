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

    let cosa = state.alpha_rad.cos();
    let sina = state.alpha_rad.sin();
    // XFOIL's XYWAKE uses the current GAM array prepared by SPECAL, including
    // the trailing-edge endpoint handling. Reuse that state here so the wake
    // streamline follows the same circulation field seen elsewhere.
    let gam: Vec<f64> = if state.gam.len() >= n {
        state.gam[..n].to_vec()
    } else {
        let mut gam: Vec<f64> = (0..n)
            .map(|j| cosa * factorized.gamu_0[j] + sina * factorized.gamu_90[j])
            .collect();
        if let Some(last) = gam.last_mut() {
            *last = 0.0;
        }
        gam
    };

    // First wake point: initial normal from TE tangents (Fortran lines 1358-1364)
    let sx = 0.5 * (geom.yp[n - 1] - geom.yp[0]);
    let sy = 0.5 * (geom.xp[0] - geom.xp[n - 1]);
    let smod = (sx * sx + sy * sy).sqrt().max(1e-12);
    let nx_init = sx / smod;
    let ny_init = sy / smod;

    let mut wake_x = Vec::with_capacity(nw);
    let mut wake_y = Vec::with_capacity(nw);
    let mut wake_nx = vec![0.0; nw];
    let mut wake_ny = vec![0.0; nw];
    let mut wake_apanel = vec![0.0; nw];

    // Place first point slightly behind TE (Fortran lines 1363-1364)
    wake_x.push(xte - 0.0001 * ny_init);
    wake_y.push(yte + 0.0001 * nx_init);
    wake_nx[0] = nx_init;
    wake_ny[0] = ny_init;

    // Calculate streamfunction gradient at first point (Fortran lines 1368-1369)
    // CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
    // CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
    let (psi_x, psi_y) = compute_psi_gradient(geom, n, wake_x[0], wake_y[0], &gam, cosa, sina);
    let grad_mag = (psi_x * psi_x + psi_y * psi_y).sqrt().max(1e-12);

    // Set unit vector normal to wake at next point (Fortran lines 1372-1373)
    let mut nx_next = -psi_x / grad_mag;
    let mut ny_next = -psi_y / grad_mag;
    if nw > 1 {
        wake_nx[1] = nx_next;
        wake_ny[1] = ny_next;
    }
    wake_apanel[0] = psi_y.atan2(psi_x);

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
            wake_nx[i + 1] = nx_next;
            wake_ny[i + 1] = ny_next;
            wake_apanel[i] = psi_y.atan2(psi_x);
        }
    }

    // Set state wake arrays
    state.wake_x = wake_x;
    state.wake_y = wake_y;
    state.wake_s = snew;
    state.wake_nx = wake_nx;
    state.wake_ny = wake_ny;
    state.wake_apanel = wake_apanel;
    state.lwake = true;
    state.awake = state.alpha_rad;
    state.lwdij = false;
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

    // Allocate per-wake-point QINVU arrays (alpha=0 and alpha=90 decomposition).
    state.wake_qinvu_0 = vec![0.0; nw];
    state.wake_qinvu_90 = vec![0.0; nw];

    // XFOIL copies the TE inviscid velocity to the first wake point.
    // The remaining wake points are computed from PSILIN.
    state.wake_qinvu_0[0] = state.qinvu_0[n - 1];
    state.wake_qinvu_90[0] = state.qinvu_90[n - 1];

    // Recover the tangential wake velocity from the streamfunction gradient
    // for the remaining wake points.
    for i in 1..nw {
        let xi = state.wake_x[i];
        let yi = state.wake_y[i];
        let nxi = state.wake_nx[i];
        let nyi = state.wake_ny[i];

        let (psi_x_0, psi_y_0) = compute_psi_gradient(
            factorized.geometry(),
            n + i,
            xi,
            yi,
            &factorized.gamu_0,
            1.0,
            0.0,
        );
        let (psi_x_90, psi_y_90) = compute_psi_gradient(
            factorized.geometry(),
            n + i,
            xi,
            yi,
            &factorized.gamu_90,
            0.0,
            1.0,
        );

        let qtan1 = psi_x_0 * nxi + psi_y_0 * nyi;
        let qtan2 = psi_x_90 * nxi + psi_y_90 * nyi;

        state.wake_qinvu_0[i] = qtan1;
        state.wake_qinvu_90[i] = qtan2;
    }

    // Combine into wake_qinv/wake_qinv_a at current alpha (XFOIL's QISET for wake).
    let cosa = state.alpha_rad.cos();
    let sina = state.alpha_rad.sin();
    state.wake_qinv = state
        .wake_qinvu_0
        .iter()
        .zip(state.wake_qinvu_90.iter())
        .map(|(&q0, &q90)| cosa * q0 + sina * q90)
        .collect();
    state.wake_qinv_a = state
        .wake_qinvu_0
        .iter()
        .zip(state.wake_qinvu_90.iter())
        .map(|(&q0, &q90)| -sina * q0 + cosa * q90)
        .collect();
}

pub fn qdcalc(state: &mut XfoilState, factorized: &FactorizedSystem) -> Result<()> {
    state.dij = factorized.build_dij_with_wake_state(
        &state.wake_x,
        &state.wake_y,
        &state.wake_nx,
        &state.wake_ny,
        &state.wake_apanel,
    )?;
    state.ladij = true;
    state.lwdij = true;
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
