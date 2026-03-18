use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use rustfoil_core::{naca, point, Body, CubicSpline, PanelingParams, Point};
use rustfoil_xfoil::oper::{solve_body_oper_point, AlphaSpec};
use rustfoil_xfoil::XfoilOptions;

struct FaithfulResult {
    alpha_deg: f64,
    cl: f64,
    cd: f64,
    cm: f64,
    converged: bool,
    iterations: usize,
    residual: f64,
    x_tr_upper: f64,
    x_tr_lower: f64,
    cd_friction: f64,
    cd_pressure: f64,
    success: bool,
    error: Option<String>,
}

fn points_from_flat(coords: &[f64]) -> Vec<Point> {
    coords.chunks(2).map(|c| point(c[0], c[1])).collect()
}

fn flat_from_points(pts: &[Point]) -> Vec<(f64, f64)> {
    pts.iter().map(|p| (p.x, p.y)).collect()
}

/// Viscous (XFOIL-faithful) analysis at a single operating point.
///
/// Returns a dict with keys: cl, cd, cm, converged, iterations, residual,
/// x_tr_upper, x_tr_lower, cd_friction, cd_pressure, alpha_deg, success, error.
#[pyfunction]
#[pyo3(signature = (coords, alpha_deg, reynolds=1.0e6, mach=0.0, ncrit=9.0, max_iterations=100))]
fn analyze_faithful(
    py: Python<'_>,
    coords: Vec<f64>,
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> PyResult<Py<PyDict>> {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        let d = PyDict::new(py);
        d.set_item("success", false)?;
        d.set_item("error", "Invalid coordinates: need at least 3 points (6 values)")?;
        return Ok(d.into());
    }

    let points = points_from_flat(&coords);
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            let d = PyDict::new(py);
            d.set_item("success", false)?;
            d.set_item("error", format!("Geometry error: {e}"))?;
            return Ok(d.into());
        }
    };

    let options = XfoilOptions {
        reynolds,
        mach,
        ncrit,
        max_iterations,
        ..Default::default()
    };

    let d = PyDict::new(py);
    match solve_body_oper_point(&body, AlphaSpec::AlphaDeg(alpha_deg), &options) {
        Ok(r) => {
            d.set_item("cl", r.cl)?;
            d.set_item("cd", r.cd)?;
            d.set_item("cm", r.cm)?;
            d.set_item("converged", r.converged)?;
            d.set_item("iterations", r.iterations)?;
            d.set_item("residual", r.residual)?;
            d.set_item("x_tr_upper", r.x_tr_upper)?;
            d.set_item("x_tr_lower", r.x_tr_lower)?;
            d.set_item("cd_friction", r.cd_friction)?;
            d.set_item("cd_pressure", r.cd_pressure)?;
            d.set_item("alpha_deg", r.alpha_deg)?;
            d.set_item("success", true)?;
            d.set_item("error", py.None())?;
        }
        Err(e) => {
            d.set_item("success", false)?;
            d.set_item("error", format!("{e}"))?;
        }
    }
    Ok(d.into())
}

/// Inviscid panel-method analysis at a single angle of attack.
///
/// Returns a dict with keys: cl, cm, cp, cp_x, success, error.
#[pyfunction]
fn analyze_inviscid(py: Python<'_>, coords: Vec<f64>, alpha_deg: f64) -> PyResult<Py<PyDict>> {
    use rustfoil_solver::inviscid::{FlowConditions, InviscidSolver};

    if coords.len() < 6 || coords.len() % 2 != 0 {
        let d = PyDict::new(py);
        d.set_item("success", false)?;
        d.set_item("error", "Invalid coordinates")?;
        return Ok(d.into());
    }

    let points = points_from_flat(&coords);
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            let d = PyDict::new(py);
            d.set_item("success", false)?;
            d.set_item("error", format!("Geometry error: {e}"))?;
            return Ok(d.into());
        }
    };

    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha_deg);

    let d = PyDict::new(py);
    match solver.solve(&[body.clone()], &flow) {
        Ok(solution) => {
            let cp_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint().x).collect();
            d.set_item("cl", solution.cl)?;
            d.set_item("cm", solution.cm)?;
            d.set_item("cp", solution.cp)?;
            d.set_item("cp_x", cp_x)?;
            d.set_item("success", true)?;
            d.set_item("error", py.None())?;
        }
        Err(e) => {
            d.set_item("success", false)?;
            d.set_item("error", format!("Solver error: {e}"))?;
        }
    }
    Ok(d.into())
}

/// Generate NACA 4-series airfoil using XFOIL's exact algorithm.
///
/// Returns a list of (x, y) tuples.
#[pyfunction]
#[pyo3(signature = (designation, n_points_per_side=None))]
fn generate_naca4(designation: u32, n_points_per_side: Option<usize>) -> Vec<(f64, f64)> {
    flat_from_points(&naca::naca4(designation, n_points_per_side))
}

/// Repanel airfoil using XFOIL's curvature-based algorithm.
///
/// Returns a list of (x, y) tuples.
#[pyfunction]
#[pyo3(signature = (coords, n_panels=160, curv_param=1.0, te_le_ratio=0.15, te_spacing_ratio=0.667))]
fn repanel_xfoil(
    coords: Vec<f64>,
    n_panels: usize,
    curv_param: f64,
    te_le_ratio: f64,
    te_spacing_ratio: f64,
) -> Vec<(f64, f64)> {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return vec![];
    }
    let points = points_from_flat(&coords);
    let spline = match CubicSpline::from_points(&points) {
        Ok(s) => s,
        Err(_) => return vec![],
    };
    let params = PanelingParams {
        curv_param,
        te_le_ratio,
        te_spacing_ratio,
    };
    flat_from_points(&spline.resample_xfoil(n_panels, &params))
}

/// Deflect a flap on an airfoil by rotating points aft of the hinge.
///
/// Returns a list of (x, y) tuples with the flap applied.
#[pyfunction]
#[pyo3(signature = (coords, hinge_x_frac, deflection_deg, hinge_y_frac=0.5))]
fn deflect_flap(
    coords: Vec<f64>,
    hinge_x_frac: f64,
    deflection_deg: f64,
    hinge_y_frac: f64,
) -> Vec<(f64, f64)> {
    if coords.len() < 8 || coords.len() % 2 != 0 {
        return vec![];
    }
    let pts = points_from_flat(&coords);
    let result = rustfoil_core::flap::xfoil_flap(&pts, hinge_x_frac, hinge_y_frac, deflection_deg);
    flat_from_points(&result)
}

fn interp_y(surface: &[Point], target_x: f64) -> f64 {
    for i in 0..surface.len().saturating_sub(1) {
        let (x0, x1) = (surface[i].x, surface[i + 1].x);
        if (target_x >= x0 && target_x <= x1) || (target_x >= x1 && target_x <= x0) {
            let dx = x1 - x0;
            if dx.abs() < 1e-15 { return surface[i].y; }
            let t = (target_x - x0) / dx;
            return surface[i].y + t * (surface[i + 1].y - surface[i].y);
        }
    }
    0.0
}

/// Parse a Selig/Lednicer .dat file and return coordinate tuples.
///
/// Returns a list of (x, y) tuples. Skips header lines automatically.
#[pyfunction]
fn parse_dat_file(path: &str) -> PyResult<Vec<(f64, f64)>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;
    let mut coords = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                coords.push((x, y));
            }
        }
    }
    Ok(coords)
}

fn solve_one_faithful(body: &Body, alpha_deg: f64, options: &XfoilOptions) -> FaithfulResult {
    match solve_body_oper_point(body, AlphaSpec::AlphaDeg(alpha_deg), options) {
        Ok(r) => FaithfulResult {
            alpha_deg: r.alpha_deg,
            cl: r.cl, cd: r.cd, cm: r.cm,
            converged: r.converged, iterations: r.iterations, residual: r.residual,
            x_tr_upper: r.x_tr_upper, x_tr_lower: r.x_tr_lower,
            cd_friction: r.cd_friction, cd_pressure: r.cd_pressure,
            success: true, error: None,
        },
        Err(e) => FaithfulResult {
            alpha_deg, cl: 0.0, cd: 0.0, cm: 0.0,
            converged: false, iterations: 0, residual: 0.0,
            x_tr_upper: 1.0, x_tr_lower: 1.0,
            cd_friction: 0.0, cd_pressure: 0.0,
            success: false, error: Some(format!("{e}")),
        },
    }
}

fn faithful_result_to_pydict(py: Python<'_>, r: &FaithfulResult) -> PyResult<Py<PyDict>> {
    let d = PyDict::new(py);
    d.set_item("alpha_deg", r.alpha_deg)?;
    d.set_item("cl", r.cl)?;
    d.set_item("cd", r.cd)?;
    d.set_item("cm", r.cm)?;
    d.set_item("converged", r.converged)?;
    d.set_item("iterations", r.iterations)?;
    d.set_item("residual", r.residual)?;
    d.set_item("x_tr_upper", r.x_tr_upper)?;
    d.set_item("x_tr_lower", r.x_tr_lower)?;
    d.set_item("cd_friction", r.cd_friction)?;
    d.set_item("cd_pressure", r.cd_pressure)?;
    d.set_item("success", r.success)?;
    d.set_item("error", r.error.as_deref())?;
    Ok(d.into())
}

/// Batch viscous analysis: solve multiple alphas in parallel via rayon.
///
/// Returns a list of dicts (same schema as analyze_faithful), one per alpha.
#[pyfunction]
#[pyo3(signature = (coords, alphas, reynolds=1.0e6, mach=0.0, ncrit=9.0, max_iterations=100))]
fn analyze_faithful_batch(
    py: Python<'_>,
    coords: Vec<f64>,
    alphas: Vec<f64>,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> PyResult<Vec<Py<PyDict>>> {
    let err_msg = if coords.len() < 6 || coords.len() % 2 != 0 {
        Some("Invalid coordinates".to_string())
    } else {
        None
    };
    if let Some(msg) = err_msg {
        return alphas.iter().map(|&a| {
            let r = FaithfulResult {
                alpha_deg: a, cl: 0.0, cd: 0.0, cm: 0.0,
                converged: false, iterations: 0, residual: 0.0,
                x_tr_upper: 1.0, x_tr_lower: 1.0,
                cd_friction: 0.0, cd_pressure: 0.0,
                success: false, error: Some(msg.clone()),
            };
            faithful_result_to_pydict(py, &r)
        }).collect();
    }

    let points = points_from_flat(&coords);
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            let msg = format!("Geometry error: {e}");
            return alphas.iter().map(|&a| {
                let r = FaithfulResult {
                    alpha_deg: a, cl: 0.0, cd: 0.0, cm: 0.0,
                    converged: false, iterations: 0, residual: 0.0,
                    x_tr_upper: 1.0, x_tr_lower: 1.0,
                    cd_friction: 0.0, cd_pressure: 0.0,
                    success: false, error: Some(msg.clone()),
                };
                faithful_result_to_pydict(py, &r)
            }).collect();
        }
    };

    let options = XfoilOptions {
        reynolds, mach, ncrit, max_iterations,
        ..Default::default()
    };

    let results: Vec<FaithfulResult> = py.allow_threads(|| {
        alphas.par_iter()
            .map(|&a| solve_one_faithful(&body, a, &options))
            .collect()
    });

    results.iter()
        .map(|r| faithful_result_to_pydict(py, r))
        .collect()
}

/// Batch inviscid analysis: solve multiple alphas in parallel via rayon.
///
/// Returns a list of dicts (same schema as analyze_inviscid), one per alpha.
#[pyfunction]
fn analyze_inviscid_batch(
    py: Python<'_>,
    coords: Vec<f64>,
    alphas: Vec<f64>,
) -> PyResult<Vec<Py<PyDict>>> {
    use rustfoil_solver::inviscid::{FlowConditions, InviscidSolver};

    if coords.len() < 6 || coords.len() % 2 != 0 {
        return alphas.iter().map(|_| {
            let d = PyDict::new(py);
            d.set_item("success", false)?;
            d.set_item("error", "Invalid coordinates")?;
            Ok(d.into())
        }).collect();
    }

    let points = points_from_flat(&coords);
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            let msg = format!("Geometry error: {e}");
            return alphas.iter().map(|_| {
                let d = PyDict::new(py);
                d.set_item("success", false)?;
                d.set_item("error", &msg)?;
                Ok(d.into())
            }).collect();
        }
    };

    let solver = InviscidSolver::new();

    struct InviscidResult {
        cl: f64, cm: f64,
        cp: Vec<f64>, cp_x: Vec<f64>,
        success: bool, error: Option<String>,
    }

    let cp_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint().x).collect();

    let results: Vec<InviscidResult> = py.allow_threads(|| {
        alphas.par_iter()
            .map(|&a| {
                let flow = FlowConditions::with_alpha_deg(a);
                match solver.solve(&[body.clone()], &flow) {
                    Ok(s) => InviscidResult {
                        cl: s.cl, cm: s.cm, cp: s.cp, cp_x: cp_x.clone(),
                        success: true, error: None,
                    },
                    Err(e) => InviscidResult {
                        cl: 0.0, cm: 0.0, cp: vec![], cp_x: vec![],
                        success: false, error: Some(format!("Solver error: {e}")),
                    },
                }
            })
            .collect()
    });

    results.iter()
        .map(|r| {
            let d = PyDict::new(py);
            d.set_item("cl", r.cl)?;
            d.set_item("cm", r.cm)?;
            d.set_item("cp", &r.cp)?;
            d.set_item("cp_x", &r.cp_x)?;
            d.set_item("success", r.success)?;
            d.set_item("error", r.error.as_deref())?;
            Ok(d.into())
        })
        .collect()
}

#[pymodule]
fn _rustfoil(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(analyze_faithful, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_inviscid, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_faithful_batch, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_inviscid_batch, m)?)?;
    m.add_function(wrap_pyfunction!(generate_naca4, m)?)?;
    m.add_function(wrap_pyfunction!(repanel_xfoil, m)?)?;
    m.add_function(wrap_pyfunction!(deflect_flap, m)?)?;
    m.add_function(wrap_pyfunction!(parse_dat_file, m)?)?;
    Ok(())
}
