use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use rustfoil_core::{naca, point, Body, CubicSpline, PanelingParams, Point};
use rustfoil_xfoil::oper::{solve_body_oper_point, AlphaSpec};
use rustfoil_xfoil::{ReType, XfoilOptions};

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
    reynolds_eff: f64,
    success: bool,
    error: Option<String>,
}

fn re_type_from_int(v: u8) -> ReType {
    match v {
        2 => ReType::Type2,
        3 => ReType::Type3,
        _ => ReType::Type1,
    }
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
#[pyo3(signature = (coords, alpha_deg, reynolds=1.0e6, mach=0.0, ncrit=9.0, max_iterations=100, re_type=1, xstrip_upper=1.0, xstrip_lower=1.0))]
fn analyze_faithful(
    py: Python<'_>,
    coords: Vec<f64>,
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    re_type: u8,
    xstrip_upper: f64,
    xstrip_lower: f64,
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
        re_type: re_type_from_int(re_type),
        xstrip_upper,
        xstrip_lower,
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
            d.set_item("reynolds_eff", r.reynolds_eff)?;
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

/// Fewest coordinates that can plausibly describe one element contour. Shorter
/// runs are treated as stray fragments, not elements (see [`fold_blocks`]).
const MIN_ELEMENT_POINTS: usize = 3;

/// XFOIL / MSES element separator sentinel. A `999.0  999.0` line marks an
/// element boundary; per XFOIL convention any coordinate pair whose both values
/// reach the sentinel is a separator, never a real point.
const ELEMENT_SEPARATOR_SENTINEL: f64 = 999.0;

/// Largest first-to-last-point distance, as a fraction of chord, still read as a
/// closed element contour. Elements in the bundled airfoil library close to
/// within 2.3% of chord, whereas a single surface run spans a full chord.
const MAX_CONTOUR_GAP_FRACTION: f64 = 0.25;

/// One lexical item of the coordinate section of a `.dat` file.
enum DatToken {
    Point((f64, f64)),
    /// Element boundary stated by the file: a `999.0 999.0` line.
    ExplicitSeparator,
    /// Blank line: an element boundary in multi-element files.
    BlankSeparator,
    /// Comment line: a boundary only when the file reads as multi-element, see
    /// [`comment_lines_delimit_elements`].
    Comment,
}

/// Split the coordinate section of a `.dat` file into tokens.
///
/// Lines that do not read as a coordinate pair (name lines, Lednicer counts,
/// prose) contribute nothing, as before: only coordinate pairs become points.
fn tokenize_dat(text: &str) -> Vec<DatToken> {
    let mut tokens = Vec::new();

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            tokens.push(DatToken::BlankSeparator);
            continue;
        }
        if trimmed.starts_with('#') {
            tokens.push(DatToken::Comment);
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                // Checked before the point is kept, so a `999.0 999.0`
                // separator can never be spliced into the contour.
                if x >= ELEMENT_SEPARATOR_SENTINEL && y >= ELEMENT_SEPARATOR_SENTINEL {
                    tokens.push(DatToken::ExplicitSeparator);
                } else {
                    tokens.push(DatToken::Point((x, y)));
                }
            }
        }
    }

    tokens
}

/// Fold a token stream into element blocks, optionally treating comment lines
/// as boundaries.
///
/// A run of fewer than [`MIN_ELEMENT_POINTS`] points cannot be an element
/// contour, and real `.dat` files do contain stray blank lines mid-contour, so
/// such a fragment is re-joined to its neighbouring block: no coordinate is lost
/// and no bogus one-point "element" appears. A fragment following an explicit
/// `999.0` separator is kept as its own block, because the file declared an
/// element there and a degenerate element should be reported, not absorbed.
fn fold_blocks(tokens: &[DatToken], comments_split: bool) -> Vec<Vec<(f64, f64)>> {
    // (points, whether the separator that opened this block was explicit)
    let mut raw: Vec<(Vec<(f64, f64)>, bool)> = Vec::new();
    let mut current: Vec<(f64, f64)> = Vec::new();
    let mut current_explicit = false;

    for token in tokens {
        let explicit = match token {
            DatToken::Point(p) => {
                current.push(*p);
                continue;
            }
            DatToken::ExplicitSeparator => true,
            DatToken::BlankSeparator => false,
            DatToken::Comment => {
                if !comments_split {
                    continue;
                }
                false
            }
        };

        if current.is_empty() {
            // Consecutive separators: an explicit one still marks the block
            // that follows.
            current_explicit = current_explicit || explicit;
        } else {
            raw.push((std::mem::take(&mut current), current_explicit));
            current_explicit = explicit;
        }
    }
    if !current.is_empty() {
        raw.push((current, current_explicit));
    }

    let mut blocks: Vec<Vec<(f64, f64)>> = Vec::new();
    // Fragment held over because there is no preceding block to join it to.
    let mut carry: Vec<(f64, f64)> = Vec::new();

    for (points, explicit) in raw {
        let mut block = std::mem::take(&mut carry);
        block.extend(points);

        if block.len() < MIN_ELEMENT_POINTS && !explicit {
            match blocks.last_mut() {
                Some(previous) => previous.extend(block),
                None => carry = block,
            }
            continue;
        }

        blocks.push(block);
    }
    if !carry.is_empty() {
        blocks.push(carry);
    }

    blocks
}

/// Whether a run of coordinates reads as one element contour on its own: enough
/// points, extent in both directions, and ends that meet near a trailing edge
/// rather than at opposite ends of the chord.
fn looks_like_element_contour(points: &[(f64, f64)]) -> bool {
    if points.len() < MIN_ELEMENT_POINTS {
        return false;
    }

    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(x, y) in points {
        if !x.is_finite() || !y.is_finite() {
            return false;
        }
        x_min = x_min.min(x);
        x_max = x_max.max(x);
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }

    let chord = x_max - x_min;
    if chord <= 0.0 || y_max - y_min <= 0.0 {
        return false;
    }

    let (first, last) = (points[0], points[points.len() - 1]);
    let gap = ((last.0 - first.0).powi(2) + (last.1 - first.1).powi(2)).sqrt();
    gap <= MAX_CONTOUR_GAP_FRACTION * chord
}

/// Whether the comment lines in a file delimit elements rather than annotate
/// one.
///
/// Conservative on purpose: splitting on comments is accepted only when it
/// yields additional blocks *and* every resulting block reads as an element
/// contour on its own. Otherwise the comments are ignored, which is how
/// annotated single-element files (a licence banner, a note between two
/// coordinates) have always been read.
///
/// `rustfoil-cli` applies its own comment rule to the same files; the two are
/// kept in single named helpers so they can be reconciled.
fn comment_lines_delimit_elements(
    candidate: &[Vec<(f64, f64)>],
    baseline: &[Vec<(f64, f64)>],
) -> bool {
    candidate.len() > baseline.len()
        && candidate
            .iter()
            .all(|block| looks_like_element_contour(block))
}

/// Parse the coordinate section of a `.dat` file into one block per element.
///
/// Blank lines and XFOIL `999.0 999.0` separator lines always split elements, so
/// a slat/main/flap file yields one block per element instead of one
/// concatenated multi-loop contour with a (999, 999) point spliced in. Comment
/// lines split only under [`comment_lines_delimit_elements`].
fn parse_dat_elements(text: &str) -> Vec<Vec<(f64, f64)>> {
    let tokens = tokenize_dat(text);
    let baseline = fold_blocks(&tokens, false);
    let candidate = fold_blocks(&tokens, true);

    if comment_lines_delimit_elements(&candidate, &baseline) {
        candidate
    } else {
        baseline
    }
}

/// Parse a Selig/Lednicer .dat file and return coordinate tuples.
///
/// Returns a list of (x, y) tuples. Skips header lines automatically. A
/// multi-element file is reported instead of being concatenated into one
/// multi-loop contour.
#[pyfunction]
fn parse_dat_file(path: &str) -> PyResult<Vec<(f64, f64)>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;
    let elements = parse_dat_elements(&text);

    if elements.len() > 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "{path}: this file contains {} elements; multi-element configurations are not yet \
             supported by parse_dat_file",
            elements.len()
        )));
    }

    // A file with no coordinates at all still returns an empty list, as before;
    // the caller reports it.
    Ok(elements.into_iter().next().unwrap_or_default())
}

fn solve_one_faithful(body: &Body, alpha_deg: f64, options: &XfoilOptions) -> FaithfulResult {
    match solve_body_oper_point(body, AlphaSpec::AlphaDeg(alpha_deg), options) {
        Ok(r) => FaithfulResult {
            alpha_deg: r.alpha_deg,
            cl: r.cl, cd: r.cd, cm: r.cm,
            converged: r.converged, iterations: r.iterations, residual: r.residual,
            x_tr_upper: r.x_tr_upper, x_tr_lower: r.x_tr_lower,
            cd_friction: r.cd_friction, cd_pressure: r.cd_pressure,
            reynolds_eff: r.reynolds_eff,
            success: true, error: None,
        },
        Err(e) => FaithfulResult {
            alpha_deg, cl: 0.0, cd: 0.0, cm: 0.0,
            converged: false, iterations: 0, residual: 0.0,
            x_tr_upper: 1.0, x_tr_lower: 1.0,
            cd_friction: 0.0, cd_pressure: 0.0,
            reynolds_eff: options.reynolds,
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
    d.set_item("reynolds_eff", r.reynolds_eff)?;
    d.set_item("success", r.success)?;
    d.set_item("error", r.error.as_deref())?;
    Ok(d.into())
}

/// Batch viscous analysis: solve multiple alphas in parallel via rayon.
///
/// Returns a list of dicts (same schema as analyze_faithful), one per alpha.
#[pyfunction]
#[pyo3(signature = (coords, alphas, reynolds=1.0e6, mach=0.0, ncrit=9.0, max_iterations=100, re_type=1, xstrip_upper=1.0, xstrip_lower=1.0))]
fn analyze_faithful_batch(
    py: Python<'_>,
    coords: Vec<f64>,
    alphas: Vec<f64>,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    re_type: u8,
    xstrip_upper: f64,
    xstrip_lower: f64,
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
                reynolds_eff: reynolds,
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
                    reynolds_eff: reynolds,
                    success: false, error: Some(msg.clone()),
                };
                faithful_result_to_pydict(py, &r)
            }).collect();
        }
    };

    let options = XfoilOptions {
        reynolds, mach, ncrit, max_iterations,
        re_type: re_type_from_int(re_type),
        xstrip_upper, xstrip_lower,
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

/// Compute boundary-layer distributions for a viscous operating point.
///
/// Returns a dict with keys: x_upper, x_lower, theta_upper, theta_lower,
/// delta_star_upper, delta_star_lower, h_upper, h_lower, cf_upper, cf_lower,
/// ue_upper, ue_lower, x_tr_upper, x_tr_lower, converged, iterations,
/// residual, success, error.
#[pyfunction]
#[pyo3(signature = (coords, alpha_deg, reynolds=1.0e6, mach=0.0, ncrit=9.0, max_iterations=100, re_type=1, xstrip_upper=1.0, xstrip_lower=1.0))]
fn get_bl_distribution(
    py: Python<'_>,
    coords: Vec<f64>,
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    re_type: u8,
    xstrip_upper: f64,
    xstrip_lower: f64,
) -> PyResult<Py<PyDict>> {
    use rustfoil_xfoil::oper::{build_state_from_coords, solve_operating_point_from_state, coords_from_body, AlphaSpec};
    use rustfoil_xfoil::XfoilOptions;

    let d = PyDict::new(py);

    if coords.len() < 6 || coords.len() % 2 != 0 {
        d.set_item("success", false)?;
        d.set_item("error", "Invalid coordinates")?;
        return Ok(d.into());
    }

    let points = points_from_flat(&coords);
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            d.set_item("success", false)?;
            d.set_item("error", format!("Geometry error: {e}"))?;
            return Ok(d.into());
        }
    };

    let body_coords = coords_from_body(&body);
    let options = XfoilOptions {
        reynolds, mach, ncrit, max_iterations,
        re_type: re_type_from_int(re_type),
        xstrip_upper, xstrip_lower,
        ..Default::default()
    };

    let (mut state, factorized) = match build_state_from_coords(
        "airfoil", &body_coords, AlphaSpec::AlphaDeg(alpha_deg), &options,
    ) {
        Ok(v) => v,
        Err(e) => {
            d.set_item("success", false)?;
            d.set_item("error", format!("{e}"))?;
            return Ok(d.into());
        }
    };

    match solve_operating_point_from_state(&mut state, &factorized, &options) {
        Ok(result) => {
            let iblte_upper = state.iblte_upper.min(state.upper_rows.len().saturating_sub(1));
            let iblte_lower = state.iblte_lower.min(state.lower_rows.len().saturating_sub(1));
            let upper = &state.upper_rows[..=iblte_upper];
            let lower = &state.lower_rows[..=iblte_lower];

            d.set_item("x_upper", upper.iter().map(|r| r.x_coord).collect::<Vec<_>>())?;
            d.set_item("x_lower", lower.iter().map(|r| r.x_coord).collect::<Vec<_>>())?;
            d.set_item("theta_upper", upper.iter().map(|r| r.theta).collect::<Vec<_>>())?;
            d.set_item("theta_lower", lower.iter().map(|r| r.theta).collect::<Vec<_>>())?;
            d.set_item("delta_star_upper", upper.iter().map(|r| r.dstr).collect::<Vec<_>>())?;
            d.set_item("delta_star_lower", lower.iter().map(|r| r.dstr).collect::<Vec<_>>())?;
            d.set_item("h_upper", upper.iter().map(|r| r.h).collect::<Vec<_>>())?;
            d.set_item("h_lower", lower.iter().map(|r| r.h).collect::<Vec<_>>())?;
            d.set_item("cf_upper", upper.iter().map(|r| r.cf).collect::<Vec<_>>())?;
            d.set_item("cf_lower", lower.iter().map(|r| r.cf).collect::<Vec<_>>())?;
            d.set_item("ue_upper", upper.iter().map(|r| r.uedg).collect::<Vec<_>>())?;
            d.set_item("ue_lower", lower.iter().map(|r| r.uedg).collect::<Vec<_>>())?;
            d.set_item("x_tr_upper", result.x_tr_upper)?;
            d.set_item("x_tr_lower", result.x_tr_lower)?;
            d.set_item("converged", result.converged)?;
            d.set_item("iterations", result.iterations)?;
            d.set_item("residual", result.residual)?;
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
    m.add_function(wrap_pyfunction!(get_bl_distribution, m)?)?;
    Ok(())
}

#[cfg(test)]
mod dat_parsing_tests {
    use super::*;
    use std::path::PathBuf;

    fn repo_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    fn read(relative: &str) -> String {
        let path = repo_root().join(relative);
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
    }

    /// The parser as it read files before element blocks existed: every line
    /// that yields a coordinate pair, in file order. Used as the reference for
    /// the corpus differential below.
    fn legacy_coords(text: &str) -> Vec<(f64, f64)> {
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
        coords
    }

    /// Chordwise stations, TE -> LE.
    const STATIONS: [f64; 7] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.05, 0.0];
    const POINTS_PER_ELEMENT: usize = STATIONS.len() * 2 - 1;

    /// A closed element contour (TE -> upper -> LE -> lower -> TE) with leading
    /// edge at `x_le` and the given chord.
    fn closed_element(x_le: f64, chord: f64) -> String {
        let half = |t: f64| chord * 0.08 * (std::f64::consts::PI * t).sin();
        let place = |t: f64| x_le + chord * t;

        let mut lines: Vec<String> = STATIONS
            .iter()
            .map(|&t| format!(" {:.6} {:.6}", place(t), half(t)))
            .collect();
        lines.extend(
            STATIONS
                .iter()
                .rev()
                .skip(1)
                .map(|&t| format!(" {:.6} {:.6}", place(t), -half(t))),
        );
        lines.join("\n")
    }

    fn three_elements(joiner: &str) -> String {
        format!(
            "MULTI\n{}{joiner}{}{joiner}{}\n",
            closed_element(-0.1, 0.15),
            closed_element(0.0, 1.0),
            closed_element(0.85, 0.45),
        )
    }

    #[test]
    fn blank_lines_split_elements() {
        let elements = parse_dat_elements(&three_elements("\n\n\n"));

        assert_eq!(elements.len(), 3);
        for element in &elements {
            assert_eq!(element.len(), POINTS_PER_ELEMENT);
        }
    }

    #[test]
    fn sentinel_line_is_a_separator_not_a_coordinate() {
        let elements = parse_dat_elements(&three_elements("\n 999.0 999.0\n"));

        assert_eq!(elements.len(), 3);
        assert_eq!(
            elements.iter().map(|e| e.len()).sum::<usize>(),
            3 * POINTS_PER_ELEMENT
        );
        assert!(
            elements
                .iter()
                .flatten()
                .all(|&(x, y)| x < ELEMENT_SEPARATOR_SENTINEL && y < ELEMENT_SEPARATOR_SENTINEL),
            "a (999, 999) point was kept as a coordinate"
        );
    }

    #[test]
    fn sentinel_tolerates_whitespace_and_format_variation() {
        for separator in [
            "\n   999.   999.  \n",
            "\n999.000000\t999.000000\n",
            "\n 1e3 1e3\n",
        ] {
            let elements = parse_dat_elements(&three_elements(separator));
            assert_eq!(
                elements.len(),
                3,
                "separator {separator:?} was not recognised"
            );
            assert!(elements
                .iter()
                .flatten()
                .all(|&(x, _)| x < ELEMENT_SEPARATOR_SENTINEL));
        }
    }

    #[test]
    fn explicitly_declared_degenerate_element_is_not_absorbed() {
        let text = format!(
            "BROKEN\n{}\n 999.0 999.0\n 0.5 0.0\n",
            closed_element(0.0, 1.0)
        );

        let elements = parse_dat_elements(&text);

        assert_eq!(elements.len(), 2);
        assert_eq!(elements[0].len(), POINTS_PER_ELEMENT);
        assert_eq!(elements[1].len(), 1);
    }

    #[test]
    fn comment_labelled_elements_are_split() {
        let text = format!(
            "# Mini 30P-30N\n# Slat\n{}\n# Main Element\n{}\n# Flap\n{}\n",
            closed_element(-0.1, 0.15),
            closed_element(0.0, 1.0),
            closed_element(0.85, 0.45),
        );

        let elements = parse_dat_elements(&text);

        assert_eq!(elements.len(), 3);
        for element in &elements {
            assert_eq!(element.len(), POINTS_PER_ELEMENT);
        }
    }

    /// A comment that annotates one element must not split it. Surface labels
    /// are the realistic case: each half spans the whole chord and its ends do
    /// not meet, so neither half reads as a contour.
    #[test]
    fn comment_between_surfaces_is_annotation_not_a_boundary() {
        let element = closed_element(0.0, 1.0);
        let lines: Vec<&str> = element.lines().collect();
        let (upper, lower) = lines.split_at(STATIONS.len());
        let text = format!(
            "ANNOTATED\n# upper surface\n{}\n# lower surface\n{}\n",
            upper.join("\n"),
            lower.join("\n"),
        );

        let elements = parse_dat_elements(&text);

        assert_eq!(elements.len(), 1, "an annotated single element was split");
        assert_eq!(elements[0].len(), POINTS_PER_ELEMENT);
        assert_eq!(elements[0], legacy_coords(&text));
    }

    #[test]
    fn stray_blank_line_inside_an_element_is_rejoined() {
        let text = format!(
            "CAP 21 (mini)\n{}\n\n\n 0.998900 -0.001006\n",
            closed_element(0.0, 1.0)
        );

        let elements = parse_dat_elements(&text);

        assert_eq!(elements.len(), 1);
        assert_eq!(elements[0].len(), POINTS_PER_ELEMENT + 1);
        assert_eq!(elements[0].last().unwrap().0, 0.9989);
    }

    #[test]
    fn short_input_is_returned_rather_than_rejected() {
        // Two points is not an airfoil, but reporting that is the caller's job:
        // the parser has always handed back whatever coordinates it found.
        let elements = parse_dat_elements("SHORT\n 1.0 0.0\n 0.0 0.0\n");

        assert_eq!(elements.len(), 1);
        assert_eq!(elements[0].len(), 2);
    }

    #[test]
    fn file_without_coordinates_yields_no_elements() {
        assert!(parse_dat_elements("JUST A NAME\n# and a note\n").is_empty());
        assert!(parse_dat_elements("").is_empty());
    }

    /// The one genuine multi-element file in the bundled airfoil library: three
    /// elements labelled with comments, previously read as a single 664-point
    /// three-loop contour.
    #[test]
    fn real_mda_30p_30n_file_is_three_elements() {
        let text = read("flexfoil-ui/public/airfoils/30p-30n.dat");

        let elements = parse_dat_elements(&text);

        assert_eq!(
            elements.iter().map(|e| e.len()).collect::<Vec<_>>(),
            vec![201, 221, 242]
        );
        // No coordinate is lost or reordered by the split.
        assert_eq!(
            elements.into_iter().flatten().collect::<Vec<_>>(),
            legacy_coords(&text)
        );
    }

    #[test]
    fn real_single_element_files_are_unchanged() {
        for relative in [
            "testdata/naca0012.dat",
            "testdata/naca2412.dat",
            // Licence banner on line 2, before any coordinate.
            "flexfoil-ui/public/airfoils/s9104.dat",
            // Closing trailing-edge point after two blank lines.
            "flexfoil-ui/public/airfoils/cap21c.dat",
        ] {
            let text = read(relative);
            let elements = parse_dat_elements(&text);

            assert_eq!(elements.len(), 1, "{relative} was split");
            assert_eq!(elements[0], legacy_coords(&text), "{relative} changed");
        }
    }

    /// Differential over the whole bundled corpus: no coordinate may be added,
    /// dropped or reordered in any file, and every file that stops reading as a
    /// single element must consist of blocks that each read as an element
    /// contour on their own.
    #[test]
    fn corpus_differential_keeps_single_element_files_intact() {
        let mut checked: Vec<String> = Vec::new();
        let mut multi_element: Vec<String> = Vec::new();

        for directory in ["flexfoil-ui/public/airfoils", "testdata"] {
            let dir = repo_root().join(directory);
            let entries = match std::fs::read_dir(&dir) {
                Ok(entries) => entries,
                // The corpus is not vendored in every checkout; the targeted
                // tests above still cover the behaviour.
                Err(_) => continue,
            };

            for entry in entries.flatten() {
                let path = entry.path();
                if path.extension().and_then(|e| e.to_str()) != Some("dat") {
                    continue;
                }
                let Ok(text) = std::fs::read_to_string(&path) else {
                    continue;
                };
                let name = entry.file_name().to_string_lossy().into_owned();

                let elements = parse_dat_elements(&text);
                let flattened: Vec<(f64, f64)> = elements.iter().flatten().copied().collect();
                assert_eq!(
                    flattened,
                    legacy_coords(&text),
                    "coordinates changed for {}",
                    path.display()
                );

                if elements.len() > 1 {
                    // The split has to justify itself: a file read as several
                    // elements must be several element contours.
                    for (i, element) in elements.iter().enumerate() {
                        assert!(
                            looks_like_element_contour(element),
                            "{}: block {i} ({} points) is not an element contour",
                            path.display(),
                            element.len()
                        );
                    }
                    multi_element.push(name.clone());
                }
                checked.push(name);
            }
        }

        assert!(!checked.is_empty(), "no .dat files were checked");
        // The bundled airfoil library is single-element apart from the MDA
        // 30P-30N high-lift section and its test fixtures.
        assert!(
            multi_element.len() * 100 < checked.len(),
            "{} of {} files were read as multi-element: {multi_element:?}",
            multi_element.len(),
            checked.len()
        );
        if checked.iter().any(|name| name == "30p-30n.dat") {
            assert!(
                multi_element.iter().any(|name| name == "30p-30n.dat"),
                "30p-30n.dat was not recognised as multi-element"
            );
        }
    }

    #[test]
    fn contour_plausibility_discriminates_surfaces_from_elements() {
        let element: Vec<(f64, f64)> = legacy_coords(&closed_element(0.0, 1.0));
        assert!(looks_like_element_contour(&element));

        // Half a contour: ends a chord apart.
        assert!(!looks_like_element_contour(&element[..STATIONS.len()]));
        // Too few points, no thickness, and no chordwise extent.
        assert!(!looks_like_element_contour(&element[..2]));
        assert!(!looks_like_element_contour(&[
            (0.0, 0.0),
            (0.5, 0.0),
            (1.0, 0.0),
            (0.0, 0.0)
        ]));
        assert!(!looks_like_element_contour(&[
            (0.5, -0.1),
            (0.5, 0.0),
            (0.5, 0.1),
            (0.5, -0.1)
        ]));
    }
}
