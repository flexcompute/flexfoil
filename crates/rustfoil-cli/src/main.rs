//! RustFoil CLI - Command-line airfoil analysis tool.
//!
//! This is primarily a development and testing tool. For production use,
//! the WASM interface is recommended.
//!
//! # Usage
//!
//! ```bash
//! # Analyze a single airfoil at one angle of attack
//! rustfoil analyze naca0012.dat --alpha 5.0
//!
//! # Generate a polar sweep
//! rustfoil polar naca0012.dat --alpha-start -5 --alpha-end 15 --alpha-step 1
//!
//! # Repanel an airfoil
//! rustfoil repanel naca0012.dat --panels 100 --output naca0012_repaneled.dat
//! ```

use clap::{Parser, Subcommand};
use rayon::prelude::*;
use rustfoil_bl;
use rustfoil_core::{point, Body, CubicSpline, GeometryError, PanelingParams, Point};
use rustfoil_solver::inviscid::{FlowConditions, InviscidSolver};
use rustfoil_solver::viscous::{
    compute_arc_from_stagnation, compute_arc_lengths, extract_surface,
    find_stagnation_by_sign_change, initialize_surface_stations, solve_viscous_two_surfaces,
    ViscousResult, ViscousSolverConfig, ViscousSetup,
};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
enum CliError {
    #[error("Failed to read file: {0}")]
    FileRead(#[from] std::io::Error),

    #[error("Geometry error: {0}")]
    Geometry(#[from] GeometryError),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Solver error: {0}")]
    Solver(String),
}

#[derive(Parser)]
#[command(name = "rustfoil")]
#[command(author = "Harry")]
#[command(version = "0.1.0")]
#[command(about = "High-performance airfoil analysis engine", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Analyze an airfoil at a single angle of attack
    Analyze {
        /// Path to airfoil coordinate file (.dat)
        #[arg(value_name = "FILE")]
        file: PathBuf,

        /// Angle of attack in degrees
        #[arg(short, long, default_value = "0.0")]
        alpha: f64,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Generate a polar (Cl vs alpha sweep)
    Polar {
        /// Path to airfoil coordinate file (.dat)
        #[arg(value_name = "FILE")]
        file: PathBuf,

        /// Starting angle of attack (degrees)
        #[arg(long, default_value = "-5.0")]
        alpha_start: f64,

        /// Ending angle of attack (degrees)
        #[arg(long, default_value = "15.0")]
        alpha_end: f64,

        /// Angle of attack step (degrees)
        #[arg(long, default_value = "1.0")]
        alpha_step: f64,
    },

    /// Repanel an airfoil using cosine spacing
    Repanel {
        /// Path to input airfoil coordinate file (.dat)
        #[arg(value_name = "FILE")]
        file: PathBuf,

        /// Number of panels
        #[arg(short = 'n', long, default_value = "100")]
        panels: usize,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Display information about an airfoil file
    Info {
        /// Path to airfoil coordinate file (.dat)
        #[arg(value_name = "FILE")]
        file: PathBuf,
    },

    /// Viscous analysis at a single angle of attack
    Viscous(ViscousCmd),

    /// Generate viscous polar (Cl, Cd vs alpha sweep)
    ViscousPolar(ViscousPolarCmd),
}

/// Viscous analysis at single angle of attack
#[derive(Parser)]
struct ViscousCmd {
    /// Input airfoil file (.dat)
    #[arg(value_name = "FILE")]
    file: PathBuf,

    /// Angle of attack (degrees)
    #[arg(short, long)]
    alpha: f64,

    /// Reynolds number
    #[arg(short, long)]
    re: f64,

    /// Mach number
    #[arg(short, long, default_value = "0.0")]
    mach: f64,

    /// Critical N factor for transition
    #[arg(short, long, default_value = "9.0")]
    ncrit: f64,

    /// Number of panels (repanels to XFOIL-style distribution)
    #[arg(short = 'p', long, default_value = "160")]
    panels: usize,

    /// Output format (json, table)
    #[arg(long, default_value = "table")]
    format: String,

    /// Debug output file for XFOIL comparison
    #[arg(long)]
    debug: Option<PathBuf>,
}

/// Viscous polar generation
#[derive(Parser)]
struct ViscousPolarCmd {
    /// Input airfoil file (.dat)
    #[arg(value_name = "FILE")]
    file: PathBuf,

    /// Alpha range (start:end:step)
    #[arg(short, long)]
    alpha: String,

    /// Reynolds number
    #[arg(short, long)]
    re: f64,

    /// Mach number
    #[arg(short, long, default_value = "0.0")]
    mach: f64,

    /// Critical N factor for transition
    #[arg(short, long, default_value = "9.0")]
    ncrit: f64,

    /// Output file (csv)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use parallel execution
    #[arg(long)]
    parallel: bool,
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Analyze {
            file,
            alpha,
            verbose,
        } => run_analyze(&file, alpha, verbose),
        Commands::Polar {
            file,
            alpha_start,
            alpha_end,
            alpha_step,
        } => run_polar(&file, alpha_start, alpha_end, alpha_step),
        Commands::Repanel {
            file,
            panels,
            output,
        } => run_repanel(&file, panels, output),
        Commands::Info { file } => run_info(&file),
        Commands::Viscous(cmd) => run_viscous(cmd),
        Commands::ViscousPolar(cmd) => run_viscous_polar(cmd),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_analyze(file: &PathBuf, alpha: f64, verbose: bool) -> Result<(), CliError> {
    let (name, points) = load_airfoil(file)?;
    let body = Body::from_points(&name, &points)?;

    if verbose {
        println!("Airfoil: {}", name);
        println!("Points: {}", points.len());
        println!("Panels: {}", body.n_panels());
        println!("Chord: {:.4}", body.chord());
        println!("Arc length: {:.4}", body.arc_length());
        println!();
    }

    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha);

    let solution = solver
        .solve(&[body], &flow)
        .map_err(|e| CliError::Solver(e.to_string()))?;

    println!("Analysis Results (α = {:.2}°):", alpha);
    println!("  Cl = {:.4}", solution.cl);
    println!("  Cm = {:.4}", solution.cm);

    if verbose {
        println!();
        println!("Cp distribution:");
        println!("  {:>8}  {:>10}", "x/c", "Cp");
        for (i, cp) in solution.cp.iter().enumerate() {
            // Approximate x/c (would be better with actual panel data)
            let x = i as f64 / solution.cp.len() as f64;
            println!("  {:>8.4}  {:>10.4}", x, cp);
        }
    }

    Ok(())
}

fn run_polar(
    file: &PathBuf,
    alpha_start: f64,
    alpha_end: f64,
    alpha_step: f64,
) -> Result<(), CliError> {
    let (name, points) = load_airfoil(file)?;
    let body = Body::from_points(&name, &points)?;

    println!("Polar for: {}", name);
    println!();
    println!("{:>8}  {:>10}  {:>10}", "Alpha", "Cl", "Cm");
    println!("{:-<8}  {:-<10}  {:-<10}", "", "", "");

    let solver = InviscidSolver::new();
    let mut alpha = alpha_start;

    while alpha <= alpha_end + 1e-9 {
        let flow = FlowConditions::with_alpha_deg(alpha);

        match solver.solve(&[body.clone()], &flow) {
            Ok(solution) => {
                println!(
                    "{:>8.2}  {:>10.4}  {:>10.4}",
                    alpha, solution.cl, solution.cm
                );
            }
            Err(e) => {
                println!("{:>8.2}  {:>10}  {:>10}", alpha, "FAIL", e);
            }
        }

        alpha += alpha_step;
    }

    Ok(())
}

fn run_repanel(file: &PathBuf, n_panels: usize, output: Option<PathBuf>) -> Result<(), CliError> {
    let (name, points) = load_airfoil(file)?;

    println!("Repaneling {} with {} panels...", name, n_panels);

    let spline = CubicSpline::from_points(&points)?;
    let new_points = spline.resample_cosine(n_panels + 1);

    // Write output
    let output_path = output.unwrap_or_else(|| {
        let stem = file.file_stem().unwrap().to_str().unwrap();
        file.with_file_name(format!("{}_repaneled.dat", stem))
    });

    let mut content = format!("{}\n", name);
    for p in &new_points {
        content.push_str(&format!("  {:.6}  {:.6}\n", p.x, p.y));
    }

    fs::write(&output_path, content)?;
    println!("Wrote {} points to {:?}", new_points.len(), output_path);

    Ok(())
}

fn run_info(file: &PathBuf) -> Result<(), CliError> {
    let (name, points) = load_airfoil(file)?;
    let body = Body::from_points(&name, &points)?;

    println!("Airfoil: {}", name);
    println!("File: {:?}", file);
    println!();
    println!("Geometry:");
    println!("  Points: {}", points.len());
    println!("  Panels: {}", body.n_panels());
    println!("  Chord: {:.6}", body.chord());
    println!("  Arc length: {:.6}", body.arc_length());
    println!();

    // Find extrema
    let x_min = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let x_max = points
        .iter()
        .map(|p| p.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let y_max = points
        .iter()
        .map(|p| p.y)
        .fold(f64::NEG_INFINITY, f64::max);

    println!("Bounds:");
    println!("  X: [{:.6}, {:.6}]", x_min, x_max);
    println!("  Y: [{:.6}, {:.6}]", y_min, y_max);
    println!("  Thickness: {:.4} (max)", y_max - y_min);

    Ok(())
}

/// Run viscous analysis at a single angle of attack.
///
/// This bridges the inviscid solver to the viscous BL solver:
/// 1. Run inviscid solution at alpha
/// 2. Extract panel geometry and edge velocities
/// 3. Split into upper/lower surfaces at stagnation
/// 4. Initialize BL stations for each surface
/// 5. Run viscous solver
fn run_viscous_analysis(
    body: &Body,
    alpha: f64,
    config: &ViscousSolverConfig,
) -> Result<ViscousResult, CliError> {
    // Step 1: Run inviscid solver
    let solver = InviscidSolver::new();
    let factorized = solver
        .factorize(&[body.clone()])
        .map_err(|e| CliError::Solver(e.to_string()))?;

    let flow = FlowConditions::with_alpha_deg(alpha);
    let inv_solution = factorized.solve_alpha(&flow);

    // Step 2: Extract node positions (p1 of each panel)
    // The inviscid solver uses node-based gamma, so we need node positions
    let panels = body.panels();
    let n_panels = panels.len();
    
    // Extract node coordinates (p1 of each panel)
    let mut node_x: Vec<f64> = panels.iter().map(|p| p.p1.x).collect();
    let mut node_y: Vec<f64> = panels.iter().map(|p| p.p1.y).collect();
    
    // Check if we need to add the closing node (blunt TE)
    let gamma_len = inv_solution.gamma.len();
    if gamma_len == n_panels + 1 {
        // Blunt TE - add the last node (p2 of last panel)
        let last_panel = panels.last().unwrap();
        node_x.push(last_panel.p2.x);
        node_y.push(last_panel.p2.y);
    }

    // In XFOIL's linear vortex panel method, gamma IS the surface velocity
    // QINV = GAMU (inviscid tangential velocity equals vortex strength)
    let ue_inviscid = inv_solution.gamma.clone();

    // Debug: check inviscid solution
    eprintln!(
        "Inviscid: CL={:.4}, CM={:.4}",
        inv_solution.cl, inv_solution.cm
    );
    
    // Check gamma/Ue values around LE
    let min_ue_idx = ue_inviscid
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    let max_ue_idx = ue_inviscid
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    eprintln!(
        "Gamma: min={:.4} at idx={} (x={:.4}), max={:.4} at idx={} (x={:.4})",
        ue_inviscid[min_ue_idx], min_ue_idx, node_x[min_ue_idx],
        ue_inviscid[max_ue_idx], max_ue_idx, node_x[max_ue_idx]
    );

    // Step 3: Find stagnation and split surfaces
    // Find the leading edge (minimum x) for reference
    let le_idx = node_x
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Use sign change for velocity-based stagnation 
    let stag_by_sign = find_stagnation_by_sign_change(&ue_inviscid);

    // Debug: show both methods
    eprintln!(
        "LE at idx={} (x={:.6}), Sign-change at idx={} (x={:.6}, Ue={:.6})",
        le_idx, node_x[le_idx], stag_by_sign, node_x[stag_by_sign], ue_inviscid[stag_by_sign]
    );

    // Split at the velocity stagnation point (where Ue changes sign)
    // This ensures neither surface includes the stagnation region
    let stag_idx = stag_by_sign;

    eprintln!(
        "Splitting at stagnation: idx={}, x={:.6}",
        stag_idx, node_x[stag_idx]
    );

    // Extract upper surface (stag to TE, going downstream)
    let (upper_x, upper_y, upper_ue) =
        extract_surface(stag_idx, &node_x, &node_y, &ue_inviscid, true);
    let upper_arc = compute_arc_from_stagnation(&upper_x, &upper_y);

    eprintln!(
        "Upper: {} pts, arc[0]={:.6}, arc[1]={:.6}, Ue[0]={:.6}, Ue[1]={:.6}",
        upper_arc.len(),
        upper_arc.get(0).unwrap_or(&0.0),
        upper_arc.get(1).unwrap_or(&0.0),
        upper_ue.get(0).unwrap_or(&0.0),
        upper_ue.get(1).unwrap_or(&0.0),
    );

    // Extract lower surface (stag to TE, reversed to go downstream)
    let (lower_x, lower_y, lower_ue) =
        extract_surface(stag_idx, &node_x, &node_y, &ue_inviscid, false);
    let lower_arc = compute_arc_from_stagnation(&lower_x, &lower_y);

    eprintln!(
        "Lower: {} pts, arc[0]={:.6}, arc[1]={:.6}, Ue[0]={:.6}, Ue[1]={:.6}",
        lower_arc.len(),
        lower_arc.get(0).unwrap_or(&0.0),
        lower_arc.get(1).unwrap_or(&0.0),
        lower_ue.get(0).unwrap_or(&0.0),
        lower_ue.get(1).unwrap_or(&0.0),
    );

    // Debug: show first few lower surface values
    eprintln!("Lower first 5: x={:?}, Ue={:?}",
        &lower_x[..5.min(lower_x.len())],
        &lower_ue[..5.min(lower_ue.len())]
    );

    // Step 4: Initialize BL stations for each surface
    // Handle edge case where a surface has fewer than 2 stations
    if upper_arc.len() < 2 || lower_arc.len() < 2 {
        return Err(CliError::Solver(format!(
            "Surface extraction failed: upper={} lower={} stations (need >= 2 each)",
            upper_arc.len(),
            lower_arc.len()
        )));
    }

    let mut upper_stations = initialize_surface_stations(&upper_arc, &upper_ue, config.reynolds);
    let mut lower_stations = initialize_surface_stations(&lower_arc, &lower_ue, config.reynolds);

    // Step 5: Build DIJ matrix for Newton iteration (full airfoil)
    let arc_lengths = compute_arc_lengths(&node_x, &node_y);
    let setup = ViscousSetup::from_raw(
        arc_lengths.clone(),
        ue_inviscid.clone(),
        node_x.clone(),
        node_y.clone(),
    );

    // Step 6: Solve viscous for both surfaces
    // For now, use the new two-surface solver that marches each surface separately
    let mut result = solve_viscous_two_surfaces(
        &mut upper_stations,
        &mut lower_stations,
        &upper_ue,
        &lower_ue,
        &setup.dij,
        config,
    )
    .map_err(|e| CliError::Solver(e.to_string()))?;

    // Set alpha in result
    result.alpha = alpha;

    Ok(result)
}

fn run_viscous(cmd: ViscousCmd) -> Result<(), CliError> {
    let (name, points) = load_airfoil(&cmd.file)?;

    // Repanel using XFOIL-style curvature-based distribution
    let spline = CubicSpline::from_points(&points)
        .map_err(|e| CliError::Geometry(e))?;
    let repaneled = spline.resample_xfoil(cmd.panels, &PanelingParams::default());
    
    eprintln!(
        "Repaneled from {} to {} points (XFOIL-style distribution)",
        points.len(),
        repaneled.len()
    );

    let body = Body::from_points(&name, &repaneled)?;

    // Initialize debug output if requested
    if let Some(ref debug_path) = cmd.debug {
        rustfoil_bl::init_debug(debug_path);
        eprintln!("Debug output enabled: {}", debug_path.display());
    }

    let config = ViscousSolverConfig {
        reynolds: cmd.re,
        mach: cmd.mach,
        ncrit: cmd.ncrit,
        ..Default::default()
    };

    let result = run_viscous_analysis(&body, cmd.alpha, &config)?;

    // Finalize debug output
    if cmd.debug.is_some() {
        rustfoil_bl::finalize_debug();
        eprintln!("Debug output written");
    }

    match cmd.format.as_str() {
        "json" => {
            // Create a JSON-serializable struct
            let json_output = serde_json::json!({
                "alpha": result.alpha,
                "cl": result.cl,
                "cd": result.cd,
                "cm": result.cm,
                "x_tr_upper": result.x_tr_upper,
                "x_tr_lower": result.x_tr_lower,
                "converged": result.converged,
                "iterations": result.iterations,
                "residual": result.residual,
                "cd_friction": result.cd_friction,
                "cd_pressure": result.cd_pressure,
                "x_separation": result.x_separation
            });
            println!("{}", serde_json::to_string_pretty(&json_output).unwrap());
        }
        _ => {
            println!("Viscous Analysis Results");
            println!("========================");
            println!("Airfoil:   {}", name);
            println!("Re:        {:.2e}", cmd.re);
            println!("Mach:      {:.2}", cmd.mach);
            println!("Ncrit:     {:.1}", cmd.ncrit);
            println!();
            println!("Alpha:     {:>8.3}°", result.alpha);
            println!("CL:        {:>8.4}", result.cl);
            println!("CD:        {:>8.5}", result.cd);
            println!("CM:        {:>8.4}", result.cm);
            println!("x_tr (U):  {:>8.4}", result.x_tr_upper);
            println!("x_tr (L):  {:>8.4}", result.x_tr_lower);
            println!("Converged: {}", result.converged);
            println!("Iters:     {}", result.iterations);
            if let Some(x_sep) = result.x_separation {
                println!("x_sep:     {:>8.4}", x_sep);
            }
        }
    }

    Ok(())
}

fn run_viscous_polar(cmd: ViscousPolarCmd) -> Result<(), CliError> {
    let (name, points) = load_airfoil(&cmd.file)?;
    let body = Body::from_points(&name, &points)?;

    // Parse alpha range "start:end:step"
    let parts: Vec<f64> = cmd
        .alpha
        .split(':')
        .map(|s| {
            s.parse().map_err(|_| CliError::Parse {
                line: 0,
                message: format!("Invalid alpha range format: {}", cmd.alpha),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    if parts.len() != 3 {
        return Err(CliError::Parse {
            line: 0,
            message: format!(
                "Alpha range must be 'start:end:step', got: {}",
                cmd.alpha
            ),
        });
    }

    let (start, end, step) = (parts[0], parts[1], parts[2]);

    let alphas: Vec<f64> = (0..)
        .map(|i| start + i as f64 * step)
        .take_while(|&a| a <= end + 0.001)
        .collect();

    let config = ViscousSolverConfig {
        reynolds: cmd.re,
        mach: cmd.mach,
        ncrit: cmd.ncrit,
        ..Default::default()
    };

    // Compute results
    let results: Vec<Result<ViscousResult, CliError>> = if cmd.parallel {
        alphas
            .par_iter()
            .map(|&alpha| run_viscous_analysis(&body, alpha, &config))
            .collect()
    } else {
        alphas
            .iter()
            .map(|&alpha| run_viscous_analysis(&body, alpha, &config))
            .collect()
    };

    // Output header
    let header = "alpha,CL,CD,CM,x_tr_u,x_tr_l,converged";

    // Build output lines
    let mut output_lines = vec![header.to_string()];
    for result in results {
        match result {
            Ok(r) => {
                output_lines.push(format!(
                    "{:.2},{:.4},{:.5},{:.4},{:.4},{:.4},{}",
                    r.alpha, r.cl, r.cd, r.cm, r.x_tr_upper, r.x_tr_lower, r.converged
                ));
            }
            Err(e) => {
                eprintln!("Error: {}", e);
            }
        }
    }

    // Write to file or stdout
    if let Some(output_path) = cmd.output {
        let mut file = File::create(&output_path)?;
        for line in &output_lines {
            writeln!(file, "{}", line)?;
        }
        println!(
            "Wrote {} polar points to {:?}",
            output_lines.len() - 1,
            output_path
        );
    } else {
        // Print header info
        println!("# Viscous Polar: {}", name);
        println!("# Re = {:.2e}, M = {:.2}, Ncrit = {:.1}", cmd.re, cmd.mach, cmd.ncrit);
        println!();
        for line in &output_lines {
            println!("{}", line);
        }
    }

    Ok(())
}

/// Load an airfoil from a Selig/XFOIL format .dat file.
///
/// Format:
/// ```text
/// NACA 0012
///   1.000000  0.000000
///   0.950000  0.011234
///   ...
/// ```
fn load_airfoil(path: &PathBuf) -> Result<(String, Vec<Point>), CliError> {
    let content = fs::read_to_string(path)?;
    let mut lines = content.lines();

    // First line is the name
    let name = lines
        .next()
        .ok_or_else(|| CliError::Parse {
            line: 1,
            message: "Empty file".to_string(),
        })?
        .trim()
        .to_string();

    let mut points = Vec::new();

    for (i, line) in lines.enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue; // Skip malformed lines
        }

        let x: f64 = parts[0].parse().map_err(|_| CliError::Parse {
            line: i + 2,
            message: format!("Invalid x coordinate: {}", parts[0]),
        })?;

        let y: f64 = parts[1].parse().map_err(|_| CliError::Parse {
            line: i + 2,
            message: format!("Invalid y coordinate: {}", parts[1]),
        })?;

        points.push(point(x, y));
    }

    if points.len() < 3 {
        return Err(CliError::Parse {
            line: 0,
            message: format!("Too few points: {} (need at least 3)", points.len()),
        });
    }

    Ok((name, points))
}
