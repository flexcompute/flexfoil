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
use rustfoil_core::{point, Body, CubicSpline, GeometryError, Point};
use rustfoil_solver::inviscid::{FlowConditions, InviscidSolver};
use std::fs;
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
