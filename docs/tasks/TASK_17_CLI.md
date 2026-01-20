# Task 17: Add Viscous CLI Commands

## Objective
Extend the CLI with viscous analysis commands.

## Prerequisites
- Task 16 (VISCAL) completed

## Deliverables

### Update rustfoil-cli/src/main.rs

Add new commands:
```rust
/// Viscous analysis at single alpha
#[derive(Parser)]
struct ViscousCmd {
    /// Input airfoil file (.dat)
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
    
    /// Critical N factor
    #[arg(short, long, default_value = "9.0")]
    ncrit: f64,
    
    /// Output format (json, table)
    #[arg(long, default_value = "table")]
    format: String,
}

/// Viscous polar generation
#[derive(Parser)]
struct ViscousPolarCmd {
    /// Input airfoil file (.dat)
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
    
    /// Critical N factor
    #[arg(short, long, default_value = "9.0")]
    ncrit: f64,
    
    /// Output file (csv)
    #[arg(short, long)]
    output: Option<PathBuf>,
    
    /// Use parallel execution
    #[arg(long)]
    parallel: bool,
}

// Add to CLI enum
#[derive(Subcommand)]
enum Commands {
    // ... existing commands ...
    
    /// Viscous analysis at single angle
    Viscous(ViscousCmd),
    
    /// Generate viscous polar
    ViscousPolar(ViscousPolarCmd),
}

// Implementation
fn run_viscous(cmd: ViscousCmd) -> Result<()> {
    let body = load_airfoil(&cmd.file)?;
    
    let config = ViscousSolverConfig {
        reynolds: cmd.re,
        mach: cmd.mach,
        ncrit: cmd.ncrit,
        ..Default::default()
    };
    
    let result = solve_viscous(&body, cmd.alpha, &config)?;
    
    match cmd.format.as_str() {
        "json" => {
            println!("{}", serde_json::to_string_pretty(&result)?);
        }
        _ => {
            println!("Viscous Analysis Results");
            println!("========================");
            println!("Alpha:     {:>8.3}°", result.alpha);
            println!("CL:        {:>8.4}", result.cl);
            println!("CD:        {:>8.5}", result.cd);
            println!("CM:        {:>8.4}", result.cm);
            println!("x_tr (U):  {:>8.4}", result.x_tr_upper);
            println!("x_tr (L):  {:>8.4}", result.x_tr_lower);
            println!("Converged: {}", result.converged);
        }
    }
    
    Ok(())
}

fn run_viscous_polar(cmd: ViscousPolarCmd) -> Result<()> {
    let body = load_airfoil(&cmd.file)?;
    
    // Parse alpha range "start:end:step"
    let parts: Vec<f64> = cmd.alpha
        .split(':')
        .map(|s| s.parse().unwrap())
        .collect();
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
    
    let results = if cmd.parallel {
        solve_viscous_polar_parallel(&body, &alphas, &config)
    } else {
        alphas.iter()
            .map(|&a| solve_viscous(&body, a, &config))
            .collect()
    };
    
    // Output
    println!("alpha,CL,CD,CM,x_tr_u,x_tr_l,converged");
    for result in results {
        match result {
            Ok(r) => {
                println!("{:.2},{:.4},{:.5},{:.4},{:.4},{:.4},{}",
                    r.alpha, r.cl, r.cd, r.cm, 
                    r.x_tr_upper, r.x_tr_lower, r.converged);
            }
            Err(e) => {
                eprintln!("Error: {}", e);
            }
        }
    }
    
    Ok(())
}
```

## Usage Examples

```bash
# Single point analysis
rustfoil viscous naca0012.dat -a 5.0 -r 3e6

# Generate polar
rustfoil viscous-polar naca0012.dat --alpha "-4:12:0.5" -r 3e6 --parallel

# With Mach number
rustfoil viscous naca0012.dat -a 5.0 -r 3e6 -m 0.3

# Custom transition
rustfoil viscous naca0012.dat -a 5.0 -r 3e6 --ncrit 7.0
```

## Next Task
After completion, proceed to TASK_18_INTEGRATION_TESTS.md
