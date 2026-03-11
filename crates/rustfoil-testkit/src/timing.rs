use std::cmp::Ordering;
use std::time::Instant;

#[derive(Debug, Clone, Copy)]
pub struct BenchmarkConfig {
    pub warmup_runs: usize,
    pub sample_runs: usize,
    pub inner_loops: usize,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            warmup_runs: 2,
            sample_runs: 7,
            inner_loops: 10_000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct BenchmarkStats {
    pub samples_seconds: Vec<f64>,
    pub median_seconds: f64,
    pub min_seconds: f64,
    pub max_seconds: f64,
    pub inner_loops: usize,
}

pub fn benchmark_closure<F>(config: BenchmarkConfig, mut f: F) -> BenchmarkStats
where
    F: FnMut(),
{
    for _ in 0..config.warmup_runs {
        for _ in 0..config.inner_loops {
            f();
        }
    }

    let mut samples = Vec::with_capacity(config.sample_runs);
    for _ in 0..config.sample_runs {
        let start = Instant::now();
        for _ in 0..config.inner_loops {
            f();
        }
        samples.push(start.elapsed().as_secs_f64());
    }

    summarize_samples(samples, config.inner_loops)
}

pub fn summarize_samples(mut samples_seconds: Vec<f64>, inner_loops: usize) -> BenchmarkStats {
    samples_seconds.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let len = samples_seconds.len();
    let median_seconds = if len == 0 {
        0.0
    } else if len % 2 == 0 {
        0.5 * (samples_seconds[len / 2 - 1] + samples_seconds[len / 2])
    } else {
        samples_seconds[len / 2]
    };
    let min_seconds = samples_seconds.first().copied().unwrap_or(0.0);
    let max_seconds = samples_seconds.last().copied().unwrap_or(0.0);

    BenchmarkStats {
        samples_seconds,
        median_seconds,
        min_seconds,
        max_seconds,
        inner_loops,
    }
}

pub fn assert_ratio_within(
    rust: &BenchmarkStats,
    fortran_median_seconds: f64,
    max_ratio: f64,
    label: &str,
) {
    assert!(
        fortran_median_seconds > 0.0,
        "{}: Fortran median must be positive, got {}",
        label,
        fortran_median_seconds
    );
    let ratio = rust.median_seconds / fortran_median_seconds;
    assert!(
        ratio <= max_ratio,
        "{}: Rust median {:.6e}s exceeds Fortran median {:.6e}s by {:.3}x (max {:.3}x)",
        label,
        rust.median_seconds,
        fortran_median_seconds,
        ratio,
        max_ratio
    );
}
