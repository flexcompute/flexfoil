# Task 15: Implement UPDATE and UESET

## Objective
Implement solution update with limiting for stability.

## Prerequisites
- Task 14 (march) completed

## Context
- UPDATE: xbl.f line 1253 - apply Newton updates with limiting
- UESET: xpanel.f line 1758 - set edge velocities from inviscid
- DSLIM: limit large changes for stability

## Deliverables

### src/update.rs
```rust
//! Solution update with limiting
//!
//! XFOIL Reference: xbl.f UPDATE (line 1253), xpanel.f UESET (line 1758)

use rustfoil_bl::state::BlStation;

/// Configuration for update limiting
pub struct UpdateConfig {
    /// Maximum relative change in θ
    pub max_theta_change: f64,
    /// Maximum relative change in δ*
    pub max_delta_star_change: f64,
    /// Relaxation factor (0-1)
    pub relaxation: f64,
}

impl Default for UpdateConfig {
    fn default() -> Self {
        Self {
            max_theta_change: 0.3,
            max_delta_star_change: 0.3,
            relaxation: 1.0,
        }
    }
}

/// Apply Newton updates to BL stations with limiting
///
/// # Arguments
/// * `stations` - BL stations to update
/// * `deltas` - Newton solution [Δδ*, Δθ, ΔN/ΔCτ] at each station
/// * `config` - Update limiting configuration
///
/// # Returns
/// Maximum relative change (for convergence check)
///
/// # Reference
/// XFOIL xbl.f UPDATE (line 1253)
pub fn update_stations(
    stations: &mut [BlStation],
    deltas: &[[f64; 3]],
    config: &UpdateConfig,
) -> f64 {
    let mut max_change = 0.0;
    
    for (station, delta) in stations.iter_mut().zip(deltas.iter()) {
        // Limit displacement thickness change
        let d_delta_star = limit_change(
            delta[0],
            station.delta_star,
            config.max_delta_star_change,
        );
        
        // Limit momentum thickness change
        let d_theta = limit_change(
            delta[1],
            station.theta,
            config.max_theta_change,
        );
        
        // Apply relaxation
        let relax = config.relaxation;
        station.delta_star += relax * d_delta_star;
        station.theta += relax * d_theta;
        
        // Third variable (N or Cτ)
        if station.is_laminar {
            station.ampl += relax * delta[2];
            // Clamp amplification factor
            station.ampl = station.ampl.max(0.0);
        } else {
            station.ctau += relax * delta[2];
            // Clamp shear stress coefficient
            station.ctau = station.ctau.max(0.0);
        }
        
        // Track maximum change
        let rel_change = (d_delta_star / station.delta_star).abs()
            .max((d_theta / station.theta).abs());
        max_change = max_change.max(rel_change);
    }
    
    max_change
}

/// Limit change to prevent instability
fn limit_change(delta: f64, current: f64, max_relative: f64) -> f64 {
    let max_abs = max_relative * current.abs();
    delta.clamp(-max_abs, max_abs)
}

/// Set edge velocities from inviscid solution
///
/// # Arguments
/// * `stations` - BL stations
/// * `ue_inviscid` - Edge velocities from inviscid solver
/// * `dij` - Mass defect influence matrix
///
/// # Reference
/// XFOIL xpanel.f UESET (line 1758)
pub fn set_edge_velocities(
    stations: &mut [BlStation],
    ue_inviscid: &[f64],
    dij: &nalgebra::DMatrix<f64>,
) {
    let n = stations.len();
    
    // Ue = Ue_inviscid + DIJ * (Ue * δ*)
    // This accounts for displacement effect on inviscid flow
    
    for i in 0..n {
        let mut ue_correction = 0.0;
        for j in 0..n {
            let mass_defect = stations[j].u * stations[j].delta_star;
            ue_correction += dij[(i, j)] * mass_defect;
        }
        stations[i].u = ue_inviscid[i] + ue_correction;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_limit_change() {
        // Large change should be limited
        let limited = limit_change(0.5, 1.0, 0.3);
        assert!((limited - 0.3).abs() < 1e-10);
        
        // Small change should pass through
        let unlimited = limit_change(0.1, 1.0, 0.3);
        assert!((unlimited - 0.1).abs() < 1e-10);
    }
}
```

## Next Task
After completion, proceed to TASK_16_VISCAL.md
