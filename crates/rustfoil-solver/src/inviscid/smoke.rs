//! Smoke particle visualization system.
//!
//! Provides a blob-based particle system for visualizing flow.
//! Uses pre-computed reference streamlines for efficiency.
//!
//! # Performance
//!
//! The system pre-computes reference streamlines when the flow field changes.
//! Per-frame updates are O(particles) - just array lookups, no velocity calculations.
//!
//! Each blob randomly samples points from the reference streamline, providing
//! visual variation while maintaining physical accuracy.

use super::velocity::{velocity_at, is_inside_airfoil};
use rustfoil_core::Point;

/// A particle in the smoke system.
/// Now just tracks which reference point it's sampling and its age.
#[derive(Debug, Clone)]
struct Particle {
    /// Index into reference streamline
    ref_index: usize,
    /// Current position (copied from reference for efficiency)
    x: f64,
    y: f64,
    /// Age of particle in seconds
    age: f64,
    /// Which blob this particle belongs to
    blob_id: usize,
    /// Time offset for this particle (creates spread within blob)
    time_offset: f64,
}

/// A blob of particles following reference streamlines.
#[derive(Debug, Clone)]
struct Blob {
    /// Which spawn point this blob originates from
    spawn_idx: usize,
    /// When this blob was created (simulation time)
    spawn_time: f64,
    /// Which streamline each particle follows (index into spawn_streamlines[spawn_idx].streamlines)
    particle_streamline_idx: Vec<usize>,
    /// Indices into the reference streamline for each particle
    particle_indices: Vec<usize>,
    /// Random time offsets for each particle (creates visual spread)
    time_offsets: Vec<f64>,
}

/// Pre-computed reference streamline for a spawn point.
#[derive(Debug, Clone)]
struct ReferenceStreamline {
    /// Points along the streamline: (x, y, time)
    /// Time is the cumulative integration time to reach this point
    points: Vec<(f64, f64, f64)>,
}

/// Multiple reference streamlines for a single spawn point (for blob spread).
#[derive(Debug, Clone)]
struct SpawnPointStreamlines {
    /// Multiple streamlines with slight y-offsets for visual variety
    streamlines: Vec<ReferenceStreamline>,
}

/// Smoke blob system for flow visualization.
/// 
/// Uses pre-computed reference streamlines for O(particles) per-frame performance.
pub struct SmokeSystem {
    /// Active blobs
    blobs: Vec<Blob>,
    /// Pre-computed reference streamlines (multiple per spawn point for blob spread)
    spawn_streamlines: Vec<SpawnPointStreamlines>,
    /// Spawn point positions
    spawn_points: Vec<(f64, f64)>,
    /// Number of particles per blob
    particles_per_blob: usize,
    /// Number of streamlines per spawn point (for blob spread)
    streamlines_per_spawn: usize,
    /// Y-offset for spread streamlines
    spread_offset: f64,
    /// Maximum particle age in seconds
    max_age: f64,
    /// Time between blob spawns in seconds
    spawn_interval: f64,
    /// Time since last spawn
    time_since_spawn: f64,
    /// Current simulation time
    current_time: f64,
    /// Next blob ID
    next_blob_id: usize,
    /// Cached flow parameters for detecting changes
    cached_alpha: f64,
    /// Number of points in reference streamlines
    ref_points_count: usize,
}

impl SmokeSystem {
    /// Create a new smoke system.
    ///
    /// # Arguments
    /// * `spawn_y_values` - Y coordinates for spawn points (all at same x)
    /// * `spawn_x` - X coordinate for all spawn points
    /// * `particles_per_blob` - Number of particles per blob
    pub fn new(spawn_y_values: &[f64], spawn_x: f64, particles_per_blob: usize) -> Self {
        let spawn_points: Vec<(f64, f64)> = spawn_y_values
            .iter()
            .map(|&y| (spawn_x, y))
            .collect();

        Self {
            blobs: Vec::new(),
            spawn_streamlines: Vec::new(),
            spawn_points,
            particles_per_blob,
            streamlines_per_spawn: 5, // 5 streamlines per spawn point for blob spread
            spread_offset: 0.015, // Y-offset between streamlines (creates ~0.06 total spread)
            max_age: 5.0,
            spawn_interval: 0.1, // Faster default
            time_since_spawn: 0.0,
            current_time: 0.0,
            next_blob_id: 0,
            cached_alpha: f64::NAN,
            ref_points_count: 1500, // Points per reference streamline
        }
    }

    /// Set spawn interval (seconds between blob spawns).
    pub fn set_spawn_interval(&mut self, interval: f64) {
        self.spawn_interval = interval.max(0.02);
    }

    /// Set maximum particle age (seconds).
    pub fn set_max_age(&mut self, max_age: f64) {
        self.max_age = max_age.max(0.1);
    }

    /// Compute reference streamlines for all spawn points.
    /// Called when flow field changes (geometry or alpha).
    /// 
    /// For each spawn point, computes multiple streamlines with small y-offsets
    /// to provide visual spread for the blob effect.
    fn compute_reference_streamlines(
        &mut self,
        nodes: &[Point],
        gamma: &[f64],
        alpha: f64,
        v_inf: f64,
    ) {
        self.spawn_streamlines.clear();
        
        // Integration parameters
        let dt_target = 0.005; // Target time step
        let x_max = 3.0; // Stop when particle exits domain
        
        // Y-offsets for multiple streamlines per spawn point
        let half_n = self.streamlines_per_spawn / 2;
        let offsets: Vec<f64> = (0..self.streamlines_per_spawn)
            .map(|i| (i as f64 - half_n as f64) * self.spread_offset)
            .collect();
        
        for &(sx, sy) in &self.spawn_points {
            let mut spawn_streamlines = SpawnPointStreamlines {
                streamlines: Vec::with_capacity(self.streamlines_per_spawn),
            };
            
            // Compute multiple streamlines with y-offsets
            for &y_offset in &offsets {
                let mut points = Vec::with_capacity(self.ref_points_count);
                let mut x = sx;
                let mut y = sy + y_offset;
                let mut t = 0.0;
                
                points.push((x, y, t));
                
                // Integrate streamline using RK4
                for _ in 0..self.ref_points_count {
                    // Check if inside airfoil
                    if is_inside_airfoil(x, y, nodes) {
                        break;
                    }
                    
                    // Check if exited domain
                    if x > x_max {
                        break;
                    }
                    
                    // Get velocity
                    let (u, vv) = velocity_at(x, y, nodes, gamma, alpha, v_inf);
                    let speed = (u * u + vv * vv).sqrt();
                    
                    if !speed.is_finite() || speed < 1e-8 {
                        break;
                    }
                    
                    // Adaptive time step based on velocity
                    let dt = (dt_target * v_inf / speed).clamp(0.001, 0.02);
                    
                    // RK4 step
                    if let Some((new_x, new_y)) = rk4_step_velocity(
                        x, y, dt, nodes, gamma, alpha, v_inf
                    ) {
                        x = new_x;
                        y = new_y;
                        t += dt;
                        points.push((x, y, t));
                    } else {
                        break;
                    }
                }
                
                spawn_streamlines.streamlines.push(ReferenceStreamline { points });
            }
            
            self.spawn_streamlines.push(spawn_streamlines);
        }
        
        self.cached_alpha = alpha;
    }

    /// Spawn new blobs at all spawn points.
    fn spawn_blobs(&mut self) {
        for (spawn_idx, spawn_streamlines) in self.spawn_streamlines.iter().enumerate() {
            // Check that we have valid streamlines
            let valid_streamlines: Vec<usize> = spawn_streamlines.streamlines.iter()
                .enumerate()
                .filter(|(_, s)| s.points.len() >= 10)
                .map(|(i, _)| i)
                .collect();
            
            if valid_streamlines.is_empty() {
                continue;
            }
            
            self.next_blob_id += 1;
            
            let mut particle_streamline_idx = Vec::with_capacity(self.particles_per_blob);
            let mut particle_indices = Vec::with_capacity(self.particles_per_blob);
            let mut time_offsets = Vec::with_capacity(self.particles_per_blob);
            
            for _ in 0..self.particles_per_blob {
                // Randomly pick which streamline this particle follows
                let streamline_idx = valid_streamlines[
                    (rand_float() * valid_streamlines.len() as f64) as usize % valid_streamlines.len()
                ];
                particle_streamline_idx.push(streamline_idx);
                
                // Get the selected streamline's length
                let max_idx = spawn_streamlines.streamlines[streamline_idx].points.len();
                
                // Random index weighted toward the front (where particles start)
                let rand_val = rand_float();
                // Use sqrt to bias toward lower indices (front of streamline)
                let idx = ((rand_val * rand_val) * max_idx as f64) as usize;
                particle_indices.push(idx.min(max_idx - 1));
                
                // Small random time offset for additional visual spread
                time_offsets.push(rand_float() * 0.08 - 0.04);
            }
            
            self.blobs.push(Blob {
                spawn_idx,
                spawn_time: self.current_time,
                particle_streamline_idx,
                particle_indices,
                time_offsets,
            });
        }
    }

    /// Update all particles using pre-computed streamlines.
    /// This is O(particles) - just time increments and array lookups.
    ///
    /// # Arguments
    /// * `nodes` - Airfoil node positions
    /// * `gamma` - Vorticity at each node
    /// * `alpha` - Angle of attack (radians)
    /// * `v_inf` - Freestream velocity
    /// * `dt` - Time step (seconds)
    pub fn update(
        &mut self,
        nodes: &[Point],
        gamma: &[f64],
        alpha: f64,
        v_inf: f64,
        dt: f64,
    ) {
        // Check if we need to recompute reference streamlines
        // (This happens when set_flow is called, which updates alpha)
        if self.spawn_streamlines.is_empty() {
            self.compute_reference_streamlines(nodes, gamma, alpha, v_inf);
        }
        
        // Update time
        self.current_time += dt;
        self.time_since_spawn += dt;
        
        // Spawn new blobs if interval elapsed
        if self.time_since_spawn >= self.spawn_interval {
            self.spawn_blobs();
            self.time_since_spawn = 0.0;
        }
        
        // Update blob positions by advancing time
        // Particles move along their pre-computed paths
        for blob in &mut self.blobs {
            let age = self.current_time - blob.spawn_time;
            
            // Update particle indices based on elapsed time
            if let Some(spawn_streamlines) = self.spawn_streamlines.get(blob.spawn_idx) {
                let time_scale = v_inf; // Particles move faster with higher v_inf
                
                for (i, idx) in blob.particle_indices.iter_mut().enumerate() {
                    let streamline_idx = blob.particle_streamline_idx[i];
                    if let Some(ref_streamline) = spawn_streamlines.streamlines.get(streamline_idx) {
                        let particle_time = age * time_scale + blob.time_offsets[i];
                        
                        // Find the reference point at this time
                        let new_idx = find_index_at_time(&ref_streamline.points, particle_time);
                        *idx = new_idx;
                    }
                }
            }
        }
        
        // Remove old blobs
        self.blobs.retain(|blob| {
            let age = self.current_time - blob.spawn_time;
            age < self.max_age
        });
    }

    /// Get all particle positions as a flat array [x0, y0, x1, y1, ...].
    pub fn get_positions(&self) -> Vec<f64> {
        let total_particles: usize = self.blobs.iter()
            .map(|b| b.particle_indices.len())
            .sum();
        
        let mut positions = Vec::with_capacity(total_particles * 2);
        
        for blob in &self.blobs {
            if let Some(spawn_streamlines) = self.spawn_streamlines.get(blob.spawn_idx) {
                for (i, &idx) in blob.particle_indices.iter().enumerate() {
                    let streamline_idx = blob.particle_streamline_idx[i];
                    if let Some(ref_streamline) = spawn_streamlines.streamlines.get(streamline_idx) {
                        if let Some(&(x, y, _)) = ref_streamline.points.get(idx) {
                            positions.push(x);
                            positions.push(y);
                        }
                    }
                }
            }
        }
        
        positions
    }

    /// Get alpha (opacity) values for each particle based on age.
    pub fn get_alphas(&self) -> Vec<f64> {
        let total_particles: usize = self.blobs.iter()
            .map(|b| b.particle_indices.len())
            .sum();
        
        let mut alphas = Vec::with_capacity(total_particles);
        
        for blob in &self.blobs {
            let age = self.current_time - blob.spawn_time;
            let t = age / self.max_age;
            
            // Fade in quickly, stay visible, fade out at end
            let alpha = if t < 0.05 {
                t / 0.05 // Quick fade in
            } else if t > 0.85 {
                (1.0 - t) / 0.15 // Fade out
            } else {
                1.0
            };
            
            // Same alpha for all particles in blob
            for _ in &blob.particle_indices {
                alphas.push(alpha.max(0.0));
            }
        }
        
        alphas
    }

    /// Get the number of active particles.
    pub fn particle_count(&self) -> usize {
        self.blobs.iter().map(|b| b.particle_indices.len()).sum()
    }

    /// Reset the system (remove all particles and cached streamlines).
    pub fn reset(&mut self) {
        self.blobs.clear();
        self.spawn_streamlines.clear();
        self.time_since_spawn = 0.0;
        self.current_time = 0.0;
        self.cached_alpha = f64::NAN;
    }
    
    /// Force recomputation of reference streamlines.
    /// Call this when flow field changes (alpha or geometry).
    pub fn invalidate_cache(&mut self) {
        self.spawn_streamlines.clear();
    }
}

/// Find the index in a reference streamline at a given time.
fn find_index_at_time(points: &[(f64, f64, f64)], target_time: f64) -> usize {
    if points.is_empty() {
        return 0;
    }
    
    // Handle negative time (particle hasn't "started" yet)
    if target_time <= 0.0 {
        return 0;
    }
    
    // Binary search for the time
    let mut lo = 0;
    let mut hi = points.len() - 1;
    
    while lo < hi {
        let mid = (lo + hi) / 2;
        if points[mid].2 < target_time {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    
    lo.min(points.len() - 1)
}

/// RK4 integration step for streamline computation.
fn rk4_step_velocity(
    x: f64,
    y: f64,
    dt: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
) -> Option<(f64, f64)> {
    let field = |px: f64, py: f64| velocity_at(px, py, nodes, gamma, alpha, v_inf);
    
    let (k1x, k1y) = field(x, y);
    let speed1 = (k1x * k1x + k1y * k1y).sqrt();
    
    if !speed1.is_finite() || speed1 < 1e-8 {
        return None;
    }
    
    let (k2x, k2y) = field(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y);
    let (k3x, k3y) = field(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y);
    let (k4x, k4y) = field(x + dt * k3x, y + dt * k3y);
    
    let new_x = x + (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    let new_y = y + (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    
    if new_x.is_finite() && new_y.is_finite() {
        Some((new_x, new_y))
    } else {
        None
    }
}

/// Simple pseudo-random number generator (no external crate needed).
/// Returns a value in [0, 1).
fn rand_float() -> f64 {
    use std::cell::Cell;
    thread_local! {
        static STATE: Cell<u64> = const { Cell::new(0x853c49e6748fea9b) };
    }
    
    STATE.with(|s| {
        let mut state = s.get();
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        s.set(state);
        (state >> 33) as f64 / (1u64 << 31) as f64
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    #[test]
    fn test_smoke_system_creation() {
        let spawn_y = vec![-0.2, 0.0, 0.2];
        let system = SmokeSystem::new(&spawn_y, -0.5, 20);
        
        assert_eq!(system.spawn_points.len(), 3);
        assert_eq!(system.particle_count(), 0);
    }

    #[test]
    fn test_find_index_at_time() {
        let points = vec![
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.1),
            (0.2, 0.0, 0.2),
            (0.3, 0.0, 0.3),
            (0.4, 0.0, 0.4),
        ];
        
        assert_eq!(find_index_at_time(&points, 0.0), 0);
        assert_eq!(find_index_at_time(&points, 0.15), 2);
        assert_eq!(find_index_at_time(&points, 0.5), 4);
        assert_eq!(find_index_at_time(&points, -0.1), 0);
    }

    #[test]
    fn test_smoke_freestream() {
        let spawn_y = vec![0.0];
        let mut system = SmokeSystem::new(&spawn_y, -0.5, 10);
        system.set_spawn_interval(0.05);
        
        // Empty nodes for freestream only
        let nodes: Vec<Point> = vec![];
        let gamma: Vec<f64> = vec![];
        
        // Update should compute streamlines and spawn
        system.update(&nodes, &gamma, 0.0, 1.0, 0.1);
        
        // Should have reference streamlines computed
        assert!(!system.reference_streamlines.is_empty());
    }
}
