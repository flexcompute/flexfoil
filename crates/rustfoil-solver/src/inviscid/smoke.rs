//! Smoke particle visualization system.
//!
//! Provides a blob-based particle system for visualizing flow.
//! Particles are spawned upstream and advected by the velocity field.

use super::velocity::velocity_at;
use rustfoil_core::Point;
use std::f64::consts::PI;

/// A particle in the smoke system.
#[derive(Debug, Clone)]
struct Particle {
    x: f64,
    y: f64,
    age: f64,
    blob_id: usize,
}

/// Smoke blob system for flow visualization.
pub struct SmokeSystem {
    particles: Vec<Particle>,
    spawn_points: Vec<(f64, f64)>,
    particles_per_blob: usize,
    blob_radius: f64,
    max_age: f64,
    spawn_interval: f64,
    time_since_spawn: f64,
    next_blob_id: usize,
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
            particles: Vec::new(),
            spawn_points,
            particles_per_blob,
            blob_radius: 0.02,
            max_age: 5.0, // 5 seconds lifetime
            spawn_interval: 0.5, // Spawn every 0.5 seconds
            time_since_spawn: 0.0,
            next_blob_id: 0,
        }
    }

    /// Set spawn interval (seconds between blob spawns).
    pub fn set_spawn_interval(&mut self, interval: f64) {
        self.spawn_interval = interval.max(0.05);
    }

    /// Set maximum particle age (seconds).
    pub fn set_max_age(&mut self, max_age: f64) {
        self.max_age = max_age.max(0.1);
    }

    /// Spawn new blobs at all spawn points.
    fn spawn_blobs(&mut self) {
        for &(sx, sy) in &self.spawn_points {
            let blob_id = self.next_blob_id;
            self.next_blob_id += 1;

            for _ in 0..self.particles_per_blob {
                // Random offset within blob radius
                let angle = rand_float() * 2.0 * PI;
                let radius = rand_float() * self.blob_radius;
                
                self.particles.push(Particle {
                    x: sx + radius * angle.cos(),
                    y: sy + radius * angle.sin(),
                    age: 0.0,
                    blob_id,
                });
            }
        }
    }

    /// Update all particles using the velocity field.
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
        // Check if we should spawn new blobs
        self.time_since_spawn += dt;
        if self.time_since_spawn >= self.spawn_interval {
            self.spawn_blobs();
            self.time_since_spawn = 0.0;
        }

        // Update each particle using RK2 integration
        for particle in &mut self.particles {
            particle.age += dt;

            // Get velocity at current position
            let (u1, v1) = velocity_at(particle.x, particle.y, nodes, gamma, alpha, v_inf);
            
            // Skip if velocity is invalid
            let speed1 = (u1 * u1 + v1 * v1).sqrt();
            if !speed1.is_finite() || speed1 < 1e-8 {
                continue;
            }

            // RK2 midpoint
            let x_mid = particle.x + 0.5 * dt * u1;
            let y_mid = particle.y + 0.5 * dt * v1;
            
            let (u2, v2) = velocity_at(x_mid, y_mid, nodes, gamma, alpha, v_inf);
            
            // Update position
            particle.x += dt * u2;
            particle.y += dt * v2;
        }

        // Remove old particles
        self.particles.retain(|p| p.age < self.max_age);
    }

    /// Get all particle positions as a flat array [x0, y0, x1, y1, ...].
    pub fn get_positions(&self) -> Vec<f64> {
        let mut positions = Vec::with_capacity(self.particles.len() * 2);
        for p in &self.particles {
            positions.push(p.x);
            positions.push(p.y);
        }
        positions
    }

    /// Get alpha (opacity) values for each particle based on age.
    pub fn get_alphas(&self) -> Vec<f64> {
        self.particles
            .iter()
            .map(|p| {
                let t = p.age / self.max_age;
                
                // Fade in quickly, stay visible, fade out at end
                if t < 0.05 {
                    t / 0.05 // Quick fade in
                } else if t > 0.9 {
                    (1.0 - t) / 0.1 // Fade out
                } else {
                    1.0
                }
            })
            .collect()
    }

    /// Get the number of active particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Reset the system (remove all particles).
    pub fn reset(&mut self) {
        self.particles.clear();
        self.time_since_spawn = 0.0;
    }
}

/// Simple pseudo-random number generator (no external crate needed).
/// Returns a value in [0, 1).
fn rand_float() -> f64 {
    // Use a simple LCG with thread-local state
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
        assert_eq!(system.particles.len(), 0);
    }

    #[test]
    fn test_smoke_spawn() {
        let spawn_y = vec![0.0];
        let mut system = SmokeSystem::new(&spawn_y, -0.5, 10);
        system.set_spawn_interval(0.1);
        
        // Empty nodes/gamma for freestream-only flow
        let nodes: Vec<Point> = vec![];
        let gamma: Vec<f64> = vec![];
        
        // Update should trigger spawn
        system.update(&nodes, &gamma, 0.0, 1.0, 0.2);
        
        assert!(system.particle_count() > 0);
    }

    #[test]
    fn test_particle_advection() {
        let spawn_y = vec![0.0];
        let mut system = SmokeSystem::new(&spawn_y, -0.5, 1);
        system.set_spawn_interval(0.01);
        system.blob_radius = 0.0; // No random offset
        
        // Empty nodes for freestream only
        let nodes: Vec<Point> = vec![];
        let gamma: Vec<f64> = vec![];
        
        // Spawn and move
        system.update(&nodes, &gamma, 0.0, 1.0, 0.1);
        
        if system.particle_count() > 0 {
            let pos = system.get_positions();
            // Particle should have moved right (freestream is +x)
            assert!(pos[0] > -0.5);
        }
    }
}
