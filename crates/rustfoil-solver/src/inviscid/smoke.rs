//! Simple smoke particle visualization system.
//!
//! Spawns blobs of particles at regular intervals. Each particle is
//! individually advected through the flow field using RK2 integration.

use super::velocity::{
    is_inside_airfoil, is_inside_polygon, psi_at, psi_at_with_sources, velocity_at,
    velocity_at_with_sources, WakePanels,
};
use rustfoil_core::Point;

/// A particle in the smoke system.
#[derive(Debug, Clone)]
struct Particle {
    x: f64,
    y: f64,
    age: f64,
}

/// Smoke particle system for flow visualization.
pub struct SmokeSystem {
    /// Active particles
    particles: Vec<Particle>,
    /// Spawn point positions (x, y)
    spawn_points: Vec<(f64, f64)>,
    /// Number of particles per blob
    particles_per_blob: usize,
    /// Radius of blob spawn area
    blob_radius: f64,
    /// Maximum particle age in seconds
    max_age: f64,
    /// Time between blob spawns in seconds
    spawn_interval: f64,
    /// Time since last spawn
    time_since_spawn: f64,
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
            max_age: 5.0,
            spawn_interval: 2.0, // Spawn every 2 seconds for clear blob separation
            time_since_spawn: 0.0,
        }
    }

    /// Set spawn points from flat array of (x, y) pairs.
    /// 
    /// # Arguments
    /// * `points` - Flat array [x0, y0, x1, y1, ...]
    pub fn set_spawn_points(&mut self, points: &[f64]) {
        self.spawn_points = points
            .chunks(2)
            .filter(|c| c.len() == 2)
            .map(|c| (c[0], c[1]))
            .collect();
    }

    /// Set spawn interval (seconds between blob spawns).
    pub fn set_spawn_interval(&mut self, interval: f64) {
        self.spawn_interval = interval.max(0.02);
    }

    /// Set maximum particle age (seconds).
    pub fn set_max_age(&mut self, max_age: f64) {
        self.max_age = max_age.max(0.1);
    }

    /// Spawn new blobs at all spawn points.
    fn spawn_blobs(&mut self) {
        for &(sx, sy) in &self.spawn_points {
            for _ in 0..self.particles_per_blob {
                // Random position within blob radius
                let angle = rand_float() * std::f64::consts::TAU;
                let r = rand_float().sqrt() * self.blob_radius;
                let x = sx + r * angle.cos();
                let y = sy + r * angle.sin();
                
                // Small random age offset for visual spread
                let age = rand_float() * 0.05;
                
                self.particles.push(Particle { x, y, age });
            }
        }
    }

    /// Update all particles.
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
        self.time_since_spawn += dt;
        
        // Spawn new blobs if interval elapsed
        if self.time_since_spawn >= self.spawn_interval {
            self.spawn_blobs();
            self.time_since_spawn = 0.0;
        }
        
        // Update each particle using RK2
        for particle in &mut self.particles {
            particle.age += dt;
            
            // Skip if too old (will be removed)
            if particle.age >= self.max_age {
                continue;
            }
            
            // Skip if inside airfoil
            if is_inside_airfoil(particle.x, particle.y, nodes) {
                particle.age = self.max_age; // Mark for removal
                continue;
            }
            
            // RK2 integration (midpoint method)
            let (u1, v1) = velocity_at(particle.x, particle.y, nodes, gamma, alpha, v_inf);
            
            let mid_x = particle.x + 0.5 * dt * u1;
            let mid_y = particle.y + 0.5 * dt * v1;
            
            let (u2, v2) = velocity_at(mid_x, mid_y, nodes, gamma, alpha, v_inf);
            
            particle.x += dt * u2;
            particle.y += dt * v2;
        }
        
        // Remove old particles
        self.particles.retain(|p| p.age < self.max_age);
    }

    /// Update particles using a source-inclusive viscous field.
    pub fn update_with_sources(
        &mut self,
        nodes: &[Point],
        gamma: &[f64],
        sigma: &[f64],
        alpha: f64,
        v_inf: f64,
        wake_panels: Option<&WakePanels>,
        effective_body: Option<&[Point]>,
        dt: f64,
    ) {
        self.time_since_spawn += dt;

        if self.time_since_spawn >= self.spawn_interval {
            self.spawn_blobs();
            self.time_since_spawn = 0.0;
        }

        for particle in &mut self.particles {
            particle.age += dt;
            if particle.age >= self.max_age {
                continue;
            }

            if is_inside_airfoil(particle.x, particle.y, nodes)
                || effective_body.is_some_and(|poly| is_inside_polygon(particle.x, particle.y, poly))
            {
                particle.age = self.max_age;
                continue;
            }

            let (u1, v1) = velocity_at_with_sources(
                particle.x,
                particle.y,
                nodes,
                gamma,
                sigma,
                alpha,
                v_inf,
                wake_panels,
            );

            let mid_x = particle.x + 0.5 * dt * u1;
            let mid_y = particle.y + 0.5 * dt * v1;

            let (u2, v2) = velocity_at_with_sources(
                mid_x,
                mid_y,
                nodes,
                gamma,
                sigma,
                alpha,
                v_inf,
                wake_panels,
            );

            let new_x = particle.x + dt * u2;
            let new_y = particle.y + dt * v2;
            if effective_body.is_some_and(|poly| is_inside_polygon(new_x, new_y, poly)) {
                particle.age = self.max_age;
                continue;
            }

            particle.x = new_x;
            particle.y = new_y;
        }

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
                    t / 0.05
                } else if t > 0.85 {
                    (1.0 - t) / 0.15
                } else {
                    1.0
                }
            })
            .collect()
    }

    /// Get stream function (psi) values for each particle.
    /// 
    /// This is used to determine which side of the dividing streamline
    /// each particle is on. Compare with psi_0 to determine above/below.
    pub fn get_psi_values(
        &self,
        nodes: &[Point],
        gamma: &[f64],
        alpha: f64,
        v_inf: f64,
    ) -> Vec<f64> {
        self.particles
            .iter()
            .map(|p| psi_at(p.x, p.y, nodes, gamma, alpha, v_inf))
            .collect()
    }

    /// Get stream function values using a source-inclusive viscous field.
    pub fn get_psi_values_with_sources(
        &self,
        nodes: &[Point],
        gamma: &[f64],
        sigma: &[f64],
        alpha: f64,
        v_inf: f64,
        wake_panels: Option<&WakePanels>,
    ) -> Vec<f64> {
        self.particles
            .iter()
            .map(|p| {
                psi_at_with_sources(p.x, p.y, nodes, gamma, sigma, alpha, v_inf, wake_panels)
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
    
    /// Invalidate cache (no-op for simple system, kept for API compatibility).
    pub fn invalidate_cache(&mut self) {
        // No cache in simple system
    }
}

/// Simple pseudo-random number generator.
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

    #[test]
    fn test_smoke_system_creation() {
        let spawn_y = vec![-0.2, 0.0, 0.2];
        let system = SmokeSystem::new(&spawn_y, -0.5, 20);
        
        assert_eq!(system.spawn_points.len(), 3);
        assert_eq!(system.particle_count(), 0);
    }
}
