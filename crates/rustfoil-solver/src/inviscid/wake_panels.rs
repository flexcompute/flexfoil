//! Wake panel geometry generation for viscous-inviscid coupling.
//!
//! This module creates wake panels extending downstream from the trailing edge.
//! Wake panels are essential for accurate drag prediction in the Semi-Direct
//! coupling method.
//!
//! # Wake Panel Properties
//!
//! Unlike body panels, wake panels:
//! - Have zero source strength (no thickness)
//! - Carry vorticity that enforces the Kutta condition
//! - Provide edge velocity distribution for wake BL marching
//!
//! # Reference
//!
//! XFOIL uses a similar wake representation with panels extending 1-2 chords
//! downstream from the trailing edge.

use rustfoil_core::{Point, point};

/// Wake panel configuration.
#[derive(Debug, Clone)]
pub struct WakeConfig {
    /// Wake length in chord lengths (default 1.0)
    pub length_chords: f64,
    /// Number of wake panels (default 10)
    pub n_panels: usize,
    /// Whether to use geometric spacing (denser near TE)
    pub geometric_spacing: bool,
    /// Geometric spacing ratio (1.2 means each panel is 1.2x longer than previous)
    pub spacing_ratio: f64,
}

impl Default for WakeConfig {
    fn default() -> Self {
        Self {
            length_chords: 1.0,
            n_panels: 10,
            geometric_spacing: true,
            spacing_ratio: 1.2,
        }
    }
}

/// Wake panels extending downstream from the trailing edge.
#[derive(Debug, Clone)]
pub struct WakePanels {
    /// Wake panel nodes (n_panels + 1 nodes for n_panels panels)
    pub nodes: Vec<Point>,
    /// Number of wake panels
    pub n_panels: usize,
    /// Total wake length
    pub total_length: f64,
    /// Arc-length coordinate at each node (0 at TE)
    pub s_coords: Vec<f64>,
    /// Trailing edge midpoint (wake origin)
    pub te_midpoint: Point,
    /// Wake direction (unit vector)
    pub direction: (f64, f64),
}

impl WakePanels {
    /// Create wake panels extending from the trailing edge.
    ///
    /// # Arguments
    /// * `te_upper` - Upper trailing edge point
    /// * `te_lower` - Lower trailing edge point
    /// * `chord` - Airfoil chord length
    /// * `alpha` - Angle of attack in radians
    /// * `config` - Wake configuration
    ///
    /// # Returns
    /// Wake panels structure with nodes and arc-length coordinates
    pub fn new(
        te_upper: Point,
        te_lower: Point,
        chord: f64,
        alpha: f64,
        config: &WakeConfig,
    ) -> Self {
        let n = config.n_panels;
        let total_length = config.length_chords * chord;
        
        // TE midpoint (Point is nalgebra Point2, accessed with [0] for x, [1] for y)
        let te_mid_x = 0.5 * (te_upper[0] + te_lower[0]);
        let te_mid_y = 0.5 * (te_upper[1] + te_lower[1]);
        let te_mid = point(te_mid_x, te_mid_y);
        
        // Wake direction: downstream in the freestream direction
        // At alpha=0, this is the +x direction
        let dir_x = alpha.cos();
        let dir_y = alpha.sin();
        
        // Generate panel spacing
        let spacing = if config.geometric_spacing {
            Self::geometric_spacing(n, config.spacing_ratio)
        } else {
            Self::uniform_spacing(n)
        };
        
        // Create nodes along wake
        let mut nodes = Vec::with_capacity(n + 1);
        let mut s_coords = Vec::with_capacity(n + 1);
        
        // First node is at TE midpoint
        nodes.push(te_mid);
        s_coords.push(0.0);
        
        let mut s = 0.0;
        for i in 0..n {
            let ds = spacing[i] * total_length;
            s += ds;
            
            let node = point(te_mid_x + s * dir_x, te_mid_y + s * dir_y);
            
            nodes.push(node);
            s_coords.push(s);
        }
        
        Self {
            nodes,
            n_panels: n,
            total_length,
            s_coords,
            te_midpoint: te_mid,
            direction: (dir_x, dir_y),
        }
    }
    
    /// Generate uniform spacing coefficients that sum to 1.
    fn uniform_spacing(n: usize) -> Vec<f64> {
        vec![1.0 / n as f64; n]
    }
    
    /// Generate geometric spacing coefficients that sum to 1.
    /// First panel is smallest, each subsequent panel is ratio times larger.
    fn geometric_spacing(n: usize, ratio: f64) -> Vec<f64> {
        if n == 0 {
            return vec![];
        }
        if (ratio - 1.0).abs() < 1e-10 {
            return Self::uniform_spacing(n);
        }
        
        // Geometric series: a, ar, ar^2, ..., ar^(n-1)
        // Sum = a * (r^n - 1) / (r - 1) = 1
        // => a = (r - 1) / (r^n - 1)
        let r_n = ratio.powi(n as i32);
        let a = (ratio - 1.0) / (r_n - 1.0);
        
        (0..n).map(|i| a * ratio.powi(i as i32)).collect()
    }
    
    /// Get the panel length at index i.
    pub fn panel_length(&self, i: usize) -> f64 {
        if i >= self.n_panels {
            return 0.0;
        }
        self.s_coords[i + 1] - self.s_coords[i]
    }
    
    /// Get the midpoint of panel i.
    pub fn panel_midpoint(&self, i: usize) -> Point {
        if i >= self.n_panels {
            return self.nodes[self.n_panels];
        }
        point(
            0.5 * (self.nodes[i][0] + self.nodes[i + 1][0]),
            0.5 * (self.nodes[i][1] + self.nodes[i + 1][1]),
        )
    }
    
    /// Get arc-length at panel midpoint.
    pub fn panel_s(&self, i: usize) -> f64 {
        if i >= self.n_panels {
            return self.total_length;
        }
        0.5 * (self.s_coords[i] + self.s_coords[i + 1])
    }
}

/// Create wake panels from body geometry.
///
/// # Arguments
/// * `body_nodes` - Nodes of the airfoil body (counterclockwise from upper TE)
/// * `chord` - Airfoil chord length
/// * `alpha` - Angle of attack in radians
/// * `config` - Wake configuration
///
/// # Returns
/// Wake panels structure
pub fn create_wake_from_body(
    body_nodes: &[Point],
    chord: f64,
    alpha: f64,
    config: &WakeConfig,
) -> WakePanels {
    if body_nodes.is_empty() {
        return WakePanels {
            nodes: vec![point(0.0, 0.0)],
            n_panels: 0,
            total_length: 0.0,
            s_coords: vec![0.0],
            te_midpoint: point(0.0, 0.0),
            direction: (1.0, 0.0),
        };
    }
    
    // XFOIL convention: node 0 is upper TE, last node is lower TE
    let te_upper = body_nodes[0];
    let te_lower = body_nodes[body_nodes.len() - 1];
    
    WakePanels::new(te_upper, te_lower, chord, alpha, config)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_wake_creation() {
        let te_upper = point(1.0, 0.01);
        let te_lower = point(1.0, -0.01);
        let chord = 1.0;
        let alpha = 0.0;
        let config = WakeConfig::default();
        
        let wake = WakePanels::new(te_upper, te_lower, chord, alpha, &config);
        
        assert_eq!(wake.n_panels, 10);
        assert_eq!(wake.nodes.len(), 11);
        assert!((wake.total_length - 1.0).abs() < 1e-10);
        
        // First node at TE midpoint
        assert!((wake.nodes[0][0] - 1.0).abs() < 1e-10);
        assert!((wake.nodes[0][1] - 0.0).abs() < 1e-10);
        
        // Last node at x = 2.0 (one chord downstream)
        assert!((wake.nodes[10][0] - 2.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_geometric_spacing() {
        let spacing = WakePanels::geometric_spacing(5, 1.5);
        
        // Sum should be 1
        let sum: f64 = spacing.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        
        // Each panel should be 1.5x the previous
        for i in 1..spacing.len() {
            let ratio = spacing[i] / spacing[i - 1];
            assert!((ratio - 1.5).abs() < 1e-10);
        }
    }
    
    #[test]
    fn test_uniform_spacing() {
        let spacing = WakePanels::uniform_spacing(4);
        
        assert_eq!(spacing.len(), 4);
        for &s in &spacing {
            assert!((s - 0.25).abs() < 1e-10);
        }
    }
    
    #[test]
    fn test_wake_at_angle() {
        let te_upper = point(1.0, 0.01);
        let te_lower = point(1.0, -0.01);
        let chord = 1.0;
        let alpha = 5.0_f64.to_radians();
        let config = WakeConfig {
            n_panels: 5,
            geometric_spacing: false,
            ..Default::default()
        };
        
        let wake = WakePanels::new(te_upper, te_lower, chord, alpha, &config);
        
        // Wake should extend at angle alpha
        let last_node = wake.nodes[5];
        let dx = last_node[0] - wake.te_midpoint[0];
        let dy = last_node[1] - wake.te_midpoint[1];
        let wake_angle = dy.atan2(dx);
        
        assert!((wake_angle - alpha).abs() < 1e-10);
    }
    
    #[test]
    fn test_panel_methods() {
        let te_upper = point(1.0, 0.0);
        let te_lower = point(1.0, 0.0);
        let chord = 1.0;
        let alpha = 0.0;
        let config = WakeConfig {
            n_panels: 4,
            geometric_spacing: false,
            length_chords: 1.0,
            spacing_ratio: 1.0,
        };
        
        let wake = WakePanels::new(te_upper, te_lower, chord, alpha, &config);
        
        // Uniform spacing: each panel is 0.25 chord
        assert!((wake.panel_length(0) - 0.25).abs() < 1e-10);
        assert!((wake.panel_length(3) - 0.25).abs() < 1e-10);
        
        // Panel midpoint at index 0 should be at x=1.125
        let mid = wake.panel_midpoint(0);
        assert!((mid[0] - 1.125).abs() < 1e-10);
        
        // Panel s at index 0 should be 0.125
        assert!((wake.panel_s(0) - 0.125).abs() < 1e-10);
    }
    
    #[test]
    fn test_create_from_body() {
        let body_nodes = vec![
            point(1.0, 0.01),   // Upper TE
            point(0.5, 0.06),
            point(0.0, 0.0),    // LE
            point(0.5, -0.06),
            point(1.0, -0.01),  // Lower TE
        ];
        
        let config = WakeConfig::default();
        let wake = create_wake_from_body(&body_nodes, 1.0, 0.0, &config);
        
        // TE midpoint should be at (1.0, 0.0)
        assert!((wake.te_midpoint[0] - 1.0).abs() < 1e-10);
        assert!(wake.te_midpoint[1].abs() < 0.01);
    }
}
