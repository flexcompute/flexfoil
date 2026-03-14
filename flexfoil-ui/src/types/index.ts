/**
 * FlexFoil UI Type Definitions
 */

/** A 2D point */
export interface Point {
  x: number;
  y: number;
}

/** An airfoil coordinate with optional metadata */
export interface AirfoilPoint extends Point {
  /** Arc-length parameter (0 to 1, TE to TE around the foil) */
  s?: number;
  /** Surface type: upper or lower */
  surface?: 'upper' | 'lower';
}

/** Control modes for airfoil manipulation */
export type ControlMode = 'parameters' | 'camber-spline' | 'thickness-spline';
export type SolverMode = 'viscous' | 'inviscid';

/** Spacing panel modes */
export type SpacingPanelMode = 'simple' | 'advanced';

/** SSP interpolation modes */
export type SSPInterpolation = 'linear' | 'spline';

/** Advanced SSP visualization modes */
export type SSPVisualization = 'plot' | 'foil';

/** A Bezier handle attached to an airfoil point */
export interface BezierHandle {
  /** Index of the airfoil point this handle belongs to */
  pointIndex: number;
  /** Position of the handle (can be off-surface) */
  position: Point;
  /** Handle type: incoming or outgoing tangent */
  type: 'in' | 'out';
}

/** A B-spline control point (typically off-surface) */
export interface BSplineControlPoint extends Point {
  /** Unique identifier */
  id: string;
  /** Weight for NURBS (default 1.0 for B-spline) */
  weight?: number;
}

/** A camber line control point */
export interface CamberControlPoint {
  /** Unique identifier */
  id: string;
  /** Chord position (0 to 1, LE to TE) */
  x: number;
  /** Camber value (typically -0.1 to 0.1) */
  y: number;
}

/** A thickness distribution control point */
export interface ThicknessControlPoint {
  /** Unique identifier */
  id: string;
  /** Chord position (0 to 1, LE to TE) */
  x: number;
  /** Half-thickness value (typically 0 to 0.15) */
  t: number;
}

/** SSP spacing knot for repaneling */
export interface SpacingKnot {
  /** Position along arc-length (0 to 1) */
  S: number;
  /** Relative spacing factor (higher = coarser spacing) */
  F: number;
}

/** NACA 4-series parameters */
export interface Naca4Params {
  /** Maximum camber as fraction of chord (0-0.09) */
  m: number;
  /** Position of max camber as fraction of chord (0-0.9) */
  p: number;
  /** Maximum thickness as fraction of chord (0-0.40) */
  t: number;
  /** Number of points to generate */
  nPoints: number;
}

/** Airfoil state for the main store */
export interface AirfoilState {
  /** Name/identifier */
  name: string;
  /** Raw coordinates (TE -> upper -> LE -> lower -> TE) */
  coordinates: AirfoilPoint[];
  /** Panelized coordinates after repaneling */
  panels: AirfoilPoint[];
  /** Current control mode */
  controlMode: ControlMode;
  /** Bezier handles (legacy, kept for compatibility) */
  bezierHandles: BezierHandle[];
  /** B-spline control points (legacy, kept for compatibility) */
  bsplineControlPoints: BSplineControlPoint[];
  /** B-spline degree */
  bsplineDegree: number;
  /** Spacing knots for SSP repaneling */
  spacingKnots: SpacingKnot[];
  /** Number of output panels */
  nPanels: number;
  /** Curvature weight for blending SSP and curvature-based spacing (0-1) */
  curvatureWeight: number;
  /** Current angle of attack for visualization (degrees) */
  displayAlpha: number;
  /** Reynolds number for faithful viscous analysis */
  reynolds: number;
  /** Freestream Mach number (compressibility correction) */
  mach: number;
  /** eN transition criterion (lower = earlier transition) */
  ncrit: number;
  /** Maximum viscous solver iterations */
  maxIterations: number;
  /** Active solver mode */
  solverMode: SolverMode;
  /** Polar sweep data */
  polarData: PolarPoint[];
  /** Spacing panel mode: simple (curvature-based) or advanced (SSP) */
  spacingPanelMode: SpacingPanelMode;
  /** SSP interpolation mode: linear (piecewise) or spline */
  sspInterpolation: SSPInterpolation;
  /** Advanced SSP visualization: plot (F vs S) or foil (normal displacement) */
  sspVisualization: SSPVisualization;
  
  // Camber/thickness editing state
  /** Camber line control points (when in camber-spline mode) */
  camberControlPoints: CamberControlPoint[];
  /** Thickness distribution control points (when in thickness-spline mode) */
  thicknessControlPoints: ThicknessControlPoint[];
  /** Thickness scale factor (1.0 = original) */
  thicknessScale: number;
  /** Camber scale factor (1.0 = original) */
  camberScale: number;
  /** Original airfoil for scaling reference */
  baseCoordinates: AirfoilPoint[];
}

/** A polar data point */
export interface PolarPoint {
  alpha: number;
  cl: number;
  cd?: number;
  cm: number;
}

/** A row from the runs SQLite table — one solver evaluation */
export interface RunRow {
  id: number;
  airfoil_name: string;
  airfoil_hash: string;
  alpha: number;
  reynolds: number;
  mach: number;
  ncrit: number;
  n_panels: number;
  max_iter: number;
  cl: number | null;
  cd: number | null;
  cm: number | null;
  converged: boolean;
  iterations: number | null;
  residual: number | null;
  x_tr_upper: number | null;
  x_tr_lower: number | null;
  success: boolean;
  error: string | null;
  created_at: string;
  session_id: string | null;
}

/** Panel configuration for FlexLayout */
export interface PanelConfig {
  id: string;
  name: string;
  component: string;
}

/** Canvas viewport state */
export interface ViewportState {
  /** Center of view in airfoil coordinates */
  center: Point;
  /** Zoom level (pixels per chord unit) */
  zoom: number;
  /** Canvas width in pixels */
  width: number;
  /** Canvas height in pixels */
  height: number;
}

/** Performance metrics for visualization */
export interface PerfMetrics {
  /** Last frame duration in milliseconds */
  frameTime: number;
  /** Rolling average frame time */
  avgFrameTime: number;
  /** Number of active smoke particles */
  particleCount: number;
  /** Frames per second */
  fps: number;
}

/** Visualization settings state */
export interface VisualizationState {
  // Display toggles
  showGrid: boolean;
  showCurve: boolean;
  showPanels: boolean;
  showPoints: boolean;
  showControls: boolean;
  showStreamlines: boolean;
  showSmoke: boolean;
  showPsiContours: boolean;
  showCp: boolean;
  showForces: boolean;
  showBoundaryLayer: boolean;
  showWake: boolean;
  
  // Animation options
  enableMorphing: boolean;
  morphDuration: number;
  
  // Streamline options
  streamlineDensity: number;
  adaptiveStreamlines: boolean;
  
  // Smoke options
  smokeDensity: number;
  smokeParticlesPerBlob: number;
  smokeSpawnInterval: number;
  smokeMaxAge: number;
  smokeWaveSpacing: number;  // Distance between smoke waves in chord lengths
  smokeResetCounter: number; // Incremented to trigger smoke reset
  
  // Flow speed multiplier
  flowSpeed: number;
  
  // Cp visualization options
  cpDisplayMode: 'color' | 'bars' | 'both';
  cpBarScale: number;
  
  // Force vector options
  forceScale: number;
  
  // Boundary layer overlay options
  blThicknessScale: number;
  
  // GPU acceleration
  useGPU: boolean;
  gpuAvailable: boolean;
  perfMetrics: PerfMetrics;
}
