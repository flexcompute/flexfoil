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
export type ControlMode = 'surface' | 'bezier' | 'bspline';

/** Spacing/paneling modes */
export type SpacingMode = 'ssp' | 'xfoil';

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
  /** Bezier handles (when in bezier mode) */
  bezierHandles: BezierHandle[];
  /** B-spline control points (when in bspline mode) */
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
  /** Polar sweep data */
  polarData: PolarPoint[];
  
  // Viscous analysis state
  /** Reynolds number */
  reynolds: number;
  /** Solver mode (inviscid or viscous) */
  solverMode: SolverMode;
  /** Turbulent boundary layer model */
  turbulentModel: TurbulentModel;
  /** Critical N-factor for transition prediction (typically 9.0) */
  nCrit: number;
  /** Cached viscous solution */
  viscousSolution: ViscousSolution | null;
  /** Boundary layer distribution data */
  blData: BLDistribution | null;
  /** Whether to auto-recompute on geometry change */
  isAutoCompute: boolean;
}

/** A polar data point */
export interface PolarPoint {
  alpha: number;
  cl: number;
  cd?: number;
  cm: number;
  reynolds?: number;
  x_tr_upper?: number;
  x_tr_lower?: number;
  converged?: boolean;
}

/** Boundary layer distribution data */
export interface BLDistribution {
  s_upper: number[];
  s_lower: number[];
  x_upper: number[];
  x_lower: number[];
  theta_upper: number[];
  theta_lower: number[];
  delta_star_upper: number[];
  delta_star_lower: number[];
  h_upper: number[];
  h_lower: number[];
  cf_upper: number[];
  cf_lower: number[];
  x_tr_upper: number;
  x_tr_lower: number;
}

/** Viscous solution result */
export interface ViscousSolution {
  cl: number;
  cd: number;
  cd_friction: number;
  cd_pressure: number;
  cm: number;
  cp: number[];
  cp_x: number[];
  x_tr_upper: number;
  x_tr_lower: number;
  converged: boolean;
  iterations: number;
  reynolds: number;
  alpha: number;
}

/** Solver mode */
export type SolverMode = 'inviscid' | 'viscous';

/** 
 * Turbulent boundary layer model selection.
 * - 0: Head's entrainment method (fast, good for attached flows)
 * - 1: XFOIL's Cτ lag-dissipation method (accurate, good for separation)
 * - 2: Green's lag-entrainment method (best history effects)
 */
export type TurbulentModel = 0 | 1 | 2;

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
