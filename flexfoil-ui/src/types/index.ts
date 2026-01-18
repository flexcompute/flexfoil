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
}

/** A polar data point */
export interface PolarPoint {
  alpha: number;
  cl: number;
  cm: number;
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
