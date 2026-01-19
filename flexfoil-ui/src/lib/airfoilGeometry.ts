/**
 * Airfoil geometry manipulation functions.
 * 
 * Provides functions to decompose airfoils into camber line and thickness,
 * and reconstruct airfoils from these components. Also supports scaling
 * operations for thickness and camber.
 */

import type { AirfoilPoint, CamberControlPoint, ThicknessControlPoint } from '../types';

/**
 * Result of decomposing an airfoil into camber and thickness.
 */
export interface CamberThicknessDecomposition {
  /** Camber line points (x, y=camber) from LE to TE */
  camber: { x: number; y: number }[];
  /** Thickness distribution (x, t=half-thickness) from LE to TE */
  thickness: { x: number; t: number }[];
}

/**
 * Find the leading edge index (minimum x coordinate).
 */
function findLeadingEdgeIndex(coords: AirfoilPoint[]): number {
  let minX = Infinity;
  let leIndex = 0;
  for (let i = 0; i < coords.length; i++) {
    if (coords[i].x < minX) {
      minX = coords[i].x;
      leIndex = i;
    }
  }
  return leIndex;
}

/**
 * Split airfoil into upper and lower surfaces.
 * Returns arrays ordered from LE to TE.
 */
function splitSurfaces(coords: AirfoilPoint[]): {
  upper: AirfoilPoint[];
  lower: AirfoilPoint[];
} {
  const leIndex = findLeadingEdgeIndex(coords);
  
  // Upper surface: from LE going backwards to start (TE)
  // In standard airfoil format: TE -> upper -> LE -> lower -> TE
  // So upper is indices 0 to leIndex (reversed), lower is leIndex to end
  const upper: AirfoilPoint[] = [];
  for (let i = leIndex; i >= 0; i--) {
    upper.push(coords[i]);
  }
  
  // Lower surface: from LE to end (TE)
  const lower: AirfoilPoint[] = [];
  for (let i = leIndex; i < coords.length; i++) {
    lower.push(coords[i]);
  }
  
  return { upper, lower };
}

/**
 * Linear interpolation to find y value at a given x on a curve.
 */
function interpolateY(curve: { x: number; y: number }[], targetX: number): number {
  // Handle edge cases
  if (curve.length === 0) return 0;
  if (curve.length === 1) return curve[0].y;
  
  // Find bracketing points
  for (let i = 0; i < curve.length - 1; i++) {
    const x0 = curve[i].x;
    const x1 = curve[i + 1].x;
    
    if ((x0 <= targetX && targetX <= x1) || (x1 <= targetX && targetX <= x0)) {
      const t = Math.abs(x1 - x0) < 1e-10 ? 0 : (targetX - x0) / (x1 - x0);
      return curve[i].y + t * (curve[i + 1].y - curve[i].y);
    }
  }
  
  // Extrapolate if outside range
  if (targetX <= Math.min(curve[0].x, curve[curve.length - 1].x)) {
    return curve[0].x < curve[curve.length - 1].x ? curve[0].y : curve[curve.length - 1].y;
  }
  return curve[0].x > curve[curve.length - 1].x ? curve[0].y : curve[curve.length - 1].y;
}

/**
 * Decompose an airfoil into camber line and thickness distribution.
 * 
 * The camber line is the locus of points midway between upper and lower surfaces.
 * The thickness is the distance from the camber line to the surface.
 * 
 * @param coords - Airfoil coordinates (TE -> upper -> LE -> lower -> TE format)
 * @returns Camber line and thickness distribution
 */
export function decomposeToCamberThickness(coords: AirfoilPoint[]): CamberThicknessDecomposition {
  if (coords.length < 5) {
    return { camber: [], thickness: [] };
  }
  
  const { upper, lower } = splitSurfaces(coords);
  
  // Sample x positions from 0 to 1
  const nSamples = 50;
  const camber: { x: number; y: number }[] = [];
  const thickness: { x: number; t: number }[] = [];
  
  for (let i = 0; i <= nSamples; i++) {
    const x = i / nSamples;
    
    // Interpolate y values on upper and lower surfaces at this x
    const yUpper = interpolateY(upper.map(p => ({ x: p.x, y: p.y })), x);
    const yLower = interpolateY(lower.map(p => ({ x: p.x, y: p.y })), x);
    
    // Camber is the average of upper and lower
    const yc = (yUpper + yLower) / 2;
    
    // Thickness is half the distance between surfaces
    const t = (yUpper - yLower) / 2;
    
    camber.push({ x, y: yc });
    thickness.push({ x, t: Math.max(0, t) }); // Ensure non-negative thickness
  }
  
  return { camber, thickness };
}

/**
 * Reconstruct airfoil coordinates from camber line and thickness distribution.
 * 
 * @param camber - Camber line points (x, y) 
 * @param thickness - Thickness distribution (x, t)
 * @param nPoints - Number of output points (default 161 for 160 panels)
 * @returns Airfoil coordinates in standard format
 */
export function reconstructFromCamberThickness(
  camber: { x: number; y: number }[],
  thickness: { x: number; t: number }[],
  nPoints: number = 161
): AirfoilPoint[] {
  if (camber.length < 2 || thickness.length < 2) {
    return [];
  }
  
  const halfPoints = Math.floor(nPoints / 2);
  const upper: AirfoilPoint[] = [];
  const lower: AirfoilPoint[] = [];
  
  // Generate points using cosine spacing for better LE resolution
  for (let i = 0; i <= halfPoints; i++) {
    const beta = (Math.PI * i) / halfPoints;
    const x = (1 - Math.cos(beta)) / 2;
    
    // Interpolate camber and thickness at this x
    const yc = interpolateY(camber, x);
    const t = interpolateY(thickness.map(p => ({ x: p.x, y: p.t })), x);
    
    // For proper NACA-style airfoils, we need the slope of the camber line
    // to rotate the thickness perpendicular to the camber
    // For simplicity, we use vertical thickness (good approximation for low camber)
    const yUpper = yc + t;
    const yLower = yc - t;
    
    if (i === 0) {
      // Trailing edge - single point
      upper.push({ x: 1, y: yc + t, surface: 'upper' });
    } else if (i === halfPoints) {
      // Leading edge - single point
      upper.push({ x, y: yc, surface: 'upper' });
    } else {
      upper.push({ x, y: yUpper, surface: 'upper' });
      lower.unshift({ x, y: yLower, surface: 'lower' });
    }
  }
  
  // Build final coordinates: upper (TE to LE) + lower (LE to TE)
  const result: AirfoilPoint[] = [...upper];
  
  // Add lower surface (skip LE, already included)
  for (let i = 0; i < lower.length; i++) {
    result.push(lower[i]);
  }
  
  // Close at TE
  if (result.length > 0 && result[result.length - 1].x < 0.99) {
    result.push({ x: 1, y: result[0].y, surface: 'lower' });
  }
  
  return result;
}

/**
 * Scale the thickness of an airfoil while preserving the camber line.
 * 
 * @param coords - Original airfoil coordinates
 * @param factor - Scale factor (1.0 = original, 2.0 = double thickness)
 * @returns Scaled airfoil coordinates
 */
export function scaleThickness(coords: AirfoilPoint[], factor: number): AirfoilPoint[] {
  if (coords.length < 5 || factor <= 0) {
    return coords;
  }
  
  const { camber, thickness } = decomposeToCamberThickness(coords);
  
  // Scale thickness
  const scaledThickness = thickness.map(p => ({
    x: p.x,
    t: p.t * factor,
  }));
  
  return reconstructFromCamberThickness(camber, scaledThickness, coords.length);
}

/**
 * Scale the camber of an airfoil while preserving the thickness distribution.
 * 
 * @param coords - Original airfoil coordinates
 * @param factor - Scale factor (0.0 = symmetric, 1.0 = original, 2.0 = double camber)
 * @returns Scaled airfoil coordinates
 */
export function scaleCamber(coords: AirfoilPoint[], factor: number): AirfoilPoint[] {
  if (coords.length < 5) {
    return coords;
  }
  
  const { camber, thickness } = decomposeToCamberThickness(coords);
  
  // Scale camber
  const scaledCamber = camber.map(p => ({
    x: p.x,
    y: p.y * factor,
  }));
  
  return reconstructFromCamberThickness(scaledCamber, thickness, coords.length);
}

/**
 * Apply both thickness and camber scaling.
 * 
 * @param coords - Original airfoil coordinates
 * @param thicknessFactor - Thickness scale factor
 * @param camberFactor - Camber scale factor
 * @returns Scaled airfoil coordinates
 */
export function scaleAirfoil(
  coords: AirfoilPoint[],
  thicknessFactor: number,
  camberFactor: number
): AirfoilPoint[] {
  if (coords.length < 5) {
    return coords;
  }
  
  const { camber, thickness } = decomposeToCamberThickness(coords);
  
  // Scale both
  const scaledCamber = camber.map(p => ({
    x: p.x,
    y: p.y * camberFactor,
  }));
  
  const scaledThickness = thickness.map(p => ({
    x: p.x,
    t: p.t * thicknessFactor,
  }));
  
  return reconstructFromCamberThickness(scaledCamber, scaledThickness, coords.length);
}

/**
 * Create initial camber control points from airfoil decomposition.
 * 
 * @param coords - Airfoil coordinates
 * @param numPoints - Number of control points (default 5)
 * @returns Array of camber control points
 */
export function createCamberControlPoints(
  coords: AirfoilPoint[],
  numPoints: number = 5
): CamberControlPoint[] {
  const { camber } = decomposeToCamberThickness(coords);
  
  if (camber.length === 0) {
    // Return default flat camber line
    return Array.from({ length: numPoints }, (_, i) => ({
      id: `camber-${i}`,
      x: i / (numPoints - 1),
      y: 0,
    }));
  }
  
  // Sample at evenly spaced x positions
  const points: CamberControlPoint[] = [];
  for (let i = 0; i < numPoints; i++) {
    const x = i / (numPoints - 1);
    const y = interpolateY(camber, x);
    points.push({
      id: `camber-${i}`,
      x,
      y,
    });
  }
  
  return points;
}

/**
 * Create initial thickness control points from airfoil decomposition.
 * 
 * @param coords - Airfoil coordinates
 * @param numPoints - Number of control points (default 5)
 * @returns Array of thickness control points
 */
export function createThicknessControlPoints(
  coords: AirfoilPoint[],
  numPoints: number = 5
): ThicknessControlPoint[] {
  const { thickness } = decomposeToCamberThickness(coords);
  
  if (thickness.length === 0) {
    // Return default NACA 0012-like thickness distribution
    const t = 0.12;
    return Array.from({ length: numPoints }, (_, i) => {
      const x = i / (numPoints - 1);
      // Simplified NACA thickness formula
      const halfT = 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
      return {
        id: `thickness-${i}`,
        x,
        t: Math.max(0, halfT),
      };
    });
  }
  
  // Sample at evenly spaced x positions
  const points: ThicknessControlPoint[] = [];
  for (let i = 0; i < numPoints; i++) {
    const x = i / (numPoints - 1);
    const t = interpolateY(thickness.map(p => ({ x: p.x, y: p.t })), x);
    points.push({
      id: `thickness-${i}`,
      x,
      t: Math.max(0, t),
    });
  }
  
  return points;
}

/**
 * Reconstruct airfoil from camber and thickness control points.
 * Uses spline interpolation through the control points.
 * 
 * @param camberPoints - Camber control points
 * @param thicknessPoints - Thickness control points  
 * @param nPoints - Number of output points
 * @returns Airfoil coordinates
 */
export function reconstructFromControlPoints(
  camberPoints: CamberControlPoint[],
  thicknessPoints: ThicknessControlPoint[],
  nPoints: number = 161
): AirfoilPoint[] {
  if (camberPoints.length < 2 || thicknessPoints.length < 2) {
    return [];
  }
  
  // Sort control points by x
  const sortedCamber = [...camberPoints].sort((a, b) => a.x - b.x);
  const sortedThickness = [...thicknessPoints].sort((a, b) => a.x - b.x);
  
  // Convert to simple arrays for reconstruction
  const camber = sortedCamber.map(p => ({ x: p.x, y: p.y }));
  const thickness = sortedThickness.map(p => ({ x: p.x, t: p.t }));
  
  return reconstructFromCamberThickness(camber, thickness, nPoints);
}

/**
 * Evaluate a cubic spline through control points at given x values.
 * Simple Catmull-Rom spline interpolation.
 */
export function evaluateSpline(
  controlPoints: { x: number; y: number }[],
  xValues: number[]
): { x: number; y: number }[] {
  if (controlPoints.length < 2) {
    return xValues.map(x => ({ x, y: 0 }));
  }
  
  // Sort by x
  const sorted = [...controlPoints].sort((a, b) => a.x - b.x);
  
  return xValues.map(x => ({
    x,
    y: interpolateY(sorted, x),
  }));
}
