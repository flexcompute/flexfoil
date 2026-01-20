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
 * Used for decomposition where we have many sample points.
 */
function interpolateYLinear(curve: { x: number; y: number }[], targetX: number): number {
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
 * Build natural cubic spline coefficients for interpolation.
 * Given points (x_i, y_i) sorted by x, computes coefficients for each interval.
 * Returns coefficients a, b, c, d where:
 *   S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
 */
interface SplineCoefficients {
  x: number[];  // x values at knots
  a: number[];  // constant coefficients (= y values)
  b: number[];  // linear coefficients
  c: number[];  // quadratic coefficients
  d: number[];  // cubic coefficients
}

function buildCubicSplineCoeffs(points: { x: number; y: number }[]): SplineCoefficients | null {
  const n = points.length;
  if (n < 2) return null;
  
  // Sort by x
  const sorted = [...points].sort((a, b) => a.x - b.x);
  
  const x = sorted.map(p => p.x);
  const y = sorted.map(p => p.y);
  
  if (n === 2) {
    // Linear interpolation for 2 points
    const h = x[1] - x[0];
    const slope = h > 0 ? (y[1] - y[0]) / h : 0;
    return {
      x: [x[0]],
      a: [y[0]],
      b: [slope],
      c: [0],
      d: [0],
    };
  }
  
  const nSeg = n - 1;
  
  // Compute h_i = x_{i+1} - x_i
  const h: number[] = [];
  for (let i = 0; i < nSeg; i++) {
    h.push(x[i + 1] - x[i]);
  }
  
  // Build tridiagonal system for second derivatives (natural spline: c[0] = c[n-1] = 0)
  const c = new Array(n).fill(0);
  
  if (n > 2) {
    const nInterior = n - 2;
    const lower = new Array(nInterior).fill(0);
    const diag = new Array(nInterior).fill(0);
    const upper = new Array(nInterior).fill(0);
    const rhs = new Array(nInterior).fill(0);
    
    for (let i = 0; i < nInterior; i++) {
      const j = i + 1; // Index in original arrays
      lower[i] = h[j - 1];
      diag[i] = 2 * (h[j - 1] + h[j]);
      upper[i] = h[j];
      rhs[i] = 3 * ((y[j + 1] - y[j]) / h[j] - (y[j] - y[j - 1]) / h[j - 1]);
    }
    
    // Thomas algorithm for tridiagonal system
    const cPrime = new Array(nInterior).fill(0);
    const dPrime = new Array(nInterior).fill(0);
    
    cPrime[0] = upper[0] / diag[0];
    dPrime[0] = rhs[0] / diag[0];
    
    for (let i = 1; i < nInterior; i++) {
      const denom = diag[i] - lower[i] * cPrime[i - 1];
      cPrime[i] = i < nInterior - 1 ? upper[i] / denom : 0;
      dPrime[i] = (rhs[i] - lower[i] * dPrime[i - 1]) / denom;
    }
    
    // Back substitution
    const solution = new Array(nInterior).fill(0);
    solution[nInterior - 1] = dPrime[nInterior - 1];
    for (let i = nInterior - 2; i >= 0; i--) {
      solution[i] = dPrime[i] - cPrime[i] * solution[i + 1];
    }
    
    // Copy to c array (interior points)
    for (let i = 0; i < nInterior; i++) {
      c[i + 1] = solution[i];
    }
  }
  
  // Compute b and d coefficients
  const a: number[] = [];
  const b: number[] = [];
  const d: number[] = [];
  
  for (let i = 0; i < nSeg; i++) {
    a.push(y[i]);
    b.push((y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3);
    d.push((c[i + 1] - c[i]) / (3 * h[i]));
  }
  
  return {
    x: x.slice(0, nSeg),
    a,
    b,
    c: c.slice(0, nSeg),
    d,
  };
}

/**
 * Evaluate cubic spline at a given x value using precomputed coefficients.
 */
function evaluateCubicSpline(coeffs: SplineCoefficients, targetX: number): number {
  const { x, a, b, c, d } = coeffs;
  const n = x.length;
  
  if (n === 0) return 0;
  
  // Find the appropriate interval
  let i = 0;
  for (let j = 0; j < n - 1; j++) {
    if (targetX >= x[j] && targetX <= x[j + 1]) {
      i = j;
      break;
    }
    if (j === n - 2) {
      i = j; // Use last segment for extrapolation
    }
  }
  
  // Handle points before first knot
  if (targetX < x[0]) {
    i = 0;
  }
  // Handle points after last knot  
  if (n > 0 && targetX > x[n - 1]) {
    i = n - 1;
  }
  
  const dx = targetX - x[i];
  return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

/**
 * Interpolate y values at given x positions using cubic spline.
 * This provides smooth C2-continuous interpolation through the control points.
 * Exported for use in other modules that need spline interpolation.
 */
export function interpolateYSpline(points: { x: number; y: number }[], targetX: number): number {
  if (points.length === 0) return 0;
  if (points.length === 1) return points[0].y;
  if (points.length === 2) {
    // Linear for 2 points
    return interpolateYLinear(points, targetX);
  }
  
  const coeffs = buildCubicSplineCoeffs(points);
  if (!coeffs) return 0;
  
  return evaluateCubicSpline(coeffs, targetX);
}

/**
 * Legacy alias for linear interpolation (used in decomposition).
 */
const interpolateY = interpolateYLinear;

/**
 * Decompose an airfoil into camber line and thickness distribution.
 * 
 * The camber line is the locus of points midway between upper and lower surfaces.
 * The thickness is the distance from the camber line to the surface.
 * 
 * IMPORTANT: Uses cosine spacing with high resolution (200+ samples) to preserve
 * the critical leading edge geometry. Uniform spacing destroys LE detail because
 * the LE region (x=0 to x=0.02) would only get ~1 sample.
 * 
 * @param coords - Airfoil coordinates (TE -> upper -> LE -> lower -> TE format)
 * @returns Camber line and thickness distribution
 */
export function decomposeToCamberThickness(coords: AirfoilPoint[]): CamberThicknessDecomposition {
  if (coords.length < 5) {
    return { camber: [], thickness: [] };
  }
  
  const { upper, lower } = splitSurfaces(coords);
  
  // Use COSINE spacing with high resolution to preserve LE geometry
  // Cosine spacing concentrates points at x=0 (LE) and x=1 (TE)
  // 200 samples ensures fine resolution even in the nose region
  const nSamples = 200;
  const camber: { x: number; y: number }[] = [];
  const thickness: { x: number; t: number }[] = [];
  
  for (let i = 0; i <= nSamples; i++) {
    // Cosine spacing: x = (1 - cos(beta)) / 2 where beta goes from 0 to PI
    // This gives x from 0 (LE) to 1 (TE) with concentration at both ends
    const beta = (Math.PI * i) / nSamples;
    const x = (1 - Math.cos(beta)) / 2;
    
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
 * Uses linear interpolation - suitable for densely sampled camber/thickness.
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
  // x goes from 1 (TE) to 0 (LE) as i goes from 0 to halfPoints
  for (let i = 0; i <= halfPoints; i++) {
    const beta = (Math.PI * i) / halfPoints;
    const x = (1 + Math.cos(beta)) / 2;  // TE (x=1) to LE (x=0)
    
    // Interpolate camber and thickness at this x (linear for dense samples)
    const yc = interpolateY(camber, x);
    const t = interpolateY(thickness.map(p => ({ x: p.x, y: p.t })), x);
    
    // For proper NACA-style airfoils, we need the slope of the camber line
    // to rotate the thickness perpendicular to the camber
    // For simplicity, we use vertical thickness (good approximation for low camber)
    const yUpper = yc + t;
    const yLower = yc - t;
    
    if (i === 0) {
      // Trailing edge - single point (x=1)
      upper.push({ x: 1, y: yUpper, surface: 'upper' });
    } else if (i === halfPoints) {
      // Leading edge - single point (x=0, thickness is ~0)
      upper.push({ x: 0, y: yc, surface: 'upper' });
    } else {
      upper.push({ x, y: yUpper, surface: 'upper' });
      lower.unshift({ x, y: yLower, surface: 'lower' });
    }
  }
  
  // Build final coordinates: upper (TE to LE) + lower (LE to TE)
  const result: AirfoilPoint[] = [...upper];
  
  // Add lower surface (goes from near-LE to near-TE)
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
 * Reconstruct airfoil coordinates from sparse camber and thickness control points.
 * Uses cubic spline interpolation for smooth C2-continuous curves.
 * 
 * Special handling for leading edge:
 * - Thickness is forced to 0 at x=0 to ensure upper/lower surfaces meet
 * - Near the LE, thickness blends smoothly to maintain roundness
 * 
 * @param camber - Camber control points (x, y) - typically 5-10 points
 * @param thickness - Thickness control points (x, t) - typically 5-10 points  
 * @param nPoints - Number of output points (default 161 for 160 panels)
 * @returns Airfoil coordinates in standard format
 */
export function reconstructFromCamberThicknessSpline(
  camber: { x: number; y: number }[],
  thickness: { x: number; t: number }[],
  nPoints: number = 161
): AirfoilPoint[] {
  if (camber.length < 2 || thickness.length < 2) {
    return [];
  }
  
  // Sort control points by x
  const sortedCamber = [...camber].sort((a, b) => a.x - b.x);
  const sortedThickness = [...thickness].sort((a, b) => a.x - b.x);
  
  // Build splines directly from user's control points
  // No anchor modification - let the user's points define the shape
  const camberCoeffs = buildCubicSplineCoeffs(sortedCamber);
  const thicknessCoeffs = buildCubicSplineCoeffs(
    sortedThickness.map(p => ({ x: p.x, y: p.t }))
  );
  
  if (!camberCoeffs || !thicknessCoeffs) {
    // Fallback to linear if spline fails
    return reconstructFromCamberThickness(camber, thickness, nPoints);
  }
  
  const halfPoints = Math.floor(nPoints / 2);
  const upper: AirfoilPoint[] = [];
  const lower: AirfoilPoint[] = [];
  
  // Pre-calculate thickness at the LE tip boundary for sqrt scaling
  const LE_TIP = 0.03;
  const tAtTip = Math.max(0, evaluateCubicSpline(thicknessCoeffs, LE_TIP));
  
  // Generate points using cosine spacing for better LE resolution
  // x goes from 1 (TE) to 0 (LE) as i goes from 0 to halfPoints
  for (let i = 0; i <= halfPoints; i++) {
    const beta = (Math.PI * i) / halfPoints;
    const x = (1 + Math.cos(beta)) / 2;  // TE (x=1) to LE (x=0)
    
    // Camber: directly from spline (no modification)
    const yc = evaluateCubicSpline(camberCoeffs, x);
    
    // Thickness: from spline, but ensure LE closure with sqrt profile
    // Only the very tip (x < 0.03) uses sqrt scaling for roundness
    let t: number;
    
    if (x < LE_TIP) {
      // Very tip: use sqrt profile scaled to match thickness at LE_TIP
      t = tAtTip * Math.sqrt(x / LE_TIP);
    } else {
      // Rest of airfoil: use spline directly
      t = Math.max(0, evaluateCubicSpline(thicknessCoeffs, x));
    }
    
    const yUpper = yc + t;
    const yLower = yc - t;
    
    if (i === 0) {
      // Trailing edge - single point (x=1)
      upper.push({ x: 1, y: yUpper, surface: 'upper' });
    } else if (i === halfPoints) {
      // Leading edge - single point (x=0, thickness is 0)
      upper.push({ x: 0, y: yc, surface: 'upper' });
    } else {
      upper.push({ x, y: yUpper, surface: 'upper' });
      lower.unshift({ x, y: yLower, surface: 'lower' });
    }
  }
  
  // Build final coordinates: upper (TE to LE) + lower (LE to TE)
  const result: AirfoilPoint[] = [...upper];
  
  // Add lower surface (goes from near-LE to near-TE)
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
 * The algorithm preserves leading edge geometry by:
 * 1. Decomposing with fine cosine-spaced resolution (200+ samples at LE)
 * 2. Applying scaling to the smooth camber/thickness distributions
 * 3. Reconstructing with cosine spacing matching the original point count
 * 
 * For best results after warping, consider repaneling with XFOIL's PANGEN
 * to get proper curvature-based panel distribution.
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
 * Note: The first point (x=0, leading edge) always has t=0 to ensure
 * the upper and lower surfaces meet properly at the LE.
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
      // Simplified NACA thickness formula (note: gives t=0 at x=0)
      const halfT = 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
      return {
        id: `thickness-${i}`,
        x,
        // Force t=0 at LE (x=0) to ensure proper closure
        t: i === 0 ? 0 : Math.max(0, halfT),
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
      // Force t=0 at LE (x=0) to ensure proper closure
      t: i === 0 ? 0 : Math.max(0, t),
    });
  }
  
  return points;
}

/**
 * Reconstruct airfoil from camber and thickness control points.
 * Uses cubic spline interpolation through the control points for smooth curves.
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
  
  // Use spline interpolation for smooth curves through sparse control points
  return reconstructFromCamberThicknessSpline(camber, thickness, nPoints);
}

/**
 * Reconstruct airfoil using original thickness distribution with modified camber.
 * This is the correct approach: the NACA thickness formula is preserved exactly,
 * only the camber line is modified via the control points.
 * 
 * @param camberPoints - Sparse camber control points defining the new camber line
 * @param originalCoords - Original airfoil coordinates (with original thickness)
 * @param nPoints - Number of output points
 * @returns Modified airfoil coordinates
 */
export function reconstructWithOriginalThickness(
  camberPoints: CamberControlPoint[],
  originalCoords: AirfoilPoint[],
  nPoints: number = 161
): AirfoilPoint[] {
  if (camberPoints.length < 2 || originalCoords.length < 5) {
    return originalCoords;
  }
  
  // Decompose original airfoil to get the EXACT thickness distribution
  const { thickness: originalThickness } = decomposeToCamberThickness(originalCoords);
  
  if (originalThickness.length === 0) {
    return originalCoords;
  }
  
  // Build camber spline from control points
  const sortedCamber = [...camberPoints].sort((a, b) => a.x - b.x);
  const camberCoeffs = buildCubicSplineCoeffs(sortedCamber);
  
  if (!camberCoeffs) {
    return originalCoords;
  }
  
  const halfPoints = Math.floor(nPoints / 2);
  const upper: AirfoilPoint[] = [];
  const lower: AirfoilPoint[] = [];
  
  // Generate points using cosine spacing
  for (let i = 0; i <= halfPoints; i++) {
    const beta = (Math.PI * i) / halfPoints;
    const x = (1 + Math.cos(beta)) / 2;  // TE (x=1) to LE (x=0)
    
    // Camber: from user's spline
    const yc = evaluateCubicSpline(camberCoeffs, x);
    
    // Thickness: from ORIGINAL airfoil (preserves exact NACA formula)
    const t = interpolateYLinear(
      originalThickness.map(p => ({ x: p.x, y: p.t })), 
      x
    );
    
    const yUpper = yc + t;
    const yLower = yc - t;
    
    if (i === 0) {
      upper.push({ x: 1, y: yUpper, surface: 'upper' });
    } else if (i === halfPoints) {
      upper.push({ x: 0, y: yc, surface: 'upper' });
    } else {
      upper.push({ x, y: yUpper, surface: 'upper' });
      lower.unshift({ x, y: yLower, surface: 'lower' });
    }
  }
  
  const result: AirfoilPoint[] = [...upper];
  for (let i = 0; i < lower.length; i++) {
    result.push(lower[i]);
  }
  
  if (result.length > 0 && result[result.length - 1].x < 0.99) {
    result.push({ x: 1, y: result[0].y, surface: 'lower' });
  }
  
  return result;
}

/**
 * Reconstruct airfoil using original camber distribution with modified thickness.
 * This preserves the original camber line exactly, only modifying thickness.
 * The nose region (x < 0.1) uses the ORIGINAL thickness to preserve the LE shape.
 * 
 * @param thicknessPoints - Sparse thickness control points defining the new thickness
 * @param originalCoords - Original airfoil coordinates (with original camber)
 * @param nPoints - Number of output points
 * @returns Modified airfoil coordinates
 */
export function reconstructWithOriginalCamber(
  thicknessPoints: ThicknessControlPoint[],
  originalCoords: AirfoilPoint[],
  nPoints: number = 161
): AirfoilPoint[] {
  if (thicknessPoints.length < 2 || originalCoords.length < 5) {
    return originalCoords;
  }
  
  // Decompose original airfoil to get the EXACT camber and thickness distributions
  const { camber: originalCamber, thickness: originalThickness } = decomposeToCamberThickness(originalCoords);
  
  if (originalCamber.length === 0 || originalThickness.length === 0) {
    return originalCoords;
  }
  
  // Build thickness spline from control points
  const sortedThickness = [...thicknessPoints].sort((a, b) => a.x - b.x);
  const thicknessCoeffs = buildCubicSplineCoeffs(
    sortedThickness.map(p => ({ x: p.x, y: p.t }))
  );
  
  if (!thicknessCoeffs) {
    return originalCoords;
  }
  
  const halfPoints = Math.floor(nPoints / 2);
  const upper: AirfoilPoint[] = [];
  const lower: AirfoilPoint[] = [];
  
  // Boundary where we transition from original thickness to user's spline
  const LE_PRESERVE = 0.1;
  
  // Get values at boundary for smooth blending
  const originalTAtBoundary = interpolateYLinear(
    originalThickness.map(p => ({ x: p.x, y: p.t })), 
    LE_PRESERVE
  );
  const splineTAtBoundary = Math.max(0, evaluateCubicSpline(thicknessCoeffs, LE_PRESERVE));
  
  // Generate points using cosine spacing
  for (let i = 0; i <= halfPoints; i++) {
    const beta = (Math.PI * i) / halfPoints;
    const x = (1 + Math.cos(beta)) / 2;  // TE (x=1) to LE (x=0)
    
    // Camber: from ORIGINAL airfoil (preserves exact shape)
    const yc = interpolateYLinear(originalCamber, x);
    
    // Thickness: blend from original (at nose) to user's spline (at body)
    let t: number;
    if (x < LE_PRESERVE) {
      // Nose region: use ORIGINAL thickness (preserves NACA LE shape)
      const originalT = interpolateYLinear(
        originalThickness.map(p => ({ x: p.x, y: p.t })), 
        x
      );
      // Scale by ratio of user's thickness to original at boundary
      // This lets user control overall thickness while preserving nose SHAPE
      const scale = splineTAtBoundary / Math.max(0.001, originalTAtBoundary);
      t = originalT * scale;
    } else {
      // Body: use user's spline directly
      t = Math.max(0, evaluateCubicSpline(thicknessCoeffs, x));
    }
    
    const yUpper = yc + t;
    const yLower = yc - t;
    
    if (i === 0) {
      upper.push({ x: 1, y: yUpper, surface: 'upper' });
    } else if (i === halfPoints) {
      upper.push({ x: 0, y: yc, surface: 'upper' });
    } else {
      upper.push({ x, y: yUpper, surface: 'upper' });
      lower.unshift({ x, y: yLower, surface: 'lower' });
    }
  }
  
  const result: AirfoilPoint[] = [...upper];
  for (let i = 0; i < lower.length; i++) {
    result.push(lower[i]);
  }
  
  if (result.length > 0 && result[result.length - 1].x < 0.99) {
    result.push({ x: 1, y: result[0].y, surface: 'lower' });
  }
  
  return result;
}

/**
 * Evaluate a cubic spline through control points at given x values.
 * Uses natural cubic spline interpolation for C2 continuity.
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
  
  // Build spline coefficients once
  const coeffs = buildCubicSplineCoeffs(sorted);
  if (!coeffs) {
    // Fallback to linear
    return xValues.map(x => ({
      x,
      y: interpolateYLinear(sorted, x),
    }));
  }
  
  return xValues.map(x => ({
    x,
    y: evaluateCubicSpline(coeffs, x),
  }));
}

/**
 * Generate a smooth spline curve through camber control points.
 * Returns densely sampled points for visualization.
 * 
 * @param controlPoints - Sparse camber control points
 * @param nSamples - Number of output samples (default 100)
 * @returns Densely sampled curve points
 */
export function generateCamberSplineCurve(
  controlPoints: CamberControlPoint[],
  nSamples: number = 100
): { x: number; y: number }[] {
  if (controlPoints.length < 2) {
    return controlPoints.map(p => ({ x: p.x, y: p.y }));
  }
  
  const sorted = [...controlPoints].sort((a, b) => a.x - b.x);
  const points = sorted.map(p => ({ x: p.x, y: p.y }));
  
  const coeffs = buildCubicSplineCoeffs(points);
  if (!coeffs) {
    return points;
  }
  
  const result: { x: number; y: number }[] = [];
  const xMin = sorted[0].x;
  const xMax = sorted[sorted.length - 1].x;
  
  for (let i = 0; i < nSamples; i++) {
    const x = xMin + (i / (nSamples - 1)) * (xMax - xMin);
    const y = evaluateCubicSpline(coeffs, x);
    result.push({ x, y });
  }
  
  return result;
}
