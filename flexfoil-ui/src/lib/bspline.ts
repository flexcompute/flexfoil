/**
 * B-Spline evaluation for airfoil shape control.
 * 
 * Implements De Boor's algorithm for evaluating B-spline curves,
 * plus least-squares fitting to find control points that approximate
 * a given set of data points.
 */

import type { BSplineControlPoint, AirfoilPoint } from '../types';

/**
 * Generate a uniform knot vector for a B-spline.
 * For n control points and degree p, we need n + p + 1 knots.
 * We use clamped (open) knot vectors where the curve passes through endpoints.
 */
function generateKnotVector(n: number, degree: number): number[] {
  const m = n + degree + 1;
  const knots: number[] = [];
  
  for (let i = 0; i < m; i++) {
    if (i <= degree) {
      knots.push(0);
    } else if (i >= m - degree - 1) {
      knots.push(1);
    } else {
      knots.push((i - degree) / (n - degree));
    }
  }
  
  return knots;
}

/**
 * Find the knot span index for parameter t.
 * Returns i such that knots[i] <= t < knots[i+1]
 */
function findSpan(n: number, degree: number, t: number, knots: number[]): number {
  // Special case for t at end
  if (t >= knots[n]) {
    return n - 1;
  }
  if (t <= knots[degree]) {
    return degree;
  }
  
  // Binary search
  let low = degree;
  let high = n;
  let mid = Math.floor((low + high) / 2);
  
  while (t < knots[mid] || t >= knots[mid + 1]) {
    if (t < knots[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = Math.floor((low + high) / 2);
  }
  
  return mid;
}

/**
 * Compute B-spline basis functions using Cox-de Boor recursion.
 */
function basisFunctions(span: number, t: number, degree: number, knots: number[]): number[] {
  const N: number[] = new Array(degree + 1).fill(0);
  const left: number[] = new Array(degree + 1).fill(0);
  const right: number[] = new Array(degree + 1).fill(0);
  
  N[0] = 1.0;
  
  for (let j = 1; j <= degree; j++) {
    left[j] = t - knots[span + 1 - j];
    right[j] = knots[span + j] - t;
    let saved = 0.0;
    
    for (let r = 0; r < j; r++) {
      const temp = N[r] / (right[r + 1] + left[j - r]);
      N[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    
    N[j] = saved;
  }
  
  return N;
}

/**
 * Evaluate a B-spline curve at parameter t using De Boor's algorithm.
 */
function evaluatePoint(
  controlPoints: BSplineControlPoint[],
  degree: number,
  t: number,
  knots: number[]
): { x: number; y: number } {
  const n = controlPoints.length;
  const span = findSpan(n, degree, t, knots);
  const N = basisFunctions(span, t, degree, knots);
  
  let x = 0;
  let y = 0;
  
  for (let i = 0; i <= degree; i++) {
    const cpIndex = span - degree + i;
    if (cpIndex >= 0 && cpIndex < n) {
      x += N[i] * controlPoints[cpIndex].x;
      y += N[i] * controlPoints[cpIndex].y;
    }
  }
  
  return { x, y };
}

/**
 * Evaluate a B-spline curve and return sampled points.
 * 
 * @param controlPoints - Array of control points
 * @param degree - Degree of the B-spline (typically 3 for cubic)
 * @param nSamples - Number of output points to generate
 * @returns Array of points along the curve
 */
export function evaluateBSpline(
  controlPoints: BSplineControlPoint[],
  degree: number,
  nSamples: number
): AirfoilPoint[] {
  if (controlPoints.length < 2) {
    return [];
  }
  
  // Clamp degree to valid range
  const p = Math.min(degree, controlPoints.length - 1);
  const n = controlPoints.length;
  
  // Generate knot vector
  const knots = generateKnotVector(n, p);
  
  // Sample the curve
  const points: AirfoilPoint[] = [];
  
  for (let i = 0; i < nSamples; i++) {
    const t = i / (nSamples - 1);
    const pt = evaluatePoint(controlPoints, p, t, knots);
    points.push({ x: pt.x, y: pt.y });
  }
  
  return points;
}

/**
 * Cubic Bezier interpolation between two points with handles.
 */
function cubicBezier(
  p0: { x: number; y: number },
  p1: { x: number; y: number },
  p2: { x: number; y: number },
  p3: { x: number; y: number },
  t: number
): { x: number; y: number } {
  const t2 = t * t;
  const t3 = t2 * t;
  const mt = 1 - t;
  const mt2 = mt * mt;
  const mt3 = mt2 * mt;
  
  return {
    x: mt3 * p0.x + 3 * mt2 * t * p1.x + 3 * mt * t2 * p2.x + t3 * p3.x,
    y: mt3 * p0.y + 3 * mt2 * t * p1.y + 3 * mt * t2 * p2.y + t3 * p3.y,
  };
}

/**
 * Evaluate a Bezier curve defined by control points and handles.
 * This creates a smooth curve that passes through the anchor points
 * with tangent directions defined by the handles.
 */
export function evaluateBezierCurve(
  coordinates: AirfoilPoint[],
  handles: { pointIndex: number; position: { x: number; y: number }; type: 'in' | 'out' }[],
  nSamples: number
): AirfoilPoint[] {
  if (coordinates.length < 2 || handles.length === 0) {
    return coordinates;
  }
  
  // Build a map of handles by point index
  const handleMap = new Map<number, { in?: { x: number; y: number }; out?: { x: number; y: number } }>();
  for (const h of handles) {
    if (!handleMap.has(h.pointIndex)) {
      handleMap.set(h.pointIndex, {});
    }
    const entry = handleMap.get(h.pointIndex)!;
    if (h.type === 'in') {
      entry.in = h.position;
    } else {
      entry.out = h.position;
    }
  }
  
  // Get sorted point indices that have handles
  const handleIndices = Array.from(handleMap.keys()).sort((a, b) => a - b);
  
  if (handleIndices.length < 2) {
    return coordinates;
  }
  
  const result: AirfoilPoint[] = [];
  const samplesPerSegment = Math.max(10, Math.floor(nSamples / (handleIndices.length - 1)));
  
  // Generate curve segments between handle points
  for (let i = 0; i < handleIndices.length - 1; i++) {
    const idx0 = handleIndices[i];
    const idx1 = handleIndices[i + 1];
    
    const p0 = coordinates[idx0];
    const p3 = coordinates[idx1];
    
    // Get control points from handles
    const h0 = handleMap.get(idx0);
    const h1 = handleMap.get(idx1);
    
    // Use handle positions as control points
    // For outgoing handle from p0, use it directly
    // For incoming handle to p3, use it directly
    // If missing, create default handles that maintain tangent continuity
    let p1: { x: number; y: number };
    let p2: { x: number; y: number };
    
    if (h0?.out) {
      p1 = h0.out;
    } else {
      // Default: 1/3 of the way from p0 to p3
      p1 = { x: p0.x + (p3.x - p0.x) / 3, y: p0.y + (p3.y - p0.y) / 3 };
    }
    
    if (h1?.in) {
      p2 = h1.in;
    } else if (h1?.out) {
      // Mirror the outgoing handle to create incoming
      const dx = p3.x - h1.out.x;
      const dy = p3.y - h1.out.y;
      p2 = { x: p3.x + dx, y: p3.y + dy };
    } else {
      // Default: 2/3 of the way from p0 to p3
      p2 = { x: p0.x + 2 * (p3.x - p0.x) / 3, y: p0.y + 2 * (p3.y - p0.y) / 3 };
    }
    
    // Sample this segment
    for (let j = 0; j < samplesPerSegment; j++) {
      const t = j / samplesPerSegment;
      const pt = cubicBezier(p0, p1, p2, p3, t);
      result.push({ x: pt.x, y: pt.y });
    }
  }
  
  // Add final point
  const lastIdx = handleIndices[handleIndices.length - 1];
  result.push({ x: coordinates[lastIdx].x, y: coordinates[lastIdx].y });
  
  return result;
}
