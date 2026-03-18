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
 * Compute all basis functions for all data points.
 * Returns a matrix N where N[i][j] = N_j,p(t_i)
 */
function computeBasisMatrix(
  dataPoints: number,
  numControlPoints: number,
  degree: number,
  knots: number[]
): number[][] {
  const N: number[][] = [];
  
  for (let i = 0; i < dataPoints; i++) {
    const t = i / (dataPoints - 1);
    const row: number[] = new Array(numControlPoints).fill(0);
    
    // Find span and compute basis functions
    const span = findSpan(numControlPoints, degree, t, knots);
    const basis = basisFunctions(span, t, degree, knots);
    
    // Fill in the non-zero basis functions
    for (let j = 0; j <= degree; j++) {
      const cpIndex = span - degree + j;
      if (cpIndex >= 0 && cpIndex < numControlPoints) {
        row[cpIndex] = basis[j];
      }
    }
    
    N.push(row);
  }
  
  return N;
}

/**
 * Solve least squares: find control points P such that N*P ≈ D
 * where D is the data points matrix and N is the basis function matrix.
 * 
 * Uses normal equations: (N^T * N) * P = N^T * D
 */
function solveLeastSquares(
  N: number[][],
  dataX: number[],
  dataY: number[]
): { x: number[]; y: number[] } {
  const m = N.length;      // number of data points
  const n = N[0].length;   // number of control points
  
  // Compute N^T * N (n x n matrix)
  const NtN: number[][] = [];
  for (let i = 0; i < n; i++) {
    NtN.push(new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < m; k++) {
        sum += N[k][i] * N[k][j];
      }
      NtN[i][j] = sum;
    }
  }
  
  // Compute N^T * D for x and y
  const NtDx: number[] = new Array(n).fill(0);
  const NtDy: number[] = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < m; k++) {
      NtDx[i] += N[k][i] * dataX[k];
      NtDy[i] += N[k][i] * dataY[k];
    }
  }
  
  // Solve using Gaussian elimination with partial pivoting
  const solveSystem = (A: number[][], b: number[]): number[] => {
    const size = A.length;
    // Create augmented matrix
    const aug: number[][] = A.map((row, i) => [...row, b[i]]);
    
    // Forward elimination with partial pivoting
    for (let col = 0; col < size; col++) {
      // Find pivot
      let maxRow = col;
      for (let row = col + 1; row < size; row++) {
        if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) {
          maxRow = row;
        }
      }
      // Swap rows
      [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
      
      // Check for singular matrix
      if (Math.abs(aug[col][col]) < 1e-12) {
        // Add small regularization
        aug[col][col] = 1e-10;
      }
      
      // Eliminate column
      for (let row = col + 1; row < size; row++) {
        const factor = aug[row][col] / aug[col][col];
        for (let j = col; j <= size; j++) {
          aug[row][j] -= factor * aug[col][j];
        }
      }
    }
    
    // Back substitution
    const x: number[] = new Array(size).fill(0);
    for (let i = size - 1; i >= 0; i--) {
      let sum = aug[i][size];
      for (let j = i + 1; j < size; j++) {
        sum -= aug[i][j] * x[j];
      }
      x[i] = sum / aug[i][i];
    }
    
    return x;
  };
  
  // Add regularization to NtN for numerical stability
  for (let i = 0; i < n; i++) {
    NtN[i][i] += 1e-8;
  }
  
  const cpX = solveSystem(NtN.map(row => [...row]), [...NtDx]);
  const cpY = solveSystem(NtN.map(row => [...row]), [...NtDy]);
  
  return { x: cpX, y: cpY };
}

/**
 * Create B-spline control points that approximate the given airfoil shape.
 * Uses least-squares fitting to find control points such that the resulting
 * B-spline curve closely matches the original airfoil.
 */
function createBSplineControlPointsFromAirfoil(
  coordinates: AirfoilPoint[],
  numControlPoints: number = 12,
  degree: number = 3
): BSplineControlPoint[] {
  if (coordinates.length < 3) {
    return [];
  }
  
  // Clamp degree
  const p = Math.min(degree, numControlPoints - 1);
  
  // Generate knot vector
  const knots = generateKnotVector(numControlPoints, p);
  
  // Extract data points
  const dataX = coordinates.map(pt => pt.x);
  const dataY = coordinates.map(pt => pt.y);
  
  // Compute basis function matrix
  const N = computeBasisMatrix(coordinates.length, numControlPoints, p, knots);
  
  // Solve least squares to find control points
  const { x: cpX, y: cpY } = solveLeastSquares(N, dataX, dataY);
  
  // Create control point objects
  const controlPoints: BSplineControlPoint[] = [];
  for (let i = 0; i < numControlPoints; i++) {
    controlPoints.push({
      id: `cp-${i}`,
      x: cpX[i],
      y: cpY[i],
      weight: 1,
    });
  }
  
  return controlPoints;
}

/**
 * Create initial Bezier handles from airfoil coordinates.
 * Creates handles at key points along the airfoil that preserve
 * the original curve shape by computing proper tangent directions
 * and handle lengths based on the local curvature.
 */
function createBezierHandlesFromAirfoil(
  coordinates: AirfoilPoint[],
  numHandles: number = 8
): { pointIndex: number; position: { x: number; y: number }; type: 'in' | 'out' }[] {
  if (coordinates.length < 3) {
    return [];
  }
  
  const handles: { pointIndex: number; position: { x: number; y: number }; type: 'in' | 'out' }[] = [];
  
  // Select key points along the airfoil
  for (let i = 0; i < numHandles; i++) {
    const t = i / (numHandles - 1);
    const pointIndex = Math.floor(t * (coordinates.length - 1));
    const p = coordinates[pointIndex];
    
    // Get neighboring points to estimate tangent (use wider stencil for smoother tangent)
    const prevIdx = Math.max(0, pointIndex - 2);
    const nextIdx = Math.min(coordinates.length - 1, pointIndex + 2);
    const prev = coordinates[prevIdx];
    const next = coordinates[nextIdx];
    
    // Tangent direction
    const dx = next.x - prev.x;
    const dy = next.y - prev.y;
    const len = Math.sqrt(dx * dx + dy * dy);
    
    if (len > 1e-6) {
      // Handle length proportional to distance to next handle point
      // This helps preserve the curve shape
      const segmentLen = len / 2;
      const handleLen = segmentLen * 0.4; // ~40% of segment for smooth Bezier
      
      const nx = dx / len;
      const ny = dy / len;
      
      // Create outgoing handle (in direction of increasing parameter)
      handles.push({
        pointIndex,
        position: {
          x: p.x + nx * handleLen,
          y: p.y + ny * handleLen,
        },
        type: 'out',
      });
      
      // Create incoming handle (opposite direction) for interior points
      if (i > 0 && i < numHandles - 1) {
        handles.push({
          pointIndex,
          position: {
            x: p.x - nx * handleLen,
            y: p.y - ny * handleLen,
          },
          type: 'in',
        });
      }
    }
  }
  
  return handles;
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
