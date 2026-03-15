/**
 * AirfoilCanvas - Canvas 2D renderer for airfoil visualization
 * 
 * Renders:
 * - Airfoil curve (smooth spline)
 * - Panel lines (linear connections between panel points)
 * - Panel points
 * - Control points (surface, bezier handles, b-spline)
 * - Grid and axes
 * 
 * Supports:
 * - Pan and zoom
 * - Point dragging
 * - Point selection
 * - Show/hide toggles for curve, panels, points, and grid
 */

import { useRef, useEffect, useCallback, useState, useMemo } from 'react';
import { useShallow } from 'zustand/react/shallow';
import { useAirfoilStore, pauseHistory, resumeHistory } from '../stores/airfoilStore';
import { useVisualizationStore } from '../stores/visualizationStore';
import { useTheme } from '../contexts/ThemeContext';
import { useLayout } from '../contexts/LayoutContext';
import type { Point, ViewportState, AirfoilPoint } from '../types';
import { computeStreamlines, computePsiGrid, createSmokeSystem, isWasmReady, WasmSmokeSystem, analyzeAirfoil, computeGamma, getBLVisualizationData, type BLVisualizationData } from '../lib/wasm';
import { useMorphingAnimation, getCpColor, computeForceVectors } from '../hooks/useMorphingAnimation';
import { generateCamberSplineCurve } from '../lib/airfoilGeometry';
import { syncToUrl, getCompleteUrlState } from '../lib/urlState';
import { WebGPURenderer, checkWebGPUSupport } from '../lib/webgpu';
// Removed d3-contour - using efficient cell-based rendering instead

/**
 * Compute cumulative arc lengths for a sequence of points.
 */
function computeArcLengths(points: AirfoilPoint[]): number[] {
  const s: number[] = [0];
  for (let i = 1; i < points.length; i++) {
    const dx = points[i].x - points[i - 1].x;
    const dy = points[i].y - points[i - 1].y;
    s.push(s[i - 1] + Math.sqrt(dx * dx + dy * dy));
  }
  return s;
}

/**
 * Build natural cubic spline coefficients.
 * Returns coefficients for evaluating the spline.
 */
function buildCubicSpline(s: number[], f: number[]): { a: number[]; b: number[]; c: number[]; d: number[] } {
  const n = s.length;
  if (n < 2) {
    return { a: [], b: [], c: [], d: [] };
  }
  
  if (n === 2) {
    const h = s[1] - s[0];
    return {
      a: [f[0]],
      b: [h > 0 ? (f[1] - f[0]) / h : 0],
      c: [0],
      d: [0],
    };
  }
  
  const nSeg = n - 1;
  const h: number[] = [];
  for (let i = 0; i < nSeg; i++) {
    h.push(s[i + 1] - s[i]);
  }
  
  // Build tridiagonal system for second derivatives
  const c = new Array(n).fill(0);
  
  if (n > 2) {
    const nInterior = n - 2;
    const lower = new Array(nInterior).fill(0);
    const diag = new Array(nInterior).fill(0);
    const upper = new Array(nInterior).fill(0);
    const rhs = new Array(nInterior).fill(0);
    
    for (let i = 0; i < nInterior; i++) {
      const j = i + 1;
      lower[i] = h[j - 1];
      diag[i] = 2 * (h[j - 1] + h[j]);
      upper[i] = h[j];
      rhs[i] = 3 * ((f[j + 1] - f[j]) / h[j] - (f[j] - f[j - 1]) / h[j - 1]);
    }
    
    // Thomas algorithm
    const cPrime = new Array(nInterior).fill(0);
    const dPrime = new Array(nInterior).fill(0);
    
    cPrime[0] = upper[0] / diag[0];
    dPrime[0] = rhs[0] / diag[0];
    
    for (let i = 1; i < nInterior; i++) {
      const denom = diag[i] - lower[i] * cPrime[i - 1];
      cPrime[i] = i < nInterior - 1 ? upper[i] / denom : 0;
      dPrime[i] = (rhs[i] - lower[i] * dPrime[i - 1]) / denom;
    }
    
    const x = new Array(nInterior).fill(0);
    x[nInterior - 1] = dPrime[nInterior - 1];
    for (let i = nInterior - 2; i >= 0; i--) {
      x[i] = dPrime[i] - cPrime[i] * x[i + 1];
    }
    
    for (let i = 0; i < nInterior; i++) {
      c[i + 1] = x[i];
    }
  }
  
  // Compute remaining coefficients
  const a: number[] = [];
  const b: number[] = [];
  const d: number[] = [];
  
  for (let i = 0; i < nSeg; i++) {
    a.push(f[i]);
    b.push((f[i + 1] - f[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3);
    d.push((c[i + 1] - c[i]) / (3 * h[i]));
  }
  
  return { a, b, c: c.slice(0, nSeg), d };
}

/**
 * Evaluate spline at parameter value.
 */
function evaluateSpline(
  sValues: number[],
  coeffs: { a: number[]; b: number[]; c: number[]; d: number[] },
  t: number
): number {
  // Find segment
  let i = 0;
  for (let j = 0; j < sValues.length - 1; j++) {
    if (t >= sValues[j] && t <= sValues[j + 1]) {
      i = j;
      break;
    }
    if (j === sValues.length - 2) {
      i = j;
    }
  }
  
  const dt = t - sValues[i];
  return coeffs.a[i] + coeffs.b[i] * dt + coeffs.c[i] * dt * dt + coeffs.d[i] * dt * dt * dt;
}

/**
 * Generate smooth spline curve points from discrete airfoil coordinates.
 */
function generateSplineCurve(coordinates: AirfoilPoint[], numPoints: number = 200): AirfoilPoint[] {
  if (coordinates.length < 3) {
    return coordinates;
  }
  
  const s = computeArcLengths(coordinates);
  const sMax = s[s.length - 1];
  
  if (sMax < 1e-10) {
    return coordinates;
  }
  
  const xVals = coordinates.map(p => p.x);
  const yVals = coordinates.map(p => p.y);
  
  const xCoeffs = buildCubicSpline(s, xVals);
  const yCoeffs = buildCubicSpline(s, yVals);
  
  const result: AirfoilPoint[] = [];
  for (let i = 0; i < numPoints; i++) {
    const t = (i / (numPoints - 1)) * sMax;
    result.push({
      x: evaluateSpline(s, xCoeffs, t),
      y: evaluateSpline(s, yCoeffs, t),
    });
  }
  
  return result;
}

/**
 * Marching squares algorithm to extract iso-contour line segments from a 2D grid.
 * Returns an array of line segments, each defined by two [x, y] points.
 */
function marchingSquares(
  grid: number[],
  nx: number,
  ny: number,
  level: number,
  x0: number,
  y0: number,
  dx: number,
  dy: number
): [[number, number], [number, number]][] {
  const segments: [[number, number], [number, number]][] = [];
  
  // Lookup table for marching squares cases (16 cases)
  // Each case maps to a list of edge pairs to connect
  // Edges: 0=bottom, 1=right, 2=top, 3=left
  
  const getValue = (ix: number, iy: number): number => {
    const idx = iy * nx + ix;
    return grid[idx];
  };
  
  const lerp = (v0: number, v1: number, t: number): number => v0 + t * (v1 - v0);
  
  const getEdgePoint = (edge: number, x: number, y: number, v: number[]): [number, number] => {
    // v = [bottom-left, bottom-right, top-right, top-left] values
    switch (edge) {
      case 0: { // bottom edge
        const t = (level - v[0]) / (v[1] - v[0]);
        return [lerp(x, x + dx, t), y];
      }
      case 1: { // right edge
        const t = (level - v[1]) / (v[2] - v[1]);
        return [x + dx, lerp(y, y + dy, t)];
      }
      case 2: { // top edge
        const t = (level - v[3]) / (v[2] - v[3]);
        return [lerp(x, x + dx, t), y + dy];
      }
      case 3: { // left edge
        const t = (level - v[0]) / (v[3] - v[0]);
        return [x, lerp(y, y + dy, t)];
      }
      default:
        return [x, y];
    }
  };
  
  // Process each cell
  for (let iy = 0; iy < ny - 1; iy++) {
    for (let ix = 0; ix < nx - 1; ix++) {
      // Get corner values (bottom-left, bottom-right, top-right, top-left)
      const v0 = getValue(ix, iy);
      const v1 = getValue(ix + 1, iy);
      const v2 = getValue(ix + 1, iy + 1);
      const v3 = getValue(ix, iy + 1);
      const v = [v0, v1, v2, v3];
      
      // Skip cells with NaN values (inside airfoil)
      if (!isFinite(v0) || !isFinite(v1) || !isFinite(v2) || !isFinite(v3)) {
        continue;
      }
      
      // Compute case index (4-bit binary: each bit is 1 if corner >= level)
      let caseIndex = 0;
      if (v0 >= level) caseIndex |= 1;
      if (v1 >= level) caseIndex |= 2;
      if (v2 >= level) caseIndex |= 4;
      if (v3 >= level) caseIndex |= 8;
      
      // Cell position
      const cellX = x0 + ix * dx;
      const cellY = y0 + iy * dy;
      
      // Generate segments based on case
      // Edges to connect for each case (saddle cases use midpoint heuristic)
      // Lookup table: 16 cases, each maps to array of edge pairs
      const edgePairs: [number, number][][] = [
        [],                   // 0: all below
        [[0, 3]],             // 1
        [[0, 1]],             // 2
        [[1, 3]],             // 3
        [[1, 2]],             // 4
        [[0, 1], [2, 3]],     // 5 (saddle - two segments)
        [[0, 2]],             // 6
        [[2, 3]],             // 7
        [[2, 3]],             // 8
        [[0, 2]],             // 9
        [[0, 3], [1, 2]],     // 10 (saddle - two segments)
        [[1, 2]],             // 11
        [[1, 3]],             // 12
        [[0, 1]],             // 13
        [[0, 3]],             // 14
        [],                   // 15: all above
      ];
      
      for (const pair of edgePairs[caseIndex]) {
        const e1 = pair[0];
        const e2 = pair[1];
        const p1 = getEdgePoint(e1, cellX, cellY, v);
        const p2 = getEdgePoint(e2, cellX, cellY, v);
        segments.push([p1, p2]);
      }
    }
  }
  
  return segments;
}

/**
 * Connect contour segments into polylines.
 * Segments that share endpoints are connected into longer polylines.
 */
function connectSegments(
  segments: [[number, number], [number, number]][]
): [number, number][][] {
  if (segments.length === 0) return [];
  
  const tolerance = 1e-8;
  const pointsClose = (a: [number, number], b: [number, number]): boolean => {
    return Math.abs(a[0] - b[0]) < tolerance && Math.abs(a[1] - b[1]) < tolerance;
  };
  
  // Track which segments have been used
  const used = new Array(segments.length).fill(false);
  const polylines: [number, number][][] = [];
  
  for (let i = 0; i < segments.length; i++) {
    if (used[i]) continue;
    
    // Start a new polyline
    const polyline: [number, number][] = [segments[i][0], segments[i][1]];
    used[i] = true;
    
    // Try to extend in both directions
    let extended = true;
    while (extended) {
      extended = false;
      const head = polyline[0];
      const tail = polyline[polyline.length - 1];
      
      for (let j = 0; j < segments.length; j++) {
        if (used[j]) continue;
        
        const [p0, p1] = segments[j];
        
        // Try to connect to tail
        if (pointsClose(tail, p0)) {
          polyline.push(p1);
          used[j] = true;
          extended = true;
        } else if (pointsClose(tail, p1)) {
          polyline.push(p0);
          used[j] = true;
          extended = true;
        }
        // Try to connect to head
        else if (pointsClose(head, p0)) {
          polyline.unshift(p1);
          used[j] = true;
          extended = true;
        } else if (pointsClose(head, p1)) {
          polyline.unshift(p0);
          used[j] = true;
          extended = true;
        }
      }
    }
    
    polylines.push(polyline);
  }
  
  return polylines;
}

/**
 * Clean up the dividing streamline (ψ = ψ₀) near the airfoil surface.
 *
 * Marching squares leaves the streamline terminating about one cell off the body.
 * We keep only the main exterior branch and snap any endpoint already near the
 * airfoil to the exact surface projection instead of free-form extrapolating it.
 */
function extrapolateDividingStreamline(
  psi0Lines: [number, number][][],
  airfoilPoints: { x: number; y: number }[]
): [number, number][][] {
  if (psi0Lines.length === 0 || airfoilPoints.length < 3) {
    return psi0Lines;
  }

  const distance = (a: [number, number], b: [number, number]): number => {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    return Math.sqrt(dx * dx + dy * dy);
  };

  const polylineLength = (line: [number, number][]): number => {
    let total = 0;
    for (let i = 1; i < line.length; i++) {
      total += distance(line[i - 1], line[i]);
    }
    return total;
  };

  // Find the closest projected point on the airfoil boundary.
  const closestPointOnAirfoil = (px: number, py: number): { x: number; y: number; dist: number } => {
    let minDist = Infinity;
    let closest = { x: airfoilPoints[0].x, y: airfoilPoints[0].y };

    for (let i = 0; i < airfoilPoints.length; i++) {
      const p1 = airfoilPoints[i];
      const p2 = airfoilPoints[(i + 1) % airfoilPoints.length];

      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      const lenSq = dx * dx + dy * dy;

      if (lenSq < 1e-10) continue;

      const t = Math.max(0, Math.min(1, ((px - p1.x) * dx + (py - p1.y) * dy) / lenSq));
      const projX = p1.x + t * dx;
      const projY = p1.y + t * dy;

      const dist = Math.sqrt((px - projX) ** 2 + (py - projY) ** 2);
      if (dist < minDist) {
        minDist = dist;
        closest = { x: projX, y: projY };
      }
    }

    return { ...closest, dist: minDist };
  };

  // The physical dividing streamline should be the dominant exterior branch.
  const sortedLines = [...psi0Lines].sort((a, b) => polylineLength(b) - polylineLength(a));
  const mainLine = sortedLines[0];
  if (!mainLine || mainLine.length < 2) {
    return psi0Lines;
  }

  const result: [number, number][] = [...mainLine];
  const start = result[0];
  const end = result[result.length - 1];
  const startSurface = closestPointOnAirfoil(start[0], start[1]);
  const endSurface = closestPointOnAirfoil(end[0], end[1]);

  // Use local contour spacing to decide whether an endpoint should snap to the body.
  const startSpacing = result.length > 1 ? distance(result[0], result[1]) : 0;
  const endSpacing = result.length > 1 ? distance(result[result.length - 2], result[result.length - 1]) : 0;
  const snapThreshold = Math.min(0.08, Math.max(0.015, 4 * Math.max(startSpacing, endSpacing)));

  if (startSurface.dist <= snapThreshold) {
    const snappedStart: [number, number] = [startSurface.x, startSurface.y];
    if (distance(snappedStart, start) > 1e-8) {
      result.unshift(snappedStart);
    } else {
      result[0] = snappedStart;
    }
  }

  if (endSurface.dist <= snapThreshold) {
    const snappedEnd: [number, number] = [endSurface.x, endSurface.y];
    if (distance(snappedEnd, end) > 1e-8) {
      result.push(snappedEnd);
    } else {
      result[result.length - 1] = snappedEnd;
    }
  }

  return [result];
}

// Constants
const PANEL_POINT_RADIUS = 2.5;
const CONTROL_RADIUS = 6;
const HIT_RADIUS = 10;

interface AirfoilCanvasProps {
  initialViewport?: { centerX: number; centerY: number; zoom: number } | null;
}

export function AirfoilCanvas({ initialViewport }: AirfoilCanvasProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const smokeCanvasRef = useRef<HTMLCanvasElement>(null);  // Overlay canvas for smoke (redraws every frame)
  const gpuCanvasRef = useRef<HTMLCanvasElement>(null);    // WebGPU canvas
  const containerRef = useRef<HTMLDivElement>(null);
  
  // WebGPU renderer
  const gpuRendererRef = useRef<WebGPURenderer | null>(null);
  const gpuAnimationRef = useRef<number | null>(null);
  
  // Theme context
  const { isDark } = useTheme();
  
  // Layout context for panel actions
  const { openPanel } = useLayout();
  
  // State from store (use shallow equality to prevent unnecessary re-renders)
  const { 
    coordinates, 
    panels, 
    controlMode,
    displayAlpha,
    reynolds,
    mach,
    ncrit,
    maxIterations,
    solverMode,
    // Camber/thickness control
    camberControlPoints,
    thicknessControlPoints,
    updateCamberControlPoint,
    updateThicknessControlPoint,
    addCamberControlPoint,
    removeCamberControlPoint,
    addThicknessControlPoint,
    removeThicknessControlPoint,
  } = useAirfoilStore(
    useShallow((state) => ({
      coordinates: state.coordinates,
      panels: state.panels,
      controlMode: state.controlMode,
      displayAlpha: state.displayAlpha,
      reynolds: state.reynolds,
      mach: state.mach,
      ncrit: state.ncrit,
      maxIterations: state.maxIterations,
      solverMode: state.solverMode,
      camberControlPoints: state.camberControlPoints,
      thicknessControlPoints: state.thicknessControlPoints,
      updateCamberControlPoint: state.updateCamberControlPoint,
      updateThicknessControlPoint: state.updateThicknessControlPoint,
      addCamberControlPoint: state.addCamberControlPoint,
      removeCamberControlPoint: state.removeCamberControlPoint,
      addThicknessControlPoint: state.addThicknessControlPoint,
      removeThicknessControlPoint: state.removeThicknessControlPoint,
    }))
  );

  // Viewport state (initialize from URL if provided)
  const [viewport, setViewport] = useState<ViewportState>(() => ({
    center: { 
      x: initialViewport?.centerX ?? 0.5, 
      y: initialViewport?.centerY ?? 0 
    },
    zoom: initialViewport?.zoom ?? 400,
    width: 800,
    height: 600,
  }));
  
  // Get airfoil state for URL sync (use shallow equality)
  const airfoilState = useAirfoilStore(
    useShallow((state) => ({
      name: state.name,
      nPanels: state.nPanels,
      spacingKnots: state.spacingKnots,
      controlMode: state.controlMode,
      displayAlpha: state.displayAlpha,
      thicknessScale: state.thicknessScale,
      camberScale: state.camberScale,
      curvatureWeight: state.curvatureWeight,
      spacingPanelMode: state.spacingPanelMode,
      sspInterpolation: state.sspInterpolation,
      sspVisualization: state.sspVisualization,
    }))
  );

  // Visualization settings from store (use shallow equality to prevent infinite loops)
  const visualizationState = useVisualizationStore(
    useShallow((state) => ({
      showGrid: state.showGrid,
      showCurve: state.showCurve,
      showPanels: state.showPanels,
      showPoints: state.showPoints,
      showControls: state.showControls,
      showStreamlines: state.showStreamlines,
      showSmoke: state.showSmoke,
      showPsiContours: state.showPsiContours,
      showCp: state.showCp,
      showForces: state.showForces,
      showBoundaryLayer: state.showBoundaryLayer,
      showWake: state.showWake,
      showDisplacementThickness: state.showDisplacementThickness,
      blThicknessScale: state.blThicknessScale,
      enableMorphing: state.enableMorphing,
      morphDuration: state.morphDuration,
      streamlineDensity: state.streamlineDensity,
      adaptiveStreamlines: state.adaptiveStreamlines,
      smokeDensity: state.smokeDensity,
      smokeParticlesPerBlob: state.smokeParticlesPerBlob,
      smokeWaveSpacing: state.smokeWaveSpacing,
      smokeResetCounter: state.smokeResetCounter,
      flowSpeed: state.flowSpeed,
      cpDisplayMode: state.cpDisplayMode,
      cpBarScale: state.cpBarScale,
      forceScale: state.forceScale,
      useGPU: state.useGPU,
      gpuAvailable: state.gpuAvailable,
    }))
  );
  
  const {
    showGrid,
    showCurve,
    showPanels,
    showPoints,
    showControls,
    showStreamlines,
    showSmoke,
    showPsiContours,
    showCp,
    showForces,
    showBoundaryLayer,
    showWake,
    showDisplacementThickness,
    blThicknessScale,
    enableMorphing,
    morphDuration,
    streamlineDensity,
    adaptiveStreamlines,
    smokeDensity,
    smokeParticlesPerBlob,
    smokeWaveSpacing,
    smokeResetCounter,
    flowSpeed,
    cpDisplayMode,
    cpBarScale,
    forceScale,
    useGPU,
    gpuAvailable: _gpuAvailable, // Used for conditional WebGPU setup
  } = visualizationState;
  void _gpuAvailable; // Suppress unused warning (used implicitly in GPU init)
  
  // Get store actions separately (not needed for URL sync)
  const setSmokeDensity = useVisualizationStore((state) => state.setSmokeDensity);
  const setGPUAvailable = useVisualizationStore((state) => state.setGPUAvailable);
  const updatePerfMetrics = useVisualizationStore((state) => state.updatePerfMetrics);
  
  // Sync all state to URL (debounced)
  useEffect(() => {
    const timeoutId = setTimeout(() => {
      const urlState = getCompleteUrlState(
        airfoilState,
        visualizationState,
        {
          centerX: viewport.center.x,
          centerY: viewport.center.y,
          zoom: viewport.zoom,
        }
      );
      syncToUrl(urlState);
    }, 500); // Debounce 500ms
    
    return () => clearTimeout(timeoutId);
  }, [airfoilState, visualizationState, viewport]);
  
  // Streamlines cache
  const [streamlines, setStreamlines] = useState<[number, number][][]>([]);
  
  // Psi contour cache - stores grid data and extracted contour lines for ψ values
  // NOTE: The dividing streamline (ψ = ψ₀) is extrapolated to intersect the airfoil surface
  // by extending the line segments until they meet the airfoil boundary. This provides a
  // complete visualization of the stagnation streamline from upstream to the body.
  const [psiContours, setPsiContours] = useState<{
    grid: number[];           // Grid values for filled contours
    bounds: [number, number, number, number];  // [xMin, xMax, yMin, yMax]
    nx: number;
    ny: number;
    psiMin: number;
    psiMax: number;
    psi0: number;             // Dividing streamline value
    lines: [number, number][][];  // All contour lines
    psi0Lines: [number, number][][];  // Dividing streamline (ψ = ψ₀), extrapolated to airfoil
  }>({ grid: [], bounds: [0, 0, 0, 0], nx: 0, ny: 0, psiMin: 0, psiMax: 0, psi0: 0, lines: [], psi0Lines: [] });
  
  // Analysis results (Cp, Cl, Cd, Cm)
  const [analysisResult, setAnalysisResult] = useState<{
    cp: number[];
    cpX: number[];
    cl: number;
    cd: number;
    cm: number;
  }>({ cp: [], cpX: [], cl: 0, cd: 0, cm: 0 });
  
  // Stable viewport for adaptive streamlines - only updates after zooming settles
  // This prevents expensive streamline recalculation on every scroll tick
  const lastZoomRef = useRef(viewport.zoom);
  const [stableViewport, setStableViewport] = useState(viewport);
  const zoomSettleTimerRef = useRef<number | null>(null);
  
  // Update stable viewport after zooming settles (debounced)
  useEffect(() => {
    // Clear any pending timer
    if (zoomSettleTimerRef.current) {
      clearTimeout(zoomSettleTimerRef.current);
    }
    
    const zoomDelta = Math.abs(viewport.zoom - lastZoomRef.current) / lastZoomRef.current;
    
    // Immediate update if zoom changed significantly (>20%)
    if (zoomDelta > 0.2) {
      lastZoomRef.current = viewport.zoom;
      setStableViewport(viewport);
    } else if (adaptiveStreamlines) {
      // For adaptive mode, also update after a settling period (2s of no changes)
      // This allows streamlines to repopulate after fine zooming
      zoomSettleTimerRef.current = window.setTimeout(() => {
        lastZoomRef.current = viewport.zoom;
        setStableViewport(viewport);
      }, 2000);
    }
    
    return () => {
      if (zoomSettleTimerRef.current) {
        clearTimeout(zoomSettleTimerRef.current);
      }
    };
  }, [viewport, adaptiveStreamlines]);
  
  // Smoke system - use refs instead of state to avoid triggering full canvas redraws
  const smokeSystemRef = useRef<WasmSmokeSystem | null>(null);
  const smokeAnimationRef = useRef<number | null>(null);
  const smokePositionsRef = useRef<Float64Array | null>(null);
  const smokeAlphasRef = useRef<Float64Array | null>(null);
  const smokePsiValuesRef = useRef<Float64Array | null>(null);
  const smokePsi0Ref = useRef<number>(0);
  const drawSmokeRef = useRef<(() => void) | null>(null);  // Ref to drawSmoke for animation loop
  
  // Compute adaptive streamline count based on stable zoom (not live zoom)
  // This avoids recreating the function on every scroll tick
  const adaptiveStreamlineCount = useMemo(() => {
    if (!adaptiveStreamlines) {
      return streamlineDensity;
    }
    // Increase count when zoomed in, decrease when zoomed out
    const zoomFactor = Math.log10(stableViewport.zoom / 100 + 1);
    return Math.min(100, Math.floor(streamlineDensity * (1 + zoomFactor * 0.5)));
  }, [adaptiveStreamlines, streamlineDensity, stableViewport.zoom]);
  
  // Track if we're actively dragging (for disabling expensive computations during drag)
  const [isDraggingPoint, setIsDraggingPoint] = useState(false);
  
  // Compute bounds for streamline domain
  // For adaptive mode: use viewport-based bounds with padding
  // For fixed mode: use large fixed bounds
  const streamlineBounds = useMemo((): [number, number, number, number] => {
    if (!adaptiveStreamlines) {
      // Fixed mode: large domain that covers most reasonable viewing areas
      return [-2.0, 4.0, -2.0, 2.0];
    }
    
    // Adaptive mode: compute bounds from stable viewport with generous padding
    // This ensures streamlines extend beyond the visible area
    const { center, zoom, width, height } = stableViewport;
    const halfWidth = (width / zoom) * 1.5;  // 50% padding on each side
    const halfHeight = (height / zoom) * 1.5;
    
    // Ensure minimum bounds for streamline quality
    const minHalfWidth = 1.5;
    const minHalfHeight = 1.0;
    
    const displayMinX = center.x - Math.max(halfWidth, minHalfWidth);
    const displayMaxX = center.x + Math.max(halfWidth, minHalfWidth);
    const displayMinY = center.y - Math.max(halfHeight, minHalfHeight);
    const displayMaxY = center.y + Math.max(halfHeight, minHalfHeight);

    // The solver operates in body coordinates, but the viewport bounds are in
    // display coordinates after the airfoil has been rotated by -alpha for
    // presentation. Rotate the visible rectangle back into body coordinates so
    // the sampled psi/streamline domain matches what is actually on screen.
    const cx = 0.25;
    const cy = 0;
    const rad = displayAlpha * Math.PI / 180;
    const cosA = Math.cos(rad);
    const sinA = Math.sin(rad);
    const rotateDisplayToBody = (x: number, y: number) => {
      const dx = x - cx;
      const dy = y - cy;
      return {
        x: cx + dx * cosA - dy * sinA,
        y: cy + dx * sinA + dy * cosA,
      };
    };

    const corners = [
      rotateDisplayToBody(displayMinX, displayMinY),
      rotateDisplayToBody(displayMinX, displayMaxY),
      rotateDisplayToBody(displayMaxX, displayMinY),
      rotateDisplayToBody(displayMaxX, displayMaxY),
    ];

    return [
      Math.min(...corners.map(p => p.x)),
      Math.max(...corners.map(p => p.x)),
      Math.min(...corners.map(p => p.y)),
      Math.max(...corners.map(p => p.y)),
    ];
  }, [adaptiveStreamlines, stableViewport, displayAlpha]);

  // Compute streamlines when enabled and alpha/panels change
  // SKIP during drag to prevent freezing
  // For adaptive mode: recomputes when viewport settles (via stableViewport -> streamlineBounds)
  // For fixed mode: uses large fixed bounds that don't change with viewport
  useEffect(() => {
    if (!showStreamlines || !isWasmReady() || panels.length < 10 || isDraggingPoint) {
      if (!isDraggingPoint) {
        setStreamlines([]);
      }
      return;
    }
    
    try {
      const result = computeStreamlines(panels, displayAlpha, reynolds, adaptiveStreamlineCount, streamlineBounds, mach, ncrit, maxIterations, solverMode);
      if (result.success) {
        setStreamlines(result.streamlines);
      }
    } catch (e) {
      console.error('Streamline computation failed:', e);
    }
  }, [showStreamlines, panels, displayAlpha, reynolds, mach, ncrit, maxIterations, solverMode, adaptiveStreamlineCount, streamlineBounds, isDraggingPoint]);
  
  // Compute stream function (ψ) contours when enabled
  // Uses marching squares to extract iso-lines from the psi grid
  // Uses adaptive bounds (same as streamlines) via stableViewport for debouncing
  // NOTE: The dividing streamline (ψ = ψ₀) is extrapolated to hit the airfoil surface
  useEffect(() => {
    if (!showPsiContours || !isWasmReady() || panels.length < 10 || isDraggingPoint) {
      if (!isDraggingPoint) {
        setPsiContours({ grid: [], bounds: [0, 0, 0, 0], nx: 0, ny: 0, psiMin: 0, psiMax: 0, psi0: 0, lines: [], psi0Lines: [] });
      }
      return;
    }
    
    try {
      // Use adaptive bounds (same as streamlines) - updates when viewport settles
      const bounds = streamlineBounds;
      // Scale resolution based on bounds size (target ~40 points per unit for performance)
      const xRange = bounds[1] - bounds[0];
      const yRange = bounds[3] - bounds[2];
      const resolution: [number, number] = [
        Math.min(200, Math.round(xRange * 40)),
        Math.min(120, Math.round(yRange * 40))
      ];
      const result = computePsiGrid(panels, displayAlpha, reynolds, bounds, resolution, mach, ncrit, maxIterations, solverMode);
      
      if (!result.success) {
        console.error('Psi grid computation failed:', result.error);
        return;
      }
      
      // Extract contours using marching squares
      const { grid, psi_0, psi_min, psi_max, nx, ny } = result;
      const dx = (bounds[1] - bounds[0]) / (nx - 1);
      const dy = (bounds[3] - bounds[2]) / (ny - 1);
      
      // Point-in-polygon test for filtering interior points
      const isInsideAirfoil = (x: number, y: number): boolean => {
        if (panels.length < 3) return false;
        let inside = false;
        for (let i = 0, j = panels.length - 1; i < panels.length; j = i++) {
          const xi = panels[i].x, yi = panels[i].y;
          const xj = panels[j].x, yj = panels[j].y;
          if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            inside = !inside;
          }
        }
        return inside;
      };
      
      // Choose contour levels split around psi0 so both sides of the dividing
      // streamline get comparable line density even when the psi range is
      // strongly asymmetric.
      const nLevels = 12;
      const levels: number[] = [];
      for (let i = 0; i <= nLevels; i++) {
        levels.push(psi_min + (psi_0 - psi_min) * (i / nLevels));
      }
      for (let i = 1; i <= nLevels; i++) {
        levels.push(psi_0 + (psi_max - psi_0) * (i / nLevels));
      }
      
      // Marching squares for regular contour lines
      const allLines: [number, number][][] = [];
      
      for (const level of levels) {
        const contourSegments = marchingSquares(grid, nx, ny, level, bounds[0], bounds[2], dx, dy);
        const polylines = connectSegments(contourSegments);
        
        for (const line of polylines) {
          if (line.length >= 2) {
            allLines.push(line);
          }
        }
      }
      
      // ALWAYS compute the dividing streamline (ψ = ψ₀) separately
      // This ensures we get exactly the right contour regardless of level spacing
      const psi0Segments = marchingSquares(grid, nx, ny, psi_0, bounds[0], bounds[2], dx, dy);
      let psi0Lines = connectSegments(psi0Segments).filter(line => line.length >= 2);
      
      // Filter out parts of dividing streamline that are inside the airfoil
      // Split lines at airfoil boundary and keep only exterior segments
      const filteredPsi0Lines: [number, number][][] = [];
      for (const line of psi0Lines) {
        let currentSegment: [number, number][] = [];
        for (const point of line) {
          const inside = isInsideAirfoil(point[0], point[1]);
          if (!inside) {
            currentSegment.push(point);
          } else if (currentSegment.length >= 2) {
            filteredPsi0Lines.push(currentSegment);
            currentSegment = [];
          } else {
            currentSegment = [];
          }
        }
        if (currentSegment.length >= 2) {
          filteredPsi0Lines.push(currentSegment);
        }
      }
      psi0Lines = filteredPsi0Lines;
      
      // Extrapolate dividing streamline (ψ₀) to hit the airfoil surface
      // This extends the line segments to the airfoil boundary for a complete visualization
      psi0Lines = extrapolateDividingStreamline(psi0Lines, panels);
      
      setPsiContours({ 
        grid, 
        bounds, 
        nx, 
        ny, 
        psiMin: psi_min, 
        psiMax: psi_max, 
        psi0: psi_0,
        lines: allLines, 
        psi0Lines 
      });
    } catch (e) {
      console.error('Psi contour computation failed:', e);
    }
  }, [showPsiContours, panels, displayAlpha, reynolds, mach, ncrit, maxIterations, solverMode, streamlineBounds, isDraggingPoint]);
  
  // Compute aerodynamic analysis (Cp, Cl, Cm)
  // SKIP during drag to prevent freezing - recalculate on drag end
  useEffect(() => {
    if (!isWasmReady() || panels.length < 10 || isDraggingPoint) {
      return;
    }
    
    try {
      const result = analyzeAirfoil(panels, displayAlpha, reynolds, mach, ncrit, maxIterations);
      if (result.success) {
        setAnalysisResult({
          cp: result.cp,
          cpX: result.cp_x,
          cl: result.cl,
          cd: result.cd,
          cm: result.cm,
        });
      }
    } catch (e) {
      console.error('Analysis failed:', e);
    }
  }, [panels, displayAlpha, reynolds, mach, ncrit, maxIterations, isDraggingPoint]);
  
  // BL visualization data (for boundary layer envelope + wake overlay)
  const [blVisData, setBlVisData] = useState<BLVisualizationData | null>(null);
  
  useEffect(() => {
    const needsBLData = showBoundaryLayer || showWake || showDisplacementThickness
      || ((showStreamlines || showPsiContours) && solverMode === 'viscous');
    if (!needsBLData || !isWasmReady() || panels.length < 10 || isDraggingPoint) {
      if (!isDraggingPoint) setBlVisData(null);
      return;
    }
    
    try {
      const data = getBLVisualizationData(panels, displayAlpha, reynolds, mach, ncrit, maxIterations);
      if (data.success) {
        setBlVisData(data);
      }
    } catch (e) {
      console.error('BL visualization data failed:', e);
    }
  }, [showBoundaryLayer, showWake, showDisplacementThickness, showStreamlines, showPsiContours, solverMode, panels, displayAlpha, reynolds, mach, ncrit, maxIterations, isDraggingPoint]);
  
  // Morphing animation integration
  const morphTarget = useMemo(() => ({
    coordinates,
    panels,
    streamlines,
    cp: analysisResult.cp,
    cpX: analysisResult.cpX,
    cl: analysisResult.cl,
    cd: analysisResult.cd,
    cm: analysisResult.cm,
  }), [coordinates, panels, streamlines, analysisResult]);
  
  const morphState = useMorphingAnimation(morphTarget, {
    duration: morphDuration,
    // Disable morphing during drag for immediate feedback
    enabled: enableMorphing && !isDraggingPoint,
  });
  
  // Memoize the spline curve to avoid recomputing on every render
  // Use morphed coordinates for smooth animation
  const splineCurve = useMemo(() => {
    return generateSplineCurve(morphState.coordinates, 300);
  }, [morphState.coordinates]);
  
  // Store panels for smoke system - only update when NOT dragging
  const smokePanelsRef = useRef(panels);
  useEffect(() => {
    if (!isDraggingPoint) {
      smokePanelsRef.current = panels;
    }
  }, [panels, isDraggingPoint]);
  
  // Compute smoke spawn bounds based on viewport (for adaptive mode)
  const smokeBounds = useMemo(() => {
    if (!adaptiveStreamlines) {
      // Fixed mode: large fixed spawn area
      return { yMin: -1.5, yMax: 1.5, spawnX: -1.5, exitX: 2.5 };
    }
    
    // Adaptive mode: spawn based on stable viewport with padding
    const { center, zoom, width, height } = stableViewport;
    const viewHalfHeight = (height / zoom) * 1.2; // 20% padding
    const viewHalfWidth = (width / zoom) * 0.5;
    
    // Ensure minimum spawn area
    const yMin = center.y - Math.max(viewHalfHeight, 1.0);
    const yMax = center.y + Math.max(viewHalfHeight, 1.0);
    const spawnX = center.x - Math.max(viewHalfWidth, 1.0) - 0.5; // Spawn to left of visible area
    const exitX = center.x + Math.max(viewHalfWidth, 1.0) + 0.5;  // Exit to right of visible area
    
    return { yMin, yMax, spawnX, exitX };
  }, [adaptiveStreamlines, stableViewport]);

  // Initialize smoke system (only when toggled on or settings change)
  // SKIP recreation during drag - use cached panels
  // For adaptive mode: recreates when viewport settles (via smokeBounds)
  useEffect(() => {
    if (!showSmoke || !isWasmReady() || smokePanelsRef.current.length < 10) {
      if (smokeAnimationRef.current) {
        cancelAnimationFrame(smokeAnimationRef.current);
        smokeAnimationRef.current = null;
      }
      smokeSystemRef.current = null;
      smokePositionsRef.current = null;
      smokeAlphasRef.current = null;
      return;
    }
    
    // Don't recreate smoke system during drag
    if (isDraggingPoint && smokeSystemRef.current) {
      return;
    }
    
    // Create smoke system with spawn points based on smokeDensity and viewport bounds
    // Spawn points are rotated by +alpha so they appear as a vertical line in display frame
    // (display applies -alpha rotation, so we pre-rotate by +alpha to cancel it out)
    const { yMin, yMax, spawnX, exitX } = smokeBounds;
    const yRange = yMax - yMin;
    const cx = 0.25; // Rotation center (quarter chord)
    const cy = 0;
    const rad = displayAlpha * Math.PI / 180; // Positive alpha to counter display's -alpha
    const cosA = Math.cos(rad);
    const sinA = Math.sin(rad);
    
    // Calculate domain width and dynamically set age/interval for 3 tranches
    const domainWidth = exitX - spawnX;
    const traversalTime = domainWidth / flowSpeed; // Time for particle to cross visible domain
    const dynamicMaxAge = traversalTime * 1.3;     // 30% margin to ensure particles exit
    const dynamicSpawnInterval = dynamicMaxAge / 3; // 3 tranches visible at once
    
    // Generate spawn points as a vertical line in display frame, then rotate to body frame
    const spawnPoints: number[] = []; // Flat array [x0, y0, x1, y1, ...]
    for (let i = 0; i < smokeDensity; i++) {
      // Desired position in display frame (vertical line at spawnX)
      const displayX = spawnX;
      const displayY = yMin + (i * yRange) / Math.max(1, smokeDensity - 1);
      
      // Rotate by +alpha around quarter chord to get body frame position
      const dx = displayX - cx;
      const dy = displayY - cy;
      const bodyX = cx + dx * cosA - dy * sinA;
      const bodyY = cy + dx * sinA + dy * cosA;
      
      spawnPoints.push(bodyX, bodyY);
    }
    
    // Create smoke system with dummy values, then set actual spawn points
    smokeSystemRef.current = createSmokeSystem([0], spawnX, smokeParticlesPerBlob);
    
    // Set the properly rotated spawn points
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const smokeSystem = smokeSystemRef.current as any;
    if (typeof smokeSystem.set_spawn_points === 'function') {
      smokeSystem.set_spawn_points(new Float64Array(spawnPoints));
    }
    
    // Use dynamically calculated values for proper tranche distribution
    smokeSystemRef.current.set_spawn_interval(dynamicSpawnInterval);
    smokeSystemRef.current.set_max_age(dynamicMaxAge);
    
    // Set initial flow using cached panels
    const smokeCoords = new Float64Array(smokePanelsRef.current.flatMap(p => [p.x, p.y]));
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const smokeSystemAny = smokeSystemRef.current as any;
    if (solverMode === 'viscous' && typeof smokeSystemAny.set_faithful_flow === 'function') {
      smokeSystemAny.set_faithful_flow(
        smokeCoords,
        displayAlpha,
        reynolds,
        mach,
        ncrit,
        maxIterations,
      );
    } else {
      smokeSystemRef.current.set_flow(smokeCoords, displayAlpha);
    }
    
    // Set initial flow speed
    if (typeof smokeSystemRef.current.set_v_inf === 'function') {
      smokeSystemRef.current.set_v_inf(flowSpeed);
    }
    
    // Performance monitoring for adaptive density
    const frameTimes: number[] = [];
    const FRAME_WINDOW = 30; // Track last 30 frames
    const MIN_FPS = 25; // Target minimum FPS
    const MIN_DENSITY = 10; // Don't go below this
    let lastPerfCheck = performance.now();
    let currentDensity = smokeDensity;
    
    // Animation loop with performance monitoring
    let lastTime = performance.now();
    const animate = (time: number) => {
      const dt = Math.min((time - lastTime) / 1000, 0.05); // Cap dt at 50ms
      const frameTime = time - lastTime;
      lastTime = time;
      
      // Track frame times
      frameTimes.push(frameTime);
      if (frameTimes.length > FRAME_WINDOW) {
        frameTimes.shift();
      }
      
      // Check performance every second
      if (time - lastPerfCheck > 1000 && frameTimes.length >= FRAME_WINDOW) {
        lastPerfCheck = time;
        const avgFrameTime = frameTimes.reduce((a, b) => a + b, 0) / frameTimes.length;
        const fps = 1000 / avgFrameTime;
        
        // If FPS is too low and we can reduce density, do so
        if (fps < MIN_FPS && currentDensity > MIN_DENSITY) {
          const newDensity = Math.max(MIN_DENSITY, Math.floor(currentDensity * 0.7));
          if (newDensity < currentDensity) {
            console.log(`Smoke: Reducing density ${currentDensity} -> ${newDensity} (${fps.toFixed(1)} FPS)`);
            currentDensity = newDensity;
            setSmokeDensity(newDensity);
          }
        }
      }
      
      if (smokeSystemRef.current) {
        smokeSystemRef.current.update(dt);
        const pos = smokeSystemRef.current.get_positions();
        const alphas = smokeSystemRef.current.get_alphas();
        // Store in refs (no state update = no full canvas redraw)
        smokePositionsRef.current = new Float64Array(pos);
        smokeAlphasRef.current = new Float64Array(alphas);
        
        // Get psi values for coloring particles by dividing streamline
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const smokeSystem = smokeSystemRef.current as any;
        if (typeof smokeSystem.get_psi_values === 'function') {
          const psiVals = smokeSystem.get_psi_values();
          const psi0Val = smokeSystem.get_psi_0();
          smokePsiValuesRef.current = new Float64Array(psiVals);
          smokePsi0Ref.current = psi0Val;
        }
        
        // Draw smoke on overlay canvas (doesn't trigger main canvas redraw)
        if (drawSmokeRef.current) {
          drawSmokeRef.current();
        }
      }
      
      smokeAnimationRef.current = requestAnimationFrame(animate);
    };
    
    smokeAnimationRef.current = requestAnimationFrame(animate);
    
    return () => {
      if (smokeAnimationRef.current) {
        cancelAnimationFrame(smokeAnimationRef.current);
      }
    };
  // Note: smokeSpawnInterval and smokeMaxAge removed - now computed dynamically from domain/flowSpeed
  }, [showSmoke, smokeDensity, smokeParticlesPerBlob, flowSpeed, setSmokeDensity, isDraggingPoint, displayAlpha, smokeBounds, solverMode, reynolds, mach, ncrit, maxIterations]);
  
  // Update flow when alpha changes (without recreating the system)
  // SKIP during drag
  useEffect(() => {
    if (smokeSystemRef.current && smokePanelsRef.current.length >= 10 && !isDraggingPoint) {
      const smokeCoords = new Float64Array(smokePanelsRef.current.flatMap(p => [p.x, p.y]));
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const smokeSystemAny = smokeSystemRef.current as any;
      if (solverMode === 'viscous' && typeof smokeSystemAny.set_faithful_flow === 'function') {
        smokeSystemAny.set_faithful_flow(
          smokeCoords,
          displayAlpha,
          reynolds,
          mach,
          ncrit,
          maxIterations,
        );
      } else {
        smokeSystemRef.current.set_flow(smokeCoords, displayAlpha);
      }
    }
  }, [displayAlpha, isDraggingPoint, solverMode, reynolds, mach, ncrit, maxIterations]);

  // Interaction state
  const [isDragging, setIsDragging] = useState(false);
  const [isPanning, setIsPanning] = useState(false);
  const [dragTarget, setDragTarget] = useState<{ type: string; index: number | string } | null>(null);
  const [hoveredPoint, setHoveredPoint] = useState<{ type: string; index: number | string } | null>(null);
  const [blHoverInfo, setBlHoverInfo] = useState<{
    surface: 'upper' | 'lower';
    index: number;
    canvasX: number;
    canvasY: number;
    x: number; y: number;
    cp: number; cf: number; dstar: number; theta: number;
    h: number; ue: number; hk: number; laminar: boolean;
  } | null>(null);
  const lastMousePos = useRef<Point>({ x: 0, y: 0 });

  // Convert airfoil coordinates to canvas coordinates
  const toCanvas = useCallback((p: Point): Point => {
    return {
      x: viewport.width / 2 + (p.x - viewport.center.x) * viewport.zoom,
      y: viewport.height / 2 - (p.y - viewport.center.y) * viewport.zoom,
    };
  }, [viewport]);

  // Convert canvas coordinates to airfoil coordinates
  const toAirfoil = useCallback((p: Point): Point => {
    return {
      x: viewport.center.x + (p.x - viewport.width / 2) / viewport.zoom,
      y: viewport.center.y - (p.y - viewport.height / 2) / viewport.zoom,
    };
  }, [viewport]);

  // Find point at canvas position
  const findPointAt = useCallback((canvasPos: Point): { type: string; index: number | string } | null => {
    // Helper to get camber at x position (same as in draw function)
    const getCamberAt = (x: number): number => {
      if (camberControlPoints.length === 0) return 0;
      const camberSpline = generateCamberSplineCurve(camberControlPoints, 50);
      for (let i = 0; i < camberSpline.length - 1; i++) {
        if (x >= camberSpline[i].x && x <= camberSpline[i + 1].x) {
          const t = (camberSpline[i + 1].x - camberSpline[i].x) > 0 
            ? (x - camberSpline[i].x) / (camberSpline[i + 1].x - camberSpline[i].x)
            : 0;
          return camberSpline[i].y + t * (camberSpline[i + 1].y - camberSpline[i].y);
        }
      }
      if (x <= camberSpline[0].x) return camberSpline[0].y;
      return camberSpline[camberSpline.length - 1].y;
    };
    
    // Check camber control points
    if (controlMode === 'camber-spline') {
      for (const cp of camberControlPoints) {
        const cpCanvas = toCanvas({ x: cp.x, y: cp.y });
        const dist = Math.hypot(canvasPos.x - cpCanvas.x, canvasPos.y - cpCanvas.y);
        if (dist < HIT_RADIUS) {
          return { type: 'camber', index: cp.id };
        }
      }
    }

    // Check thickness control points  
    if (controlMode === 'thickness-spline') {
      for (const cp of thicknessControlPoints) {
        // Thickness points are shown on the upper surface at y = camber + t
        const camber = getCamberAt(cp.x);
        const cpCanvas = toCanvas({ x: cp.x, y: camber + cp.t });
        const dist = Math.hypot(canvasPos.x - cpCanvas.x, canvasPos.y - cpCanvas.y);
        if (dist < HIT_RADIUS) {
          return { type: 'thickness', index: cp.id };
        }
      }
    }

    return null;
  }, [controlMode, camberControlPoints, thicknessControlPoints, toCanvas]);

  // Draw smoke particles on overlay canvas (called every animation frame)
  // This is separate from main draw() to avoid expensive redraws of psi contours etc.
  const drawSmoke = useCallback(() => {
    const canvas = smokeCanvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    
    // Clear the overlay canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    const smokePositions = smokePositionsRef.current;
    const smokeAlphas = smokeAlphasRef.current;
    const smokePsiValues = smokePsiValuesRef.current;
    const smokePsi0 = smokePsi0Ref.current;
    
    if (!showSmoke || !smokePositions || !smokeAlphas) return;
    
    const count = smokePositions.length / 2;
    const hasPsiData = smokePsiValues && smokePsiValues.length === count;
    
    // Colors for above/below dividing streamline
    const colorAbove = isDark ? '#ff6b6b' : '#e63946';  // Red for upper surface
    const colorBelow = isDark ? '#4dabf7' : '#1971c2';  // Blue for lower surface
    const defaultColor = isDark ? '#61dafb' : '#228be6';
    
    // Helper to rotate a point around quarter chord by alpha
    const rotatePoint = (p: Point): Point => {
      if (displayAlpha === 0) return p;
      const cx = 0.25;
      const cy = 0;
      const rad = -displayAlpha * Math.PI / 180;
      const cos = Math.cos(rad);
      const sin = Math.sin(rad);
      const dx = p.x - cx;
      const dy = p.y - cy;
      return {
        x: cx + dx * cos - dy * sin,
        y: cy + dx * sin + dy * cos,
      };
    };
    
    for (let i = 0; i < count; i++) {
      const x = smokePositions[i * 2];
      const y = smokePositions[i * 2 + 1];
      const alpha = smokeAlphas[i] || 0;
      
      if (alpha > 0.01) {
        const pCanvas = toCanvas(rotatePoint({ x, y }));
        ctx.beginPath();
        
        // Choose color based on psi value relative to psi0
        if (hasPsiData) {
          const psi = smokePsiValues[i];
          ctx.fillStyle = psi > smokePsi0 ? colorAbove : colorBelow;
        } else {
          ctx.fillStyle = defaultColor;
        }
        
        ctx.globalAlpha = alpha * 0.6;
        ctx.arc(pCanvas.x, pCanvas.y, 2, 0, Math.PI * 2);
        ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
  }, [showSmoke, displayAlpha, toCanvas, isDark]);
  
  // Keep ref updated so animation loop can call it
  useEffect(() => {
    drawSmokeRef.current = drawSmoke;
  }, [drawSmoke]);
  
  // Clear smoke canvas when smoke is disabled or when viewport changes
  useEffect(() => {
    const canvas = smokeCanvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    
    if (!showSmoke) {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }
  }, [showSmoke, viewport]);
  
  // Draw the canvas
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const { width, height } = viewport;

    // Theme-based colors
    const colors = {
      bgPrimary: isDark ? '#0f0f0f' : '#ffffff',
      foilLine: isDark ? '#00d4aa' : '#00a88a',
      foilGrid: isDark ? '#333333' : '#dee2e6',
      foilPoint: isDark ? '#ffffff' : '#212529',
      foilPointSelected: isDark ? '#ffaa00' : '#e69500',
      foilControl: isDark ? '#0099ff' : '#0077cc',
      foilHandle: isDark ? '#ff6666' : '#dc3545',
      accentSecondary: isDark ? '#0099ff' : '#0077cc',
      accentWarning: isDark ? '#ffaa00' : '#e69500',
      textSecondary: isDark ? '#a0a0a0' : '#495057',
      textPrimary: isDark ? '#e0e0e0' : '#212529',
    };

    // Clear canvas
    ctx.fillStyle = colors.bgPrimary;
    ctx.fillRect(0, 0, width, height);

    // Draw grid (if enabled)
    if (showGrid) {
      drawGrid(ctx);
      drawAxes(ctx);
    }
    
    // Helper to rotate a point around quarter chord by alpha
    const rotatePoint = (p: Point): Point => {
      if (displayAlpha === 0) return p;
      const cx = 0.25; // Rotation center at quarter chord
      const cy = 0;
      const rad = -displayAlpha * Math.PI / 180; // Negative for visual convention
      const cos = Math.cos(rad);
      const sin = Math.sin(rad);
      const dx = p.x - cx;
      const dy = p.y - cy;
      return {
        x: cx + dx * cos - dy * sin,
        y: cy + dx * sin + dy * cos,
      };
    };
    
    // Draw stream function visualization using filled contour bands
    // Build "effective body" polygon for clipping contour lines and streamlines
    // at the BL envelope. Shared by psi contours and streamline rendering.
    let effectiveBodyPoly: { x: number; y: number }[] | null = null;

    if ((showBoundaryLayer || showWake) && blVisData && blVisData.success) {
      const sc = blThicknessScale;
      const poly: { x: number; y: number }[] = [];

      const envNormal = (pts: { x: number[]; y: number[] }, i: number, sign: number) => {
        let tx: number, ty: number;
        if (pts.x.length < 2) return { nx: 0, ny: sign };
        if (i === 0) { tx = pts.x[1] - pts.x[0]; ty = pts.y[1] - pts.y[0]; }
        else if (i === pts.x.length - 1) { tx = pts.x[i] - pts.x[i - 1]; ty = pts.y[i] - pts.y[i - 1]; }
        else { tx = pts.x[i + 1] - pts.x[i - 1]; ty = pts.y[i + 1] - pts.y[i - 1]; }
        const len = Math.sqrt(tx * tx + ty * ty) || 1e-10;
        return { nx: -ty / len * sign, ny: tx / len * sign };
      };

      const uData = blVisData.upper;
      for (let i = 0; i < uData.x.length; i++) {
        const n = envNormal(uData, i, 1);
        const ds = uData.delta_star[i] * sc;
        poly.push({ x: uData.x[i] + n.nx * ds, y: uData.y[i] + n.ny * ds });
      }

      const wk = blVisData.wake;
      const fU = blVisData.wake_upper_fraction;
      const fL = 1 - fU;
      if (showWake && wk.x.length > 0) {
        const blendN = Math.min(5, wk.x.length);
        const teNU = uData.x.length > 1 ? envNormal(uData, uData.x.length - 1, 1) : { nx: 0, ny: 1 };
        for (let i = 0; i < wk.x.length; i++) {
          const wn = envNormal(wk, i, 1);
          const t = i < blendN ? i / blendN : 1;
          const nx = teNU.nx * (1 - t) + wn.nx * t;
          const ny = teNU.ny * (1 - t) + wn.ny * t;
          const len = Math.sqrt(nx * nx + ny * ny) || 1;
          const ds = wk.delta_star[i] * sc * fU;
          poly.push({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds });
        }
        const lData2 = blVisData.lower;
        const teNL = lData2.x.length > 1 ? envNormal(lData2, lData2.x.length - 1, -1) : { nx: 0, ny: -1 };
        for (let i = wk.x.length - 1; i >= 0; i--) {
          const wn = envNormal(wk, i, 1);
          const t = i < blendN ? i / blendN : 1;
          const nx = teNL.nx * (1 - t) + (-wn.nx) * t;
          const ny = teNL.ny * (1 - t) + (-wn.ny) * t;
          const len = Math.sqrt(nx * nx + ny * ny) || 1;
          const ds = wk.delta_star[i] * sc * fL;
          poly.push({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds });
        }
      }

      const lData = blVisData.lower;
      for (let i = lData.x.length - 1; i >= 0; i--) {
        const n = envNormal(lData, i, -1);
        const ds = lData.delta_star[i] * sc;
        poly.push({ x: lData.x[i] + n.nx * ds, y: lData.y[i] + n.ny * ds });
      }

      if (poly.length > 4) effectiveBodyPoly = poly;
    }

    const insideEffectiveBody = effectiveBodyPoly ? (px: number, py: number): boolean => {
      const poly = effectiveBodyPoly!;
      let inside = false;
      for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
        const yi = poly[i].y, yj = poly[j].y;
        if ((yi > py) !== (yj > py)) {
          const xi = poly[i].x, xj = poly[j].x;
          const xIntersect = xi + (py - yi) / (yj - yi) * (xj - xi);
          if (px < xIntersect) inside = !inside;
        }
      }
      return inside;
    } : null;

    // Uses marching squares to generate iso-contour polygons at each threshold
    // Open polylines are closed by tracing along the grid boundary
    // - Blue tones: flow going under the airfoil (ψ < ψ₀)
    // - Red tones: flow going over the airfoil (ψ > ψ₀)
    if (showPsiContours && psiContours.grid.length > 0) {
      const { grid, bounds, nx, ny, psiMin, psiMax, psi0 } = psiContours;
      
      const [xMin, xMax, yMin, yMax] = bounds;
      const dx = (xMax - xMin) / (nx - 1);
      const dy = (yMax - yMin) / (ny - 1);
      
      // Generate all thresholds including psi0
      const nLevels = 12;
      const allThresholds: number[] = [];
      
      // Levels below ψ₀ (from psiMin to psi0)
      for (let i = 0; i <= nLevels; i++) {
        allThresholds.push(psiMin + (psi0 - psiMin) * (i / nLevels));
      }
      // Levels above ψ₀ (from psi0 to psiMax, skip duplicate psi0)
      for (let i = 1; i <= nLevels; i++) {
        allThresholds.push(psi0 + (psiMax - psi0) * (i / nLevels));
      }
      allThresholds.sort((a, b) => a - b);
      
      const alpha = isDark ? 0.65 : 0.75;
      
      // Get color based on stream function strength (distance from psi0)
      // Color intensity increases with distance from dividing streamline
      const getColor = (psiValue: number): [number, number, number] => {
        if (psiValue < psi0) {
          // Blue gradient: stronger blue for more negative (further from psi0)
          const t = Math.min(1, Math.max(0, (psi0 - psiValue) / (psi0 - psiMin + 1e-10)));
          const intensity = Math.pow(t, 0.5); // Square root for more variation in lighter tones
          return [
            Math.round(220 - intensity * 180),  // 220 -> 40
            Math.round(230 - intensity * 160),  // 230 -> 70
            Math.round(255 - intensity * 55)    // 255 -> 200
          ];
        } else {
          // Red gradient: stronger red for more positive (further from psi0)
          const t = Math.min(1, Math.max(0, (psiValue - psi0) / (psiMax - psi0 + 1e-10)));
          const intensity = Math.pow(t, 0.5); // Square root for more variation
          return [
            Math.round(255 - intensity * 55),   // 255 -> 200
            Math.round(220 - intensity * 160),  // 220 -> 60
            Math.round(210 - intensity * 150)   // 210 -> 60
          ];
        }
      };
      
      // Orient polyline so it goes left to right (by x-coordinate of start)
      const orientLeftToRight = (polyline: [number, number][]): [number, number][] => {
        if (polyline.length < 2) return polyline;
        const startX = polyline[0][0];
        const endX = polyline[polyline.length - 1][0];
        if (startX > endX) {
          return [...polyline].reverse();
        }
        return polyline;
      };
      
      // Find closest point on airfoil boundary
      const closestPointOnAirfoil = (x: number, y: number): { x: number; y: number; idx: number } => {
        let minDist = Infinity;
        let closest = { x: panels[0]?.x || 0, y: panels[0]?.y || 0, idx: 0 };
        for (let i = 0; i < panels.length; i++) {
          const dx = panels[i].x - x;
          const dy = panels[i].y - y;
          const dist = dx * dx + dy * dy;
          if (dist < minDist) {
            minDist = dist;
            closest = { x: panels[i].x, y: panels[i].y, idx: i };
          }
        }
        return closest;
      };
      
      // Trace along airfoil from index i1 to i2 (shortest path)
      const traceAirfoil = (i1: number, i2: number): [number, number][] => {
        const n = panels.length;
        const points: [number, number][] = [];
        
        // Determine direction (shorter path around the airfoil)
        const fwdDist = (i2 - i1 + n) % n;
        const bwdDist = (i1 - i2 + n) % n;
        
        if (fwdDist <= bwdDist) {
          // Go forward
          let i = i1;
          while (i !== i2) {
            points.push([panels[i].x, panels[i].y]);
            i = (i + 1) % n;
          }
          points.push([panels[i2].x, panels[i2].y]);
        } else {
          // Go backward
          let i = i1;
          while (i !== i2) {
            points.push([panels[i].x, panels[i].y]);
            i = (i - 1 + n) % n;
          }
          points.push([panels[i2].x, panels[i2].y]);
        }
        return points;
      };
      
      // Check if a point is on the grid boundary
      const isOnGridBoundary = (x: number, y: number): boolean => {
        const tol = dx * 2;
        return Math.abs(x - xMin) < tol || Math.abs(x - xMax) < tol ||
               Math.abs(y - yMin) < tol || Math.abs(y - yMax) < tol;
      };
      
      // Check if a point is in the interior (not on grid boundary) - likely near airfoil
      const isInInterior = (x: number, y: number): boolean => {
        return !isOnGridBoundary(x, y);
      };
      
      // Get contour for a threshold, combining and sorting all segments
      // For contours with interior endpoints (like dividing streamline), connect via airfoil trace
      const getPrimaryContour = (threshold: number, connectInterior: boolean = false): [number, number][] | null => {
        const segments = marchingSquares(grid, nx, ny, threshold, xMin, yMin, dx, dy);
        const polylines = connectSegments(segments);
        if (polylines.length === 0) return null;
        
        // Orient all polylines L→R and sort them by their leftmost x
        const oriented = polylines.map(pl => orientLeftToRight(pl));
        oriented.sort((a, b) => a[0][0] - b[0][0]);
        
        // Combine all polylines into one, connecting gaps
        let result: [number, number][] = [];
        for (const pl of oriented) {
          if (result.length === 0) {
            result = [...pl];
          } else {
            // Connect previous end to this start
            result.push(...pl);
          }
        }
        
        if (result.length < 2) return null;
        
        // Final orientation check
        result = orientLeftToRight(result);
        
        // If both endpoints are in interior, connect them by tracing along airfoil
        if (connectInterior) {
          const start = result[0];
          const end = result[result.length - 1];
          const startInInterior = isInInterior(start[0], start[1]);
          const endInInterior = isInInterior(end[0], end[1]);
          
          if (startInInterior && endInInterior) {
            // Both endpoints near airfoil - trace along airfoil surface to close the loop
            const p1 = closestPointOnAirfoil(start[0], start[1]);
            const p2 = closestPointOnAirfoil(end[0], end[1]);
            const airfoilPath = traceAirfoil(p2.idx, p1.idx);
            result = [...result, ...airfoilPath];
          } else if (startInInterior || endInInterior) {
            // One end in interior - find closest point on airfoil and add it
            if (startInInterior) {
              const closest = closestPointOnAirfoil(start[0], start[1]);
              result = [[closest.x, closest.y], ...result];
            }
            if (endInInterior) {
              const closest = closestPointOnAirfoil(end[0], end[1]);
              result = [...result, [closest.x, closest.y]];
            }
          }
        }
        
        return result;
      };
      
      
      // Build filled bands between adjacent contours
      for (let i = 0; i < allThresholds.length - 1; i++) {
        const t1 = allThresholds[i];
        const t2 = allThresholds[i + 1];
        const tMid = (t1 + t2) / 2;
        
        // Connect interior endpoints to complete all contours
        const contour1 = getPrimaryContour(t1, true);
        const contour2 = getPrimaryContour(t2, true);
        
        if (!contour1 || !contour2 || contour1.length < 2 || contour2.length < 2) {
          continue;
        }
        
        const polygon: [number, number][] = [];
        
        // Add contour1 (already L→R)
        for (const pt of contour1) {
          polygon.push(pt);
        }
        
        // Check if we need to trace along airfoil between contour1 end and contour2 end
        const c1End = contour1[contour1.length - 1];
        const c2End = contour2[contour2.length - 1];
        
        // If either end is in interior (near airfoil), trace along airfoil
        const c1EndInterior = isInInterior(c1End[0], c1End[1]);
        const c2EndInterior = isInInterior(c2End[0], c2End[1]);
        
        if (c1EndInterior || c2EndInterior) {
          const p1 = closestPointOnAirfoil(c1End[0], c1End[1]);
          const p2 = closestPointOnAirfoil(c2End[0], c2End[1]);
          if (c1EndInterior && c2EndInterior) {
            // Both in interior - trace along airfoil
            const airfoilPath = traceAirfoil(p1.idx, p2.idx);
            for (const pt of airfoilPath) {
              polygon.push(pt);
            }
          } else if (c1EndInterior) {
            // Only c1 end in interior - add its closest airfoil point
            polygon.push([p1.x, p1.y]);
          } else {
            // Only c2 end in interior - add its closest airfoil point  
            polygon.push([p2.x, p2.y]);
          }
        }
        
        // Add contour2 reversed (R→L)
        for (let j = contour2.length - 1; j >= 0; j--) {
          polygon.push(contour2[j]);
        }
        
        // Check if we need to trace along airfoil between contour2 start and contour1 start
        const c2Start = contour2[0];
        const c1Start = contour1[0];
        
        // If either start is in interior (near airfoil), trace along airfoil
        const c2StartInterior = isInInterior(c2Start[0], c2Start[1]);
        const c1StartInterior = isInInterior(c1Start[0], c1Start[1]);
        
        if (c2StartInterior || c1StartInterior) {
          const p2 = closestPointOnAirfoil(c2Start[0], c2Start[1]);
          const p1 = closestPointOnAirfoil(c1Start[0], c1Start[1]);
          if (c2StartInterior && c1StartInterior) {
            // Both in interior - trace along airfoil
            const airfoilPath = traceAirfoil(p2.idx, p1.idx);
            for (const pt of airfoilPath) {
              polygon.push(pt);
            }
          } else if (c2StartInterior) {
            // Only c2 start in interior - add its closest airfoil point
            polygon.push([p2.x, p2.y]);
          } else {
            // Only c1 start in interior - add its closest airfoil point
            polygon.push([p1.x, p1.y]);
          }
        }
        
        if (polygon.length < 3) continue;
        
        // Color based on stream function strength (distance from psi0)
        const [r, g, b] = getColor(tMid);
        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${alpha})`;
        
        ctx.beginPath();
        const first = toCanvas(rotatePoint({ x: polygon[0][0], y: polygon[0][1] }));
        ctx.moveTo(first.x, first.y);
        for (let j = 1; j < polygon.length; j++) {
          const p = toCanvas(rotatePoint({ x: polygon[j][0], y: polygon[j][1] }));
          ctx.lineTo(p.x, p.y);
        }
        ctx.closePath();
        ctx.fill();
      }
      
      // Draw all contour branches, clipping inside the BL envelope.
      ctx.strokeStyle = isDark ? 'rgba(255, 255, 255, 0.3)' : 'rgba(0, 0, 0, 0.2)';
      ctx.lineWidth = 0.8;
      for (const line of psiContours.lines) {
        if (line.length < 2) continue;
        let drawing = false;
        for (let j = 0; j < line.length; j++) {
          const px = line[j][0];
          const py = line[j][1];
          if (insideEffectiveBody && insideEffectiveBody(px, py)) {
            if (drawing) { ctx.stroke(); drawing = false; }
            continue;
          }
          const p = toCanvas(rotatePoint({ x: px, y: py }));
          if (!drawing) {
            ctx.beginPath();
            ctx.moveTo(p.x, p.y);
            drawing = true;
          } else {
            ctx.lineTo(p.x, p.y);
          }
        }
        if (drawing) ctx.stroke();
      }
      
      // Mask out the airfoil interior
      if (panels.length > 2) {
        ctx.fillStyle = isDark ? '#1a1a2e' : '#ffffff';
        ctx.beginPath();
        const firstPanel = toCanvas(rotatePoint(panels[0]));
        ctx.moveTo(firstPanel.x, firstPanel.y);
        for (let i = 1; i < panels.length; i++) {
          const p = toCanvas(rotatePoint(panels[i]));
          ctx.lineTo(p.x, p.y);
        }
        ctx.closePath();
        ctx.fill();
      }
      
      // Draw dividing streamline, clipping inside BL envelope
      if (psiContours.psi0Lines.length > 0) {
        ctx.strokeStyle = isDark ? 'rgba(255, 255, 255, 0.95)' : 'rgba(0, 0, 0, 0.9)';
        ctx.lineWidth = 2.5;
        ctx.setLineDash([]);
        
        for (const line of psiContours.psi0Lines) {
          if (line.length < 2) continue;
          let drawing = false;
          for (let i = 0; i < line.length; i++) {
            const px = line[i][0];
            const py = line[i][1];
            if (insideEffectiveBody && insideEffectiveBody(px, py)) {
              if (drawing) { ctx.stroke(); drawing = false; }
              continue;
            }
            const p = toCanvas(rotatePoint({ x: px, y: py }));
            if (!drawing) {
              ctx.beginPath(); ctx.moveTo(p.x, p.y); drawing = true;
            } else {
              ctx.lineTo(p.x, p.y);
            }
          }
          if (drawing) ctx.stroke();
        }
      }
    }
    
    // Draw streamlines (behind airfoil) - rotated by alpha to match airfoil
    // Uses morphed streamlines for smooth animation
    // Skip CPU streamlines when GPU mode is enabled - GPU uses same data
    if (showStreamlines && morphState.streamlines.length > 0 && !useGPU) {
      ctx.strokeStyle = colors.accentSecondary;
      ctx.globalAlpha = 0.5;
      ctx.lineWidth = 1;

      for (const line of morphState.streamlines) {
        if (line.length < 2) continue;
        let drawing = false;
        
        for (let i = 0; i < line.length; i++) {
          const px = line[i][0];
          const py = line[i][1];
          const inside = insideEffectiveBody ? insideEffectiveBody(px, py) : false;

          if (inside) {
            if (drawing) {
              ctx.stroke();
              drawing = false;
            }
            continue;
          }

          const p = toCanvas(rotatePoint({ x: px, y: py }));
          if (!drawing) {
            ctx.beginPath();
            ctx.moveTo(p.x, p.y);
            drawing = true;
          } else {
            ctx.lineTo(p.x, p.y);
          }
        }
        if (drawing) ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }
    
    // NOTE: Smoke is drawn on a separate overlay canvas (smokeCanvasRef) for performance
    // This avoids redrawing expensive psi contours on every animation frame

    // Draw smooth spline curve (if enabled) - rotated by alpha
    if (showCurve && splineCurve.length > 1) {
      ctx.beginPath();
      ctx.strokeStyle = colors.foilLine;
      ctx.lineWidth = 2;
      
      const first = toCanvas(rotatePoint(splineCurve[0]));
      ctx.moveTo(first.x, first.y);
      
      for (let i = 1; i < splineCurve.length; i++) {
        const p = toCanvas(rotatePoint(splineCurve[i]));
        ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
    }
    
    // Draw panel lines (linear connections between panel points) - rotated by alpha
    // Uses morphed panels for smooth animation
    const morphedPanels = morphState.panels;
    if (showPanels && morphedPanels.length > 1) {
      ctx.beginPath();
      ctx.strokeStyle = colors.accentWarning;
      ctx.globalAlpha = 0.6;
      ctx.lineWidth = 1;
      
      const first = toCanvas(rotatePoint(morphedPanels[0]));
      ctx.moveTo(first.x, first.y);
      
      for (let i = 1; i < morphedPanels.length; i++) {
        const p = toCanvas(rotatePoint(morphedPanels[i]));
        ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
      ctx.globalAlpha = 1; // Reset alpha
    }


    // Draw camber line and control points (if in camber-spline mode)
    if (showControls && controlMode === 'camber-spline' && camberControlPoints.length > 0) {
      // Generate smooth spline curve through control points
      const splineCurve = generateCamberSplineCurve(camberControlPoints, 100);
      
      ctx.beginPath();
      ctx.strokeStyle = colors.foilControl;
      ctx.globalAlpha = 0.7;
      ctx.lineWidth = 2;
      ctx.setLineDash([6, 4]);
      
      if (splineCurve.length > 0) {
        const firstCamber = toCanvas(rotatePoint({ x: splineCurve[0].x, y: splineCurve[0].y }));
        ctx.moveTo(firstCamber.x, firstCamber.y);
        
        for (let i = 1; i < splineCurve.length; i++) {
          const p = toCanvas(rotatePoint({ x: splineCurve[i].x, y: splineCurve[i].y }));
          ctx.lineTo(p.x, p.y);
        }
      }
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      
      // Draw camber control points
      for (const cp of camberControlPoints) {
        const cpCanvas = toCanvas(rotatePoint({ x: cp.x, y: cp.y }));
        const isHovered = hoveredPoint?.type === 'camber' && hoveredPoint.index === cp.id;
        
        ctx.beginPath();
        ctx.fillStyle = isHovered ? colors.foilPointSelected : colors.foilControl;
        ctx.strokeStyle = colors.foilPoint;
        ctx.lineWidth = 2;
        ctx.arc(cpCanvas.x, cpCanvas.y, CONTROL_RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

    // Draw thickness envelope and control points (if in thickness-spline mode)
    if (showControls && controlMode === 'thickness-spline' && thicknessControlPoints.length > 0) {
      // Generate smooth spline curves for camber and thickness
      const camberSpline = generateCamberSplineCurve(camberControlPoints, 100);
      
      // Get camber value at x using spline (with fallback)
      const getCamberAt = (x: number): number => {
        if (camberSpline.length === 0) return 0;
        // Find closest points and interpolate
        for (let i = 0; i < camberSpline.length - 1; i++) {
          if (x >= camberSpline[i].x && x <= camberSpline[i + 1].x) {
            const t = (camberSpline[i + 1].x - camberSpline[i].x) > 0 
              ? (x - camberSpline[i].x) / (camberSpline[i + 1].x - camberSpline[i].x)
              : 0;
            return camberSpline[i].y + t * (camberSpline[i + 1].y - camberSpline[i].y);
          }
        }
        if (x <= camberSpline[0].x) return camberSpline[0].y;
        return camberSpline[camberSpline.length - 1].y;
      };
      
      // Generate smooth thickness envelope using spline through thickness control points
      const sortedThickness = [...thicknessControlPoints].sort((a, b) => a.x - b.x);
      const thicknessPoints = sortedThickness.map(p => ({ x: p.x, y: p.t }));
      const thicknessSpline = generateCamberSplineCurve(
        thicknessPoints.map((p, i) => ({ id: `t-${i}`, x: p.x, y: p.y })),
        100
      );
      
      // Draw upper thickness envelope as smooth curve
      ctx.beginPath();
      ctx.strokeStyle = colors.accentWarning;
      ctx.globalAlpha = 0.5;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      
      if (thicknessSpline.length > 0) {
        const firstUpper = toCanvas(rotatePoint({ 
          x: thicknessSpline[0].x, 
          y: getCamberAt(thicknessSpline[0].x) + Math.max(0, thicknessSpline[0].y)
        }));
        ctx.moveTo(firstUpper.x, firstUpper.y);
        
        for (let i = 1; i < thicknessSpline.length; i++) {
          const camber = getCamberAt(thicknessSpline[i].x);
          const p = toCanvas(rotatePoint({ 
            x: thicknessSpline[i].x, 
            y: camber + Math.max(0, thicknessSpline[i].y)
          }));
          ctx.lineTo(p.x, p.y);
        }
      }
      ctx.stroke();
      
      // Draw lower thickness envelope as smooth curve
      ctx.beginPath();
      if (thicknessSpline.length > 0) {
        const firstLower = toCanvas(rotatePoint({ 
          x: thicknessSpline[0].x, 
          y: getCamberAt(thicknessSpline[0].x) - Math.max(0, thicknessSpline[0].y)
        }));
        ctx.moveTo(firstLower.x, firstLower.y);
        
        for (let i = 1; i < thicknessSpline.length; i++) {
          const camber = getCamberAt(thicknessSpline[i].x);
          const p = toCanvas(rotatePoint({ 
            x: thicknessSpline[i].x, 
            y: camber - Math.max(0, thicknessSpline[i].y)
          }));
          ctx.lineTo(p.x, p.y);
        }
      }
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      
      // Draw thickness control points (shown at upper surface)
      for (const cp of thicknessControlPoints) {
        const camber = getCamberAt(cp.x);
        const cpCanvas = toCanvas(rotatePoint({ x: cp.x, y: camber + cp.t }));
        const isHovered = hoveredPoint?.type === 'thickness' && hoveredPoint.index === cp.id;
        
        ctx.beginPath();
        ctx.fillStyle = isHovered ? colors.foilPointSelected : colors.accentWarning;
        ctx.strokeStyle = colors.foilPoint;
        ctx.lineWidth = 2;
        ctx.arc(cpCanvas.x, cpCanvas.y, CONTROL_RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

    // Draw panel points (small white dots to show spacing distribution) - rotated by alpha
    if (showPoints) {
      ctx.fillStyle = colors.foilPoint;
      ctx.strokeStyle = colors.textPrimary;
      ctx.globalAlpha = 0.6;
      ctx.lineWidth = 1;
      for (const pt of morphedPanels) {
        const pCanvas = toCanvas(rotatePoint(pt));
        ctx.beginPath();
        ctx.arc(pCanvas.x, pCanvas.y, PANEL_POINT_RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
      ctx.globalAlpha = 1; // Reset alpha
    }
    
    // Draw Cp (pressure coefficient) visualization
    // Cp values are ordered by panel index, so use panel points directly
    if (showCp && morphState.cp.length > 0) {
      const cpLen = Math.min(morphState.cp.length, morphedPanels.length);
      
      for (let i = 0; i < cpLen; i++) {
        const cpVal = morphState.cp[i];
        const closestPt = morphedPanels[i];
        
        if (!closestPt) continue;
        
        const cpColor = getCpColor(cpVal, isDark);
        
        // Draw Cp bar (always pointing outward, length = |Cp|, color indicates sign)
        if (cpDisplayMode === 'bars' || cpDisplayMode === 'both') {
          // Bar length proportional to |Cp|, always pointing outward
          const barLength = Math.abs(cpVal) * cpBarScale * viewport.zoom;
          const pCanvas = toCanvas(rotatePoint(closestPt));
          
          // Use panel index directly (i is already the panel index)
          const prevIdx = Math.max(i - 1, 0);
          const nextIdx = Math.min(i + 1, morphedPanels.length - 1);
          const dx = morphedPanels[nextIdx].x - morphedPanels[prevIdx].x;
          const dy = morphedPanels[nextIdx].y - morphedPanels[prevIdx].y;
          const len = Math.sqrt(dx * dx + dy * dy) || 1;
          
          // Normal perpendicular to tangent
          let nx = -dy / len;
          let ny = dx / len;
          
          // Compute airfoil centroid (approximate center)
          let cx = 0, cy = 0;
          for (const pt of morphedPanels) {
            cx += pt.x;
            cy += pt.y;
          }
          cx /= morphedPanels.length;
          cy /= morphedPanels.length;
          
          // Vector from centroid to surface point
          const toCentroidX = cx - closestPt.x;
          const toCentroidY = cy - closestPt.y;
          
          // If normal points toward centroid, flip it
          const dotProduct = nx * toCentroidX + ny * toCentroidY;
          if (dotProduct > 0) {
            nx = -nx;
            ny = -ny;
          }
          
          ctx.beginPath();
          ctx.strokeStyle = cpColor;
          ctx.lineWidth = 2;
          ctx.globalAlpha = 0.8;
          ctx.moveTo(pCanvas.x, pCanvas.y);
          ctx.lineTo(pCanvas.x + nx * barLength, pCanvas.y - ny * barLength);
          ctx.stroke();
        }
        
        // Draw colored point on surface
        if (cpDisplayMode === 'color' || cpDisplayMode === 'both') {
          const pCanvas = toCanvas(rotatePoint(closestPt));
          ctx.beginPath();
          ctx.fillStyle = cpColor;
          ctx.globalAlpha = 0.9;
          ctx.arc(pCanvas.x, pCanvas.y, 3, 0, Math.PI * 2);
          ctx.fill();
        }
      }
      ctx.globalAlpha = 1;
    }
    
    // Draw force vectors (lift and drag)
    if (showForces && Math.abs(morphState.cl) > 0.001) {
      const forces = computeForceVectors(morphState.cl, displayAlpha, forceScale);
      
      // Draw at quarter chord
      const quarterChord = toCanvas(rotatePoint({ x: 0.25, y: 0 }));
      
      // Lift vector (perpendicular to flow, usually up)
      const liftEnd = {
        x: quarterChord.x + forces.lift.x * viewport.zoom,
        y: quarterChord.y - forces.lift.y * viewport.zoom, // Canvas y is inverted
      };
      
      ctx.beginPath();
      ctx.strokeStyle = isDark ? '#00ff88' : '#00aa55';
      ctx.lineWidth = 3;
      ctx.globalAlpha = 0.9;
      ctx.moveTo(quarterChord.x, quarterChord.y);
      ctx.lineTo(liftEnd.x, liftEnd.y);
      ctx.stroke();
      
      // Draw arrowhead for lift
      const liftAngle = Math.atan2(quarterChord.y - liftEnd.y, liftEnd.x - quarterChord.x);
      const arrowLen = 10;
      ctx.beginPath();
      ctx.moveTo(liftEnd.x, liftEnd.y);
      ctx.lineTo(
        liftEnd.x - arrowLen * Math.cos(liftAngle - Math.PI / 6),
        liftEnd.y + arrowLen * Math.sin(liftAngle - Math.PI / 6)
      );
      ctx.moveTo(liftEnd.x, liftEnd.y);
      ctx.lineTo(
        liftEnd.x - arrowLen * Math.cos(liftAngle + Math.PI / 6),
        liftEnd.y + arrowLen * Math.sin(liftAngle + Math.PI / 6)
      );
      ctx.stroke();
      
      // Label for Cl
      ctx.fillStyle = isDark ? '#00ff88' : '#00aa55';
      ctx.font = '12px sans-serif';
      ctx.fillText(`Cl = ${morphState.cl.toFixed(3)}`, liftEnd.x + 10, liftEnd.y);
      
      // Drag would go here if we had Cd (inviscid method doesn't compute drag)
      
      ctx.globalAlpha = 1;
    }

    // Draw boundary layer displacement thickness envelope + wake as one continuous shape
    if ((showBoundaryLayer || showWake) && blVisData && blVisData.success) {
      const scale = blThicknessScale;
      const laminarColor = isDark ? 'rgba(0, 180, 255, 0.25)' : 'rgba(0, 120, 220, 0.2)';
      const turbulentColor = isDark ? 'rgba(255, 100, 60, 0.25)' : 'rgba(220, 60, 30, 0.2)';
      const laminarStroke = isDark ? 'rgba(0, 180, 255, 0.6)' : 'rgba(0, 120, 220, 0.5)';
      const turbulentStroke = isDark ? 'rgba(255, 100, 60, 0.6)' : 'rgba(220, 60, 30, 0.5)';
      const wakeColor = isDark ? 'rgba(255, 200, 50, 0.2)' : 'rgba(180, 140, 30, 0.15)';
      const wakeStroke = isDark ? 'rgba(255, 200, 50, 0.5)' : 'rgba(180, 140, 30, 0.4)';

      // Helper: compute outward normal for a point in a polyline
      const computeNormal = (pts: { x: number[]; y: number[] }, i: number, sign: number) => {
        let tx: number, ty: number;
        if (pts.x.length < 2) return { nx: 0, ny: sign };
        if (i === 0) {
          tx = pts.x[1] - pts.x[0]; ty = pts.y[1] - pts.y[0];
        } else if (i === pts.x.length - 1) {
          tx = pts.x[i] - pts.x[i - 1]; ty = pts.y[i] - pts.y[i - 1];
        } else {
          tx = pts.x[i + 1] - pts.x[i - 1]; ty = pts.y[i + 1] - pts.y[i - 1];
        }
        const len = Math.sqrt(tx * tx + ty * ty) || 1e-10;
        return { nx: -ty / len * sign, ny: tx / len * sign };
      };

      // Build surface envelope + surface points for each side
      if (showBoundaryLayer) {
        for (const surface of ['upper', 'lower'] as const) {
          const data = blVisData[surface];
          if (data.x.length < 2) continue;
          const sign = surface === 'upper' ? 1 : -1;

          const surfacePoints: Point[] = [];
          const envelopePoints: Point[] = [];
          for (let i = 0; i < data.x.length; i++) {
            const n = computeNormal(data, i, sign);
            const ds = data.delta_star[i] * scale;
            surfacePoints.push({ x: data.x[i], y: data.y[i] });
            envelopePoints.push({ x: data.x[i] + n.nx * ds, y: data.y[i] + n.ny * ds });
          }

          // Draw laminar/turbulent segments as filled regions
          let segStart = 0;
          for (let i = 0; i <= data.x.length; i++) {
            const isLam = i < data.x.length ? data.is_laminar[i] : !data.is_laminar[data.x.length - 1];
            const prevLam = i > 0 ? data.is_laminar[i - 1] : data.is_laminar[0];
            
            if (i === data.x.length || isLam !== prevLam) {
              const segEnd = Math.min(i, data.x.length - 1);
              if (segEnd > segStart) {
                const segLaminar = data.is_laminar[segStart];
                ctx.beginPath();
                ctx.fillStyle = segLaminar ? laminarColor : turbulentColor;
                const firstSurf = toCanvas(rotatePoint(surfacePoints[segStart]));
                ctx.moveTo(firstSurf.x, firstSurf.y);
                for (let j = segStart + 1; j <= segEnd; j++) {
                  const p = toCanvas(rotatePoint(surfacePoints[j]));
                  ctx.lineTo(p.x, p.y);
                }
                for (let j = segEnd; j >= segStart; j--) {
                  const p = toCanvas(rotatePoint(envelopePoints[j]));
                  ctx.lineTo(p.x, p.y);
                }
                ctx.closePath();
                ctx.fill();

                ctx.beginPath();
                ctx.strokeStyle = segLaminar ? laminarStroke : turbulentStroke;
                ctx.lineWidth = 1.5;
                const firstEnv = toCanvas(rotatePoint(envelopePoints[segStart]));
                ctx.moveTo(firstEnv.x, firstEnv.y);
                for (let j = segStart + 1; j <= segEnd; j++) {
                  const p = toCanvas(rotatePoint(envelopePoints[j]));
                  ctx.lineTo(p.x, p.y);
                }
                ctx.stroke();
              }
              segStart = i;
            }
          }
        }

        // Transition markers
        for (const [xTr, label, color] of [
          [blVisData.x_tr_upper, 'Tr', isDark ? '#00b4ff' : '#0078dc'],
          [blVisData.x_tr_lower, 'Tr', isDark ? '#ff643c' : '#dc3c1e'],
        ] as [number, string, string][]) {
          if (xTr < 1 && xTr > 0) {
            const p = toCanvas(rotatePoint({ x: xTr, y: 0 }));
            ctx.beginPath();
            ctx.fillStyle = color;
            ctx.arc(p.x, p.y, 5, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = isDark ? '#fff' : '#000';
            ctx.font = '10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(label, p.x, p.y - 8);
          }
        }
      }

      // Draw wake: continuous from upper TE → downstream → lower TE
      if (showWake) {
        const wk = blVisData.wake;
        const fU = blVisData.wake_upper_fraction;
        const fL = 1 - fU;
        const wakeGeoX = blVisData.wake_geometry_x;
        const wakeGeoY = blVisData.wake_geometry_y;

        if (wakeGeoX.length > 1) {
          ctx.beginPath();
          ctx.strokeStyle = wakeStroke;
          ctx.lineWidth = 1;
          ctx.setLineDash([6, 4]);
          const first = toCanvas(rotatePoint({ x: wakeGeoX[0], y: wakeGeoY[0] }));
          ctx.moveTo(first.x, first.y);
          for (let i = 1; i < wakeGeoX.length; i++) {
            const p = toCanvas(rotatePoint({ x: wakeGeoX[i], y: wakeGeoY[i] }));
            ctx.lineTo(p.x, p.y);
          }
          ctx.stroke();
          ctx.setLineDash([]);
        }

        if (wk.x.length > 1) {
          const blendStations = Math.min(5, wk.x.length);
          let teNormalUpper = { nx: 0, ny: 1 };
          let teNormalLower = { nx: 0, ny: -1 };
          if (showBoundaryLayer && blVisData.upper.x.length > 1) {
            const uD = blVisData.upper;
            teNormalUpper = computeNormal(uD, uD.x.length - 1, 1);
          }
          if (showBoundaryLayer && blVisData.lower.x.length > 1) {
            const lD = blVisData.lower;
            teNormalLower = computeNormal(lD, lD.x.length - 1, -1);
          }

          // Upper edge: fraction fU of total wake δ*, blended normal
          const upperEnv: Point[] = [];
          for (let i = 0; i < wk.x.length; i++) {
            const wakeN = computeNormal(wk, i, 1);
            const t = i < blendStations ? i / blendStations : 1;
            const nx = teNormalUpper.nx * (1 - t) + wakeN.nx * t;
            const ny = teNormalUpper.ny * (1 - t) + wakeN.ny * t;
            const len = Math.sqrt(nx * nx + ny * ny) || 1;
            const ds = wk.delta_star[i] * scale * fU;
            upperEnv.push({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds });
          }

          // Lower edge: fraction fL of total wake δ*, blended normal
          const lowerEnv: Point[] = [];
          for (let i = 0; i < wk.x.length; i++) {
            const wakeN = computeNormal(wk, i, 1);
            const t = i < blendStations ? i / blendStations : 1;
            const nx = teNormalLower.nx * (1 - t) + (-wakeN.nx) * t;
            const ny = teNormalLower.ny * (1 - t) + (-wakeN.ny) * t;
            const len = Math.sqrt(nx * nx + ny * ny) || 1;
            const ds = wk.delta_star[i] * scale * fL;
            lowerEnv.push({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds });
          }

          // Fill: bridge from upper BL TE → wake upper → wake lower → lower BL TE
          ctx.beginPath();
          ctx.fillStyle = wakeColor;

          if (showBoundaryLayer && blVisData.upper.x.length > 0) {
            const uData = blVisData.upper;
            const li = uData.x.length - 1;
            const nU = computeNormal(uData, li, 1);
            const dsU = uData.delta_star[li] * scale;
            const p = toCanvas(rotatePoint({ x: uData.x[li] + nU.nx * dsU, y: uData.y[li] + nU.ny * dsU }));
            ctx.moveTo(p.x, p.y);
          } else {
            const p = toCanvas(rotatePoint(upperEnv[0]));
            ctx.moveTo(p.x, p.y);
          }

          for (let i = 0; i < upperEnv.length; i++) {
            const p = toCanvas(rotatePoint(upperEnv[i]));
            ctx.lineTo(p.x, p.y);
          }
          for (let i = lowerEnv.length - 1; i >= 0; i--) {
            const p = toCanvas(rotatePoint(lowerEnv[i]));
            ctx.lineTo(p.x, p.y);
          }

          if (showBoundaryLayer && blVisData.lower.x.length > 0) {
            const lData = blVisData.lower;
            const li = lData.x.length - 1;
            const nL = computeNormal(lData, li, -1);
            const dsL = lData.delta_star[li] * scale;
            const p = toCanvas(rotatePoint({ x: lData.x[li] + nL.nx * dsL, y: lData.y[li] + nL.ny * dsL }));
            ctx.lineTo(p.x, p.y);
          }

          ctx.closePath();
          ctx.fill();

          // Stroke edges, bridging from BL TE points
          ctx.strokeStyle = wakeStroke;
          ctx.lineWidth = 1;

          ctx.beginPath();
          if (showBoundaryLayer && blVisData.upper.x.length > 0) {
            const uData = blVisData.upper;
            const li = uData.x.length - 1;
            const nU = computeNormal(uData, li, 1);
            const dsU = uData.delta_star[li] * scale;
            const p = toCanvas(rotatePoint({ x: uData.x[li] + nU.nx * dsU, y: uData.y[li] + nU.ny * dsU }));
            ctx.moveTo(p.x, p.y);
          } else {
            const p = toCanvas(rotatePoint(upperEnv[0]));
            ctx.moveTo(p.x, p.y);
          }
          for (const pt of upperEnv) {
            const p = toCanvas(rotatePoint(pt));
            ctx.lineTo(p.x, p.y);
          }
          ctx.stroke();

          ctx.beginPath();
          if (showBoundaryLayer && blVisData.lower.x.length > 0) {
            const lData = blVisData.lower;
            const li = lData.x.length - 1;
            const nL = computeNormal(lData, li, -1);
            const dsL = lData.delta_star[li] * scale;
            const p = toCanvas(rotatePoint({ x: lData.x[li] + nL.nx * dsL, y: lData.y[li] + nL.ny * dsL }));
            ctx.moveTo(p.x, p.y);
          } else {
            const p = toCanvas(rotatePoint(lowerEnv[0]));
            ctx.moveTo(p.x, p.y);
          }
          for (const pt of lowerEnv) {
            const p = toCanvas(rotatePoint(pt));
            ctx.lineTo(p.x, p.y);
          }
          ctx.stroke();
        }
      }

      // Highlight hovered BL station with a crosshair + ring
      if (blHoverInfo) {
        const hp = toCanvas(rotatePoint({ x: blHoverInfo.x, y: blHoverInfo.y }));
        const r = 6;
        const isUpper = blHoverInfo.surface === 'upper';
        const markerColor = isUpper
          ? (isDark ? '#00e5ff' : '#0090d0')
          : (isDark ? '#ff6eb4' : '#d0306a');

        // Filled ring
        ctx.beginPath();
        ctx.arc(hp.x, hp.y, r, 0, Math.PI * 2);
        ctx.fillStyle = markerColor;
        ctx.globalAlpha = 0.35;
        ctx.fill();
        ctx.globalAlpha = 1;
        ctx.strokeStyle = markerColor;
        ctx.lineWidth = 2;
        ctx.stroke();

        // Crosshair lines
        const arm = 12;
        ctx.beginPath();
        ctx.moveTo(hp.x - arm, hp.y); ctx.lineTo(hp.x + arm, hp.y);
        ctx.moveTo(hp.x, hp.y - arm); ctx.lineTo(hp.x, hp.y + arm);
        ctx.strokeStyle = markerColor;
        ctx.lineWidth = 1;
        ctx.globalAlpha = 0.6;
        ctx.stroke();
        ctx.globalAlpha = 1;
      }
    }

    // Draw displacement thickness (δ*) offset curve -- the "effective inviscid body"
    if (showDisplacementThickness && blVisData && blVisData.success) {
      const scale = blThicknessScale;
      const curveColor = isDark ? 'rgba(0, 255, 200, 0.85)' : 'rgba(0, 160, 120, 0.85)';

      const dsNormal = (pts: { x: number[]; y: number[] }, i: number, sign: number) => {
        let tx: number, ty: number;
        if (pts.x.length < 2) return { nx: 0, ny: sign };
        if (i === 0) { tx = pts.x[1] - pts.x[0]; ty = pts.y[1] - pts.y[0]; }
        else if (i === pts.x.length - 1) { tx = pts.x[i] - pts.x[i - 1]; ty = pts.y[i] - pts.y[i - 1]; }
        else { tx = pts.x[i + 1] - pts.x[i - 1]; ty = pts.y[i + 1] - pts.y[i - 1]; }
        const len = Math.sqrt(tx * tx + ty * ty) || 1e-10;
        return { nx: -ty / len * sign, ny: tx / len * sign };
      };

      ctx.strokeStyle = curveColor;
      ctx.lineWidth = 2;
      ctx.setLineDash([6, 3]);

      // Upper surface δ* curve
      const uData = blVisData.upper;
      if (uData.x.length > 1) {
        ctx.beginPath();
        for (let i = 0; i < uData.x.length; i++) {
          const n = dsNormal(uData, i, 1);
          const ds = uData.delta_star[i] * scale;
          const p = toCanvas(rotatePoint({ x: uData.x[i] + n.nx * ds, y: uData.y[i] + n.ny * ds }));
          if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
        }
        ctx.stroke();
      }

      // Lower surface δ* curve
      const lData = blVisData.lower;
      if (lData.x.length > 1) {
        ctx.beginPath();
        for (let i = 0; i < lData.x.length; i++) {
          const n = dsNormal(lData, i, -1);
          const ds = lData.delta_star[i] * scale;
          const p = toCanvas(rotatePoint({ x: lData.x[i] + n.nx * ds, y: lData.y[i] + n.ny * ds }));
          if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
        }
        ctx.stroke();
      }

      // Wake δ* curves (upper and lower edges)
      if (showWake) {
        const wk = blVisData.wake;
        const fU = blVisData.wake_upper_fraction;
        const fL = 1 - fU;
        if (wk.x.length > 1) {
          const blendN = Math.min(5, wk.x.length);
          const teNU = uData.x.length > 1 ? dsNormal(uData, uData.x.length - 1, 1) : { nx: 0, ny: 1 };
          const teNL = lData.x.length > 1 ? dsNormal(lData, lData.x.length - 1, -1) : { nx: 0, ny: -1 };

          // Wake upper edge
          ctx.beginPath();
          for (let i = 0; i < wk.x.length; i++) {
            const wn = dsNormal(wk, i, 1);
            const t = i < blendN ? i / blendN : 1;
            const nx = teNU.nx * (1 - t) + wn.nx * t;
            const ny = teNU.ny * (1 - t) + wn.ny * t;
            const len = Math.sqrt(nx * nx + ny * ny) || 1;
            const ds = wk.delta_star[i] * scale * fU;
            const p = toCanvas(rotatePoint({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds }));
            if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
          }
          ctx.stroke();

          // Wake lower edge
          ctx.beginPath();
          for (let i = 0; i < wk.x.length; i++) {
            const wn = dsNormal(wk, i, 1);
            const t = i < blendN ? i / blendN : 1;
            const nx = teNL.nx * (1 - t) + (-wn.nx) * t;
            const ny = teNL.ny * (1 - t) + (-wn.ny) * t;
            const len = Math.sqrt(nx * nx + ny * ny) || 1;
            const ds = wk.delta_star[i] * scale * fL;
            const p = toCanvas(rotatePoint({ x: wk.x[i] + (nx / len) * ds, y: wk.y[i] + (ny / len) * ds }));
            if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
          }
          ctx.stroke();
        }
      }

      ctx.setLineDash([]);
    }

  // NOTE: Smoke state (smokePositions, etc.) intentionally NOT in dependencies
  // Smoke is drawn on separate overlay canvas to avoid expensive redraws
  }, [viewport, morphState, splineCurve, controlMode, camberControlPoints, thicknessControlPoints, hoveredPoint, showGrid, showCurve, showPanels, showPoints, showControls, showStreamlines, showPsiContours, psiContours, displayAlpha, toCanvas, isDark, showCp, showForces, cpDisplayMode, cpBarScale, forceScale, showBoundaryLayer, showWake, showDisplacementThickness, blVisData, blThicknessScale, blHoverInfo]);

  // Draw grid
  const drawGrid = useCallback((ctx: CanvasRenderingContext2D) => {
    const { width, height, zoom, center } = viewport;
    
    const gridColor = isDark ? '#333333' : '#dee2e6';
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 0.5;

    // Determine grid spacing based on zoom
    let gridSpacing = 0.1;
    if (zoom < 100) gridSpacing = 0.5;
    if (zoom < 50) gridSpacing = 1.0;
    if (zoom > 500) gridSpacing = 0.05;
    if (zoom > 1000) gridSpacing = 0.01;
    if (zoom > 5000) gridSpacing = 0.002;
    if (zoom > 20000) gridSpacing = 0.0005;
    if (zoom > 100000) gridSpacing = 0.0001;

    // Calculate visible range
    const left = center.x - width / (2 * zoom);
    const right = center.x + width / (2 * zoom);
    const top = center.y + height / (2 * zoom);
    const bottom = center.y - height / (2 * zoom);

    // Draw vertical lines
    const startX = Math.floor(left / gridSpacing) * gridSpacing;
    for (let x = startX; x <= right; x += gridSpacing) {
      const canvasX = toCanvas({ x, y: 0 }).x;
      ctx.beginPath();
      ctx.moveTo(canvasX, 0);
      ctx.lineTo(canvasX, height);
      ctx.stroke();
    }

    // Draw horizontal lines
    const startY = Math.floor(bottom / gridSpacing) * gridSpacing;
    for (let y = startY; y <= top; y += gridSpacing) {
      const canvasY = toCanvas({ x: 0, y }).y;
      ctx.beginPath();
      ctx.moveTo(0, canvasY);
      ctx.lineTo(width, canvasY);
      ctx.stroke();
    }
  }, [viewport, toCanvas, isDark]);

  // Draw axes
  const drawAxes = useCallback((ctx: CanvasRenderingContext2D) => {
    const { width, height } = viewport;
    const origin = toCanvas({ x: 0, y: 0 });

    const axesColor = isDark ? '#a0a0a0' : '#495057';
    ctx.strokeStyle = axesColor;
    ctx.globalAlpha = 0.3;
    ctx.lineWidth = 1;

    // X axis
    ctx.beginPath();
    ctx.moveTo(0, origin.y);
    ctx.lineTo(width, origin.y);
    ctx.stroke();

    // Y axis
    ctx.beginPath();
    ctx.moveTo(origin.x, 0);
    ctx.lineTo(origin.x, height);
    ctx.stroke();
    
    ctx.globalAlpha = 1; // Reset alpha
  }, [viewport, toCanvas, isDark]);

  // Handle resize
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const resizeObserver = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        setViewport((v) => ({ ...v, width, height }));
      }
    });

    resizeObserver.observe(container);
    
    // Initial size
    const rect = container.getBoundingClientRect();
    setViewport((v) => ({ ...v, width: rect.width, height: rect.height }));

    return () => resizeObserver.disconnect();
  }, []);

  // Initialize WebGPU renderer
  useEffect(() => {
    const initGPU = async () => {
      // Check WebGPU support first
      const support = await checkWebGPUSupport();
      setGPUAvailable(support.available);
      
      if (!support.available) {
        console.log('WebGPU not available:', support.info);
        return;
      }
      
      console.log('WebGPU available:', support.info);
      
      // Wait for canvas to be ready
      const canvas = gpuCanvasRef.current;
      if (!canvas) return;
      
      // Initialize renderer
      const renderer = await WebGPURenderer.create(canvas);
      if (renderer) {
        gpuRendererRef.current = renderer;
        renderer.setTheme(isDark);
        console.log('WebGPU renderer initialized');
      }
    };
    
    initGPU();
    
    return () => {
      if (gpuAnimationRef.current) {
        cancelAnimationFrame(gpuAnimationRef.current);
      }
      gpuRendererRef.current?.destroy();
      gpuRendererRef.current = null;
    };
  }, [setGPUAvailable, isDark]);

  // Update WebGPU renderer theme
  useEffect(() => {
    gpuRendererRef.current?.setTheme(isDark);
  }, [isDark]);

  // Update WebGPU geometry when panels change
  useEffect(() => {
    if (!useGPU || !gpuRendererRef.current || panels.length < 10 || !isWasmReady()) return;
    
    try {
      // Compute gamma values from WASM solver (real values, not approximations)
      const gammaResult = computeGamma(panels, displayAlpha, reynolds);
      if (!gammaResult.success || !gammaResult.gamma) {
        console.warn('GPU gamma computation failed:', gammaResult.error);
        return;
      }
      
      const gamma = Array.from(gammaResult.gamma);
      
      // Use psi_0 from solver (dividing streamline value)
      // Fall back to psiContours values if needed
      const psi0 = gammaResult.psi_0 || psiContours.psi0 || 0;
      const psiMin = psiContours.psiMin || -1;
      const psiMax = psiContours.psiMax || 1;
      
      // Update GPU renderer geometry with real solver values
      gpuRendererRef.current.updateGeometry(
        panels,
        gamma,
        displayAlpha,
        streamlineBounds,
        psi0,
        psiMin,
        psiMax
      );
      
    } catch (e) {
      console.error('GPU geometry update failed:', e);
    }
  }, [useGPU, panels, displayAlpha, reynolds, streamlineBounds, psiContours.psi0, psiContours.psiMin, psiContours.psiMax]);

  // Update GPU smoke spawn points (separate from geometry since spawn doesn't need gamma)
  useEffect(() => {
    if (!useGPU || !gpuRendererRef.current) return;
    
    // Create spawn points as vertical line in DISPLAY frame, then rotate to BODY frame
    // This matches CPU smoke - spawn line appears vertical when viewed
    const { yMin, yMax, spawnX } = smokeBounds;
    const yRange = yMax - yMin;
    
    // Rotation center (quarter chord) and angle
    const cx = 0.25;
    const cy = 0;
    const rad = displayAlpha * Math.PI / 180;  // +alpha to rotate from display to body
    const cosA = Math.cos(rad);
    const sinA = Math.sin(rad);
    
    const spawnPoints: { x: number; y: number }[] = [];
    for (let i = 0; i < smokeDensity; i++) {
      // Desired position in DISPLAY frame (vertical line at spawnX)
      const displayX = spawnX;
      const displayY = yMin + (i * yRange) / Math.max(1, smokeDensity - 1);
      
      // Rotate by +alpha around quarter chord to get BODY frame position
      const dx = displayX - cx;
      const dy = displayY - cy;
      const bodyX = cx + dx * cosA - dy * sinA;
      const bodyY = cy + dx * sinA + dy * cosA;
      
      spawnPoints.push({ x: bodyX, y: bodyY });
    }
    gpuRendererRef.current.setSpawnPoints(spawnPoints);
  }, [useGPU, smokeBounds, smokeDensity, displayAlpha]);
  
  // Update GPU smoke wave spacing and timing parameters
  useEffect(() => {
    if (!useGPU || !gpuRendererRef.current) return;
    
    // Calculate spawn interval from wave spacing (distance / speed = time)
    // waveSpacing is in chord lengths, flowSpeed is chord lengths per second
    const spawnInterval = smokeWaveSpacing / flowSpeed;
    
    // Calculate max age based on domain traversal time
    // Particles should live long enough to cross the visible domain
    const { spawnX, exitX } = smokeBounds;
    const domainWidth = exitX - spawnX;
    const traversalTime = domainWidth / flowSpeed;
    const maxAge = traversalTime * 1.3; // 30% margin
    
    gpuRendererRef.current.setSmokeParams(spawnInterval, maxAge);
    gpuRendererRef.current.resetSmoke(); // Reset to redistribute particles
  }, [useGPU, smokeWaveSpacing, flowSpeed, smokeBounds, displayAlpha]);

  // Handle manual smoke reset request (from UI button)
  useEffect(() => {
    if (smokeResetCounter === 0) return; // Skip initial mount
    
    // Reset CPU smoke
    if (smokeSystemRef.current) {
      smokeSystemRef.current.reset();
    }
    
    // Reset GPU smoke
    if (gpuRendererRef.current) {
      gpuRendererRef.current.resetSmoke();
    }
  }, [smokeResetCounter]);

  // WebGPU render loop
  useEffect(() => {
    if (!useGPU || !gpuRendererRef.current) {
      // Clean up any existing GPU animation
      if (gpuAnimationRef.current) {
        cancelAnimationFrame(gpuAnimationRef.current);
        gpuAnimationRef.current = null;
      }
      return;
    }
    
    let lastTime = performance.now();
    
    const animate = (time: number) => {
      const dt = Math.min((time - lastTime) / 1000, 0.05);
      lastTime = time;
      
      if (!gpuRendererRef.current) return;
      
      // Update smoke if enabled
      if (showSmoke) {
        gpuRendererRef.current.updateSmoke(dt * flowSpeed);
      }
      
      // Render
      gpuRendererRef.current.render(
        {
          centerX: viewport.center.x,
          centerY: viewport.center.y,
          zoom: viewport.zoom,
          width: viewport.width,
          height: viewport.height,
        },
        {
          smoke: showSmoke,
          contours: showPsiContours,
          streamlines: showStreamlines,
          alpha: displayAlpha,
          isDark,
        }
      );
      
      // Update performance metrics
      const metrics = gpuRendererRef.current.getPerformanceMetrics();
      updatePerfMetrics(metrics);
      
      gpuAnimationRef.current = requestAnimationFrame(animate);
    };
    
    gpuAnimationRef.current = requestAnimationFrame(animate);
    
    return () => {
      if (gpuAnimationRef.current) {
        cancelAnimationFrame(gpuAnimationRef.current);
        gpuAnimationRef.current = null;
      }
    };
  }, [useGPU, showSmoke, showPsiContours, showStreamlines, flowSpeed, viewport, displayAlpha, isDark, updatePerfMetrics]);

  // Update GPU streamlines when streamlines state changes
  // Use direct streamlines state, not morphState, to avoid timing issues
  useEffect(() => {
    if (!useGPU || !gpuRendererRef.current) return;
    if (!showStreamlines || streamlines.length === 0) return;
    
    // Transform streamlines to display coordinates (apply rotation)
    const rad = displayAlpha * Math.PI / 180;
    const cosA = Math.cos(rad);
    const sinA = Math.sin(rad);
    const cx = 0.25;
    const cy = 0;
    
    const transformedStreamlines: [number, number][][] = streamlines.map(line => 
      line.map(([x, y]) => {
        // Rotate around quarter chord
        const dx = x - cx;
        const dy = y - cy;
        const rx = cx + dx * cosA + dy * sinA;
        const ry = cy - dx * sinA + dy * cosA;
        return [rx, ry] as [number, number];
      })
    );
    
    gpuRendererRef.current.setStreamlines(transformedStreamlines);
  }, [useGPU, showStreamlines, streamlines, displayAlpha]);

  // Redraw when state changes
  useEffect(() => {
    draw();
  }, [draw]);

  // Mouse event handlers
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    lastMousePos.current = canvasPos;

    // Check for point hit
    const hit = findPointAt(canvasPos);
    if (hit) {
      pauseHistory(); // Pause history tracking during drag
      setIsDragging(true);
      setIsDraggingPoint(true); // Disable morph animation during drag
      setDragTarget(hit);
    } else {
      // Start panning
      setIsPanning(true);
    }
  }, [findPointAt]);

  const handleMouseMove = useCallback((e: React.MouseEvent) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilPos = toAirfoil(canvasPos);

    if (isDragging && dragTarget) {
      // Update point position
      if (dragTarget.type === 'camber' && typeof dragTarget.index === 'string') {
        // For camber points, only allow y to change (x is fixed to chord position)
        const cp = camberControlPoints.find(p => p.id === dragTarget.index);
        if (cp) {
          updateCamberControlPoint(dragTarget.index, { y: airfoilPos.y });
        }
      } else if (dragTarget.type === 'thickness' && typeof dragTarget.index === 'string') {
        // For thickness points, only allow t to change (derived from y position)
        const cp = thicknessControlPoints.find(p => p.id === dragTarget.index);
        if (cp) {
          // Get camber at this x position using spline interpolation
          const getCamberAt = (x: number): number => {
            if (camberControlPoints.length === 0) return 0;
            const camberSpline = generateCamberSplineCurve(camberControlPoints, 50);
            // Find closest points and interpolate
            for (let i = 0; i < camberSpline.length - 1; i++) {
              if (x >= camberSpline[i].x && x <= camberSpline[i + 1].x) {
                const t = (camberSpline[i + 1].x - camberSpline[i].x) > 0 
                  ? (x - camberSpline[i].x) / (camberSpline[i + 1].x - camberSpline[i].x)
                  : 0;
                return camberSpline[i].y + t * (camberSpline[i + 1].y - camberSpline[i].y);
              }
            }
            if (x <= camberSpline[0].x) return camberSpline[0].y;
            return camberSpline[camberSpline.length - 1].y;
          };
          const camber = getCamberAt(cp.x);
          // Thickness is the distance from camber line (positive only)
          const newT = Math.max(0, airfoilPos.y - camber);
          updateThicknessControlPoint(dragTarget.index, { t: newT });
        }
      }
      
      // Auto-zoom when dragging near edge of viewport
      // Define margin (in pixels) where we start zooming out
      const edgeMargin = 50;
      const { width, height } = viewport;
      
      // Check if canvas position is near edge
      const nearLeftEdge = canvasPos.x < edgeMargin;
      const nearRightEdge = canvasPos.x > width - edgeMargin;
      const nearTopEdge = canvasPos.y < edgeMargin;
      const nearBottomEdge = canvasPos.y > height - edgeMargin;
      
      if (nearLeftEdge || nearRightEdge || nearTopEdge || nearBottomEdge) {
        // Calculate how far into the margin we are (0 = at edge, 1 = at margin boundary)
        let maxOverlap = 0;
        
        if (nearLeftEdge) maxOverlap = Math.max(maxOverlap, 1 - canvasPos.x / edgeMargin);
        if (nearRightEdge) maxOverlap = Math.max(maxOverlap, 1 - (width - canvasPos.x) / edgeMargin);
        if (nearTopEdge) maxOverlap = Math.max(maxOverlap, 1 - canvasPos.y / edgeMargin);
        if (nearBottomEdge) maxOverlap = Math.max(maxOverlap, 1 - (height - canvasPos.y) / edgeMargin);
        
        // Zoom out proportionally (more aggressive near edge)
        // Zoom factor ranges from 0.99 (at margin) to 0.95 (at edge)
        const zoomFactor = 1 - (0.01 + maxOverlap * 0.04);
        const newZoom = Math.max(50, viewport.zoom * zoomFactor);
        
        // Also pan slightly toward the dragged point to keep it visible
        const panSpeed = 0.002;
        const panX = nearRightEdge ? panSpeed : (nearLeftEdge ? -panSpeed : 0);
        const panY = nearBottomEdge ? -panSpeed : (nearTopEdge ? panSpeed : 0);
        
        setViewport((v) => ({
          ...v,
          zoom: newZoom,
          center: {
            x: v.center.x + panX,
            y: v.center.y + panY,
          },
        }));
      }
    } else if (isPanning) {
      // Pan viewport
      const dx = (canvasPos.x - lastMousePos.current.x) / viewport.zoom;
      const dy = (canvasPos.y - lastMousePos.current.y) / viewport.zoom;
      setViewport((v) => ({
        ...v,
        center: { x: v.center.x - dx, y: v.center.y + dy },
      }));
    } else {
      // Update hover state
      const hit = findPointAt(canvasPos);
      setHoveredPoint(hit);

      // BL station hover: find nearest surface point when BL envelope is visible.
      // Compare in canvas-pixel space so the hit zone is zoom-independent.
      // BL data is in body coords; we rotate + project to canvas, same as rendering.
      if (showBoundaryLayer && blVisData?.success && !hit) {
        const HIT_PX = 20; // pixel radius
        const hitPx2 = HIT_PX * HIT_PX;
        let best: typeof blHoverInfo = null;
        let bestDist = hitPx2;

        // Rotate body → visual (same transform used for drawing)
        const rad = -displayAlpha * Math.PI / 180;
        const cosA = Math.cos(rad);
        const sinA = Math.sin(rad);
        const rcx = 0.25, rcy = 0;

        for (const surface of ['upper', 'lower'] as const) {
          const data = blVisData[surface];
          for (let i = 0; i < data.x.length; i++) {
            // Rotate BL station to visual space, then project to canvas
            let vx = data.x[i], vy = data.y[i];
            if (displayAlpha !== 0) {
              const dx2 = vx - rcx, dy2 = vy - rcy;
              vx = rcx + dx2 * cosA - dy2 * sinA;
              vy = rcy + dx2 * sinA + dy2 * cosA;
            }
            const cx2 = viewport.width / 2 + (vx - viewport.center.x) * viewport.zoom;
            const cy2 = viewport.height / 2 - (vy - viewport.center.y) * viewport.zoom;
            const dpx = cx2 - canvasPos.x;
            const dpy = cy2 - canvasPos.y;
            const d2 = dpx * dpx + dpy * dpy;
            if (d2 < bestDist) {
              bestDist = d2;
              const ue = data.ue[i];
              const dstar = data.delta_star[i];
              const theta = data.theta[i];
              best = {
                surface, index: i,
                canvasX: canvasPos.x, canvasY: canvasPos.y,
                x: data.x[i], y: data.y[i],
                cp: 1 - ue * ue,
                cf: data.cf[i],
                dstar, theta,
                h: theta > 1e-15 ? dstar / theta : 0,
                ue, hk: data.hk[i],
                laminar: data.is_laminar[i],
              };
            }
          }
        }
        setBlHoverInfo(best);
      } else {
        setBlHoverInfo(null);
      }
    }

    lastMousePos.current = canvasPos;
  }, [isDragging, isPanning, dragTarget, viewport, toAirfoil, findPointAt, camberControlPoints, thicknessControlPoints, updateCamberControlPoint, updateThicknessControlPoint, showBoundaryLayer, blVisData, displayAlpha]);

  const handleMouseLeave = useCallback(() => {
    setBlHoverInfo(null);
  }, []);

  const handleMouseUp = useCallback(() => {
    if (isDragging) {
      resumeHistory(); // Resume history tracking after drag
    }
    setIsDragging(false);
    setIsDraggingPoint(false); // Re-enable morph animation
    setIsPanning(false);
    setDragTarget(null);
  }, [isDragging]);

  // Handle wheel events with native listener for proper preventDefault
  const handleWheel = useCallback((e: WheelEvent) => {
    e.preventDefault();
    
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilBefore = toAirfoil(canvasPos);

    // Zoom
    const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
    const newZoom = Math.max(50, Math.min(500000, viewport.zoom * zoomFactor));

    // Adjust center to zoom toward mouse position
    setViewport((v) => {
      const newViewport = { ...v, zoom: newZoom };
      const airfoilAfter = {
        x: v.center.x + (canvasPos.x - v.width / 2) / newZoom,
        y: v.center.y - (canvasPos.y - v.height / 2) / newZoom,
      };
      return {
        ...newViewport,
        center: {
          x: v.center.x + (airfoilBefore.x - airfoilAfter.x),
          y: v.center.y + (airfoilBefore.y - airfoilAfter.y),
        },
      };
    });
  }, [viewport, toAirfoil]);

  // Attach wheel listener with passive: false to allow preventDefault
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    canvas.addEventListener('wheel', handleWheel, { passive: false });
    return () => canvas.removeEventListener('wheel', handleWheel);
  }, [handleWheel]);

  // Double-click to add control point
  const handleDoubleClick = useCallback((e: React.MouseEvent) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilPos = toAirfoil(canvasPos);
    
    // Only add points in camber-spline or thickness-spline modes
    if (controlMode === 'camber-spline') {
      // Check if there's already a point too close (within 0.02 in x)
      const MIN_GAP = 0.02;
      const tooClose = camberControlPoints.some(p => Math.abs(p.x - airfoilPos.x) < MIN_GAP);
      if (tooClose) return;
      
      // Don't allow points outside [0, 1] range
      if (airfoilPos.x < 0 || airfoilPos.x > 1) return;
      
      // Add new camber control point
      const newId = `camber-${Date.now()}`;
      addCamberControlPoint({
        id: newId,
        x: airfoilPos.x,
        y: airfoilPos.y,
      });
    } else if (controlMode === 'thickness-spline') {
      // Check if there's already a point too close (within 0.02 in x)
      const MIN_GAP = 0.02;
      const tooClose = thicknessControlPoints.some(p => Math.abs(p.x - airfoilPos.x) < MIN_GAP);
      if (tooClose) return;
      
      // Don't allow points outside [0, 1] range
      if (airfoilPos.x < 0 || airfoilPos.x > 1) return;
      
      // Get camber at this x position to compute thickness
      const getCamberAt = (x: number): number => {
        if (camberControlPoints.length === 0) return 0;
        const camberSpline = generateCamberSplineCurve(camberControlPoints, 50);
        for (let i = 0; i < camberSpline.length - 1; i++) {
          if (x >= camberSpline[i].x && x <= camberSpline[i + 1].x) {
            const t = (camberSpline[i + 1].x - camberSpline[i].x) > 0 
              ? (x - camberSpline[i].x) / (camberSpline[i + 1].x - camberSpline[i].x)
              : 0;
            return camberSpline[i].y + t * (camberSpline[i + 1].y - camberSpline[i].y);
          }
        }
        if (x <= camberSpline[0].x) return camberSpline[0].y;
        return camberSpline[camberSpline.length - 1].y;
      };
      
      const camber = getCamberAt(airfoilPos.x);
      // Thickness is the distance from camber line (positive only)
      const newT = Math.max(0, airfoilPos.y - camber);
      
      // Add new thickness control point
      const newId = `thickness-${Date.now()}`;
      addThicknessControlPoint({
        id: newId,
        x: airfoilPos.x,
        t: newT,
      });
    }
  }, [controlMode, camberControlPoints, thicknessControlPoints, addCamberControlPoint, addThicknessControlPoint, toAirfoil]);

  // Right-click to remove control point
  const handleContextMenu = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    
    // Check if we clicked on a point
    const hit = findPointAt(canvasPos);
    if (!hit) return;
    
    if (hit.type === 'camber' && typeof hit.index === 'string') {
      // Don't allow removing if only 2 points left
      if (camberControlPoints.length <= 2) return;
      // Don't allow removing endpoints (x=0 or x=1)
      const point = camberControlPoints.find(p => p.id === hit.index);
      if (point && (point.x <= 0.01 || point.x >= 0.99)) return;
      
      removeCamberControlPoint(hit.index);
    } else if (hit.type === 'thickness' && typeof hit.index === 'string') {
      // Don't allow removing if only 2 points left
      if (thicknessControlPoints.length <= 2) return;
      // Don't allow removing endpoints (x=0 or x=1)
      const point = thicknessControlPoints.find(p => p.id === hit.index);
      if (point && (point.x <= 0.01 || point.x >= 0.99)) return;
      
      removeThicknessControlPoint(hit.index);
    }
  }, [findPointAt, camberControlPoints, thicknessControlPoints, removeCamberControlPoint, removeThicknessControlPoint]);

  // Reset view
  const handleResetView = useCallback(() => {
    setViewport((v) => ({
      ...v,
      center: { x: 0.5, y: 0 },
      zoom: Math.min(v.width, v.height) * 0.8,
    }));
  }, []);

  return (
    <div 
      ref={containerRef} 
      className="canvas-container"
      style={{ position: 'relative' }}
    >
      <canvas
        ref={canvasRef}
        width={viewport.width}
        height={viewport.height}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={(e) => { handleMouseUp(); handleMouseLeave(); }}
        onDoubleClick={handleDoubleClick}
        onContextMenu={handleContextMenu}
        style={{ cursor: isDragging ? 'grabbing' : isPanning ? 'grabbing' : (hoveredPoint || blHoverInfo) ? 'pointer' : 'crosshair' }}
      />
      
      {/* Smoke overlay canvas - renders smoke particles without triggering expensive main canvas redraws */}
      <canvas
        ref={smokeCanvasRef}
        width={viewport.width}
        height={viewport.height}
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          pointerEvents: 'none',  // Let clicks pass through to main canvas
          display: useGPU ? 'none' : 'block',  // Hide when using GPU
        }}
      />
      
      {/* WebGPU canvas - GPU-accelerated rendering for smoke, contours, streamlines */}
      <canvas
        ref={gpuCanvasRef}
        width={viewport.width}
        height={viewport.height}
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          pointerEvents: 'none',
          display: useGPU ? 'block' : 'none',  // Show only when using GPU
        }}
      />
      
      {/* BL hover tooltip - positioned away from the airfoil surface */}
      {blHoverInfo && (() => {
        const tooltipW = 195, tooltipH = 185;
        const pad = 10;
        // Place tooltip on the outward side of the surface:
        // upper surface → above the cursor, lower surface → below
        const above = blHoverInfo.surface === 'upper';
        let tx = blHoverInfo.canvasX - tooltipW / 2;
        let ty = above
          ? blHoverInfo.canvasY - tooltipH - 30
          : blHoverInfo.canvasY + 30;
        // Clamp to viewport
        tx = Math.max(pad, Math.min(viewport.width - tooltipW - pad, tx));
        ty = Math.max(pad, Math.min(viewport.height - tooltipH - pad, ty));
        return (
        <div
          style={{
            position: 'absolute',
            left: tx,
            top: ty,
            pointerEvents: 'none',
            zIndex: 20,
            background: 'var(--bg-secondary, rgba(15,20,35,0.92))',
            backdropFilter: 'blur(12px)',
            border: `1px solid ${blHoverInfo.surface === 'upper' ? 'rgba(0,229,255,0.3)' : 'rgba(255,110,180,0.3)'}`,
            borderRadius: '8px',
            padding: '8px 10px',
            fontSize: '11px',
            lineHeight: '1.55',
            color: 'var(--text-primary, #edf2ff)',
            boxShadow: '0 8px 24px rgba(0,0,0,0.35)',
            width: tooltipW,
            fontFamily: 'ui-monospace, SFMono-Regular, Menlo, monospace',
          }}
        >
          <div style={{ fontWeight: 700, marginBottom: 4, letterSpacing: '0.03em', color: blHoverInfo.surface === 'upper' ? 'var(--accent-primary, #61dafb)' : 'var(--accent-secondary, #f472b6)' }}>
            {blHoverInfo.surface === 'upper' ? '▲ Upper' : '▼ Lower'} &nbsp;
            <span style={{ fontWeight: 400, opacity: 0.7 }}>{blHoverInfo.laminar ? 'laminar' : 'turbulent'}</span>
          </div>
          <div>x/c&nbsp;&nbsp;&nbsp;&nbsp;= {blHoverInfo.x.toFixed(4)}</div>
          <div>Cp&nbsp;&nbsp;&nbsp;&nbsp; = {blHoverInfo.cp.toFixed(4)}</div>
          <div>Cf&nbsp;&nbsp;&nbsp;&nbsp; = {(blHoverInfo.cf * 1e3).toFixed(3)} <span style={{ opacity: 0.5 }}>×10⁻³</span></div>
          <div>δ*&nbsp;&nbsp;&nbsp;&nbsp; = {(blHoverInfo.dstar * 1e3).toFixed(3)} <span style={{ opacity: 0.5 }}>×10⁻³</span></div>
          <div>θ&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = {(blHoverInfo.theta * 1e3).toFixed(3)} <span style={{ opacity: 0.5 }}>×10⁻³</span></div>
          <div>H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = {blHoverInfo.h.toFixed(2)}</div>
          <div>Hk&nbsp;&nbsp;&nbsp;&nbsp;= {blHoverInfo.hk.toFixed(2)}</div>
          <div>Ue/V∞ = {blHoverInfo.ue.toFixed(4)}</div>
        </div>
        );
      })()}

      {/* Overlay controls - simplified, full controls in Visualization panel */}
      <div
        style={{
          position: 'absolute',
          bottom: '12px',
          right: '12px',
          display: 'flex',
          gap: '8px',
          alignItems: 'center',
          zIndex: 10,
          pointerEvents: 'auto',
          background: 'var(--bg-secondary)',
          backdropFilter: 'blur(8px)',
          padding: '6px 10px',
          borderRadius: '6px',
          border: '1px solid var(--border-color)',
        }}
      >
        <button 
          onClick={() => openPanel('visualization')} 
          style={{ 
            padding: '4px 8px', 
            fontSize: '11px',
            background: 'var(--accent-primary)',
            color: 'var(--bg-primary)',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
          }}
          title="Open Visualization Options panel"
        >
          Vis Options
        </button>
        <button onClick={handleResetView} style={{ padding: '4px 8px', fontSize: '11px' }}>
          Reset View
        </button>
        <span style={{ 
          background: 'var(--bg-tertiary)', 
          padding: '4px 8px', 
          borderRadius: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
        }}>
          Zoom: {viewport.zoom.toFixed(0)}x
        </span>
      </div>
    </div>
  );
}
