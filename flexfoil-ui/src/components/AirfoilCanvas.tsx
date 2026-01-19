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
import { useAirfoilStore, pauseHistory, resumeHistory } from '../stores/airfoilStore';
import { useVisualizationStore } from '../stores/visualizationStore';
import { useTheme } from '../contexts/ThemeContext';
import { useLayout } from '../contexts/LayoutContext';
import type { Point, ViewportState, AirfoilPoint } from '../types';
import { computeStreamlines, createSmokeSystem, isWasmReady, WasmSmokeSystem, analyzeAirfoil } from '../lib/wasm';
import { useMorphingAnimation, getCpColor, computeForceVectors } from '../hooks/useMorphingAnimation';
import { generateCamberSplineCurve } from '../lib/airfoilGeometry';

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

// Constants
const PANEL_POINT_RADIUS = 2.5;
const CONTROL_RADIUS = 6;
const HIT_RADIUS = 10;

export function AirfoilCanvas() {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  // Theme context
  const { isDark } = useTheme();
  
  // Layout context for panel actions
  const { openPanel } = useLayout();
  
  // State from store
  const { 
    coordinates, 
    panels, 
    controlMode,
    displayAlpha,
    // Camber/thickness control
    camberControlPoints,
    thicknessControlPoints,
    updateCamberControlPoint,
    updateThicknessControlPoint,
    addCamberControlPoint,
    removeCamberControlPoint,
    addThicknessControlPoint,
    removeThicknessControlPoint,
  } = useAirfoilStore();

  // Viewport state
  const [viewport, setViewport] = useState<ViewportState>({
    center: { x: 0.5, y: 0 },
    zoom: 400,
    width: 800,
    height: 600,
  });

  // Visualization settings from store
  const {
    showGrid,
    showCurve,
    showPanels,
    showPoints,
    showControls,
    showStreamlines,
    showSmoke,
    showCp,
    showForces,
    enableMorphing,
    morphDuration,
    streamlineDensity,
    adaptiveStreamlines,
    smokeDensity,
    smokeParticlesPerBlob,
    smokeSpawnInterval,
    smokeMaxAge,
    flowSpeed,
    cpDisplayMode,
    cpBarScale,
    forceScale,
  } = useVisualizationStore();
  
  // Streamlines cache
  const [streamlines, setStreamlines] = useState<[number, number][][]>([]);
  
  // Analysis results (Cp, Cl, Cm)
  const [analysisResult, setAnalysisResult] = useState<{
    cp: number[];
    cpX: number[];
    cl: number;
    cm: number;
  }>({ cp: [], cpX: [], cl: 0, cm: 0 });
  
  // Track last zoom for adaptive streamlines
  const lastZoomRef = useRef(viewport.zoom);
  
  // Smoke system
  const smokeSystemRef = useRef<WasmSmokeSystem | null>(null);
  const smokeAnimationRef = useRef<number | null>(null);
  const [smokePositions, setSmokePositions] = useState<Float64Array | null>(null);
  const [smokeAlphas, setSmokeAlphas] = useState<Float64Array | null>(null);
  
  // Compute adaptive streamline count based on zoom
  const getAdaptiveStreamlineCount = useCallback(() => {
    if (!adaptiveStreamlines) {
      return streamlineDensity;
    }
    // Increase count when zoomed in, decrease when zoomed out
    const zoomFactor = Math.log10(viewport.zoom / 100 + 1);
    return Math.min(100, Math.floor(streamlineDensity * (1 + zoomFactor * 0.5)));
  }, [adaptiveStreamlines, streamlineDensity, viewport.zoom]);
  
  // Track if we're actively dragging (for disabling expensive computations during drag)
  const [isDraggingPoint, setIsDraggingPoint] = useState(false);
  
  // Compute bounds for streamline domain - use large fixed bounds since Rust is fast
  const getVisibleBounds = useCallback((): [number, number, number, number] => {
    // Large domain that covers most reasonable viewing areas
    // Streamlines will fill the entire screen regardless of zoom/pan
    return [-2.0, 4.0, -2.0, 2.0];
  }, []);

  // Compute streamlines when enabled and alpha/panels/zoom change
  // SKIP during drag to prevent freezing
  useEffect(() => {
    if (!showStreamlines || !isWasmReady() || panels.length < 10 || isDraggingPoint) {
      if (!isDraggingPoint) {
        setStreamlines([]);
      }
      return;
    }
    
    // Check if zoom changed significantly for adaptive mode
    const zoomChanged = adaptiveStreamlines && 
      Math.abs(viewport.zoom - lastZoomRef.current) / lastZoomRef.current > 0.2;
    
    if (zoomChanged) {
      lastZoomRef.current = viewport.zoom;
    }
    
    try {
      const seedCount = getAdaptiveStreamlineCount();
      const bounds = getVisibleBounds();
      const result = computeStreamlines(panels, displayAlpha, seedCount, bounds);
      if (result.success) {
        setStreamlines(result.streamlines);
      }
    } catch (e) {
      console.error('Streamline computation failed:', e);
    }
  }, [showStreamlines, panels, displayAlpha, streamlineDensity, adaptiveStreamlines, viewport.zoom, getAdaptiveStreamlineCount, getVisibleBounds, isDraggingPoint]);
  
  // Compute aerodynamic analysis (Cp, Cl, Cm)
  // SKIP during drag to prevent freezing - recalculate on drag end
  useEffect(() => {
    if (!isWasmReady() || panels.length < 10 || isDraggingPoint) {
      return;
    }
    
    try {
      const result = analyzeAirfoil(panels, displayAlpha);
      if (result.success) {
        setAnalysisResult({
          cp: result.cp,
          cpX: result.cp_x,
          cl: result.cl,
          cm: result.cm,
        });
      }
    } catch (e) {
      console.error('Analysis failed:', e);
    }
  }, [panels, displayAlpha, isDraggingPoint]);
  
  // Morphing animation integration
  const morphTarget = useMemo(() => ({
    coordinates,
    panels,
    streamlines,
    cp: analysisResult.cp,
    cpX: analysisResult.cpX,
    cl: analysisResult.cl,
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
  
  // Get setSmokeDensity for adaptive performance
  const { setSmokeDensity } = useVisualizationStore();
  
  // Store panels for smoke system - only update when NOT dragging
  const smokePanelsRef = useRef(panels);
  useEffect(() => {
    if (!isDraggingPoint) {
      smokePanelsRef.current = panels;
    }
  }, [panels, isDraggingPoint]);
  
  // Initialize smoke system (only when toggled on or settings change)
  // SKIP recreation during drag - use cached panels
  useEffect(() => {
    if (!showSmoke || !isWasmReady() || smokePanelsRef.current.length < 10) {
      if (smokeAnimationRef.current) {
        cancelAnimationFrame(smokeAnimationRef.current);
        smokeAnimationRef.current = null;
      }
      smokeSystemRef.current = null;
      setSmokePositions(null);
      setSmokeAlphas(null);
      return;
    }
    
    // Don't recreate smoke system during drag
    if (isDraggingPoint && smokeSystemRef.current) {
      return;
    }
    
    // Create smoke system with spawn points based on smokeDensity
    // Use large Y range to fill screen (-1.5 to 1.5) and spawn further back (-1.5)
    const spawnY = Array.from({ length: smokeDensity }, (_, i) => 
      -1.5 + (i * 3.0) / Math.max(1, smokeDensity - 1)
    );
    smokeSystemRef.current = createSmokeSystem(spawnY, -1.5, smokeParticlesPerBlob);
    smokeSystemRef.current.set_spawn_interval(smokeSpawnInterval);
    smokeSystemRef.current.set_max_age(smokeMaxAge);
    
    // Set initial flow using cached panels
    smokeSystemRef.current.set_flow(
      new Float64Array(smokePanelsRef.current.flatMap(p => [p.x, p.y])),
      displayAlpha
    );
    
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
        setSmokePositions(new Float64Array(pos));
        setSmokeAlphas(new Float64Array(alphas));
      }
      
      smokeAnimationRef.current = requestAnimationFrame(animate);
    };
    
    smokeAnimationRef.current = requestAnimationFrame(animate);
    
    return () => {
      if (smokeAnimationRef.current) {
        cancelAnimationFrame(smokeAnimationRef.current);
      }
    };
  }, [showSmoke, smokeDensity, smokeParticlesPerBlob, smokeSpawnInterval, smokeMaxAge, flowSpeed, setSmokeDensity, isDraggingPoint, displayAlpha]);
  
  // Update flow when alpha changes (without recreating the system)
  // SKIP during drag
  useEffect(() => {
    if (smokeSystemRef.current && smokePanelsRef.current.length >= 10 && !isDraggingPoint) {
      smokeSystemRef.current.set_flow(
        new Float64Array(smokePanelsRef.current.flatMap(p => [p.x, p.y])),
        displayAlpha
      );
    }
  }, [displayAlpha, isDraggingPoint]);

  // Interaction state
  const [isDragging, setIsDragging] = useState(false);
  const [isPanning, setIsPanning] = useState(false);
  const [dragTarget, setDragTarget] = useState<{ type: string; index: number | string } | null>(null);
  const [hoveredPoint, setHoveredPoint] = useState<{ type: string; index: number | string } | null>(null);
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
    
    // Draw streamlines (behind airfoil) - rotated by alpha to match airfoil
    // Uses morphed streamlines for smooth animation
    if (showStreamlines && morphState.streamlines.length > 0) {
      ctx.strokeStyle = colors.accentSecondary;
      ctx.globalAlpha = 0.5;
      ctx.lineWidth = 1;
      
      for (const line of morphState.streamlines) {
        if (line.length < 2) continue;
        ctx.beginPath();
        const first = toCanvas(rotatePoint({ x: line[0][0], y: line[0][1] }));
        ctx.moveTo(first.x, first.y);
        
        for (let i = 1; i < line.length; i++) {
          const p = toCanvas(rotatePoint({ x: line[i][0], y: line[i][1] }));
          ctx.lineTo(p.x, p.y);
        }
        ctx.stroke();
      }
      ctx.globalAlpha = 1; // Reset alpha
    }
    
    // Draw smoke particles (behind airfoil) - rotated by alpha to match airfoil
    if (showSmoke && smokePositions && smokeAlphas) {
      const count = smokePositions.length / 2;
      for (let i = 0; i < count; i++) {
        const x = smokePositions[i * 2];
        const y = smokePositions[i * 2 + 1];
        const alpha = smokeAlphas[i] || 0;
        
        if (alpha > 0.01) {
          const pCanvas = toCanvas(rotatePoint({ x, y }));
          ctx.beginPath();
          ctx.fillStyle = colors.accentSecondary;
          ctx.globalAlpha = alpha * 0.6;
          ctx.arc(pCanvas.x, pCanvas.y, 2, 0, Math.PI * 2);
          ctx.fill();
        }
      }
      ctx.globalAlpha = 1; // Reset alpha
    }


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
    if (showCp && morphState.cp.length > 0 && morphState.cpX.length > 0) {
      const cpLen = Math.min(morphState.cp.length, morphState.cpX.length);
      
      // Find corresponding y values on the airfoil for each Cp point
      // cpX are x positions along the chord
      for (let i = 0; i < cpLen; i++) {
        const cpX = morphState.cpX[i];
        const cpVal = morphState.cp[i];
        
        // Find the closest panel point to get y
        let closestPt = morphedPanels[0] || { x: 0, y: 0 };
        let minDist = Infinity;
        for (const pt of morphedPanels) {
          const dist = Math.abs(pt.x - cpX);
          if (dist < minDist) {
            minDist = dist;
            closestPt = pt;
          }
        }
        
        const cpColor = getCpColor(cpVal, isDark);
        
        // Draw Cp bar (perpendicular to surface)
        if (cpDisplayMode === 'bars' || cpDisplayMode === 'both') {
          const barLength = -cpVal * cpBarScale * viewport.zoom; // Negative Cp = suction = outward
          const pCanvas = toCanvas(rotatePoint(closestPt));
          
          // Compute normal direction (approximate)
          const nextIdx = Math.min(i + 1, morphedPanels.length - 1);
          const prevIdx = Math.max(i - 1, 0);
          const dx = morphedPanels[nextIdx].x - morphedPanels[prevIdx].x;
          const dy = morphedPanels[nextIdx].y - morphedPanels[prevIdx].y;
          const len = Math.sqrt(dx * dx + dy * dy) || 1;
          const nx = -dy / len; // Normal points outward
          const ny = dx / len;
          
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

  }, [viewport, morphState, splineCurve, controlMode, camberControlPoints, thicknessControlPoints, hoveredPoint, showGrid, showCurve, showPanels, showPoints, showControls, showStreamlines, showSmoke, smokePositions, smokeAlphas, displayAlpha, toCanvas, isDark, showCp, showForces, cpDisplayMode, cpBarScale, forceScale]);

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
    }

    lastMousePos.current = canvasPos;
  }, [isDragging, isPanning, dragTarget, viewport, toAirfoil, findPointAt, camberControlPoints, thicknessControlPoints, updateCamberControlPoint, updateThicknessControlPoint]);

  const handleMouseUp = useCallback(() => {
    if (isDragging) {
      resumeHistory(); // Resume history tracking after drag
    }
    setIsDragging(false);
    setIsDraggingPoint(false); // Re-enable morph animation
    setIsPanning(false);
    setDragTarget(null);
  }, [isDragging]);

  const handleWheel = useCallback((e: React.WheelEvent) => {
    e.preventDefault();
    
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilBefore = toAirfoil(canvasPos);

    // Zoom
    const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
    const newZoom = Math.max(50, Math.min(5000, viewport.zoom * zoomFactor));

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
        onMouseLeave={handleMouseUp}
        onDoubleClick={handleDoubleClick}
        onContextMenu={handleContextMenu}
        onWheel={handleWheel}
        style={{ cursor: isDragging ? 'grabbing' : isPanning ? 'grabbing' : hoveredPoint ? 'pointer' : 'crosshair' }}
      />
      
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
