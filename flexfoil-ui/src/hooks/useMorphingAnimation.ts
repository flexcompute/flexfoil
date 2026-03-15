/**
 * Morphing Animation Hook (WASM-Accelerated)
 * 
 * Provides smooth interpolation between airfoil states for animated transitions.
 * All visualizations (airfoil shape, streamlines, Cp, forces) morph together.
 * 
 * Performance optimizations:
 * - WASM-based interpolation for coordinates, panels, Cp, cpX
 * - Pre-allocated buffers for streamlines (still JS for encoding overhead)
 * - Ref-based animation state (avoids React renders during animation)
 * - Fast paths for same-length arrays
 * - Batched state updates
 */

import { useRef, useEffect, useState, useCallback, useMemo } from 'react';
import type { AirfoilPoint } from '../types';
import { 
    isWasmReady, 
    wasmLerpMorphState, 
    wasmLerpStreamlines,
} from '../lib/wasm';

export interface MorphState {
  /** Interpolated airfoil coordinates */
  coordinates: AirfoilPoint[];
  /** Interpolated panel coordinates */
  panels: AirfoilPoint[];
  /** Interpolated streamlines */
  streamlines: [number, number][][];
  /** Interpolated Cp values */
  cp: number[];
  /** Interpolated Cp x positions */
  cpX: number[];
  /** Interpolated lift coefficient */
  cl: number;
  /** Interpolated drag coefficient */
  cd: number;
  /** Interpolated moment coefficient */
  cm: number;
  /** Whether animation is in progress */
  isAnimating: boolean;
  /** Animation progress (0 to 1) */
  progress: number;
}

interface MorphTarget {
  coordinates: AirfoilPoint[];
  panels: AirfoilPoint[];
  streamlines: [number, number][][];
  cp: number[];
  cpX: number[];
  cl: number;
  cd: number;
  cm: number;
}

// Easing function for smooth animation
function easeInOutCubic(t: number): number {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
}

// Interpolate between two numbers (inline for speed)
const lerp = (a: number, b: number, t: number): number => a + (b - a) * t;

/**
 * Fast path: interpolate same-length point arrays in-place
 * Avoids object creation by reusing result array
 */
function lerpPointsSameLength(
  from: AirfoilPoint[],
  to: AirfoilPoint[],
  t: number,
  result: AirfoilPoint[]
): void {
  const len = to.length;
  // Ensure result array is correct size
  if (result.length !== len) {
    result.length = len;
    for (let i = 0; i < len; i++) {
      result[i] = { x: 0, y: 0 };
    }
  }
  for (let i = 0; i < len; i++) {
    result[i].x = from[i].x + (to[i].x - from[i].x) * t;
    result[i].y = from[i].y + (to[i].y - from[i].y) * t;
  }
}

/**
 * General case: interpolate arrays of different lengths
 * Uses index mapping with linear interpolation
 */
function lerpPointsDiffLength(
  from: AirfoilPoint[],
  to: AirfoilPoint[],
  t: number,
  result: AirfoilPoint[]
): void {
  const fromLen = from.length;
  const toLen = to.length;
  
  if (fromLen === 0) {
    result.length = toLen;
    for (let i = 0; i < toLen; i++) {
      result[i] = { x: to[i].x, y: to[i].y };
    }
    return;
  }
  if (toLen === 0) {
    result.length = fromLen;
    for (let i = 0; i < fromLen; i++) {
      result[i] = { x: from[i].x, y: from[i].y };
    }
    return;
  }
  
  // Ensure result array is correct size
  if (result.length !== toLen) {
    result.length = toLen;
    for (let i = 0; i < toLen; i++) {
      result[i] = { x: 0, y: 0 };
    }
  }
  
  const scale = (fromLen - 1) / (toLen - 1);
  for (let i = 0; i < toLen; i++) {
    const srcIdx = i * scale;
    const srcIdxLow = srcIdx | 0; // Fast floor
    const srcIdxHigh = Math.min(srcIdxLow + 1, fromLen - 1);
    const srcT = srcIdx - srcIdxLow;
    
    // Interpolate source point
    const srcX = from[srcIdxLow].x + (from[srcIdxHigh].x - from[srcIdxLow].x) * srcT;
    const srcY = from[srcIdxLow].y + (from[srcIdxHigh].y - from[srcIdxLow].y) * srcT;
    
    // Interpolate between source and target
    result[i].x = srcX + (to[i].x - srcX) * t;
    result[i].y = srcY + (to[i].y - srcY) * t;
  }
}

/**
 * Fast path: interpolate same-length number arrays
 */
function lerpArraySameLength(
  from: number[],
  to: number[],
  t: number,
  result: number[]
): void {
  const len = to.length;
  if (result.length !== len) result.length = len;
  for (let i = 0; i < len; i++) {
    result[i] = from[i] + (to[i] - from[i]) * t;
  }
}

/**
 * General case: interpolate number arrays of different lengths
 */
function lerpArrayDiffLength(
  from: number[],
  to: number[],
  t: number,
  result: number[]
): void {
  const fromLen = from.length;
  const toLen = to.length;
  
  if (fromLen === 0) {
    result.length = toLen;
    for (let i = 0; i < toLen; i++) result[i] = to[i];
    return;
  }
  if (toLen === 0) {
    result.length = fromLen;
    for (let i = 0; i < fromLen; i++) result[i] = from[i];
    return;
  }
  
  if (result.length !== toLen) result.length = toLen;
  
  const scale = (fromLen - 1) / (toLen - 1);
  for (let i = 0; i < toLen; i++) {
    const srcIdx = i * scale;
    const srcIdxLow = srcIdx | 0;
    const srcIdxHigh = Math.min(srcIdxLow + 1, fromLen - 1);
    const srcT = srcIdx - srcIdxLow;
    
    const srcVal = from[srcIdxLow] + (from[srcIdxHigh] - from[srcIdxLow]) * srcT;
    result[i] = srcVal + (to[i] - srcVal) * t;
  }
}

/**
 * Interpolate streamlines - the most expensive operation (JS fallback)
 * Optimized with in-place updates and fast paths
 */
function lerpStreamlinesJS(
  from: [number, number][][],
  to: [number, number][][],
  t: number,
  result: [number, number][][]
): void {
  const fromLen = from.length;
  const toLen = to.length;
  
  if (fromLen === 0) {
    result.length = toLen;
    for (let i = 0; i < toLen; i++) {
      result[i] = to[i].slice(); // Shallow copy
    }
    return;
  }
  if (toLen === 0) {
    result.length = fromLen;
    for (let i = 0; i < fromLen; i++) {
      result[i] = from[i].slice();
    }
    return;
  }
  
  const maxLen = Math.max(fromLen, toLen);
  if (result.length !== maxLen) result.length = maxLen;
  
  for (let i = 0; i < maxLen; i++) {
    const fromLine = from[i % fromLen];
    const toLine = to[i % toLen];
    
    if (!fromLine || fromLine.length === 0) {
      result[i] = toLine ? toLine.slice() : [];
      continue;
    }
    if (!toLine || toLine.length === 0) {
      result[i] = fromLine.slice();
      continue;
    }
    
    const fromPts = fromLine.length;
    const toPts = toLine.length;
    
    // Ensure result line exists and has correct length
    if (!result[i]) result[i] = [];
    const lineResult = result[i];
    if (lineResult.length !== toPts) {
      lineResult.length = toPts;
      for (let j = 0; j < toPts; j++) {
        lineResult[j] = [0, 0];
      }
    }
    
    if (fromPts === toPts) {
      // Fast path: same length
      for (let j = 0; j < toPts; j++) {
        lineResult[j][0] = fromLine[j][0] + (toLine[j][0] - fromLine[j][0]) * t;
        lineResult[j][1] = fromLine[j][1] + (toLine[j][1] - fromLine[j][1]) * t;
      }
    } else {
      // Different lengths: need index mapping
      const scale = (fromPts - 1) / (toPts - 1);
      for (let j = 0; j < toPts; j++) {
        const srcIdx = j * scale;
        const srcIdxLow = srcIdx | 0;
        const srcIdxHigh = Math.min(srcIdxLow + 1, fromPts - 1);
        const srcT = srcIdx - srcIdxLow;
        
        const srcX = fromLine[srcIdxLow][0] + (fromLine[srcIdxHigh][0] - fromLine[srcIdxLow][0]) * srcT;
        const srcY = fromLine[srcIdxLow][1] + (fromLine[srcIdxHigh][1] - fromLine[srcIdxLow][1]) * srcT;
        
        lineResult[j][0] = srcX + (toLine[j][0] - srcX) * t;
        lineResult[j][1] = srcY + (toLine[j][1] - srcY) * t;
      }
    }
  }
}

export interface UseMorphingAnimationOptions {
  /** Animation duration in milliseconds */
  duration?: number;
  /** Whether morphing is enabled */
  enabled?: boolean;
}

/**
 * Hook for morphing animation between airfoil states.
 * Optimized for 60fps rendering with minimal GC pressure.
 */
export function useMorphingAnimation(
  target: MorphTarget,
  options: UseMorphingAnimationOptions = {}
): MorphState {
  const { duration = 300, enabled = true } = options;
  
  // Store previous state for interpolation
  const prevStateRef = useRef<MorphTarget | null>(null);
  const animationRef = useRef<number | null>(null);
  const startTimeRef = useRef<number>(0);
  
  // Pre-allocated buffers for interpolation (avoids GC during animation)
  const buffersRef = useRef({
    coordinates: [] as AirfoilPoint[],
    panels: [] as AirfoilPoint[],
    streamlines: [] as [number, number][][],
    cp: [] as number[],
    cpX: [] as number[],
  });
  
  // Frame throttling - only update React state every N frames
  // This dramatically reduces React reconciliation overhead
  const frameCountRef = useRef(0);
  const FRAMES_PER_UPDATE = 2; // Update every 2 frames (~30fps state updates)
  
  // Current interpolated state (only updated at key frames or end)
  const [morphState, setMorphState] = useState<MorphState>(() => ({
    coordinates: target.coordinates,
    panels: target.panels,
    streamlines: target.streamlines,
    cp: target.cp,
    cpX: target.cpX,
    cl: target.cl,
    cd: target.cd,
    cm: target.cm,
    isAnimating: false,
    progress: 1,
  }));
  
  // Memoize target identity for change detection
  const targetId = useMemo(() => ({
    coordLen: target.coordinates.length,
    panelLen: target.panels.length,
    streamLen: target.streamlines.length,
    cpLen: target.cp.length,
    cl: target.cl,
    cd: target.cd,
    cm: target.cm,
    // Sample a few points for change detection
    firstCoord: target.coordinates[0],
    lastCoord: target.coordinates[target.coordinates.length - 1],
  }), [target.coordinates, target.panels, target.streamlines, target.cp, target.cl, target.cd, target.cm]);
  
  // Check if target has changed significantly
  const hasTargetChanged = useCallback((prev: MorphTarget | null, curr: MorphTarget): boolean => {
    if (!prev) return false;
    
    // Check array length changes (cheap)
    if (prev.coordinates.length !== curr.coordinates.length) return true;
    if (prev.panels.length !== curr.panels.length) return true;
    if (prev.streamlines.length !== curr.streamlines.length) return true;
    if (prev.cp.length !== curr.cp.length) return true;
    
    // Check Cl/Cm change (indicates alpha changed)
    if (Math.abs(prev.cl - curr.cl) > 0.0005) return true;
    if (Math.abs(prev.cd - curr.cd) > 0.0005) return true;
    if (Math.abs(prev.cm - curr.cm) > 0.0005) return true;
    
    // Check sample points for significant change
    if (prev.coordinates.length > 0 && curr.coordinates.length > 0) {
      const threshold = 0.0005;
      const firstDiff = Math.abs(prev.coordinates[0].x - curr.coordinates[0].x) +
                       Math.abs(prev.coordinates[0].y - curr.coordinates[0].y);
      const lastIdx = prev.coordinates.length - 1;
      const lastDiff = Math.abs(prev.coordinates[lastIdx].x - curr.coordinates[lastIdx].x) +
                      Math.abs(prev.coordinates[lastIdx].y - curr.coordinates[lastIdx].y);
      
      if (firstDiff > threshold || lastDiff > threshold) return true;
    }
    
    return false;
  }, []);
  
  // Animation frame callback - WASM-accelerated interpolation
  const animate = useCallback((timestamp: number) => {
    const elapsed = timestamp - startTimeRef.current;
    const rawProgress = Math.min(elapsed / duration, 1);
    const progress = easeInOutCubic(rawProgress);
    const isComplete = rawProgress >= 1;
    
    const prev = prevStateRef.current;
    if (!prev) {
      setMorphState({
        ...target,
        isAnimating: false,
        progress: 1,
      });
      return;
    }
    
    // Frame throttling - only update React state periodically
    // Always update on completion for final state accuracy
    frameCountRef.current++;
    const shouldUpdateState = isComplete || (frameCountRef.current % FRAMES_PER_UPDATE === 0);
    
    if (shouldUpdateState) {
      let interpolated: MorphState;
      
      // Use WASM for interpolation if available (much faster)
      if (isWasmReady()) {
        // Batch interpolation of coordinates, panels, cp, cpX in single WASM call
        const morphed = wasmLerpMorphState(
          prev.coordinates,
          target.coordinates,
          prev.panels,
          target.panels,
          prev.cp,
          target.cp,
          prev.cpX,
          target.cpX,
          progress
        );
        
        // Streamlines - also via WASM
        const streamlines = wasmLerpStreamlines(
          prev.streamlines,
          target.streamlines,
          progress
        );
        
        interpolated = {
          coordinates: morphed.coordinates,
          panels: morphed.panels,
          streamlines,
          cp: morphed.cp,
          cpX: morphed.cpX,
          cl: lerp(prev.cl, target.cl, progress),
          cd: lerp(prev.cd, target.cd, progress),
          cm: lerp(prev.cm, target.cm, progress),
          isAnimating: !isComplete,
          progress,
        };
      } else {
        // Fallback to JS implementation if WASM not ready
        const buffers = buffersRef.current;
        
        // Use fast paths when array lengths match (common case)
        if (prev.coordinates.length === target.coordinates.length) {
          lerpPointsSameLength(prev.coordinates, target.coordinates, progress, buffers.coordinates);
        } else {
          lerpPointsDiffLength(prev.coordinates, target.coordinates, progress, buffers.coordinates);
        }
        
        if (prev.panels.length === target.panels.length) {
          lerpPointsSameLength(prev.panels, target.panels, progress, buffers.panels);
        } else {
          lerpPointsDiffLength(prev.panels, target.panels, progress, buffers.panels);
        }
        
        if (prev.cp.length === target.cp.length) {
          lerpArraySameLength(prev.cp, target.cp, progress, buffers.cp);
        } else {
          lerpArrayDiffLength(prev.cp, target.cp, progress, buffers.cp);
        }
        
        if (prev.cpX.length === target.cpX.length) {
          lerpArraySameLength(prev.cpX, target.cpX, progress, buffers.cpX);
        } else {
          lerpArrayDiffLength(prev.cpX, target.cpX, progress, buffers.cpX);
        }
        
        // Streamlines - most expensive, but optimized
        lerpStreamlinesJS(prev.streamlines, target.streamlines, progress, buffers.streamlines);
        
        interpolated = {
          coordinates: buffers.coordinates.map(p => ({ x: p.x, y: p.y })),
          panels: buffers.panels.map(p => ({ x: p.x, y: p.y })),
          streamlines: buffers.streamlines.map(line => line.map(pt => [pt[0], pt[1]] as [number, number])),
          cp: [...buffers.cp],
          cpX: [...buffers.cpX],
          cl: lerp(prev.cl, target.cl, progress),
          cd: lerp(prev.cd, target.cd, progress),
          cm: lerp(prev.cm, target.cm, progress),
          isAnimating: !isComplete,
          progress,
        };
      }
      
      setMorphState(interpolated);
    }
    
    if (!isComplete) {
      animationRef.current = requestAnimationFrame(animate);
    } else {
      // Animation complete, update prev state
      prevStateRef.current = target;
      frameCountRef.current = 0;
    }
  }, [target, duration]);
  
  // Start animation when target changes
  useEffect(() => {
    // Cancel any existing animation
    if (animationRef.current) {
      cancelAnimationFrame(animationRef.current);
      animationRef.current = null;
    }
    
    if (!enabled) {
      // No morphing, just set state directly
      setMorphState({
        ...target,
        isAnimating: false,
        progress: 1,
      });
      prevStateRef.current = target;
      return;
    }
    
    // First render - set initial state without animation
    if (!prevStateRef.current) {
      prevStateRef.current = target;
      setMorphState({
        ...target,
        isAnimating: false,
        progress: 1,
      });
      return;
    }
    
    // Check if we should animate (significant change)
    if (hasTargetChanged(prevStateRef.current, target)) {
      // Start new animation from current state
      startTimeRef.current = performance.now();
      animationRef.current = requestAnimationFrame(animate);
    } else {
      // Minor change or no change - still update state but don't animate
      setMorphState({
        ...target,
        isAnimating: false,
        progress: 1,
      });
      prevStateRef.current = target;
    }
    
    return () => {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [targetId, enabled, animate, hasTargetChanged]);
  
  return morphState;
}

/**
 * Get color for Cp value visualization.
 * Blue = negative Cp (suction) - upper surface typically
 * Red = positive Cp (pressure) - stagnation point, lower surface
 * Color intensity increases with magnitude (similar to stream function)
 */
export function getCpColor(cp: number, isDark: boolean): string {
  // Typical range: Cp from -3 (strong suction) to +1 (stagnation)
  const cpClamped = Math.max(-4, Math.min(1.5, cp));
  
  if (cpClamped < 0) {
    // Negative Cp (suction) -> Blue (upper surface typically)
    const t = Math.min(1, Math.abs(cpClamped) / 3); // Normalize to ~0-1 for Cp in [-3, 0]
    const intensity = Math.pow(t, 0.5);
    
    if (isDark) {
      return `rgb(${Math.round(180 - intensity * 140)}, ${Math.round(200 - intensity * 130)}, ${Math.round(255 - intensity * 25)})`;
    } else {
      return `rgb(${Math.round(150 - intensity * 120)}, ${Math.round(180 - intensity * 130)}, ${Math.round(255 - intensity * 55)})`;
    }
  } else {
    // Positive Cp (pressure/stagnation) -> Red (stagnation, lower surface)
    const t = Math.min(1, cpClamped / 1.0); // Normalize to ~0-1 for Cp in [0, 1]
    const intensity = Math.pow(t, 0.5);
    
    if (isDark) {
      return `rgb(${Math.round(255 - intensity * 25)}, ${Math.round(180 - intensity * 120)}, ${Math.round(170 - intensity * 110)})`;
    } else {
      return `rgb(${Math.round(255 - intensity * 35)}, ${Math.round(160 - intensity * 120)}, ${Math.round(150 - intensity * 110)})`;
    }
  }
}

/**
 * Compute lift and drag vectors from Cl.
 * Note: Panel method gives Cl but not Cd (inviscid).
 * Returns vectors in display coordinates (where freestream is horizontal).
 */
export function computeForceVectors(
  cl: number,
  _alphaDeg: number,
  scale: number = 0.1
): { lift: { x: number; y: number }; drag: { x: number; y: number } } {
  const liftMag = cl * scale;
  
  // Lift is perpendicular to freestream.
  // In the display, the airfoil is rotated so freestream appears horizontal.
  // Therefore lift should point straight up (positive y) for positive Cl.
  const lift = {
    x: 0,
    y: liftMag,
  };
  
  // Drag parallel to freestream (inviscid = 0)
  // For panel method, we don't have real drag
  const drag = {
    x: 0,
    y: 0,
  };
  
  return { lift, drag };
}
