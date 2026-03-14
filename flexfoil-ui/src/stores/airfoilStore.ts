/**
 * Zustand store for airfoil state management
 * 
 * Uses zundo for undo/redo functionality.
 */

import { create, useStore } from 'zustand';
import { temporal, type TemporalState } from 'zundo';
import type { 
  AirfoilState, 
  AirfoilPoint, 
  ControlMode, 
  SolverMode,
  BezierHandle, 
  BSplineControlPoint,
  SpacingKnot,
  Naca4Params,
  PolarPoint,
  SpacingPanelMode,
  SSPInterpolation,
  SSPVisualization,
  CamberControlPoint,
  ThicknessControlPoint,
} from '../types';
import {
  generateNaca4 as wasmGenerateNaca4,
  generateNaca4Xfoil,
  repanelWithSpacingAndCurvature,
  repanelXfoil,
  isWasmReady
} from '../lib/wasm';
import { evaluateBSpline, evaluateBezierCurve } from '../lib/bspline';
import {
  scaleAirfoil,
  createCamberControlPoints,
  createThicknessControlPoints,
  reconstructWithOriginalThickness,
  reconstructWithOriginalCamber,
} from '../lib/airfoilGeometry';
import { 
  syncToUrl, 
  loadFromUrl, 
  parseNacaFromName,
  type UrlState 
} from '../lib/urlState';

/**
 * State that is tracked for undo/redo.
 * We only track meaningful edits, not UI state like displayAlpha.
 */
type TrackedState = Pick<
  AirfoilState,
  | 'name'
  | 'coordinates'
  | 'panels'
  | 'bezierHandles'
  | 'bsplineControlPoints'
  | 'bsplineDegree'
  | 'spacingKnots'
  | 'nPanels'
  | 'curvatureWeight'
  | 'spacingPanelMode'
  | 'sspInterpolation'
  | 'camberControlPoints'
  | 'thicknessControlPoints'
  | 'thicknessScale'
  | 'camberScale'
  | 'baseCoordinates'
>;

// Default NACA 0012 coordinates (simplified)
const DEFAULT_NACA0012: AirfoilPoint[] = [
  { x: 1.0, y: 0.0 },
  { x: 0.9, y: 0.0185 },
  { x: 0.8, y: 0.0352 },
  { x: 0.7, y: 0.0498 },
  { x: 0.6, y: 0.0618 },
  { x: 0.5, y: 0.0704 },
  { x: 0.4, y: 0.0747 },
  { x: 0.3, y: 0.0736 },
  { x: 0.2, y: 0.0654 },
  { x: 0.1, y: 0.0470 },
  { x: 0.05, y: 0.0345 },
  { x: 0.025, y: 0.0252 },
  { x: 0.0, y: 0.0 },
  { x: 0.025, y: -0.0252 },
  { x: 0.05, y: -0.0345 },
  { x: 0.1, y: -0.0470 },
  { x: 0.2, y: -0.0654 },
  { x: 0.3, y: -0.0736 },
  { x: 0.4, y: -0.0747 },
  { x: 0.5, y: -0.0704 },
  { x: 0.6, y: -0.0618 },
  { x: 0.7, y: -0.0498 },
  { x: 0.8, y: -0.0352 },
  { x: 0.9, y: -0.0185 },
  { x: 1.0, y: 0.0 },
];

// Default spacing: cosine-like with fine spacing at LE and TE
// In SSP: higher F = more points (denser), lower F = fewer points (sparser)
// S=0 is TE (start), S=0.5 is LE, S=1 is TE (end)
const DEFAULT_SPACING_KNOTS: SpacingKnot[] = [
  { S: 0, F: 1.5 },    // Dense at TE
  { S: 0.25, F: 0.4 }, // Sparse between TE and LE
  { S: 0.5, F: 1.5 },  // Dense at LE (high curvature region)
  { S: 0.75, F: 0.4 }, // Sparse between LE and TE
  { S: 1, F: 1.5 },    // Dense at TE
];

interface AirfoilStore extends AirfoilState {
  // Actions
  setCoordinates: (coords: AirfoilPoint[]) => void;
  setPanels: (panels: AirfoilPoint[]) => void;
  setControlMode: (mode: ControlMode) => void;
  setName: (name: string) => void;
  setDisplayAlpha: (alpha: number) => void;
  setReynolds: (reynolds: number) => void;
  setMach: (mach: number) => void;
  setNcrit: (ncrit: number) => void;
  setMaxIterations: (maxIterations: number) => void;
  setSolverMode: (mode: SolverMode) => void;
  
  // Point manipulation (legacy, kept for compatibility)
  updatePoint: (index: number, point: AirfoilPoint) => void;
  addPoint: (index: number, point: AirfoilPoint) => void;
  removePoint: (index: number) => void;
  
  // Bezier handles (legacy, kept for compatibility)
  setBezierHandles: (handles: BezierHandle[]) => void;
  updateBezierHandle: (index: number, handle: BezierHandle) => void;
  
  // B-spline control points (legacy, kept for compatibility)
  setBSplineControlPoints: (points: BSplineControlPoint[]) => void;
  updateBSplineControlPoint: (id: string, point: Partial<BSplineControlPoint>) => void;
  addBSplineControlPoint: (point: BSplineControlPoint) => void;
  removeBSplineControlPoint: (id: string) => void;
  setBSplineDegree: (degree: number) => void;
  
  // Camber/thickness scaling (parameters mode)
  setThicknessScale: (scale: number, skipRepanel?: boolean) => void;
  setCamberScale: (scale: number, skipRepanel?: boolean) => void;
  applyScaling: () => void;
  repanelCoordinates: () => void;
  
  // Camber control points (camber-spline mode)
  setCamberControlPoints: (points: CamberControlPoint[]) => void;
  updateCamberControlPoint: (id: string, point: Partial<CamberControlPoint>) => void;
  addCamberControlPoint: (point: CamberControlPoint) => void;
  removeCamberControlPoint: (id: string) => void;
  initializeCamberControlPoints: () => void;
  
  // Thickness control points (thickness-spline mode)
  setThicknessControlPoints: (points: ThicknessControlPoint[]) => void;
  updateThicknessControlPoint: (id: string, point: Partial<ThicknessControlPoint>) => void;
  addThicknessControlPoint: (point: ThicknessControlPoint) => void;
  removeThicknessControlPoint: (id: string) => void;
  initializeThicknessControlPoints: () => void;
  
  // Spacing
  setSpacingKnots: (knots: SpacingKnot[]) => void;
  updateSpacingKnot: (index: number, knot: SpacingKnot) => void;
  addSpacingKnot: (knot: SpacingKnot) => void;
  removeSpacingKnot: (index: number) => void;
  setNPanels: (n: number) => void;
  setCurvatureWeight: (weight: number) => void;
  setSpacingPanelMode: (mode: SpacingPanelMode) => void;
  setSSPInterpolation: (interp: SSPInterpolation) => void;
  setSSPVisualization: (viz: SSPVisualization) => void;
  
  // NACA generation
  generateNaca4: (params: Naca4Params) => void;
  
  // Repaneling
  repanel: () => void;
  repanelWithXfoil: () => void;
  
  // Polar data
  setPolarData: (data: PolarPoint[]) => void;
  addPolarPoint: (point: PolarPoint) => void;
  clearPolar: () => void;
  
  // Reset
  reset: () => void;
  
  // Initialize default airfoil (call after WASM ready)
  initializeDefaultAirfoil: () => void;
}

/**
 * Deep equality check for arrays of objects (used by temporal middleware)
 */
function deepEqual(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  if (typeof a !== typeof b) return false;
  if (typeof a !== 'object' || a === null || b === null) return false;
  
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return false;
    return a.every((item, i) => deepEqual(item, b[i]));
  }
  
  const keysA = Object.keys(a as object);
  const keysB = Object.keys(b as object);
  if (keysA.length !== keysB.length) return false;
  
  return keysA.every(key => 
    deepEqual((a as Record<string, unknown>)[key], (b as Record<string, unknown>)[key])
  );
}

export const useAirfoilStore = create<AirfoilStore>()(
  temporal(
    (set) => ({
      // Initial state
      name: 'NACA 0012',
      coordinates: DEFAULT_NACA0012,
      panels: DEFAULT_NACA0012,
      controlMode: 'parameters',  // Default to parameters mode
      bezierHandles: [],
      bsplineControlPoints: [],
      bsplineDegree: 3,
      spacingKnots: DEFAULT_SPACING_KNOTS,
      nPanels: 160,  // XFOIL's default NPAN
      curvatureWeight: 0,
      displayAlpha: 0,
      reynolds: 1e6,
      mach: 0,
      ncrit: 9,
      maxIterations: 100,
      solverMode: 'viscous',
      polarData: [],
      spacingPanelMode: 'simple',  // Default to simple curvature-based
      sspInterpolation: 'linear',  // Default to linear (Mark Drela's original)
      sspVisualization: 'plot',    // Default to S-F plot view
      
      // Camber/thickness editing state
      camberControlPoints: [],
      thicknessControlPoints: [],
      thicknessScale: 1.0,
      camberScale: 1.0,
      baseCoordinates: DEFAULT_NACA0012,

      // Actions
      setCoordinates: (coords) => set({ coordinates: coords }),
      setPanels: (panels) => set({ panels }),
      setControlMode: (mode) => set({ controlMode: mode }),
      setName: (name) => set({ name }),
      setDisplayAlpha: (alpha) => set({ displayAlpha: alpha }),
      setReynolds: (reynolds) => set({ reynolds }),
      setMach: (mach) => set({ mach: Math.max(0, Math.min(0.8, mach)) }),
      setNcrit: (ncrit) => set({ ncrit: Math.max(1, Math.min(14, ncrit)) }),
      setMaxIterations: (maxIterations) => set({ maxIterations: Math.max(10, Math.min(500, maxIterations)) }),
      setSolverMode: (solverMode) => set({ solverMode }),

      updatePoint: (index, point) => set((state) => {
        const newCoords = [...state.coordinates];
        newCoords[index] = point;
        // In surface mode, also update panels to show immediate feedback
        return { coordinates: newCoords, panels: newCoords };
      }),

      addPoint: (index, point) => set((state) => {
        const newCoords = [...state.coordinates];
        newCoords.splice(index, 0, point);
        return { coordinates: newCoords };
      }),

      removePoint: (index) => set((state) => {
        if (state.coordinates.length <= 3) return state;
        const newCoords = state.coordinates.filter((_, i) => i !== index);
        return { coordinates: newCoords };
      }),

      setBezierHandles: (handles) => set(() => {
        // When initially setting handles, don't re-evaluate - keep original shape
        // The curve will only change when handles are explicitly moved
        return { bezierHandles: handles };
      }),
      
      updateBezierHandle: (index, handle) => set((state) => {
        const newHandles = [...state.bezierHandles];
        newHandles[index] = handle;
        // Re-evaluate Bezier curve when handles change
        if (state.coordinates.length > 0) {
          const newPanels = evaluateBezierCurve(state.coordinates, newHandles, state.nPanels + 1);
          return { bezierHandles: newHandles, panels: newPanels };
        }
        return { bezierHandles: newHandles };
      }),

      setBSplineControlPoints: (points) => set((state) => {
        // When setting control points, also evaluate the B-spline curve
        if (points.length >= 2) {
          const newPanels = evaluateBSpline(points, state.bsplineDegree, state.nPanels + 1);
          return { bsplineControlPoints: points, panels: newPanels, coordinates: newPanels };
        }
        return { bsplineControlPoints: points };
      }),
      
      updateBSplineControlPoint: (id, point) => set((state) => {
        const newPoints = state.bsplineControlPoints.map((p) =>
          p.id === id ? { ...p, ...point } : p
        );
        // Re-evaluate B-spline curve when control points change
        if (newPoints.length >= 2) {
          const newPanels = evaluateBSpline(newPoints, state.bsplineDegree, state.nPanels + 1);
          return { bsplineControlPoints: newPoints, panels: newPanels, coordinates: newPanels };
        }
        return { bsplineControlPoints: newPoints };
      }),

      addBSplineControlPoint: (point) => set((state) => ({
        bsplineControlPoints: [...state.bsplineControlPoints, point],
      })),

      removeBSplineControlPoint: (id) => set((state) => ({
        bsplineControlPoints: state.bsplineControlPoints.filter((p) => p.id !== id),
      })),

      setBSplineDegree: (degree) => set((state) => {
        const newDegree = Math.max(1, Math.min(5, degree));
        // Re-evaluate B-spline curve when degree changes
        if (state.bsplineControlPoints.length >= 2) {
          const newPanels = evaluateBSpline(state.bsplineControlPoints, newDegree, state.nPanels + 1);
          return { bsplineDegree: newDegree, panels: newPanels, coordinates: newPanels };
        }
        return { bsplineDegree: newDegree };
      }),

      setSpacingKnots: (knots) => set({ spacingKnots: knots }),
      
      updateSpacingKnot: (index, knot) => set((state) => {
        const newKnots = [...state.spacingKnots];
        newKnots[index] = knot;
        return { spacingKnots: newKnots };
      }),

      addSpacingKnot: (knot) => set((state) => {
        const newKnots = [...state.spacingKnots, knot].sort((a, b) => a.S - b.S);
        return { spacingKnots: newKnots };
      }),

      removeSpacingKnot: (index) => set((state) => {
        if (state.spacingKnots.length <= 2) return state;
        if (index === 0 || index === state.spacingKnots.length - 1) return state;
        return { spacingKnots: state.spacingKnots.filter((_, i) => i !== index) };
      }),

      setNPanels: (n) => set({ nPanels: Math.max(10, Math.min(500, n)) }),

      setCurvatureWeight: (weight) => set({ curvatureWeight: Math.max(0, Math.min(1, weight)) }),

      setSpacingPanelMode: (mode) => set({ spacingPanelMode: mode }),
      
      setSSPInterpolation: (interp) => set({ sspInterpolation: interp }),
      
      setSSPVisualization: (viz) => set({ sspVisualization: viz }),

      // Thickness/camber scaling actions
      // Pattern: immediate coordinate update for responsiveness, debounced repanel for flow
      setThicknessScale: (scale, skipRepanel = false) => set((state) => {
        const newScale = Math.max(0.1, Math.min(3.0, scale));
        // Apply scaling to base coordinates
        if (state.baseCoordinates.length > 0) {
          // scaleAirfoil uses high-resolution cosine decomposition to preserve LE
          const scaledCoords = scaleAirfoil(state.baseCoordinates, newScale, state.camberScale);
          
          // Skip repanel if requested (for real-time slider dragging)
          if (skipRepanel) {
            return { 
              thicknessScale: newScale, 
              coordinates: scaledCoords,
              // Keep existing panels - will be updated by debounced repanel
            };
          }
          
          // Auto-repanel after warping for proper curvature-based distribution
          let panels = scaledCoords;
          if (isWasmReady()) {
            try {
              const repaneled = repanelXfoil(scaledCoords, state.nPanels);
              if (repaneled.length > 0) {
                panels = repaneled.map(pt => ({ x: pt.x, y: pt.y }));
              }
            } catch (e) {
              console.warn('Auto-repanel after thickness scale failed:', e);
            }
          }
          
          return { 
            thicknessScale: newScale, 
            coordinates: scaledCoords,
            panels,
          };
        }
        return { thicknessScale: newScale };
      }),
      
      setCamberScale: (scale, skipRepanel = false) => set((state) => {
        const newScale = Math.max(0, Math.min(3.0, scale));
        // Apply scaling to base coordinates
        if (state.baseCoordinates.length > 0) {
          // scaleAirfoil uses high-resolution cosine decomposition to preserve LE
          const scaledCoords = scaleAirfoil(state.baseCoordinates, state.thicknessScale, newScale);
          
          // Skip repanel if requested (for real-time slider dragging)
          if (skipRepanel) {
            return { 
              camberScale: newScale, 
              coordinates: scaledCoords,
              // Keep existing panels - will be updated by debounced repanel
            };
          }
          
          // Auto-repanel after warping for proper curvature-based distribution
          let panels = scaledCoords;
          if (isWasmReady()) {
            try {
              const repaneled = repanelXfoil(scaledCoords, state.nPanels);
              if (repaneled.length > 0) {
                panels = repaneled.map(pt => ({ x: pt.x, y: pt.y }));
              }
            } catch (e) {
              console.warn('Auto-repanel after camber scale failed:', e);
            }
          }
          
          return { 
            camberScale: newScale, 
            coordinates: scaledCoords,
            panels,
          };
        }
        return { camberScale: newScale };
      }),
      
      // Trigger repaneling of current coordinates (for debounced use after slider changes)
      repanelCoordinates: () => set((state) => {
        if (state.coordinates.length < 5 || !isWasmReady()) {
          return { panels: state.coordinates };
        }
        try {
          const repaneled = repanelXfoil(state.coordinates, state.nPanels);
          if (repaneled.length > 0) {
            return { panels: repaneled.map(pt => ({ x: pt.x, y: pt.y })) };
          }
        } catch (e) {
          console.warn('Repanel failed:', e);
        }
        return { panels: state.coordinates };
      }),
      
      applyScaling: () => set((state) => {
        // Save current scaled airfoil as the new base
        return {
          baseCoordinates: [...state.coordinates],
          thicknessScale: 1.0,
          camberScale: 1.0,
        };
      }),
      
      // Camber control point actions
      setCamberControlPoints: (points) => set({ camberControlPoints: points }),
      
      updateCamberControlPoint: (id, point) => set((state) => {
        const newPoints = state.camberControlPoints.map((p) =>
          p.id === id ? { ...p, ...point } : p
        );
        // Reconstruct airfoil using the new camber line with CURRENT thickness.
        // Using state.coordinates (not baseCoordinates) preserves any thickness
        // modifications the user has made in thickness-spline mode.
        // The thickness is extracted via dense decomposition (200+ cosine samples),
        // which accurately captures any modifications including the LE profile.
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          const newCoords = reconstructWithOriginalThickness(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          
          // Auto-repanel after warping for proper curvature-based distribution
          let panels: AirfoilPoint[] = newCoords;
          if (isWasmReady()) {
            try {
              const repaneled = repanelXfoil(newCoords, state.nPanels);
              if (repaneled.length > 0) {
                panels = repaneled.map(pt => ({ x: pt.x, y: pt.y }));
              }
            } catch (e) {
              // Silent fail - use unrepaneled coords
            }
          }
          
          return { 
            camberControlPoints: newPoints, 
            coordinates: newCoords,
            panels,
          };
        }
        return { camberControlPoints: newPoints };
      }),
      
      addCamberControlPoint: (point) => set((state) => {
        const newPoints = [...state.camberControlPoints, point].sort((a, b) => a.x - b.x);
        // Reconstruct using current thickness (preserves any thickness modifications)
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          const newCoords = reconstructWithOriginalThickness(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          return { 
            camberControlPoints: newPoints, 
            coordinates: newCoords,
            panels: newCoords,
          };
        }
        return { camberControlPoints: newPoints };
      }),
      
      removeCamberControlPoint: (id) => set((state) => {
        if (state.camberControlPoints.length <= 2) return state;
        const newPoints = state.camberControlPoints.filter((p) => p.id !== id);
        // Reconstruct using current thickness (preserves any thickness modifications)
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          const newCoords = reconstructWithOriginalThickness(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          return { 
            camberControlPoints: newPoints, 
            coordinates: newCoords,
            panels: newCoords,
          };
        }
        return { camberControlPoints: newPoints };
      }),
      
      initializeCamberControlPoints: () => set((state) => {
        const points = createCamberControlPoints(state.coordinates, 7);
        return { camberControlPoints: points };
      }),
      
      // Thickness control point actions
      setThicknessControlPoints: (points) => set({ thicknessControlPoints: points }),
      
      updateThicknessControlPoint: (id, point) => set((state) => {
        const newPoints = state.thicknessControlPoints.map((p) =>
          p.id === id ? { ...p, ...point } : p
        );
        // Reconstruct airfoil with new thickness, preserving CURRENT camber.
        // Using current coordinates (not baseCoordinates) preserves any camber
        // modifications the user has made in camber-spline mode.
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          // Use camber from current coordinates (may have been modified)
          const newCoords = reconstructWithOriginalCamber(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          
          // Auto-repanel after warping for proper curvature-based distribution
          let panels: AirfoilPoint[] = newCoords;
          if (isWasmReady()) {
            try {
              const repaneled = repanelXfoil(newCoords, state.nPanels);
              if (repaneled.length > 0) {
                panels = repaneled.map(pt => ({ x: pt.x, y: pt.y }));
              }
            } catch (e) {
              // Silent fail - use unrepaneled coords
            }
          }
          
          return { 
            thicknessControlPoints: newPoints, 
            coordinates: newCoords,
            panels,
          };
        }
        return { thicknessControlPoints: newPoints };
      }),
      
      addThicknessControlPoint: (point) => set((state) => {
        const newPoints = [...state.thicknessControlPoints, point].sort((a, b) => a.x - b.x);
        // Reconstruct using current camber (preserves any camber modifications)
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          const newCoords = reconstructWithOriginalCamber(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          return { 
            thicknessControlPoints: newPoints, 
            coordinates: newCoords,
            panels: newCoords,
          };
        }
        return { thicknessControlPoints: newPoints };
      }),
      
      removeThicknessControlPoint: (id) => set((state) => {
        if (state.thicknessControlPoints.length <= 2) return state;
        const newPoints = state.thicknessControlPoints.filter((p) => p.id !== id);
        // Reconstruct using current camber (preserves any camber modifications)
        if (newPoints.length >= 2 && state.coordinates.length >= 5) {
          const newCoords = reconstructWithOriginalCamber(
            newPoints,
            state.coordinates,
            state.nPanels + 1
          );
          return { 
            thicknessControlPoints: newPoints, 
            coordinates: newCoords,
            panels: newCoords,
          };
        }
        return { thicknessControlPoints: newPoints };
      }),
      
      initializeThicknessControlPoints: () => set((state) => {
        const points = createThicknessControlPoints(state.coordinates, 7);
        return { thicknessControlPoints: points };
      }),

      generateNaca4: (params) => {
        const { m, p, t, nPoints } = params;
        
        // Generate name (params are already fractions: m=0.02, p=0.4, t=0.12)
        const mInt = Math.round(m * 100);
        const pInt = Math.round(p * 10);
        const tInt = Math.round(t * 100);
        const name = `NACA ${mInt}${pInt}${tInt.toString().padStart(2, '0')}`;
        
        // Use WASM if available, otherwise fallback to JS
        if (isWasmReady()) {
          try {
            // XFOIL procedure:
            // 1. NACA4 generates buffer coordinates (XB/YB)
            // 2. SCALC + SEGSPL fit splines to buffer
            // 3. PANGEN creates working coordinates (X/Y) from splined buffer
            //
            // We replicate this exactly:
            // 1. Generate buffer with XFOIL's exact NACA generator
            const designation = mInt * 1000 + pInt * 100 + tInt;  // e.g., 0012 -> 12, 2412 -> 2412
            const bufferCoords = generateNaca4Xfoil(designation, 123);  // 245 buffer points
            const buffer: AirfoilPoint[] = bufferCoords.map(pt => ({ x: pt.x, y: pt.y }));
            
            // 2. Apply XFOIL paneling (PANGEN) - fits spline and distributes panels
            // Use nPanels from current state, defaulting to 160 (XFOIL's NPAN default)
            const currentNPanels = useAirfoilStore.getState().nPanels || 160;
            const paneledCoords = repanelXfoil(buffer, currentNPanels);
            const panels: AirfoilPoint[] = paneledCoords.map(pt => ({ x: pt.x, y: pt.y }));
            
            set({ 
              coordinates: buffer,  // Store original buffer as "coordinates"
              panels: panels,       // Store paneled result as "panels" for analysis
              name,
              bezierHandles: [],
              bsplineControlPoints: [],
              baseCoordinates: buffer,  // Store as base for scaling
              thicknessScale: 1.0,
              camberScale: 1.0,
              camberControlPoints: [],
              thicknessControlPoints: [],
            });
            return;
          } catch (e) {
            console.warn('WASM XFOIL NACA generation failed, trying generic:', e);
            try {
              // Fallback to generic NACA generator (no auto-paneling)
              const wasmCoords = wasmGenerateNaca4(mInt, pInt, tInt, nPoints);
              const coords: AirfoilPoint[] = wasmCoords.map(pt => ({ x: pt.x, y: pt.y }));
              
              set({ 
                coordinates: coords, 
                panels: coords,
                name,
                bezierHandles: [],
                bsplineControlPoints: [],
                baseCoordinates: coords,
                thicknessScale: 1.0,
                camberScale: 1.0,
                camberControlPoints: [],
                thicknessControlPoints: [],
              });
              return;
            } catch (e2) {
              console.warn('WASM NACA generation failed, using JS fallback:', e2);
            }
          }
        }
        
        // JavaScript fallback
        const coords: AirfoilPoint[] = [];
        const halfPoints = Math.floor(nPoints / 2);
        
        for (let i = 0; i <= halfPoints; i++) {
          const beta = (Math.PI * i) / halfPoints;
          const x = (1 - Math.cos(beta)) / 2;
          
          const yt = 5 * t * (
            0.2969 * Math.sqrt(x) 
            - 0.1260 * x 
            - 0.3516 * x * x 
            + 0.2843 * x * x * x 
            - 0.1015 * x * x * x * x
          );
          
          let yc = 0;
          let dyc = 0;
          if (m > 0 && p > 0) {
            if (x < p) {
              yc = (m / (p * p)) * (2 * p * x - x * x);
              dyc = (2 * m / (p * p)) * (p - x);
            } else {
              yc = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
              dyc = (2 * m / ((1 - p) * (1 - p))) * (p - x);
            }
          }
          
          const theta = Math.atan(dyc);
          
          if (i < halfPoints) {
            const xu = x - yt * Math.sin(theta);
            const yu = yc + yt * Math.cos(theta);
            coords.unshift({ x: xu, y: yu, surface: 'upper' });
          }
          
          if (i === halfPoints) {
            coords.unshift({ x: 0, y: yc });
          }
          
          if (i > 0) {
            const xl = x + yt * Math.sin(theta);
            const yl = yc - yt * Math.cos(theta);
            coords.push({ x: xl, y: yl, surface: 'lower' });
          }
        }
        
        set({ 
          coordinates: coords, 
          panels: coords,
          name,
          bezierHandles: [],
          bsplineControlPoints: [],
          baseCoordinates: coords,
          thicknessScale: 1.0,
          camberScale: 1.0,
          camberControlPoints: [],
          thicknessControlPoints: [],
        });
      },

      repanel: () => set((state) => {
        if (!isWasmReady()) {
          return state;
        }
        
        try {
          // Convert spacing knots to wasm format
          const spacingKnots = state.spacingKnots.map(k => ({ s: k.S, f: k.F }));
          const newPanels = repanelWithSpacingAndCurvature(
            state.coordinates,
            spacingKnots,
            state.nPanels,
            state.curvatureWeight
          );
          
          if (newPanels.length === 0) {
            return state;
          }
          
          return { panels: newPanels.map(pt => ({ x: pt.x, y: pt.y })) };
        } catch (e) {
          console.error('Repaneling failed:', e);
          return state;
        }
      }),

      repanelWithXfoil: () => set((state) => {
        if (!isWasmReady()) {
          return state;
        }
        
        try {
          // Use XFOIL's PANGEN algorithm with curvature-based paneling
          const newPanels = repanelXfoil(state.coordinates, state.nPanels);
          
          if (newPanels.length === 0) {
            return state;
          }
          
          return { panels: newPanels.map(pt => ({ x: pt.x, y: pt.y })) };
        } catch (e) {
          console.error('XFOIL repaneling failed:', e);
          return state;
        }
      }),
      
      // Polar data actions
      setPolarData: (data) => set({ polarData: data }),
      addPolarPoint: (point) => set((state) => ({
        polarData: [...state.polarData, point]
      })),
      clearPolar: () => set({ polarData: [] }),

      reset: () => set({
        name: 'NACA 0012',
        coordinates: DEFAULT_NACA0012,
        panels: DEFAULT_NACA0012,
        controlMode: 'parameters',
        bezierHandles: [],
        bsplineControlPoints: [],
        bsplineDegree: 3,
        spacingKnots: DEFAULT_SPACING_KNOTS,
        nPanels: 160,  // XFOIL's default NPAN
        curvatureWeight: 0,
        displayAlpha: 0,
        mach: 0,
        ncrit: 9,
        maxIterations: 100,
        polarData: [],
        spacingPanelMode: 'simple',
        sspInterpolation: 'linear',
        sspVisualization: 'plot',
        camberControlPoints: [],
        thicknessControlPoints: [],
        thicknessScale: 1.0,
        camberScale: 1.0,
        baseCoordinates: DEFAULT_NACA0012,
      }),
      
      initializeDefaultAirfoil: () => {
        if (!isWasmReady()) {
          console.warn('WASM not ready, cannot initialize default airfoil');
          return;
        }
        
        try {
          // Generate NACA 0012 using XFOIL's exact algorithm
          const bufferCoords = generateNaca4Xfoil(12, 123);
          const buffer: AirfoilPoint[] = bufferCoords.map(pt => ({ x: pt.x, y: pt.y }));
          
          // Apply XFOIL paneling for proper analysis
          const currentNPanels = useAirfoilStore.getState().nPanels || 160;
          const paneledCoords = repanelXfoil(buffer, currentNPanels);
          const panels: AirfoilPoint[] = paneledCoords.map(pt => ({ x: pt.x, y: pt.y }));
          
          set({
            name: 'NACA 0012',
            coordinates: buffer,
            panels: panels,
            bezierHandles: [],
            bsplineControlPoints: [],
            baseCoordinates: buffer,
            thicknessScale: 1.0,
            camberScale: 1.0,
            camberControlPoints: [],
            thicknessControlPoints: [],
          });
        } catch (e) {
          console.error('Failed to initialize default airfoil:', e);
        }
      },
    }),
    {
      // Temporal middleware options
      limit: 100, // Maximum history size
      
      // Only track meaningful state changes (not UI state like displayAlpha, controlMode, polarData)
      partialize: (state): TrackedState => ({
        name: state.name,
        coordinates: state.coordinates,
        panels: state.panels,
        bezierHandles: state.bezierHandles,
        bsplineControlPoints: state.bsplineControlPoints,
        bsplineDegree: state.bsplineDegree,
        spacingKnots: state.spacingKnots,
        nPanels: state.nPanels,
        curvatureWeight: state.curvatureWeight,
        spacingPanelMode: state.spacingPanelMode,
        sspInterpolation: state.sspInterpolation,
        camberControlPoints: state.camberControlPoints,
        thicknessControlPoints: state.thicknessControlPoints,
        thicknessScale: state.thicknessScale,
        camberScale: state.camberScale,
        baseCoordinates: state.baseCoordinates,
      }),
      
      // Deep equality check to avoid duplicate history entries
      equality: (pastState, currentState) => deepEqual(pastState, currentState),
    }
  )
);

/**
 * Access the temporal store for undo/redo operations.
 * 
 * Usage:
 *   const canUndo = useTemporalStore((state) => state.pastStates.length > 0);
 *   const canRedo = useTemporalStore((state) => state.futureStates.length > 0);
 */
export const useTemporalStore = <T>(
  selector: (state: TemporalState<TrackedState>) => T
): T => {
  return useStore(useAirfoilStore.temporal, selector);
};

/**
 * Pause history tracking (useful during drag operations).
 * Call this at the START of a drag.
 */
export function pauseHistory(): void {
  useAirfoilStore.temporal.getState().pause();
}

/**
 * Resume history tracking and optionally create a single history entry.
 * Call this at the END of a drag.
 */
export function resumeHistory(): void {
  useAirfoilStore.temporal.getState().resume();
}

/**
 * Perform undo operation.
 */
export function undo(): void {
  useAirfoilStore.temporal.getState().undo();
}

/**
 * Perform redo operation.
 */
export function redo(): void {
  useAirfoilStore.temporal.getState().redo();
}

/**
 * Clear all history (past and future states).
 */
export function clearHistory(): void {
  useAirfoilStore.temporal.getState().clear();
}

/**
 * Get current state for URL encoding
 * Note: Viewport and visualization state are added separately by the components
 */
export function getUrlState(): UrlState {
  const state = useAirfoilStore.getState();
  const naca = parseNacaFromName(state.name);
  
  return {
    naca: naca || undefined,
    custom: !naca,
    nPanels: state.nPanels,
    spacing: state.spacingKnots,
    mode: state.controlMode,
    alpha: state.displayAlpha,
  };
}

/**
 * Hydrate store from URL state
 */
export function hydrateFromUrl(): boolean {
  const urlState = loadFromUrl();
  if (!urlState) return false;
  
  const store = useAirfoilStore.getState();
  
  // Apply NACA airfoil if specified
  if (urlState.naca && !urlState.custom) {
    const designation = parseInt(urlState.naca, 10);
    const m = Math.floor(designation / 1000) / 100;
    const p = Math.floor((designation % 1000) / 100) / 10;
    const t = (designation % 100) / 100;
    
    // Generate after WASM is ready
    setTimeout(() => {
      if (isWasmReady()) {
        store.generateNaca4({ m, p, t, nPoints: 100 });
      }
    }, 100);
  }
  
  // Apply paneling settings
  if (urlState.nPanels) {
    store.setNPanels(urlState.nPanels);
  }
  if (urlState.spacing) {
    store.setSpacingKnots(urlState.spacing);
  }
  
  // Apply control mode
  if (urlState.mode) {
    store.setControlMode(urlState.mode);
  }
  
  return true;
}

/**
 * Sync store state to URL (debounced)
 */
let syncTimeout: ReturnType<typeof setTimeout> | null = null;

export function syncStoreToUrl(): void {
  if (syncTimeout) {
    clearTimeout(syncTimeout);
  }
  
  syncTimeout = setTimeout(() => {
    const urlState = getUrlState();
    syncToUrl(urlState);
  }, 500); // Debounce 500ms
}

/**
 * Subscribe to store changes and sync to URL
 */
export function subscribeToUrlSync(): () => void {
  return useAirfoilStore.subscribe((state, prevState) => {
    // Only sync on meaningful changes
    if (
      state.name !== prevState.name ||
      state.nPanels !== prevState.nPanels ||
      state.controlMode !== prevState.controlMode ||
      state.spacingKnots !== prevState.spacingKnots
    ) {
      syncStoreToUrl();
    }
  });
}
