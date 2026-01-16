/**
 * Zustand store for airfoil state management
 */

import { create } from 'zustand';
import type { 
  AirfoilState, 
  AirfoilPoint, 
  ControlMode, 
  BezierHandle, 
  BSplineControlPoint,
  SpacingKnot,
  Naca4Params 
} from '../types';
import {
  generateNaca4 as wasmGenerateNaca4,
  generateNaca4Xfoil,
  repanelWithSpacingAndCurvature,
  repanelXfoil,
  isWasmReady
} from '../lib/wasm';
import { evaluateBSpline, evaluateBezierCurve } from '../lib/bspline';

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
  
  // Point manipulation
  updatePoint: (index: number, point: AirfoilPoint) => void;
  addPoint: (index: number, point: AirfoilPoint) => void;
  removePoint: (index: number) => void;
  
  // Bezier handles
  setBezierHandles: (handles: BezierHandle[]) => void;
  updateBezierHandle: (index: number, handle: BezierHandle) => void;
  
  // B-spline control points
  setBSplineControlPoints: (points: BSplineControlPoint[]) => void;
  updateBSplineControlPoint: (id: string, point: Partial<BSplineControlPoint>) => void;
  addBSplineControlPoint: (point: BSplineControlPoint) => void;
  removeBSplineControlPoint: (id: string) => void;
  setBSplineDegree: (degree: number) => void;
  
  // Spacing
  setSpacingKnots: (knots: SpacingKnot[]) => void;
  updateSpacingKnot: (index: number, knot: SpacingKnot) => void;
  addSpacingKnot: (knot: SpacingKnot) => void;
  removeSpacingKnot: (index: number) => void;
  setNPanels: (n: number) => void;
  setCurvatureWeight: (weight: number) => void;
  
  // NACA generation
  generateNaca4: (params: Naca4Params) => void;
  
  // Repaneling
  repanel: () => void;
  repanelWithXfoil: () => void;
  
  // Reset
  reset: () => void;
}

export const useAirfoilStore = create<AirfoilStore>((set) => ({
  // Initial state
  name: 'NACA 0012',
  coordinates: DEFAULT_NACA0012,
  panels: DEFAULT_NACA0012,
  controlMode: 'surface',
  bezierHandles: [],
  bsplineControlPoints: [],
  bsplineDegree: 3,
  spacingKnots: DEFAULT_SPACING_KNOTS,
  nPanels: 160,  // XFOIL's default NPAN
  curvatureWeight: 0,

  // Actions
  setCoordinates: (coords) => set({ coordinates: coords }),
  setPanels: (panels) => set({ panels }),
  setControlMode: (mode) => set({ controlMode: mode }),
  setName: (name) => set({ name }),

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

  reset: () => set({
    name: 'NACA 0012',
    coordinates: DEFAULT_NACA0012,
    panels: DEFAULT_NACA0012,
    controlMode: 'surface',
    bezierHandles: [],
    bsplineControlPoints: [],
    bsplineDegree: 3,
    spacingKnots: DEFAULT_SPACING_KNOTS,
    nPanels: 160,  // XFOIL's default NPAN
    curvatureWeight: 0,
  }),
}));
