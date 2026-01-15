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
  repanelWithSpacing,
  isWasmReady 
} from '../lib/wasm';

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

const DEFAULT_SPACING_KNOTS: SpacingKnot[] = [
  { S: 0, F: 0.5 },   // Dense at TE
  { S: 0.5, F: 1.5 }, // Coarser at mid
  { S: 1, F: 0.5 },   // Dense at TE again
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
  
  // NACA generation
  generateNaca4: (params: Naca4Params) => void;
  
  // Repaneling
  repanel: () => void;
  
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
  nPanels: 50,

  // Actions
  setCoordinates: (coords) => set({ coordinates: coords }),
  setPanels: (panels) => set({ panels }),
  setControlMode: (mode) => set({ controlMode: mode }),
  setName: (name) => set({ name }),

  updatePoint: (index, point) => set((state) => {
    const newCoords = [...state.coordinates];
    newCoords[index] = point;
    return { coordinates: newCoords };
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

  setBezierHandles: (handles) => set({ bezierHandles: handles }),
  
  updateBezierHandle: (index, handle) => set((state) => {
    const newHandles = [...state.bezierHandles];
    newHandles[index] = handle;
    return { bezierHandles: newHandles };
  }),

  setBSplineControlPoints: (points) => set({ bsplineControlPoints: points }),
  
  updateBSplineControlPoint: (id, point) => set((state) => {
    const newPoints = state.bsplineControlPoints.map((p) =>
      p.id === id ? { ...p, ...point } : p
    );
    return { bsplineControlPoints: newPoints };
  }),

  addBSplineControlPoint: (point) => set((state) => ({
    bsplineControlPoints: [...state.bsplineControlPoints, point],
  })),

  removeBSplineControlPoint: (id) => set((state) => ({
    bsplineControlPoints: state.bsplineControlPoints.filter((p) => p.id !== id),
  })),

  setBSplineDegree: (degree) => set({ bsplineDegree: Math.max(1, Math.min(5, degree)) }),

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
        // wasmGenerateNaca4 expects integer digits (0-9 for m, 0-9 for p, 00-99 for t)
        // and converts them internally to fractions
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
      } catch (e) {
        console.warn('WASM NACA generation failed, using JS fallback:', e);
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
      const newPanels = repanelWithSpacing(
        state.coordinates,
        spacingKnots,
        state.nPanels
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

  reset: () => set({
    name: 'NACA 0012',
    coordinates: DEFAULT_NACA0012,
    panels: DEFAULT_NACA0012,
    controlMode: 'surface',
    bezierHandles: [],
    bsplineControlPoints: [],
    bsplineDegree: 3,
    spacingKnots: DEFAULT_SPACING_KNOTS,
    nPanels: 50,
  }),
}));
