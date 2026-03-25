/**
 * Zustand store for CFD solver state management.
 *
 * Manages grid settings, physics mode, solver runtime state,
 * convergence history, and force coefficient history.
 */

import { create } from 'zustand';
import type {
  CfdConfig,
  PhysicsMode,
  ReconstructionMode,
  TimeSteppingMode,
} from '../lib/webgpu/cfd/CfdConfig';
import { DEFAULT_CFD_CONFIG } from '../lib/webgpu/cfd/CfdConfig';

export interface ConvergencePoint {
  iteration: number;
  residualL2: number;
}

export interface ForcePoint {
  iteration: number;
  cl: number;
  cd: number;
  cm: number;
}

interface CfdState {
  // Configuration
  config: CfdConfig;

  // Runtime state
  isRunning: boolean;
  iteration: number;
  converged: boolean;
  meshGenerated: boolean;

  // History for plots
  convergenceHistory: ConvergencePoint[];
  forceHistory: ForcePoint[];

  // Latest values
  cl: number;
  cd: number;
  cm: number;
  residualL2: number;
}

interface CfdActions {
  // Config setters
  setNi: (ni: number) => void;
  setNj: (nj: number) => void;
  setMach: (mach: number) => void;
  setAlphaDeg: (alpha: number) => void;
  setReynolds: (re: number) => void;
  setCfl: (cfl: number) => void;
  setPhysics: (mode: PhysicsMode) => void;
  setReconstruction: (mode: ReconstructionMode) => void;
  setTimeStepping: (mode: TimeSteppingMode) => void;
  setFarField: (dist: number) => void;
  setDs0: (ds0: number) => void;

  // Runtime actions
  setRunning: (running: boolean) => void;
  setMeshGenerated: (generated: boolean) => void;
  appendResult: (result: {
    iteration: number;
    cl: number;
    cd: number;
    cm: number;
    residualL2: number;
  }) => void;
  setConverged: (converged: boolean) => void;
  reset: () => void;
}

type CfdStore = CfdState & CfdActions;

const initialState: CfdState = {
  config: { ...DEFAULT_CFD_CONFIG },
  isRunning: false,
  iteration: 0,
  converged: false,
  meshGenerated: false,
  convergenceHistory: [],
  forceHistory: [],
  cl: 0,
  cd: 0,
  cm: 0,
  residualL2: 1,
};

export const useCfdStore = create<CfdStore>((set) => ({
  ...initialState,

  setNi: (ni) => set((s) => ({ config: { ...s.config, ni } })),
  setNj: (nj) => set((s) => ({ config: { ...s.config, nj } })),
  setMach: (machInf) => set((s) => ({ config: { ...s.config, machInf } })),
  setAlphaDeg: (alphaDeg) => set((s) => ({ config: { ...s.config, alphaDeg } })),
  setReynolds: (reynolds) => set((s) => ({ config: { ...s.config, reynolds } })),
  setCfl: (cfl) => set((s) => ({ config: { ...s.config, cfl } })),
  setPhysics: (physics) => set((s) => ({ config: { ...s.config, physics } })),
  setReconstruction: (reconstruction) => set((s) => ({ config: { ...s.config, reconstruction } })),
  setTimeStepping: (timeStepping) => set((s) => ({ config: { ...s.config, timeStepping } })),
  setFarField: (farField) => set((s) => ({ config: { ...s.config, farField } })),
  setDs0: (ds0) => set((s) => ({ config: { ...s.config, ds0 } })),

  setRunning: (isRunning) => set({ isRunning }),
  setMeshGenerated: (meshGenerated) => set({ meshGenerated }),
  setConverged: (converged) => set({ converged }),

  appendResult: (result) =>
    set((s) => ({
      iteration: result.iteration,
      cl: result.cl,
      cd: result.cd,
      cm: result.cm,
      residualL2: result.residualL2,
      convergenceHistory: [
        ...s.convergenceHistory,
        { iteration: result.iteration, residualL2: result.residualL2 },
      ],
      forceHistory: [
        ...s.forceHistory,
        {
          iteration: result.iteration,
          cl: result.cl,
          cd: result.cd,
          cm: result.cm,
        },
      ],
    })),

  reset: () => set(initialState),
}));
