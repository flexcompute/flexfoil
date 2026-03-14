/**
 * Zustand store for visualization state management
 * 
 * Manages display toggles, streamline options, smoke options, flow speed,
 * and GPU acceleration settings.
 */

import { create } from 'zustand';
import type { VisualizationState, PerfMetrics } from '../types';

interface VisualizationStore extends VisualizationState {
  // Display toggle actions
  setShowGrid: (show: boolean) => void;
  setShowCurve: (show: boolean) => void;
  setShowPanels: (show: boolean) => void;
  setShowPoints: (show: boolean) => void;
  setShowControls: (show: boolean) => void;
  setShowStreamlines: (show: boolean) => void;
  setShowSmoke: (show: boolean) => void;
  setShowPsiContours: (show: boolean) => void;
  setShowCp: (show: boolean) => void;
  setShowForces: (show: boolean) => void;
  setShowBoundaryLayer: (show: boolean) => void;
  setShowWake: (show: boolean) => void;
  
  // Animation options actions
  setEnableMorphing: (enable: boolean) => void;
  setMorphDuration: (duration: number) => void;
  
  // Streamline options actions
  setStreamlineDensity: (density: number) => void;
  setAdaptiveStreamlines: (adaptive: boolean) => void;
  
  // Smoke options actions
  setSmokeDensity: (density: number) => void;
  setSmokeParticlesPerBlob: (count: number) => void;
  setSmokeSpawnInterval: (interval: number) => void;
  setSmokeMaxAge: (age: number) => void;
  setSmokeWaveSpacing: (spacing: number) => void;
  requestSmokeReset: () => void;
  
  // Flow speed action
  setFlowSpeed: (speed: number) => void;
  
  // Cp visualization options actions
  setCpDisplayMode: (mode: 'color' | 'bars' | 'both') => void;
  setCpBarScale: (scale: number) => void;
  
  // Force vector options actions
  setForceScale: (scale: number) => void;
  
  // Boundary layer overlay options
  setBlThicknessScale: (scale: number) => void;
  
  // GPU acceleration actions
  setUseGPU: (use: boolean) => void;
  setGPUAvailable: (available: boolean) => void;
  updatePerfMetrics: (metrics: Partial<PerfMetrics>) => void;
  
  // Reset
  resetVisualization: () => void;
}

const DEFAULT_PERF_METRICS: PerfMetrics = {
  frameTime: 0,
  avgFrameTime: 0,
  particleCount: 0,
  fps: 0,
};

const DEFAULT_STATE: VisualizationState = {
  // Display toggles
  showGrid: false,
  showCurve: true,
  showPanels: false,
  showPoints: false,
  showControls: false,
  showStreamlines: false,
  showSmoke: false,
  showPsiContours: false,
  showCp: false,
  showForces: false,
  showBoundaryLayer: false,
  showWake: false,
  
  // Animation options
  enableMorphing: true,
  morphDuration: 300,
  
  // Streamline options
  streamlineDensity: 50,
  adaptiveStreamlines: true,  // Adaptive bounds based on viewport (recommended)
  
  // Smoke options
  smokeDensity: 30,          // Number of spawn points
  smokeParticlesPerBlob: 15, // Particles per blob
  smokeSpawnInterval: 2.0,   // Seconds between blob spawns (for clear separation)
  smokeMaxAge: 6.0,          // Particle lifetime in seconds
  smokeWaveSpacing: 0.5,     // Distance between smoke waves in chord lengths
  smokeResetCounter: 0,      // Incremented to trigger smoke reset
  
  // Flow speed
  flowSpeed: 1.0,
  
  // Cp visualization options
  cpDisplayMode: 'both',
  cpBarScale: 0.1,
  
  // Force vector options
  forceScale: 0.15,
  
  // Boundary layer overlay options
  blThicknessScale: 1,
  
  // GPU acceleration
  useGPU: false,             // Default to CPU, enable when detected
  gpuAvailable: false,       // Set by feature detection
  perfMetrics: DEFAULT_PERF_METRICS,
};

export const useVisualizationStore = create<VisualizationStore>((set) => ({
  ...DEFAULT_STATE,

  // Display toggle actions
  setShowGrid: (show) => set({ showGrid: show }),
  setShowCurve: (show) => set({ showCurve: show }),
  setShowPanels: (show) => set({ showPanels: show }),
  setShowPoints: (show) => set({ showPoints: show }),
  setShowControls: (show) => set({ showControls: show }),
  setShowStreamlines: (show) => set({ showStreamlines: show }),
  setShowSmoke: (show) => set({ showSmoke: show }),
  setShowPsiContours: (show) => set({ showPsiContours: show }),
  setShowCp: (show) => set({ showCp: show }),
  setShowForces: (show) => set({ showForces: show }),
  setShowBoundaryLayer: (show) => set({ showBoundaryLayer: show }),
  setShowWake: (show) => set({ showWake: show }),
  
  // Animation options actions
  setEnableMorphing: (enable) => set({ enableMorphing: enable }),
  setMorphDuration: (duration) => set({ 
    morphDuration: Math.max(50, Math.min(1000, duration)) 
  }),
  
  // Streamline options actions
  setStreamlineDensity: (density) => set({ 
    streamlineDensity: Math.max(10, Math.min(150, density)) 
  }),
  setAdaptiveStreamlines: (adaptive) => set({ adaptiveStreamlines: adaptive }),
  
  // Smoke options actions
  setSmokeDensity: (density) => set({ 
    smokeDensity: Math.max(10, Math.min(80, density)) 
  }),
  setSmokeParticlesPerBlob: (count) => set({ 
    smokeParticlesPerBlob: Math.max(3, Math.min(30, count)) 
  }),
  setSmokeSpawnInterval: (interval) => set({ 
    smokeSpawnInterval: Math.max(0.1, Math.min(20.0, interval)) 
  }),
  setSmokeMaxAge: (age) => set({ 
    smokeMaxAge: Math.max(0.5, Math.min(30.0, age)) 
  }),
  setSmokeWaveSpacing: (spacing) => set({ 
    smokeWaveSpacing: Math.max(0.2, Math.min(2.0, spacing)) 
  }),
  requestSmokeReset: () => set((state) => ({ 
    smokeResetCounter: state.smokeResetCounter + 1 
  })),
  
  // Flow speed action
  setFlowSpeed: (speed) => set({ 
    flowSpeed: Math.max(0.1, Math.min(5.0, speed)) 
  }),
  
  // Cp visualization options actions
  setCpDisplayMode: (mode) => set({ cpDisplayMode: mode }),
  setCpBarScale: (scale) => set({ 
    cpBarScale: Math.max(0.01, Math.min(0.5, scale)) 
  }),
  
  // Force vector options actions
  setForceScale: (scale) => set({ 
    forceScale: Math.max(0.05, Math.min(0.5, scale)) 
  }),
  
  // Boundary layer overlay options
  setBlThicknessScale: (scale) => set({
    blThicknessScale: Math.max(1, Math.min(100, scale))
  }),
  
  // GPU acceleration actions
  setUseGPU: (use) => set((state) => ({ 
    useGPU: use && state.gpuAvailable  // Can only enable if GPU is available
  })),
  setGPUAvailable: (available) => set({ gpuAvailable: available }),
  updatePerfMetrics: (metrics) => set((state) => ({
    perfMetrics: { ...state.perfMetrics, ...metrics }
  })),
  
  // Reset
  resetVisualization: () => set(DEFAULT_STATE),
}));
