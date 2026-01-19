/**
 * Zustand store for visualization state management
 * 
 * Manages display toggles, streamline options, smoke options, and flow speed.
 */

import { create } from 'zustand';
import type { VisualizationState } from '../types';

interface VisualizationStore extends VisualizationState {
  // Display toggle actions
  setShowGrid: (show: boolean) => void;
  setShowCurve: (show: boolean) => void;
  setShowPanels: (show: boolean) => void;
  setShowPoints: (show: boolean) => void;
  setShowControls: (show: boolean) => void;
  setShowStreamlines: (show: boolean) => void;
  setShowSmoke: (show: boolean) => void;
  setShowCp: (show: boolean) => void;
  setShowForces: (show: boolean) => void;
  
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
  
  // Flow speed action
  setFlowSpeed: (speed: number) => void;
  
  // Cp visualization options actions
  setCpDisplayMode: (mode: 'color' | 'bars' | 'both') => void;
  setCpBarScale: (scale: number) => void;
  
  // Force vector options actions
  setForceScale: (scale: number) => void;
  
  // Reset
  resetVisualization: () => void;
}

const DEFAULT_STATE: VisualizationState = {
  // Display toggles
  showGrid: true,
  showCurve: true,
  showPanels: false,
  showPoints: true,
  showControls: true,
  showStreamlines: false,
  showSmoke: false,
  showCp: false,
  showForces: false,
  
  // Animation options
  enableMorphing: true,
  morphDuration: 300,
  
  // Streamline options
  streamlineDensity: 50,
  adaptiveStreamlines: false,  // Use fixed density since we have large bounds
  
  // Smoke options
  smokeDensity: 30,          // Number of spawn points
  smokeParticlesPerBlob: 15, // Particles per blob (increased for better visuals)
  smokeSpawnInterval: 0.05,  // Seconds between spawns (much faster now that it's cheap)
  smokeMaxAge: 6.0,          // Particle lifetime in seconds (longer for larger domain)
  
  // Flow speed
  flowSpeed: 1.0,
  
  // Cp visualization options
  cpDisplayMode: 'both',
  cpBarScale: 0.1,
  
  // Force vector options
  forceScale: 0.15,
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
  setShowCp: (show) => set({ showCp: show }),
  setShowForces: (show) => set({ showForces: show }),
  
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
    smokeSpawnInterval: Math.max(0.05, Math.min(1.0, interval)) 
  }),
  setSmokeMaxAge: (age) => set({ 
    smokeMaxAge: Math.max(0.5, Math.min(10.0, age)) 
  }),
  
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
  
  // Reset
  resetVisualization: () => set(DEFAULT_STATE),
}));
