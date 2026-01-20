/**
 * URL State Encoding/Decoding
 * 
 * Encodes application state into compact URL hash for sharing.
 * Uses lz-string for compression.
 */

import { compressToEncodedURIComponent, decompressFromEncodedURIComponent } from 'lz-string';
import type { ControlMode, SpacingKnot } from '../types';

// Default values (don't encode if unchanged)
const DEFAULTS = {
  naca: '0012',
  nPanels: 160,
  mode: 'parameters' as ControlMode,
  alpha: 0,
  polarStart: -5,
  polarEnd: 15,
  polarStep: 1,
  theme: 'dark' as 'dark' | 'light',
  // Viewport defaults
  viewCenterX: 0.5,
  viewCenterY: 0,
  viewZoom: 400,
  // Display toggles defaults
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
  // Animation defaults
  enableMorphing: true,
  morphDuration: 300,
  // Streamline defaults
  streamlineDensity: 50,
  adaptiveStreamlines: true,
  // Smoke defaults
  smokeDensity: 30,
  smokeParticlesPerBlob: 15,
  // Flow speed default
  flowSpeed: 1.0,
  // Cp defaults
  cpDisplayMode: 'both' as 'color' | 'bars' | 'both',
  cpBarScale: 0.1,
  // Force defaults
  forceScale: 0.15,
  // Parameters mode defaults
  thicknessScale: 1.0,
  camberScale: 1.0,
  // Spacing panel defaults
  curvatureWeight: 0,
  spacingPanelMode: 'simple' as 'simple' | 'advanced',
  sspInterpolation: 'linear' as 'linear' | 'spline',
  sspVisualization: 'plot' as 'plot' | 'foil',
};

const DEFAULT_SPACING: SpacingKnot[] = [
  { S: 0, F: 1.5 },
  { S: 0.25, F: 0.4 },
  { S: 0.5, F: 1.5 },
  { S: 0.75, F: 0.4 },
  { S: 1, F: 1.5 },
];

/**
 * State that can be encoded in URL
 */
export interface UrlState {
  // Airfoil
  naca?: string;           // NACA designation (e.g., "0012", "2412")
  custom?: boolean;        // True if airfoil is not NACA
  
  // Paneling
  nPanels?: number;
  spacing?: SpacingKnot[];
  
  // Control mode
  mode?: ControlMode;
  
  // Solver settings
  alpha?: number;
  polarStart?: number;
  polarEnd?: number;
  polarStep?: number;
  
  // Theme
  theme?: 'dark' | 'light';
  
  // Panel layout (only if modified from default)
  layout?: object;
  
  // Viewport (zoom/pan state)
  viewCenterX?: number;
  viewCenterY?: number;
  viewZoom?: number;
  
  // Display toggles
  showGrid?: boolean;
  showCurve?: boolean;
  showPanels?: boolean;
  showPoints?: boolean;
  showControls?: boolean;
  showStreamlines?: boolean;
  showSmoke?: boolean;
  showPsiContours?: boolean;
  showCp?: boolean;
  showForces?: boolean;
  
  // Animation options
  enableMorphing?: boolean;
  morphDuration?: number;
  
  // Streamline options
  streamlineDensity?: number;
  adaptiveStreamlines?: boolean;
  
  // Smoke options
  smokeDensity?: number;
  smokeParticlesPerBlob?: number;
  
  // Flow speed
  flowSpeed?: number;
  
  // Cp visualization options
  cpDisplayMode?: 'color' | 'bars' | 'both';
  cpBarScale?: number;
  
  // Force vector options
  forceScale?: number;
  
  // Parameters mode (thickness/camber scaling)
  thicknessScale?: number;
  camberScale?: number;
  
  // Spacing panel settings
  curvatureWeight?: number;
  spacingPanelMode?: 'simple' | 'advanced';
  sspInterpolation?: 'linear' | 'spline';
  sspVisualization?: 'plot' | 'foil';
}

/**
 * Compact state for encoding (short keys, omit defaults)
 */
interface CompactState {
  n?: string;      // naca
  c?: 1;           // custom (flag)
  p?: number;      // nPanels
  s?: number[][];  // spacing [[S,F], ...]
  m?: string;      // mode (first char: p/c/t)
  a?: number;      // alpha
  ps?: number;     // polar start
  pe?: number;     // polar end
  pd?: number;     // polar delta
  t?: string;      // theme (d/l)
  l?: object;      // layout
  // Viewport
  vcx?: number;    // view center x
  vcy?: number;    // view center y
  vz?: number;     // view zoom
  // Display toggles (only store if true to save space)
  dg?: 1;          // showGrid
  dc?: 1;          // showCurve
  dp?: 1;          // showPanels
  dpt?: 1;         // showPoints
  dct?: 1;         // showControls
  dst?: 1;         // showStreamlines
  dsm?: 1;         // showSmoke
  dps?: 1;         // showPsiContours
  dcp?: 1;         // showCp
  df?: 1;          // showForces
  // Animation
  em?: 1;          // enableMorphing
  md?: number;     // morphDuration
  // Streamlines
  sld?: number;    // streamlineDensity
  sla?: 1;         // adaptiveStreamlines
  // Smoke
  smd?: number;    // smokeDensity
  smp?: number;    // smokeParticlesPerBlob
  // Flow
  fs?: number;     // flowSpeed
  // Cp
  cpm?: string;    // cpDisplayMode (c/b/x for color/bars/both)
  cpb?: number;    // cpBarScale
  // Force
  fsc?: number;    // forceScale
  // Parameters mode
  ts?: number;     // thicknessScale
  cs?: number;     // camberScale
  // Spacing panel
  cw?: number;     // curvatureWeight
  spm?: string;    // spacingPanelMode (s/a for simple/advanced)
  sspi?: string;   // sspInterpolation (l/s for linear/spline)
  sspv?: string;   // sspVisualization (p/f for plot/foil)
}

/**
 * Encode state to URL-safe string
 */
export function encodeUrlState(state: UrlState): string {
  const compact: CompactState = {};
  
  // NACA designation
  if (state.naca && state.naca !== DEFAULTS.naca) {
    compact.n = state.naca;
  }
  if (state.custom) {
    compact.c = 1;
  }
  
  // Paneling
  if (state.nPanels && state.nPanels !== DEFAULTS.nPanels) {
    compact.p = state.nPanels;
  }
  if (state.spacing && !isDefaultSpacing(state.spacing)) {
    compact.s = state.spacing.map(k => [k.S, k.F]);
  }
  
  // Control mode
  if (state.mode && state.mode !== DEFAULTS.mode) {
    // p = parameters, c = camber-spline, t = thickness-spline
    switch (state.mode) {
      case 'parameters': compact.m = 'p'; break;
      case 'camber-spline': compact.m = 'c'; break;
      case 'thickness-spline': compact.m = 't'; break;
    }
  }
  
  // Solver settings
  if (state.alpha !== undefined && state.alpha !== DEFAULTS.alpha) {
    compact.a = state.alpha;
  }
  if (state.polarStart !== undefined && state.polarStart !== DEFAULTS.polarStart) {
    compact.ps = state.polarStart;
  }
  if (state.polarEnd !== undefined && state.polarEnd !== DEFAULTS.polarEnd) {
    compact.pe = state.polarEnd;
  }
  if (state.polarStep !== undefined && state.polarStep !== DEFAULTS.polarStep) {
    compact.pd = state.polarStep;
  }
  
  // Theme
  if (state.theme && state.theme !== DEFAULTS.theme) {
    compact.t = state.theme[0]; // 'd' or 'l'
  }
  
  // Layout (only if provided)
  if (state.layout) {
    compact.l = state.layout;
  }
  
  // Viewport
  if (state.viewCenterX !== undefined && state.viewCenterX !== DEFAULTS.viewCenterX) {
    compact.vcx = state.viewCenterX;
  }
  if (state.viewCenterY !== undefined && state.viewCenterY !== DEFAULTS.viewCenterY) {
    compact.vcy = state.viewCenterY;
  }
  if (state.viewZoom !== undefined && state.viewZoom !== DEFAULTS.viewZoom) {
    compact.vz = state.viewZoom;
  }
  
  // Display toggles (only encode if different from default)
  if (state.showGrid !== undefined && state.showGrid !== DEFAULTS.showGrid) {
    compact.dg = state.showGrid ? 1 : undefined;
  }
  if (state.showCurve !== undefined && state.showCurve !== DEFAULTS.showCurve) {
    compact.dc = state.showCurve ? 1 : undefined;
  }
  if (state.showPanels !== undefined && state.showPanels !== DEFAULTS.showPanels) {
    compact.dp = state.showPanels ? 1 : undefined;
  }
  if (state.showPoints !== undefined && state.showPoints !== DEFAULTS.showPoints) {
    compact.dpt = state.showPoints ? 1 : undefined;
  }
  if (state.showControls !== undefined && state.showControls !== DEFAULTS.showControls) {
    compact.dct = state.showControls ? 1 : undefined;
  }
  if (state.showStreamlines !== undefined && state.showStreamlines !== DEFAULTS.showStreamlines) {
    compact.dst = state.showStreamlines ? 1 : undefined;
  }
  if (state.showSmoke !== undefined && state.showSmoke !== DEFAULTS.showSmoke) {
    compact.dsm = state.showSmoke ? 1 : undefined;
  }
  if (state.showPsiContours !== undefined && state.showPsiContours !== DEFAULTS.showPsiContours) {
    compact.dps = state.showPsiContours ? 1 : undefined;
  }
  if (state.showCp !== undefined && state.showCp !== DEFAULTS.showCp) {
    compact.dcp = state.showCp ? 1 : undefined;
  }
  if (state.showForces !== undefined && state.showForces !== DEFAULTS.showForces) {
    compact.df = state.showForces ? 1 : undefined;
  }
  
  // Animation
  if (state.enableMorphing !== undefined && state.enableMorphing !== DEFAULTS.enableMorphing) {
    compact.em = state.enableMorphing ? 1 : undefined;
  }
  if (state.morphDuration !== undefined && state.morphDuration !== DEFAULTS.morphDuration) {
    compact.md = state.morphDuration;
  }
  
  // Streamlines
  if (state.streamlineDensity !== undefined && state.streamlineDensity !== DEFAULTS.streamlineDensity) {
    compact.sld = state.streamlineDensity;
  }
  if (state.adaptiveStreamlines !== undefined && state.adaptiveStreamlines !== DEFAULTS.adaptiveStreamlines) {
    compact.sla = state.adaptiveStreamlines ? 1 : undefined;
  }
  
  // Smoke
  if (state.smokeDensity !== undefined && state.smokeDensity !== DEFAULTS.smokeDensity) {
    compact.smd = state.smokeDensity;
  }
  if (state.smokeParticlesPerBlob !== undefined && state.smokeParticlesPerBlob !== DEFAULTS.smokeParticlesPerBlob) {
    compact.smp = state.smokeParticlesPerBlob;
  }
  
  // Flow speed
  if (state.flowSpeed !== undefined && state.flowSpeed !== DEFAULTS.flowSpeed) {
    compact.fs = state.flowSpeed;
  }
  
  // Cp
  if (state.cpDisplayMode && state.cpDisplayMode !== DEFAULTS.cpDisplayMode) {
    switch (state.cpDisplayMode) {
      case 'color': compact.cpm = 'c'; break;
      case 'bars': compact.cpm = 'b'; break;
      case 'both': compact.cpm = 'x'; break;
    }
  }
  if (state.cpBarScale !== undefined && state.cpBarScale !== DEFAULTS.cpBarScale) {
    compact.cpb = state.cpBarScale;
  }
  
  // Force
  if (state.forceScale !== undefined && state.forceScale !== DEFAULTS.forceScale) {
    compact.fsc = state.forceScale;
  }
  
  // Parameters mode
  if (state.thicknessScale !== undefined && state.thicknessScale !== DEFAULTS.thicknessScale) {
    compact.ts = state.thicknessScale;
  }
  if (state.camberScale !== undefined && state.camberScale !== DEFAULTS.camberScale) {
    compact.cs = state.camberScale;
  }
  
  // Spacing panel
  if (state.curvatureWeight !== undefined && state.curvatureWeight !== DEFAULTS.curvatureWeight) {
    compact.cw = state.curvatureWeight;
  }
  if (state.spacingPanelMode && state.spacingPanelMode !== DEFAULTS.spacingPanelMode) {
    compact.spm = state.spacingPanelMode === 'simple' ? 's' : 'a';
  }
  if (state.sspInterpolation && state.sspInterpolation !== DEFAULTS.sspInterpolation) {
    compact.sspi = state.sspInterpolation === 'linear' ? 'l' : 's';
  }
  if (state.sspVisualization && state.sspVisualization !== DEFAULTS.sspVisualization) {
    compact.sspv = state.sspVisualization === 'plot' ? 'p' : 'f';
  }
  
  // If nothing to encode, return empty
  if (Object.keys(compact).length === 0) {
    return '';
  }
  
  // Compress and return
  const json = JSON.stringify(compact);
  return compressToEncodedURIComponent(json);
}

/**
 * Decode URL string to state
 */
export function decodeUrlState(encoded: string): UrlState | null {
  if (!encoded) return null;
  
  try {
    const json = decompressFromEncodedURIComponent(encoded);
    if (!json) return null;
    
    const compact: CompactState = JSON.parse(json);
    const state: UrlState = {};
    
    // NACA designation
    if (compact.n) {
      state.naca = compact.n;
    }
    if (compact.c) {
      state.custom = true;
    }
    
    // Paneling
    if (compact.p) {
      state.nPanels = compact.p;
    }
    if (compact.s) {
      state.spacing = compact.s.map(([S, F]) => ({ S, F }));
    }
    
    // Control mode
    if (compact.m) {
      switch (compact.m) {
        case 'p': state.mode = 'parameters'; break;
        case 'c': state.mode = 'camber-spline'; break;
        case 't': state.mode = 'thickness-spline'; break;
      }
    }
    
    // Solver settings
    if (compact.a !== undefined) state.alpha = compact.a;
    if (compact.ps !== undefined) state.polarStart = compact.ps;
    if (compact.pe !== undefined) state.polarEnd = compact.pe;
    if (compact.pd !== undefined) state.polarStep = compact.pd;
    
    // Theme
    if (compact.t) {
      state.theme = compact.t === 'l' ? 'light' : 'dark';
    }
    
    // Layout
    if (compact.l) {
      state.layout = compact.l;
    }
    
    // Viewport
    if (compact.vcx !== undefined) state.viewCenterX = compact.vcx;
    if (compact.vcy !== undefined) state.viewCenterY = compact.vcy;
    if (compact.vz !== undefined) state.viewZoom = compact.vz;
    
    // Display toggles
    if (compact.dg !== undefined) state.showGrid = !!compact.dg;
    if (compact.dc !== undefined) state.showCurve = !!compact.dc;
    if (compact.dp !== undefined) state.showPanels = !!compact.dp;
    if (compact.dpt !== undefined) state.showPoints = !!compact.dpt;
    if (compact.dct !== undefined) state.showControls = !!compact.dct;
    if (compact.dst !== undefined) state.showStreamlines = !!compact.dst;
    if (compact.dsm !== undefined) state.showSmoke = !!compact.dsm;
    if (compact.dps !== undefined) state.showPsiContours = !!compact.dps;
    if (compact.dcp !== undefined) state.showCp = !!compact.dcp;
    if (compact.df !== undefined) state.showForces = !!compact.df;
    
    // Animation
    if (compact.em !== undefined) state.enableMorphing = !!compact.em;
    if (compact.md !== undefined) state.morphDuration = compact.md;
    
    // Streamlines
    if (compact.sld !== undefined) state.streamlineDensity = compact.sld;
    if (compact.sla !== undefined) state.adaptiveStreamlines = !!compact.sla;
    
    // Smoke
    if (compact.smd !== undefined) state.smokeDensity = compact.smd;
    if (compact.smp !== undefined) state.smokeParticlesPerBlob = compact.smp;
    
    // Flow speed
    if (compact.fs !== undefined) state.flowSpeed = compact.fs;
    
    // Cp
    if (compact.cpm) {
      switch (compact.cpm) {
        case 'c': state.cpDisplayMode = 'color'; break;
        case 'b': state.cpDisplayMode = 'bars'; break;
        case 'x': state.cpDisplayMode = 'both'; break;
      }
    }
    if (compact.cpb !== undefined) state.cpBarScale = compact.cpb;
    
    // Force
    if (compact.fsc !== undefined) state.forceScale = compact.fsc;
    
    // Parameters mode
    if (compact.ts !== undefined) state.thicknessScale = compact.ts;
    if (compact.cs !== undefined) state.camberScale = compact.cs;
    
    // Spacing panel
    if (compact.cw !== undefined) state.curvatureWeight = compact.cw;
    if (compact.spm) {
      state.spacingPanelMode = compact.spm === 's' ? 'simple' : 'advanced';
    }
    if (compact.sspi) {
      state.sspInterpolation = compact.sspi === 'l' ? 'linear' : 'spline';
    }
    if (compact.sspv) {
      state.sspVisualization = compact.sspv === 'p' ? 'plot' : 'foil';
    }
    
    return state;
  } catch (e) {
    console.warn('Failed to decode URL state:', e);
    return null;
  }
}

/**
 * Update URL hash without triggering navigation
 */
export function syncToUrl(state: UrlState): void {
  const encoded = encodeUrlState(state);
  
  if (encoded) {
    window.history.replaceState(null, '', `#s=${encoded}`);
  } else {
    // Clear hash if no state to encode
    window.history.replaceState(null, '', window.location.pathname);
  }
}

/**
 * Load state from current URL hash
 */
export function loadFromUrl(): UrlState | null {
  const hash = window.location.hash;
  if (!hash || !hash.startsWith('#s=')) {
    return null;
  }
  
  const encoded = hash.slice(3); // Remove '#s='
  return decodeUrlState(encoded);
}

/**
 * Generate shareable URL with current state
 */
export function getShareableUrl(state: UrlState): string {
  const encoded = encodeUrlState(state);
  const base = window.location.origin + window.location.pathname;
  return encoded ? `${base}#s=${encoded}` : base;
}

/**
 * Parse NACA designation from airfoil name
 */
export function parseNacaFromName(name: string): string | null {
  const match = name.match(/NACA\s*(\d{4})/i);
  return match ? match[1] : null;
}

/**
 * Get complete URL state from all stores
 * This should be called by components that need to sync to URL
 */
export function getCompleteUrlState(
  airfoilState: {
    name: string;
    nPanels: number;
    spacingKnots: any[];
    controlMode: any;
    displayAlpha: number;
    thicknessScale: number;
    camberScale: number;
    curvatureWeight: number;
    spacingPanelMode: 'simple' | 'advanced';
    sspInterpolation: 'linear' | 'spline';
    sspVisualization: 'plot' | 'foil';
  },
  visualizationState: {
    showGrid: boolean;
    showCurve: boolean;
    showPanels: boolean;
    showPoints: boolean;
    showControls: boolean;
    showStreamlines: boolean;
    showSmoke: boolean;
    showPsiContours: boolean;
    showCp: boolean;
    showForces: boolean;
    enableMorphing: boolean;
    morphDuration: number;
    streamlineDensity: number;
    adaptiveStreamlines: boolean;
    smokeDensity: number;
    smokeParticlesPerBlob: number;
    flowSpeed: number;
    cpDisplayMode: 'color' | 'bars' | 'both';
    cpBarScale: number;
    forceScale: number;
  },
  viewportState?: {
    centerX: number;
    centerY: number;
    zoom: number;
  }
): UrlState {
  const naca = parseNacaFromName(airfoilState.name);
  
  return {
    // Airfoil
    naca: naca || undefined,
    custom: !naca,
    nPanels: airfoilState.nPanels,
    spacing: airfoilState.spacingKnots,
    mode: airfoilState.controlMode,
    alpha: airfoilState.displayAlpha,
    
    // Viewport
    viewCenterX: viewportState?.centerX,
    viewCenterY: viewportState?.centerY,
    viewZoom: viewportState?.zoom,
    
    // Display toggles
    showGrid: visualizationState.showGrid,
    showCurve: visualizationState.showCurve,
    showPanels: visualizationState.showPanels,
    showPoints: visualizationState.showPoints,
    showControls: visualizationState.showControls,
    showStreamlines: visualizationState.showStreamlines,
    showSmoke: visualizationState.showSmoke,
    showPsiContours: visualizationState.showPsiContours,
    showCp: visualizationState.showCp,
    showForces: visualizationState.showForces,
    
    // Animation
    enableMorphing: visualizationState.enableMorphing,
    morphDuration: visualizationState.morphDuration,
    
    // Streamlines
    streamlineDensity: visualizationState.streamlineDensity,
    adaptiveStreamlines: visualizationState.adaptiveStreamlines,
    
    // Smoke
    smokeDensity: visualizationState.smokeDensity,
    smokeParticlesPerBlob: visualizationState.smokeParticlesPerBlob,
    
    // Flow
    flowSpeed: visualizationState.flowSpeed,
    
    // Cp
    cpDisplayMode: visualizationState.cpDisplayMode,
    cpBarScale: visualizationState.cpBarScale,
    
    // Force
    forceScale: visualizationState.forceScale,
    
    // Parameters mode
    thicknessScale: airfoilState.thicknessScale,
    camberScale: airfoilState.camberScale,
    
    // Spacing panel
    curvatureWeight: airfoilState.curvatureWeight,
    spacingPanelMode: airfoilState.spacingPanelMode,
    sspInterpolation: airfoilState.sspInterpolation,
    sspVisualization: airfoilState.sspVisualization,
  };
}

/**
 * Check if spacing matches default
 */
function isDefaultSpacing(spacing: SpacingKnot[]): boolean {
  if (spacing.length !== DEFAULT_SPACING.length) return false;
  
  for (let i = 0; i < spacing.length; i++) {
    if (Math.abs(spacing[i].S - DEFAULT_SPACING[i].S) > 0.001) return false;
    if (Math.abs(spacing[i].F - DEFAULT_SPACING[i].F) > 0.001) return false;
  }
  
  return true;
}
