/**
 * URL State Encoding/Decoding
 * 
 * Encodes application state into compact URL hash for sharing.
 * Uses lz-string for compression.
 */

import { compressToEncodedURIComponent, decompressFromEncodedURIComponent } from 'lz-string';
import type { ControlMode, FlapDefinition, SpacingKnot } from '../types';

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

  // Flap definitions
  flaps?: FlapDefinition[];
}

/**
 * Compact state for encoding (single-char keys where possible).
 *
 * Display/animation booleans are packed into a single bitmask `v`:
 *   bit 0: showGrid          bit 6:  showSmoke
 *   bit 1: showCurve         bit 7:  showPsiContours
 *   bit 2: showPanels        bit 8:  showCp
 *   bit 3: showPoints        bit 9:  showForces
 *   bit 4: showControls      bit 10: enableMorphing
 *   bit 5: showStreamlines   bit 11: adaptiveStreamlines
 */
interface CompactState {
  n?: string;      // naca
  c?: 1;           // custom
  p?: number;      // nPanels
  s?: number[][];  // spacing
  m?: string;      // mode
  a?: number;      // alpha
  b?: number;      // polar start (begin)
  e?: number;      // polar end
  d?: number;      // polar step (delta)
  t?: string;      // theme
  l?: object;      // layout
  x?: number;      // view center x
  y?: number;      // view center y
  z?: number;      // view zoom
  v?: number;      // vis-flag bitmask (see above)
  h?: number;      // morph duration
  i?: number;      // streamline density
  j?: number;      // smoke density
  k?: number;      // smoke particles per blob
  f?: number;      // flow speed
  o?: string;      // cp display mode (c/b/x)
  q?: number;      // cp bar scale
  r?: number;      // force scale
  u?: number;      // thickness scale
  w?: number;      // camber scale
  g?: number;      // curvature weight
  // Spacing panel: mode + interpolation + vis packed as "sap" etc.
  sp?: string;
  // Flaps: [name, hingeX, hingeYFrac, deflection] (id regenerated on load)
  F?: Array<[string, number, number, number]>;

  // --- Legacy keys (accepted on decode, never emitted) ---
  ps?: number; pe?: number; pd?: number;
  vcx?: number; vcy?: number; vz?: number;
  dg?: 1; dc?: 1; dp?: 1; dpt?: 1; dct?: 1;
  dst?: 1; dsm?: 1; dps?: 1; dcp?: 1; df?: 1;
  em?: 1; md?: number; sld?: number; sla?: 1;
  smd?: number; smp?: number; fs?: number;
  cpm?: string; cpb?: number; fsc?: number;
  ts?: number; cs?: number; cw?: number;
  spm?: string; sspi?: string; sspv?: string;
  fl?: Array<[string, string, number, number, number]>;
}

/**
 * Encode state to URL-safe string
 */
export function encodeUrlState(state: UrlState): string {
  const compact: CompactState = {};

  if (state.naca && state.naca !== DEFAULTS.naca) compact.n = state.naca;
  if (state.custom) compact.c = 1;
  if (state.nPanels && state.nPanels !== DEFAULTS.nPanels) compact.p = state.nPanels;
  if (state.spacing && !isDefaultSpacing(state.spacing)) compact.s = state.spacing.map(k => [k.S, k.F]);

  if (state.mode && state.mode !== DEFAULTS.mode) {
    switch (state.mode) {
      case 'parameters': compact.m = 'p'; break;
      case 'camber-spline': compact.m = 'c'; break;
      case 'thickness-spline': compact.m = 't'; break;
      case 'inverse-design': compact.m = 'i'; break;
      case 'geometry-design': compact.m = 'g'; break;
    }
  }

  if (state.alpha !== undefined && state.alpha !== DEFAULTS.alpha) compact.a = state.alpha;
  if (state.polarStart !== undefined && state.polarStart !== DEFAULTS.polarStart) compact.b = state.polarStart;
  if (state.polarEnd !== undefined && state.polarEnd !== DEFAULTS.polarEnd) compact.e = state.polarEnd;
  if (state.polarStep !== undefined && state.polarStep !== DEFAULTS.polarStep) compact.d = state.polarStep;
  if (state.theme && state.theme !== DEFAULTS.theme) compact.t = state.theme[0];
  if (state.layout) compact.l = state.layout;

  // Viewport
  if (state.viewCenterX !== undefined && state.viewCenterX !== DEFAULTS.viewCenterX) compact.x = state.viewCenterX;
  if (state.viewCenterY !== undefined && state.viewCenterY !== DEFAULTS.viewCenterY) compact.y = state.viewCenterY;
  if (state.viewZoom !== undefined && state.viewZoom !== DEFAULTS.viewZoom) compact.z = state.viewZoom;

  // Pack 12 boolean flags into one bitmask
  const flags = [
    state.showGrid, state.showCurve, state.showPanels, state.showPoints,
    state.showControls, state.showStreamlines, state.showSmoke, state.showPsiContours,
    state.showCp, state.showForces, state.enableMorphing, state.adaptiveStreamlines,
  ];
  const defaults = [
    DEFAULTS.showGrid, DEFAULTS.showCurve, DEFAULTS.showPanels, DEFAULTS.showPoints,
    DEFAULTS.showControls, DEFAULTS.showStreamlines, DEFAULTS.showSmoke, DEFAULTS.showPsiContours,
    DEFAULTS.showCp, DEFAULTS.showForces, DEFAULTS.enableMorphing, DEFAULTS.adaptiveStreamlines,
  ];
  const anyDifferent = flags.some((f, i) => f !== undefined && f !== defaults[i]);
  if (anyDifferent) {
    let bits = 0;
    for (let i = 0; i < flags.length; i++) {
      if (flags[i] ?? defaults[i]) bits |= (1 << i);
    }
    compact.v = bits;
  }

  if (state.morphDuration !== undefined && state.morphDuration !== DEFAULTS.morphDuration) compact.h = state.morphDuration;
  if (state.streamlineDensity !== undefined && state.streamlineDensity !== DEFAULTS.streamlineDensity) compact.i = state.streamlineDensity;
  if (state.smokeDensity !== undefined && state.smokeDensity !== DEFAULTS.smokeDensity) compact.j = state.smokeDensity;
  if (state.smokeParticlesPerBlob !== undefined && state.smokeParticlesPerBlob !== DEFAULTS.smokeParticlesPerBlob) compact.k = state.smokeParticlesPerBlob;
  if (state.flowSpeed !== undefined && state.flowSpeed !== DEFAULTS.flowSpeed) compact.f = state.flowSpeed;

  if (state.cpDisplayMode && state.cpDisplayMode !== DEFAULTS.cpDisplayMode) {
    compact.o = state.cpDisplayMode === 'color' ? 'c' : state.cpDisplayMode === 'bars' ? 'b' : 'x';
  }
  if (state.cpBarScale !== undefined && state.cpBarScale !== DEFAULTS.cpBarScale) compact.q = state.cpBarScale;
  if (state.forceScale !== undefined && state.forceScale !== DEFAULTS.forceScale) compact.r = state.forceScale;
  if (state.thicknessScale !== undefined && state.thicknessScale !== DEFAULTS.thicknessScale) compact.u = state.thicknessScale;
  if (state.camberScale !== undefined && state.camberScale !== DEFAULTS.camberScale) compact.w = state.camberScale;
  if (state.curvatureWeight !== undefined && state.curvatureWeight !== DEFAULTS.curvatureWeight) compact.g = state.curvatureWeight;

  // Spacing panel: pack mode+interp+vis into 3-char string
  const spMode = state.spacingPanelMode ?? DEFAULTS.spacingPanelMode;
  const spInterp = state.sspInterpolation ?? DEFAULTS.sspInterpolation;
  const spVis = state.sspVisualization ?? DEFAULTS.sspVisualization;
  if (spMode !== DEFAULTS.spacingPanelMode || spInterp !== DEFAULTS.sspInterpolation || spVis !== DEFAULTS.sspVisualization) {
    compact.sp = (spMode === 'simple' ? 's' : 'a') + (spInterp === 'linear' ? 'l' : 's') + (spVis === 'plot' ? 'p' : 'f');
  }

  // Flaps (omit id — regenerated on decode)
  if (state.flaps && state.flaps.length > 0) {
    compact.F = state.flaps.map(f => [f.name, f.hingeX, f.hingeYFrac, f.deflection]);
  }

  if (Object.keys(compact).length === 0) return '';
  return compressToEncodedURIComponent(JSON.stringify(compact));
}

/**
 * Decode URL string to state
 */
export function decodeUrlState(encoded: string): UrlState | null {
  if (!encoded) return null;
  
  try {
    const json = decompressFromEncodedURIComponent(encoded);
    if (!json) return null;
    
    const C: CompactState = JSON.parse(json);
    const state: UrlState = {};
    
    if (C.n) state.naca = C.n;
    if (C.c) state.custom = true;
    if (C.p) state.nPanels = C.p;
    if (C.s) state.spacing = C.s.map(([S, F]) => ({ S, F }));

    if (C.m) {
      switch (C.m) {
        case 'p': state.mode = 'parameters'; break;
        case 'c': state.mode = 'camber-spline'; break;
        case 't': state.mode = 'thickness-spline'; break;
        case 'i': state.mode = 'inverse-design'; break;
        case 'g': state.mode = 'geometry-design'; break;
      }
    }

    if (C.a !== undefined) state.alpha = C.a;
    state.polarStart = C.b ?? C.ps;
    state.polarEnd = C.e ?? C.pe;
    state.polarStep = C.d ?? C.pd;
    if (C.t) state.theme = C.t === 'l' ? 'light' : 'dark';
    if (C.l) state.layout = C.l;

    // Viewport (new keys x/y/z, fallback to legacy vcx/vcy/vz)
    state.viewCenterX = C.x ?? C.vcx;
    state.viewCenterY = C.y ?? C.vcy;
    state.viewZoom = C.z ?? C.vz;

    // Display bitmask (new key `v`) or legacy individual flags
    if (C.v !== undefined) {
      const bit = (i: number) => !!(C.v! & (1 << i));
      state.showGrid = bit(0);
      state.showCurve = bit(1);
      state.showPanels = bit(2);
      state.showPoints = bit(3);
      state.showControls = bit(4);
      state.showStreamlines = bit(5);
      state.showSmoke = bit(6);
      state.showPsiContours = bit(7);
      state.showCp = bit(8);
      state.showForces = bit(9);
      state.enableMorphing = bit(10);
      state.adaptiveStreamlines = bit(11);
    } else {
      if (C.dg !== undefined) state.showGrid = !!C.dg;
      if (C.dc !== undefined) state.showCurve = !!C.dc;
      if (C.dp !== undefined) state.showPanels = !!C.dp;
      if (C.dpt !== undefined) state.showPoints = !!C.dpt;
      if (C.dct !== undefined) state.showControls = !!C.dct;
      if (C.dst !== undefined) state.showStreamlines = !!C.dst;
      if (C.dsm !== undefined) state.showSmoke = !!C.dsm;
      if (C.dps !== undefined) state.showPsiContours = !!C.dps;
      if (C.dcp !== undefined) state.showCp = !!C.dcp;
      if (C.df !== undefined) state.showForces = !!C.df;
      if (C.em !== undefined) state.enableMorphing = !!C.em;
      if (C.sla !== undefined) state.adaptiveStreamlines = !!C.sla;
    }

    state.morphDuration = C.h ?? C.md;
    state.streamlineDensity = C.i ?? C.sld;
    state.smokeDensity = C.j ?? C.smd;
    state.smokeParticlesPerBlob = C.k ?? C.smp;
    state.flowSpeed = C.f ?? C.fs;

    // Cp (new key o/q, legacy cpm/cpb)
    const cpMode = C.o ?? C.cpm;
    if (cpMode) {
      switch (cpMode) {
        case 'c': state.cpDisplayMode = 'color'; break;
        case 'b': state.cpDisplayMode = 'bars'; break;
        case 'x': state.cpDisplayMode = 'both'; break;
      }
    }
    state.cpBarScale = C.q ?? C.cpb;
    state.forceScale = C.r ?? C.fsc;
    state.thicknessScale = C.u ?? C.ts;
    state.camberScale = C.w ?? C.cs;
    state.curvatureWeight = C.g ?? C.cw;

    // Spacing panel (new packed `sp`, legacy spm/sspi/sspv)
    if (C.sp && C.sp.length === 3) {
      state.spacingPanelMode = C.sp[0] === 's' ? 'simple' : 'advanced';
      state.sspInterpolation = C.sp[1] === 'l' ? 'linear' : 'spline';
      state.sspVisualization = C.sp[2] === 'p' ? 'plot' : 'foil';
    } else {
      if (C.spm) state.spacingPanelMode = C.spm === 's' ? 'simple' : 'advanced';
      if (C.sspi) state.sspInterpolation = C.sspi === 'l' ? 'linear' : 'spline';
      if (C.sspv) state.sspVisualization = C.sspv === 'p' ? 'plot' : 'foil';
    }

    // Flaps (new compact `F` without id, legacy `fl` with id)
    if (C.F && Array.isArray(C.F)) {
      state.flaps = C.F.map(([name, hingeX, hingeYFrac, deflection], i) => ({
        id: `flap_${Date.now().toString(36)}_${i}`,
        name: String(name),
        hingeX: Number(hingeX),
        hingeYFrac: Number(hingeYFrac),
        deflection: Number(deflection),
      }));
    } else if (C.fl && Array.isArray(C.fl)) {
      state.flaps = C.fl.map(([id, name, hingeX, hingeYFrac, deflection]) => ({
        id: String(id),
        name: String(name),
        hingeX: Number(hingeX),
        hingeYFrac: Number(hingeYFrac),
        deflection: Number(deflection),
      }));
    }

    // Strip undefined values
    for (const key of Object.keys(state) as (keyof UrlState)[]) {
      if (state[key] === undefined) delete state[key];
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
function getShareableUrl(state: UrlState): string {
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
function getCompleteUrlState(
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
