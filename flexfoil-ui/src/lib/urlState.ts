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
  mode: 'surface' as ControlMode,
  alpha: 0,
  polarStart: -5,
  polarEnd: 15,
  polarStep: 1,
  theme: 'dark' as 'dark' | 'light',
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
}

/**
 * Compact state for encoding (short keys, omit defaults)
 */
interface CompactState {
  n?: string;      // naca
  c?: 1;           // custom (flag)
  p?: number;      // nPanels
  s?: number[][];  // spacing [[S,F], ...]
  m?: string;      // mode (first char: s/b/b)
  a?: number;      // alpha
  ps?: number;     // polar start
  pe?: number;     // polar end
  pd?: number;     // polar delta
  t?: string;      // theme (d/l)
  l?: object;      // layout
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
    compact.m = state.mode[0]; // 's', 'b', or 'b' (bezier/bspline distinguished later)
    if (state.mode === 'bspline') compact.m = 'p'; // use 'p' for bspline
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
        case 's': state.mode = 'surface'; break;
        case 'b': state.mode = 'bezier'; break;
        case 'p': state.mode = 'bspline'; break;
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
