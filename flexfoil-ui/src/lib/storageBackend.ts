/**
 * Storage backend abstraction.
 *
 * Lets the app run in two modes:
 * - **BrowserBackend** (default): sql.js + IndexedDB, fully offline.
 * - **LocalAPIBackend**: REST calls to the local Python server (`flexfoil serve`).
 *
 * The active backend is chosen once at startup and stays for the session.
 */

import type { RunRow, AirfoilPoint, RunGeometrySnapshot, SolverMode } from '../types';

export interface RunInsert {
  airfoil_name: string;
  airfoil_hash: string;
  alpha: number;
  reynolds: number;
  mach: number;
  ncrit: number;
  n_panels: number;
  max_iter: number;
  cl: number | null;
  cd: number | null;
  cm: number | null;
  converged: boolean;
  iterations: number | null;
  residual: number | null;
  x_tr_upper: number | null;
  x_tr_lower: number | null;
  success: boolean;
  error: string | null;
  solver_mode: SolverMode;
  coordinates_json: string | null;
  panels_json: string | null;
  flaps_json: string | null;
}

export interface StorageBackend {
  init(): Promise<void>;
  insertRun(run: RunInsert): Promise<number>;
  /** Batch-insert multiple runs in a single transaction / persist cycle. */
  insertRunBatch(runs: RunInsert[]): Promise<number>;
  lookupCache(
    airfoilHash: string,
    alpha: number,
    reynolds: number,
    mach: number,
    ncrit: number,
    nPanels: number,
    maxIter: number,
  ): RunRow | null;
  queryAllRuns(): RunRow[];
  updateRunAirfoilName(id: number, newName: string): Promise<void>;
  clearAllRuns(): Promise<void>;
  pruneOldRuns(): Promise<number>;
  exportDatabase(): Uint8Array;
  importDatabase(data: Uint8Array): Promise<void>;
  onExternalChange?(callback: () => void): () => void;
}

// ---------------------------------------------------------------------------
// Detection: is the local Python server serving us?
// ---------------------------------------------------------------------------

let _isLocal: boolean | null = null;

export function isLocalMode(): boolean {
  if (_isLocal !== null) return _isLocal;
  if (typeof window === 'undefined') {
    _isLocal = false;
    return false;
  }
  // Explicit opt-in via query param or data attribute
  if (
    window.location.search.includes('local=1') ||
    document.documentElement.dataset.flexfoilLocal === 'true'
  ) {
    _isLocal = true;
    return true;
  }
  // The Python server injects a meta tag into the served index.html
  const meta = document.querySelector('meta[name="flexfoil-local"]');
  if (meta) {
    _isLocal = true;
    return true;
  }
  // Heuristic: localhost/127.0.0.1 with a non-standard port (not 5173/4173
  // which are Vite dev/preview) is very likely the Python server.
  const { hostname, port } = window.location;
  const isLocalhost = hostname === 'localhost' || hostname === '127.0.0.1';
  const isDevServer = port === '5173' || port === '4173' || port === '4174' || port === '4175' || port === '4176' || port === '';
  _isLocal = isLocalhost && !isDevServer;
  return _isLocal;
}

export function getApiBase(): string {
  return `${window.location.protocol}//${window.location.host}`;
}
