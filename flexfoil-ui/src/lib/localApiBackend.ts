/**
 * LocalAPIBackend — talks to the Python server's REST API.
 *
 * Used when the app is served by `flexfoil serve` (localhost:8420).
 * All data lives in a real SQLite file on disk.
 */

import type { RunRow, SolverMode, AirfoilPoint, FlapDefinition, RunGeometrySnapshot } from '../types';
import type { StorageBackend, RunInsert } from './storageBackend';
import { getApiBase } from './storageBackend';

function parsePointArray(json: unknown): AirfoilPoint[] | null {
  if (typeof json !== 'string' || json.trim() === '') return null;
  try {
    const parsed = JSON.parse(json);
    if (!Array.isArray(parsed)) return null;
    return parsed.filter(
      (p: any): p is AirfoilPoint =>
        typeof p === 'object' && p !== null && typeof p.x === 'number' && typeof p.y === 'number',
    );
  } catch {
    return null;
  }
}

function parseGeometrySnapshot(
  coordinatesJson: unknown,
  panelsJson: unknown,
): RunGeometrySnapshot | null {
  const coordinates = parsePointArray(coordinatesJson);
  const panels = parsePointArray(panelsJson);
  if (!coordinates || !panels) return null;
  return { coordinates, panels };
}

function parseFlapsJson(json: unknown): FlapDefinition[] | null {
  if (typeof json !== 'string' || json.trim() === '') return null;
  try {
    const parsed = JSON.parse(json);
    if (!Array.isArray(parsed)) return null;
    return parsed.filter(
      (f: any): f is FlapDefinition =>
        typeof f === 'object' && f !== null &&
        typeof f.hingeX === 'number' &&
        typeof f.deflection === 'number',
    );
  } catch {
    return null;
  }
}

function apiRowToRunRow(obj: Record<string, unknown>): RunRow {
  const flaps = parseFlapsJson(obj.flaps_json);
  return {
    id: obj.id as number,
    airfoil_name: obj.airfoil_name as string,
    airfoil_hash: obj.airfoil_hash as string,
    alpha: obj.alpha as number,
    reynolds: obj.reynolds as number,
    mach: obj.mach as number,
    ncrit: obj.ncrit as number,
    n_panels: obj.n_panels as number,
    max_iter: obj.max_iter as number,
    cl: obj.cl as number | null,
    cd: obj.cd as number | null,
    cm: obj.cm as number | null,
    converged: obj.converged === 1 || obj.converged === true,
    iterations: obj.iterations as number | null,
    residual: obj.residual as number | null,
    x_tr_upper: obj.x_tr_upper as number | null,
    x_tr_lower: obj.x_tr_lower as number | null,
    solver_mode: (obj.solver_mode as SolverMode | null) ?? 'viscous',
    success: obj.success === 1 || obj.success === true,
    error: obj.error as string | null,
    created_at: obj.created_at as string,
    session_id: obj.session_id as string | null,
    geometry_snapshot: parseGeometrySnapshot(obj.coordinates_json, obj.panels_json),
    flaps,
    ld:
      obj.cl != null && obj.cd != null && Math.abs(obj.cd as number) > 1e-10
        ? (obj.cl as number) / (obj.cd as number)
        : null,
    flap_deflection: flaps && flaps.length > 0 ? flaps[0].deflection : null,
    flap_hinge_x: flaps && flaps.length > 0 ? flaps[0].hingeX : null,
  };
}

function runInsertPayload(run: RunInsert): Record<string, unknown> {
  return {
    airfoil_name: run.airfoil_name,
    airfoil_hash: run.airfoil_hash,
    alpha: run.alpha,
    reynolds: run.reynolds,
    mach: run.mach,
    ncrit: run.ncrit,
    n_panels: run.n_panels,
    max_iter: run.max_iter,
    cl: run.cl,
    cd: run.cd,
    cm: run.cm,
    converged: run.converged,
    iterations: run.iterations,
    residual: run.residual,
    x_tr_upper: run.x_tr_upper,
    x_tr_lower: run.x_tr_lower,
    solver_mode: run.solver_mode,
    success: run.success,
    error: run.error,
    coordinates_json: run.coordinates_json,
    panels_json: run.panels_json,
    flaps_json: run.flaps_json,
  };
}

let _cachedRuns: RunRow[] = [];
let _changeCallbacks: Array<() => void> = [];
let _eventSource: EventSource | null = null;

function startSSE(): void {
  if (_eventSource) return;
  const base = getApiBase();
  _eventSource = new EventSource(`${base}/api/events`);

  const handleUpdate = () => {
    fetchAllRuns().then((rows) => {
      _cachedRuns = rows;
      for (const cb of _changeCallbacks) cb();
    });
  };

  _eventSource.addEventListener('run_added', handleUpdate);
  _eventSource.addEventListener('runs_cleared', handleUpdate);
  _eventSource.addEventListener('db_imported', handleUpdate);
}

async function fetchAllRuns(): Promise<RunRow[]> {
  const base = getApiBase();
  const resp = await fetch(`${base}/api/runs?limit=50000`);
  const data: Record<string, unknown>[] = await resp.json();
  return data.map(apiRowToRunRow);
}

export const localApiBackend: StorageBackend = {
  async init() {
    _cachedRuns = await fetchAllRuns();
    startSSE();
  },

  async insertRun(run: RunInsert): Promise<number> {
    const base = getApiBase();
    const resp = await fetch(`${base}/api/runs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(runInsertPayload(run)),
    });
    const { id } = await resp.json();
    _cachedRuns = await fetchAllRuns();
    return id;
  },

  async insertRunBatch(runs: RunInsert[]): Promise<number> {
    if (runs.length === 0) return 0;
    const base = getApiBase();
    // Try the batch endpoint; fall back to sequential POSTs.
    try {
      const resp = await fetch(`${base}/api/runs/batch`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(runs.map(runInsertPayload)),
      });
      if (resp.ok) {
        _cachedRuns = await fetchAllRuns();
        return runs.length;
      }
    } catch { /* batch endpoint not available */ }

    for (const run of runs) {
      await fetch(`${base}/api/runs`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(runInsertPayload(run)),
      });
    }
    _cachedRuns = await fetchAllRuns();
    return runs.length;
  },

  lookupCache(airfoilHash, alpha, reynolds, mach, ncrit, nPanels, maxIter): RunRow | null {
    return (
      _cachedRuns.find(
        (r) =>
          r.airfoil_hash === airfoilHash &&
          r.alpha === alpha &&
          r.reynolds === reynolds &&
          r.mach === mach &&
          r.ncrit === ncrit &&
          r.n_panels === nPanels &&
          r.max_iter === maxIter,
      ) ?? null
    );
  },

  queryAllRuns(): RunRow[] {
    return _cachedRuns;
  },

  async updateRunAirfoilName(_id: number, _newName: string): Promise<void> {
    // TODO: implement PATCH endpoint
    _cachedRuns = await fetchAllRuns();
  },

  async clearAllRuns(): Promise<void> {
    const base = getApiBase();
    await fetch(`${base}/api/runs`, { method: 'DELETE' });
    _cachedRuns = [];
  },

  async pruneOldRuns(): Promise<number> {
    return 0; // server handles pruning
  },

  exportDatabase(): Uint8Array {
    // For local mode, we direct the user to the /api/db/export endpoint
    throw new Error('Use /api/db/export endpoint directly for local mode export');
  },

  async importDatabase(_data: Uint8Array): Promise<void> {
    const base = getApiBase();
    await fetch(`${base}/api/db/import`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/octet-stream' },
      body: _data,
    });
    _cachedRuns = await fetchAllRuns();
  },

  onExternalChange(callback: () => void): () => void {
    _changeCallbacks.push(callback);
    return () => {
      _changeCallbacks = _changeCallbacks.filter((cb) => cb !== callback);
    };
  },
};
