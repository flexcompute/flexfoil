/**
 * Zustand store for run data management.
 *
 * Bridges the storage backend and React UI — exposes rows,
 * filtered rows (from AG Grid), and mutation helpers.
 *
 * The backend is selected at startup: sql.js+IndexedDB (browser) or
 * REST API to the local Python server (local mode).
 */

import { create } from 'zustand';
import type { RunRow } from '../types';
import type { RunInsert } from '../lib/storageBackend';
export type { RunInsert } from '../lib/storageBackend';
import { getBackend } from '../lib/activeBackend';
import { computeAirfoilHash } from '../lib/airfoilHash';
import { pauseHistory, resumeHistory, useAirfoilStore } from './airfoilStore';

interface RunStoreState {
  /** All rows from the database, newest first */
  allRuns: RunRow[];
  /** Subset filtered by the AG Grid (defaults to allRuns) */
  filteredRuns: RunRow[];
  /** Most recently restored historical run */
  selectedRunId: number | null;
  /** Increments when a historical run is restored */
  restorationRevision: number;
  /** Whether the DB has been initialised */
  ready: boolean;

  /** Initialise sql.js + load persisted DB from IndexedDB */
  init: () => Promise<void>;
  /** Refresh allRuns from the DB */
  refresh: () => void;
  /** Update the filtered subset (called by DataExplorer on filter change) */
  setFilteredRuns: (rows: RunRow[]) => void;

  /** Insert a new run (returns the cached-or-new row) */
  addRun: (run: RunInsert) => Promise<void>;
  /** Batch-insert runs in one transaction, refresh once at the end */
  addRunBatch: (runs: RunInsert[]) => Promise<void>;
  /** Look up a cached result; returns null on miss */
  lookup: (
    airfoilHash: string,
    alpha: number,
    reynolds: number,
    mach: number,
    ncrit: number,
    nPanels: number,
    maxIter: number,
  ) => RunRow | null;
  /** Compute a deterministic hash for a set of panels */
  hashPanels: (panels: { x: number; y: number }[]) => Promise<string>;
  /** Restore a persisted run snapshot into the shared airfoil state */
  restoreRunById: (id: number) => boolean;
  /** Rename the airfoil_name for a single run in the DB */
  renameRun: (id: number, newName: string) => Promise<void>;

  /** Delete every row and refresh */
  clearAll: () => Promise<void>;
  /** Export the raw .sqlite bytes */
  exportDb: () => Uint8Array;
  /** Import a .sqlite file, replacing current DB */
  importDb: (data: Uint8Array) => Promise<void>;
}

export const useRunStore = create<RunStoreState>()((set, get) => ({
  allRuns: [],
  filteredRuns: [],
  selectedRunId: null,
  restorationRevision: 0,
  ready: false,

  init: async () => {
    if (get().ready) return;
    const backend = getBackend();
    await backend.init();
    const rows = backend.queryAllRuns();
    set({ allRuns: rows, filteredRuns: rows, ready: true });

    // Subscribe to external changes (e.g. SSE from local Python server)
    backend.onExternalChange?.(() => get().refresh());
  },

  refresh: () => {
    const rows = getBackend().queryAllRuns();
    set({ allRuns: rows, filteredRuns: rows });
  },

  setFilteredRuns: (rows) => set({ filteredRuns: rows }),

  addRun: async (run) => {
    const backend = getBackend();
    if (!get().ready) {
      await backend.init();
      set({ ready: true });
    }
    await backend.insertRun(run);
    await backend.pruneOldRuns();
    get().refresh();
  },

  addRunBatch: async (runs) => {
    if (runs.length === 0) return;
    const backend = getBackend();
    if (!get().ready) {
      await backend.init();
      set({ ready: true });
    }
    await backend.insertRunBatch(runs);
    await backend.pruneOldRuns();
    get().refresh();
  },

  lookup: (hash, alpha, re, mach, ncrit, nPanels, maxIter) => {
    try {
      return getBackend().lookupCache(hash, alpha, re, mach, ncrit, nPanels, maxIter);
    } catch {
      return null;
    }
  },

  hashPanels: (panels) => computeAirfoilHash(panels),
  restoreRunById: (id) => {
    const run = get().allRuns.find((row) => row.id === id);
    if (!run?.geometry_snapshot) return false;
    pauseHistory();
    try {
      useAirfoilStore.getState().restoreRunSnapshot(run);
    } finally {
      resumeHistory();
    }
    set((state) => ({
      selectedRunId: id,
      restorationRevision: state.restorationRevision + 1,
    }));
    return true;
  },

  renameRun: async (id, newName) => {
    await getBackend().updateRunAirfoilName(id, newName);
    get().refresh();
  },

  clearAll: async () => {
    await getBackend().clearAllRuns();
    set({ allRuns: [], filteredRuns: [], selectedRunId: null, restorationRevision: 0 });
  },

  exportDb: () => getBackend().exportDatabase(),

  importDb: async (data) => {
    await getBackend().importDatabase(data);
    get().refresh();
  },
}));
