/**
 * Zustand store for run data management.
 *
 * Bridges the SQLite database and React UI — exposes rows,
 * filtered rows (from AG Grid), and mutation helpers.
 */

import { create } from 'zustand';
import type { RunRow } from '../types';
import {
  initRunDatabase,
  ensureRunDatabase,
  insertRun,
  lookupCache,
  queryAllRuns,
  clearAllRuns,
  pruneOldRuns,
  exportDatabase,
  importDatabase,
  type RunInsert,
} from '../lib/runDatabase';
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
    await initRunDatabase();
    const rows = queryAllRuns();
    set({ allRuns: rows, filteredRuns: rows, ready: true });
  },

  refresh: () => {
    const rows = queryAllRuns();
    set({ allRuns: rows, filteredRuns: rows });
  },

  setFilteredRuns: (rows) => set({ filteredRuns: rows }),

  addRun: async (run) => {
    await ensureRunDatabase();
    if (!get().ready) {
      const rows = queryAllRuns();
      set({ allRuns: rows, filteredRuns: rows, ready: true });
    }
    await insertRun(run);
    await pruneOldRuns();
    get().refresh();
  },

  lookup: (hash, alpha, re, mach, ncrit, nPanels, maxIter) => {
    try {
      return lookupCache(hash, alpha, re, mach, ncrit, nPanels, maxIter);
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

  clearAll: async () => {
    await clearAllRuns();
    set({ allRuns: [], filteredRuns: [], selectedRunId: null, restorationRevision: 0 });
  },

  exportDb: () => exportDatabase(),

  importDb: async (data) => {
    await importDatabase(data);
    get().refresh();
  },
}));
