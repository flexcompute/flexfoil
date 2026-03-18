/**
 * BrowserBackend — wraps the existing sql.js + IndexedDB storage.
 *
 * This is the default backend when the app is hosted on the web (foil.flexcompute.com).
 * Zero changes to existing behavior.
 */

import type { RunRow } from '../types';
import type { StorageBackend, RunInsert } from './storageBackend';
import {
  initRunDatabase,
  insertRun as dbInsertRun,
  insertRunBatch as dbInsertRunBatch,
  lookupCache as dbLookupCache,
  queryAllRuns as dbQueryAllRuns,
  clearAllRuns as dbClearAllRuns,
  pruneOldRuns as dbPruneOldRuns,
  exportDatabase as dbExportDatabase,
  importDatabase as dbImportDatabase,
  updateRunAirfoilName as dbUpdateRunAirfoilName,
  ensureRunDatabase,
} from './runDatabase';

export const browserBackend: StorageBackend = {
  async init() {
    await initRunDatabase();
  },

  async insertRun(run: RunInsert): Promise<number> {
    await ensureRunDatabase();
    return dbInsertRun(run);
  },

  async insertRunBatch(runs: RunInsert[]): Promise<number> {
    await ensureRunDatabase();
    return dbInsertRunBatch(runs);
  },

  lookupCache(airfoilHash, alpha, reynolds, mach, ncrit, nPanels, maxIter): RunRow | null {
    return dbLookupCache(airfoilHash, alpha, reynolds, mach, ncrit, nPanels, maxIter);
  },

  queryAllRuns(): RunRow[] {
    return dbQueryAllRuns();
  },

  async updateRunAirfoilName(id: number, newName: string): Promise<void> {
    await dbUpdateRunAirfoilName(id, newName);
  },

  async clearAllRuns(): Promise<void> {
    await dbClearAllRuns();
  },

  async pruneOldRuns(): Promise<number> {
    return dbPruneOldRuns();
  },

  exportDatabase(): Uint8Array {
    return dbExportDatabase();
  },

  async importDatabase(data: Uint8Array): Promise<void> {
    await dbImportDatabase(data);
  },
};
