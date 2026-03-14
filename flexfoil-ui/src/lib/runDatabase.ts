/**
 * Browser-local SQLite database for solver run history & caching.
 *
 * Uses sql.js (SQLite compiled to WASM) with IndexedDB for persistence.
 */

import initSqlJs, { type Database } from 'sql.js';
import { openDB, type IDBPDatabase } from 'idb';
import type { RunRow } from '../types';
import sqlWasmUrl from 'sql.js/dist/sql-wasm.wasm?url';

const DB_NAME = 'flexfoil-runs';
const IDB_STORE = 'sqlite';
const IDB_KEY = 'db';
const MAX_ROWS = 50_000;

let db: Database | null = null;
let idb: IDBPDatabase | null = null;
let sessionId: string | null = null;

function getSessionId(): string {
  if (!sessionId) {
    sessionId = `s_${Date.now().toString(36)}_${Math.random().toString(36).slice(2, 6)}`;
  }
  return sessionId;
}

async function getIdb(): Promise<IDBPDatabase> {
  if (!idb) {
    idb = await openDB(DB_NAME, 1, {
      upgrade(db) {
        if (!db.objectStoreNames.contains(IDB_STORE)) {
          db.createObjectStore(IDB_STORE);
        }
      },
    });
  }
  return idb;
}

async function persistToIdb(): Promise<void> {
  if (!db) return;
  const store = await getIdb();
  const data = db.export();
  await store.put(IDB_STORE, data, IDB_KEY);
}

const SCHEMA = `
CREATE TABLE IF NOT EXISTS runs (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  airfoil_name  TEXT NOT NULL,
  airfoil_hash  TEXT NOT NULL,
  alpha         REAL NOT NULL,
  reynolds      REAL NOT NULL,
  mach          REAL NOT NULL,
  ncrit         REAL NOT NULL,
  n_panels      INTEGER NOT NULL,
  max_iter      INTEGER NOT NULL,
  cl            REAL,
  cd            REAL,
  cm            REAL,
  converged     INTEGER NOT NULL DEFAULT 0,
  iterations    INTEGER,
  residual      REAL,
  x_tr_upper    REAL,
  x_tr_lower    REAL,
  success       INTEGER NOT NULL DEFAULT 0,
  error         TEXT,
  created_at    TEXT DEFAULT (datetime('now')),
  session_id    TEXT
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_cache_key
  ON runs(airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter);
`;

export async function initRunDatabase(): Promise<void> {
  if (db) return;

  const SQL = await initSqlJs({
    locateFile: () => sqlWasmUrl,
  });

  const store = await getIdb();
  const saved = await store.get(IDB_STORE, IDB_KEY);

  if (saved) {
    db = new SQL.Database(new Uint8Array(saved as ArrayBuffer));
    db.run(SCHEMA);
  } else {
    db = new SQL.Database();
    db.run(SCHEMA);
  }
}

let initPromise: Promise<void> | null = null;

/** Ensure the DB is ready, coalescing concurrent callers. */
export async function ensureRunDatabase(): Promise<void> {
  if (db) return;
  if (!initPromise) initPromise = initRunDatabase();
  await initPromise;
}

function getDb(): Database {
  if (!db) throw new Error('RunDatabase not initialised — call ensureRunDatabase() first');
  return db;
}

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
}

export async function insertRun(run: RunInsert): Promise<number> {
  const d = getDb();
  d.run(
    `INSERT OR IGNORE INTO runs
       (airfoil_name, airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter,
        cl, cd, cm, converged, iterations, residual, x_tr_upper, x_tr_lower,
        success, error, session_id)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)`,
    [
      run.airfoil_name,
      run.airfoil_hash,
      run.alpha,
      run.reynolds,
      run.mach,
      run.ncrit,
      run.n_panels,
      run.max_iter,
      run.cl,
      run.cd,
      run.cm,
      run.converged ? 1 : 0,
      run.iterations,
      run.residual,
      run.x_tr_upper,
      run.x_tr_lower,
      run.success ? 1 : 0,
      run.error,
      getSessionId(),
    ]
  );
  await persistToIdb();

  const result = d.exec('SELECT last_insert_rowid() as id');
  return (result[0]?.values[0]?.[0] as number) ?? -1;
}

export function lookupCache(
  airfoilHash: string,
  alpha: number,
  reynolds: number,
  mach: number,
  ncrit: number,
  nPanels: number,
  maxIter: number,
): RunRow | null {
  const d = getDb();
  const rows = d.exec(
    `SELECT * FROM runs
     WHERE airfoil_hash = ? AND alpha = ? AND reynolds = ? AND mach = ?
       AND ncrit = ? AND n_panels = ? AND max_iter = ?
     LIMIT 1`,
    [airfoilHash, alpha, reynolds, mach, ncrit, nPanels, maxIter]
  );
  if (!rows.length || !rows[0].values.length) return null;
  return rowToRunRow(rows[0].columns, rows[0].values[0]);
}

export function queryAllRuns(): RunRow[] {
  const d = getDb();
  const rows = d.exec('SELECT * FROM runs ORDER BY id DESC');
  if (!rows.length) return [];
  return rows[0].values.map((v: (string | number | null | Uint8Array)[]) => rowToRunRow(rows[0].columns, v));
}

export function getRowCount(): number {
  const d = getDb();
  const r = d.exec('SELECT COUNT(*) FROM runs');
  return (r[0]?.values[0]?.[0] as number) ?? 0;
}

export async function pruneOldRuns(): Promise<number> {
  const count = getRowCount();
  if (count <= MAX_ROWS) return 0;
  const excess = count - MAX_ROWS;
  const d = getDb();
  d.run(
    `DELETE FROM runs WHERE id IN (SELECT id FROM runs ORDER BY id ASC LIMIT ?)`,
    [excess]
  );
  await persistToIdb();
  return excess;
}

export async function clearAllRuns(): Promise<void> {
  const d = getDb();
  d.run('DELETE FROM runs');
  await persistToIdb();
}

export function exportDatabase(): Uint8Array {
  return getDb().export();
}

export async function importDatabase(data: Uint8Array): Promise<void> {
  const SQL = await initSqlJs({
    locateFile: () => sqlWasmUrl,
  });
  if (db) db.close();
  db = new SQL.Database(data);
  db.run(SCHEMA);
  await persistToIdb();
}

function rowToRunRow(columns: string[], values: (string | number | null | Uint8Array)[]): RunRow {
  const obj: Record<string, unknown> = {};
  columns.forEach((col, i) => {
    obj[col] = values[i];
  });
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
    converged: (obj.converged as number) === 1,
    iterations: obj.iterations as number | null,
    residual: obj.residual as number | null,
    x_tr_upper: obj.x_tr_upper as number | null,
    x_tr_lower: obj.x_tr_lower as number | null,
    success: (obj.success as number) === 1,
    error: obj.error as string | null,
    created_at: obj.created_at as string,
    session_id: obj.session_id as string | null,
  };
}
