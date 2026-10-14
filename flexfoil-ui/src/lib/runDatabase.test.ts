import { beforeAll, expect, test, vi } from 'vitest';
import type { RunInsert } from './runDatabase';
vi.mock('idb', () => ({openDB: async () => ({get: async () => undefined, put: async () => undefined})}));
vi.mock('sql.js', async (original) => {
  const sql = await original<typeof import('sql.js')>();
  return {...sql, default: () => sql.default({locateFile: () => 'node_modules/sql.js/dist/sql-wasm.wasm'})};
});
import initSqlJs from 'sql.js';
import { initRunDatabase, solverCacheKey, insertRun, lookupCache, queryAllRuns, clearAllRuns, exportDatabase, importDatabase } from './runDatabase';
const key = solverCacheKey('v1:content-hashed.wasm', 'viscous', 1, 1, 1)!;
const run: RunInsert = {airfoil_name:'foil', airfoil_hash:'geometry', alpha:4, reynolds:1e6, mach:.25, ncrit:9, n_panels:160, max_iter:100, cl:.15, cd:.007, cm:0, converged:false, iterations:100, residual:1, x_tr_upper:1, x_tr_lower:1, success:true, error:null, solver_mode:'viscous', coordinates_json:null, panels_json:null, flaps_json:null, cache_key:key};
const lookup = (cacheKey: string | null = key) => lookupCache('geometry',4,1e6,.25,9,160,100,cacheKey);
beforeAll(initRunDatabase);

test('legacy history and failed attempts survive a corrected solve and roundtrip', async () => {
  await clearAllRuns();
  await insertRun({...run,cache_key:null});
  const SQL = await initSqlJs();
  const old = new SQL.Database(exportDatabase());
  old.run('DROP INDEX idx_cache_lookup; ALTER TABLE runs DROP COLUMN cache_key; CREATE UNIQUE INDEX idx_cache_key ON runs(airfoil_hash,alpha,reynolds,mach,ncrit,n_panels,max_iter)');
  await importDatabase(old.export());
  old.close();
  expect(queryAllRuns()[0].cl).toBe(.15);
  expect(lookup()).toBeNull();
  await insertRun(run);
  expect(lookup()).toBeNull();
  const fresh = await insertRun({...run,converged:true,cl:.4432});
  expect(lookup()?.id).toBe(fresh);
  expect(lookup()?.cl).toBe(.4432);
  expect(queryAllRuns()).toHaveLength(3);
  await importDatabase(exportDatabase());
  expect(queryAllRuns()).toHaveLength(3);
  expect(lookup()?.cl).toBe(.4432);
  expect(lookup(null)).toBeNull();
  for (const different of [solverCacheKey('new-engine','viscous',1,1,1), solverCacheKey('v1:content-hashed.wasm','viscous',2,1,1), solverCacheKey('v1:content-hashed.wasm','viscous',1,.3,1), solverCacheKey('v1:content-hashed.wasm','viscous',1,1,.4), solverCacheKey('v1:content-hashed.wasm','inviscid',1,1,1)]) expect(lookup(different)).toBeNull();
  expect(solverCacheKey(null,'viscous',1,1,1)).toBeNull();
});

test('invalid imports cannot replace existing history', async () => {
  const before = queryAllRuns();
  await expect(importDatabase(new TextEncoder().encode('not SQLite'))).rejects.toThrow();
  expect(queryAllRuns()).toEqual(before);
});
