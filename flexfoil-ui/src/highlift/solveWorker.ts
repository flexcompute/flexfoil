/**
 * Web Worker: runs the coupled multi-element solver off the main thread so the
 * UI stays responsive. Handles two request kinds — 'cp' (per-element Cp at one α)
 * and 'polar' (α-sweep of per-element Cl). Newer requests of the same kind
 * supersede older in-flight ones (the polar loop yields between α steps and bails
 * when a newer request has arrived).
 */

import { initWasm } from '../lib/wasm';
import { analyzeConfiguration } from './solve';
import type { HighLiftAirfoil } from './estolConfig';

const ready = initWasm();
let currentCp = 0;
let currentPolar = 0;

interface CpMsg { id: number; kind: 'cp'; cfg: HighLiftAirfoil; alphaDeg: number; }
interface PolarMsg { id: number; kind: 'polar'; cfg: HighLiftAirfoil; }

self.onmessage = async (e: MessageEvent<CpMsg | PolarMsg>) => {
  const msg = e.data;
  await ready;

  if (msg.kind === 'cp') {
    currentCp = msg.id;
    const res = analyzeConfiguration(msg.cfg, msg.alphaDeg);
    if (msg.id === currentCp) self.postMessage({ id: msg.id, kind: 'cp', res });
    return;
  }

  currentPolar = msg.id;
  const rows: Record<string, { alpha: number; cl: number }[]> = { main: [], vane: [], flaperon: [] };
  for (let a = -5; a <= 15.001; a += 2) {
    await new Promise((r) => setTimeout(r)); // let a newer request preempt this one
    if (msg.id !== currentPolar) return;
    for (const r of analyzeConfiguration(msg.cfg, a)) {
      rows[r.name].push({ alpha: a, cl: r.cl });
    }
  }
  if (msg.id === currentPolar) self.postMessage({ id: msg.id, kind: 'polar', rows });
};
