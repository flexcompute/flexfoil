/**
 * Coupled multi-element aero seam.
 *
 * Builds the three deployed element contours (in the global frame) and runs them
 * through one coupled inviscid solve — element interaction is fully modelled.
 * Inviscid for now (Cd = 0; per-element Cl/Cp); the viscous/drag phase will swap
 * the WASM call behind this same signature.
 */

import { analyzeMultiElement } from '../lib/wasm';
import { buildConfiguration, toPoints } from './geometry';
import type { HighLiftAirfoil } from './estolConfig';

export type ElementName = 'main' | 'vane' | 'flaperon';

export interface ElementResult {
  name: ElementName;
  cl: number;
  cm: number;
  /** Pressure coefficient at each node. */
  cp: number[];
  /** Node x-coordinates (global frame), same length as cp. */
  cpX: number[];
}

const ORDER: ElementName[] = ['main', 'vane', 'flaperon'];

/** Coupled solve of the whole configuration at freestream angle `alphaDeg`.
 *  The operating point (deployment + flaperon angle) is carried in `cfg`. */
export function analyzeConfiguration(cfg: HighLiftAirfoil, alphaDeg: number): ElementResult[] {
  const { main, vane, flaperon } = buildConfiguration(cfg);
  const r = analyzeMultiElement([toPoints(main), toPoints(vane), toPoints(flaperon)], alphaDeg);
  return r.elements.map((e, i) => ({ name: ORDER[i], cl: e.cl, cm: e.cm, cp: e.cp, cpX: e.cpX }));
}
