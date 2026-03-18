/**
 * Shared flap-application pipeline.
 *
 * Used by both the reactive GeometryDesignPanel (live preview) and the
 * sweep engine (computing geometry per sweep step).
 */

import {
  deflectFlap,
  isWasmReady,
  repanelXfoil,
} from './wasm';
import type { AirfoilPoint, FlapDefinition } from '../types';

const HIRES_POINTS = 1000;

/**
 * Compute the signed area of a (possibly open) polygon contour.
 * Negative = clockwise (Selig/XFOIL convention), positive = counterclockwise.
 */
function signedArea(pts: { x: number; y: number }[]): number {
  let area = 0;
  const n = pts.length;
  for (let i = 0; i < n - 1; i++) {
    area += (pts[i + 1].x - pts[i].x) * (pts[i + 1].y + pts[i].y);
  }
  area += (pts[0].x - pts[n - 1].x) * (pts[0].y + pts[n - 1].y);
  return area * 0.5;
}

/**
 * Ensure `output` has the same winding direction as `reference`.
 * If they differ (one CW, the other CCW), reverse the output.
 * This preserves whatever convention the upstream pipeline uses,
 * without forcing a specific winding.
 */
function matchWinding<T extends { x: number; y: number }>(
  output: T[],
  reference: { x: number; y: number }[],
): T[] {
  const refSign = signedArea(reference);
  const outSign = signedArea(output);
  // Same sign → same winding → keep as-is.  Different sign → reverse.
  if ((refSign >= 0) !== (outSign >= 0)) {
    return output.slice().reverse();
  }
  return output;
}

/**
 * Apply an ordered list of flaps to a base airfoil and return both
 * high-res visual coordinates and solver-resolution panels.
 *
 * Pipeline:
 * 1. Repanel base to 1000 points (clean spline, no kink)
 * 2. Apply each flap rotation sequentially
 * 3. Return high-res result as `coordinates` and repaneled-to-nPanels as `panels`
 * 4. Ensure repaneling preserved the winding of the flapped geometry
 */
export function applyFlapsToBase(
  base: { x: number; y: number }[],
  flaps: FlapDefinition[],
  nPanels: number,
): { coordinates: AirfoilPoint[]; panels: AirfoilPoint[] } | null {
  if (!isWasmReady() || base.length < 4) return null;

  const hasDeflection = flaps.some(f => Math.abs(f.deflection) >= 0.01);

  if (flaps.length === 0 || !hasDeflection) {
    return repanelBoth(base, nPanels);
  }

  const hiRes = repanelXfoil(base, HIRES_POINTS);
  if (hiRes.length === 0) return null;

  let pts: { x: number; y: number }[] = hiRes;
  for (const flap of flaps) {
    if (Math.abs(flap.deflection) < 0.01) continue;
    const result = deflectFlap(pts, flap.hingeX, flap.hingeYFrac, flap.deflection);
    if (result.length > 0) pts = result;
  }

  // `pts` is the ground-truth flapped geometry from flap_impl (winding preserved
  // from the base). Use it as the reference for any repaneled output.
  const coordinates: AirfoilPoint[] = pts.map(p => ({ x: p.x, y: p.y }));

  const solverPanels = repanelXfoil(pts, nPanels);
  const panels: AirfoilPoint[] = solverPanels.length > 0
    ? matchWinding(solverPanels.map(pt => ({ x: pt.x, y: pt.y })), pts)
    : coordinates.slice();

  return { coordinates, panels };
}

/**
 * Repanel coordinates and return the same smooth result for both
 * `coordinates` (visual) and `panels` (solver).
 */
export function repanelBoth(
  coords: { x: number; y: number }[],
  nPanels: number,
): { coordinates: AirfoilPoint[]; panels: AirfoilPoint[] } {
  const raw: AirfoilPoint[] = coords.map(p => ({ x: p.x, y: p.y }));
  if (isWasmReady()) {
    try {
      const repaneled = repanelXfoil(coords, nPanels);
      if (repaneled.length > 0) {
        const smooth: AirfoilPoint[] = matchWinding(
          repaneled.map(pt => ({ x: pt.x, y: pt.y })),
          coords,
        );
        return { coordinates: smooth, panels: smooth };
      }
    } catch { /* keep raw */ }
  }
  return { coordinates: raw, panels: raw };
}
