/**
 * Geometry engine for the 3-element STOL high-lift airfoil.
 *
 * Faithful TypeScript port of the geometric primitives in
 * ~/himalaya/geometry/airfoils/build_estol_geometry.py. Each builder returns a
 * Selig-ordered closed polyline (the `flat` contour the Python writes via
 * write_selig), validated against the reference .dat files in
 * ./__fixtures__/ by geometry.test.ts to ~1e-6.
 *
 * Conventions match the reference: points are [x, y] tuples; main contour runs
 * upper-lip-tip → LE → lower-cut-pt → cove; NACA elements run upper TE corner →
 * LE → lower TE corner. Use toPoints() to convert to the {x, y} form the UI and
 * WASM glue expect.
 */

import type {
  AirfoilCoords,
  CoveCutouts,
  HighLiftAirfoil,
  NacaElementCfg,
  TrackConfig,
  V2,
} from './estolConfig';

// ---------------------------------------------------------------------------
// Vector / array helpers (numpy-equivalent semantics)
// ---------------------------------------------------------------------------

const sub = (a: V2, b: V2): V2 => [a[0] - b[0], a[1] - b[1]];
const add = (a: V2, b: V2): V2 => [a[0] + b[0], a[1] + b[1]];
const scale = (a: V2, s: number): V2 => [a[0] * s, a[1] * s];
const dot = (a: V2, b: V2): number => a[0] * b[0] + a[1] * b[1];
const norm = (v: V2): number => Math.hypot(v[0], v[1]);
const unit = (v: V2): V2 => scale(v, 1 / norm(v));

/** Rotate points by `angleRad` (CCW) about `center`. Matches numpy `(p-c)@R.T+c`. */
function rotate(points: V2[], angleRad: number, center: V2 = [0, 0]): V2[] {
  const c = Math.cos(angleRad);
  const s = Math.sin(angleRad);
  return points.map(([x, y]) => {
    const dx = x - center[0];
    const dy = y - center[1];
    return [dx * c - dy * s + center[0], dx * s + dy * c + center[1]] as V2;
  });
}

/** Inclusive linspace of `n` samples from `a` to `b` (numpy.linspace). */
function linspace(a: number, b: number, n: number): number[] {
  if (n === 1) return [a];
  const step = (b - a) / (n - 1);
  const out = new Array<number>(n);
  for (let i = 0; i < n; i++) out[i] = a + step * i;
  out[n - 1] = b; // exact endpoint, like numpy
  return out;
}

/**
 * Scalar linear interpolation matching numpy.interp: `xp` must be ascending;
 * values outside the range clamp to the endpoints.
 */
function interp(x: number, xp: number[], fp: number[]): number {
  const n = xp.length;
  if (x <= xp[0]) return fp[0];
  if (x >= xp[n - 1]) return fp[n - 1];
  let i = 1;
  while (i < n && xp[i] < x) i++;
  const x0 = xp[i - 1];
  const x1 = xp[i];
  const t = (x - x0) / (x1 - x0);
  return fp[i - 1] + (fp[i] - fp[i - 1]) * t;
}

/** Convert internal [x, y] tuples to the {x, y} form used by the UI / WASM glue. */
export function toPoints(contour: V2[]): { x: number; y: number }[] {
  return contour.map(([x, y]) => ({ x, y }));
}

// ---------------------------------------------------------------------------
// NACA 4-digit airfoil (unit chord, LE at origin, TE near (1,0))
// ---------------------------------------------------------------------------

interface NacaSurfaces {
  xu: number[];
  yu: number[];
  xl: number[];
  yl: number[];
}

function nacaSurfaces(xs: number[], m: number, p: number, t: number): NacaSurfaces {
  const xu: number[] = [];
  const yu: number[] = [];
  const xl: number[] = [];
  const yl: number[] = [];
  for (const x of xs) {
    const yt =
      5 * t *
      (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1015 * x ** 4);
    let yc: number;
    let dyc: number;
    if (x < p) {
      yc = (m / (p * p)) * (2 * p * x - x * x);
      dyc = (2 * m / (p * p)) * (p - x);
    } else {
      const q = (1 - p) * (1 - p);
      yc = (m / q) * (1 - 2 * p + 2 * p * x - x * x);
      dyc = (2 * m / q) * (p - x);
    }
    const theta = Math.atan(dyc);
    const st = Math.sin(theta);
    const ct = Math.cos(theta);
    xu.push(x - yt * st);
    yu.push(yc + yt * ct);
    xl.push(x + yt * st);
    yl.push(yc - yt * ct);
  }
  return { xu, yu, xl, yl };
}

/**
 * Selig-ordered NACA 4-digit coordinates (unit chord) with an optional blunt TE
 * of vertical thickness `teThickness`. Port of `naca4()`.
 */
export function naca4(code: string, teThickness = 0.0, nPanels = 160): V2[] {
  const m = parseInt(code[0], 10) / 100.0;
  const p = parseInt(code[1], 10) / 10.0;
  const t = parseInt(code.slice(2, 4), 10) / 100.0;

  // Find the chordwise cut where the surface-to-surface gap equals teThickness.
  const xScan = linspace(0.5, 1.0, 5001);
  const s = nacaSurfaces(xScan, m, p, t);
  const gap = xScan.map((_, i) => Math.hypot(s.xu[i] - s.xl[i], s.yu[i] - s.yl[i]));
  // numpy interp requires ascending xp: gap decreases toward TE, so reverse both.
  const gapRev = gap.slice().reverse();
  const xScanRev = xScan.slice().reverse();
  const xCut = interp(teThickness, gapRev, xScanRev);

  // Cosine-clustered chord stations scaled to the cut.
  const beta = linspace(0.0, Math.PI, Math.floor(nPanels / 2) + 1);
  const x = beta.map((b) => 0.5 * (1.0 - Math.cos(b)) * xCut);
  const surf = nacaSurfaces(x, m, p, t);

  const out: V2[] = [];
  // upper TE corner → LE (reversed)
  for (let i = x.length - 1; i >= 0; i--) out.push([surf.xu[i], surf.yu[i]]);
  // LE → lower TE corner (skip duplicated LE point)
  for (let i = 1; i < x.length; i++) out.push([surf.xl[i], surf.yl[i]]);
  return out;
}

// ---------------------------------------------------------------------------
// Anchor a unit-chord airfoil between explicit LE/TE coordinates
// ---------------------------------------------------------------------------

/** Port of `anchor_airfoil()`. */
export function anchorAirfoil(raw: V2[], targetLe: V2, targetTe: V2): V2[] {
  // raw LE = min-x point; raw TE = midpoint of first/last (the two TE corners).
  let leIdx = 0;
  for (let i = 1; i < raw.length; i++) if (raw[i][0] < raw[leIdx][0]) leIdx = i;
  const rawLe = raw[leIdx];
  const rawTe: V2 = [0.5 * (raw[0][0] + raw[raw.length - 1][0]), 0.5 * (raw[0][1] + raw[raw.length - 1][1])];
  const rawVec = sub(rawTe, rawLe);
  const tgtVec = sub(targetTe, targetLe);
  const sc = norm(tgtVec) / norm(rawVec);
  const angle = Math.atan2(tgtVec[1], tgtVec[0]) - Math.atan2(rawVec[1], rawVec[0]);
  const translated = raw.map((pt) => sub(pt, rawLe));
  const rotated = rotate(translated, angle);
  return rotated.map((pt) => add(scale(pt, sc), targetLe));
}

// ---------------------------------------------------------------------------
// Cove fillet
// ---------------------------------------------------------------------------

export interface FilletResult {
  /** Tangent point on the leg toward `pIn`. */
  tIn: V2;
  /** Arc samples from the `pIn` tangent to the `pOut` tangent. */
  arc: V2[];
  /** Tangent point on the leg toward `pOut`. */
  tOut: V2;
}

/** Port of `fillet_corner()`: round the corner at `vertex` with `radius`. */
export function filletCorner(pIn: V2, vertex: V2, pOut: V2, radius: number, nArc = 40): FilletResult {
  const v1 = unit(sub(pIn, vertex));
  const v2 = unit(sub(pOut, vertex));
  const alpha = Math.acos(Math.max(-1.0, Math.min(1.0, dot(v1, v2))));
  const half = alpha / 2.0;
  const tIn = add(vertex, scale(v1, radius / Math.tan(half)));
  const tOut = add(vertex, scale(v2, radius / Math.tan(half)));
  const center = add(vertex, scale(unit(add(v1, v2)), radius / Math.sin(half)));

  const aIn = Math.atan2(tIn[1] - center[1], tIn[0] - center[0]);
  const aOut = Math.atan2(tOut[1] - center[1], tOut[0] - center[0]);
  // Wrap the swept angle into (-pi, pi], matching the Python modulo expression.
  const delta = mod(aOut - aIn + Math.PI, 2 * Math.PI) - Math.PI;
  const arc = linspace(0.0, delta, nArc).map((d) => {
    const a = aIn + d;
    return [center[0] + radius * Math.cos(a), center[1] + radius * Math.sin(a)] as V2;
  });
  return { tIn, arc, tOut };
}

/** Python-style modulo (result has the sign of the divisor). */
function mod(a: number, b: number): number {
  return ((a % b) + b) % b;
}

// ---------------------------------------------------------------------------
// Element builders
// ---------------------------------------------------------------------------

/** Coved main wing. Port of `build_main_wing()` (returns the `flat` contour). */
export function buildMainWing(airfoil: AirfoilCoords, cutouts: CoveCutouts, bluntThickness: number): V2[] {
  const upper = airfoil.upper;
  const lower = airfoil.lower;
  const xLoCut = cutouts.lowerCutX;
  const xUpCut = cutouts.upperLipCutX;
  const vertex: V2 = [cutouts.coveVertexX, cutouts.coveVertexY];
  const rFillet = cutouts.coveFilletRadius;

  const upX = upper.map((q) => q[0]);
  const upY = upper.map((q) => q[1]);
  const loX = lower.map((q) => q[0]);
  const loY = lower.map((q) => q[1]);
  const yUpAtCut = interp(xUpCut, upX, upY);
  const yLoAtCut = interp(xLoCut, loX, loY);

  const upKept: V2[] = [...upper.filter((q) => q[0] < xUpCut), [xUpCut, yUpAtCut]];
  const loKept: V2[] = [...lower.filter((q) => q[0] < xLoCut), [xLoCut, yLoAtCut]];
  const lipBottom: V2 = [xUpCut, yUpAtCut - bluntThickness];
  const lowerCutPt: V2 = [xLoCut, yLoAtCut];

  // fillet_corner returns tangent points on each leg: tIn toward lip_bottom
  // (upper side), tOut toward lower_cut_pt (lower side). The arc traces
  // lip-side → lower-side.
  const { arc, tOut: tLower } = filletCorner(lipBottom, vertex, lowerCutPt, rFillet);
  // main_contour: upper-lip-tip → LE → lower-cut-pt
  const mainContour: V2[] = [...upKept.slice().reverse(), ...loKept.slice(1)];
  // Reverse so the flat traversal runs lower-side → upper-side.
  const arcLowerToUpper = arc.slice().reverse();

  return [...mainContour, tLower, ...arcLowerToUpper.slice(1), lipBottom];
}

/** Vane / aft flap: anchored NACA airfoil. Port of `build_naca_element()` (`flat`). */
export function buildNacaElement(
  code: string,
  teThickness: number,
  targetLe: V2,
  targetTe: V2,
  nPanels = 160,
): V2[] {
  const raw = naca4(code, teThickness, nPanels);
  return anchorAirfoil(raw, targetLe, targetTe);
}

// ---------------------------------------------------------------------------
// Flap-track kinematics (1-DOF: lead-axel arc-length `sLead`)
// ---------------------------------------------------------------------------

/** Point on the track at arc-length `s`: linear segment, then tangent arc. */
export function trackPoint(track: TrackConfig, s: number): V2 {
  const phi = (track.angleDeg * Math.PI) / 180.0;
  const t: V2 = [Math.cos(phi), Math.sin(phi)];
  if (s <= track.linearLength) return add(track.anchor, scale(t, s));
  const n: V2 = [-Math.sin(phi), Math.cos(phi)];
  const j = add(track.anchor, scale(t, track.linearLength));
  const c = add(j, scale(n, track.arcRadius));
  const ang0 = Math.atan2(j[1] - c[1], j[0] - c[0]);
  const ang = ang0 + (s - track.linearLength) / track.arcRadius;
  const r = Math.abs(track.arcRadius);
  return [c[0] + r * Math.cos(ang), c[1] + r * Math.sin(ang)];
}

/**
 * Arc-length of the trailing axel: the point on the track at chord distance `d`
 * behind the lead axel. Chord distance grows monotonically as we move back from
 * the lead, so a bisection on [0, sLead - d] converges to the single root.
 */
export function trailingAxel(track: TrackConfig, d: number, sLead: number): number {
  const lead = trackPoint(track, sLead);
  let lo = 0.0;
  let hi = sLead - d;
  for (let i = 0; i < 60; i++) {
    const m = 0.5 * (lo + hi);
    if (norm(sub(trackPoint(track, m), lead)) > d) lo = m;
    else hi = m;
  }
  return 0.5 * (lo + hi);
}

/**
 * Rigid transform carrying the stowed axel segment to the deployed one. Stowed
 * places the lead axel at s=d and the trailing axel at s=0, so deployment is the
 * identity at `sLead = d`.
 */
function trackMotion(track: TrackConfig, d: number, sLead: number): (p: V2) => V2 {
  const a0 = trackPoint(track, d);
  const b0 = trackPoint(track, 0);
  const a = trackPoint(track, sLead);
  const b = trackPoint(track, trailingAxel(track, d, sLead));
  const ang = Math.atan2(b[1] - a[1], b[0] - a[0]) - Math.atan2(b0[1] - a0[1], b0[0] - a0[0]);
  const c = Math.cos(ang);
  const s = Math.sin(ang);
  return ([x, y]) => {
    const dx = x - a0[0];
    const dy = y - a0[1];
    return [a[0] + dx * c - dy * s, a[1] + dx * s + dy * c];
  };
}

// ---------------------------------------------------------------------------
// Full configuration assembly
// ---------------------------------------------------------------------------

export interface Configuration {
  /** Main element contour (fixed; cove cut). */
  main: V2[];
  /** Vane contour at the current deployment. */
  vane: V2[];
  /** Flaperon contour at the current deployment. */
  flaperon: V2[];
  /** Flaperon hinge pivot, carried to the current deployment. */
  flaperonPivot: V2;
}

const chord = (cfg: NacaElementCfg): number => norm(sub(cfg.stowedTe, cfg.stowedLe));

/**
 * Build the element contours for the configuration's operating point. The lead-axel
 * arc-length is `sLead = design.axelSpacing + operation.deploy` (deploy = 0 is
 * stowed). The main element is fixed; the vane and flaperon deploy together as one
 * rigid assembly riding the track. The flaperon is first rotated about its hinge
 * (relative to the assembly), then carried by the track motion along with its pivot.
 * The blunt TE thickness (overall-chord units) is expressed in each element's own
 * chord units before truncation, matching the reference build.
 */
export function buildConfiguration(cfg: HighLiftAirfoil): Configuration {
  const { main: m, design: d, operation: op } = cfg;
  const sLead = d.axelSpacing + op.deploy;
  const main = buildMainWing(m.coords, d.cove, m.bluntThickness);
  const vaneFlat = buildNacaElement(d.vane.naca, m.bluntThickness / chord(d.vane), d.vane.stowedLe, d.vane.stowedTe);
  const flaperonFlat = buildNacaElement(d.flaperon.naca, m.bluntThickness / chord(d.flaperon), d.flaperon.stowedLe, d.flaperon.stowedTe);
  const move = trackMotion(d.track, d.axelSpacing, sLead);
  const pivot = d.flaperonHinge.pivot;
  const flaperonHinged = rotate(flaperonFlat, (op.flaperonAngleDeg * Math.PI) / 180.0, pivot);
  return {
    main,
    vane: vaneFlat.map(move),
    flaperon: flaperonHinged.map(move),
    flaperonPivot: move(pivot),
  };
}
