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
  FlaperonShape,
  HighLiftAirfoil,
  TrackConfig,
  V2,
  VaneShape,
} from './estolConfig';
import { evaluateBSpline } from '../lib/bspline';

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

/**
 * Coved main wing. The cove follows the airfoil's own surfaces offset inward by
 * `bluntThickness`: the ceiling parallels the upper surface from `upperCutX`
 * forward to `coveX`; the floor parallels the lower surface from `coveX` out to
 * `lowerCutX`; a vertical back wall at `coveX` joins them with sharp corners. Both
 * cuts carry a blunt-TE edge of `bluntThickness` (the lip dropped below the upper
 * surface, the lower cut raised above the lower surface). Returns the closed `flat`
 * contour: upper-lip-tip → LE → lower-surface → lower-cut → (blunt) → floor → back
 * wall → ceiling → lip-bottom, closing the upper blunt TE back to the lip tip.
 */
export function buildMainWing(airfoil: AirfoilCoords, cutouts: CoveCutouts, bluntThickness: number): V2[] {
  const upper = airfoil.upper;
  const lower = airfoil.lower;
  const { upperCutX: xUpCut, lowerCutX: xLoCut, coveX } = cutouts;

  const upX = upper.map((q) => q[0]);
  const upY = upper.map((q) => q[1]);
  const loX = lower.map((q) => q[0]);
  const loY = lower.map((q) => q[1]);
  const yUpAtCut = interp(xUpCut, upX, upY);
  const yLoAtCut = interp(xLoCut, loX, loY);
  const yUpCove = interp(coveX, upX, upY);
  const yLoCove = interp(coveX, loX, loY);

  const upKept: V2[] = [...upper.filter((q) => q[0] < xUpCut), [xUpCut, yUpAtCut]];
  const loKept: V2[] = [...lower.filter((q) => q[0] < xLoCut), [xLoCut, yLoAtCut]];
  const lipBottom: V2 = [xUpCut, yUpAtCut - bluntThickness];     // upper blunt-TE inner point
  const lowerCutTop: V2 = [xLoCut, yLoAtCut + bluntThickness];   // lower blunt-TE inner point

  // Cove floor: lower surface offset +bt, from the lower cut forward to coveX.
  const floorMid: V2[] = lower
    .filter((q) => q[0] > coveX && q[0] < xLoCut)
    .sort((a, b) => b[0] - a[0])                                 // decreasing x
    .map((q) => [q[0], q[1] + bluntThickness]);
  const floorCove: V2 = [coveX, yLoCove + bluntThickness];       // back-wall bottom corner (sharp)
  // Cove ceiling: upper surface offset -bt, from coveX aft to the upper cut.
  const ceilCove: V2 = [coveX, yUpCove - bluntThickness];        // back-wall top corner (sharp)
  const ceilMid: V2[] = upper
    .filter((q) => q[0] > coveX && q[0] < xUpCut)
    .sort((a, b) => a[0] - b[0])                                 // increasing x
    .map((q) => [q[0], q[1] - bluntThickness]);

  // main_contour: upper-lip-tip → LE → lower-cut-pt (the kept outer surfaces)
  const mainContour: V2[] = [...upKept.slice().reverse(), ...loKept.slice(1)];

  return [
    ...mainContour,        // lip tip → LE → lower cut point
    lowerCutTop,           // lower blunt-TE edge
    ...floorMid,           // floor: lowerCut → coveX (offset +bt)
    floorCove,             // back-wall bottom corner
    ceilCove,              // back-wall top corner
    ...ceilMid,            // ceiling: coveX → upperCut (offset -bt)
    lipBottom,             // upper blunt-TE inner point; closes to the lip tip
  ];
}

/** Vane: anchored NACA airfoil. Port of `build_naca_element()` (`flat`). */
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

/** Forward (decreasing-x) unit tangent of a tabulated surface at chordwise `xc`. */
function forwardTangent(xs: number[], ys: number[], xc: number): V2 {
  const h = 1e-3;
  const x0 = Math.max(xs[0], xc - h);
  const x1 = Math.min(xs[xs.length - 1], xc + h);
  const dx = x1 - x0;
  const dy = interp(x1, xs, ys) - interp(x0, xs, ys);
  return unit([-dx, -dy]);
}

/**
 * Flaperon contour, derived from the main airfoil (not an independent foil). Aft
 * of the cuts the upper/lower surfaces are exactly the airfoil, blunt-truncated at
 * the TE by `bluntThickness` (matching the other elements). Forward of the cuts the
 * nose is a clamped cubic B-spline through 5 control points: the two cut points,
 * two tangent-handle points along the airfoil slope at each cut (lengths
 * `noseTangent{Upper,Lower}`), and the free middle `noseTip` — so the nose meets the
 * airfoil with matching slope (C1) at each cut. Built in the airfoil (global) frame,
 * Selig-ordered: upper-TE-corner → upper-cut → nose → lower-cut → lower-TE-corner.
 */
/**
 * The five B-spline nose control points (stowed/airfoil frame): the two cut points
 * (P0, P4), the two slope-locked tangent handles (P1, P3), and the free middle
 * `noseTip` (P2). Exposed so the UI can draw them.
 */
export function flaperonControlPoints(airfoil: AirfoilCoords, shape: FlaperonShape): [V2, V2, V2, V2, V2] {
  const { upperCutX, lowerCutX, noseTangentUpper, noseTangentLower, noseTip } = shape;
  const upX = airfoil.upper.map((q) => q[0]);
  const upY = airfoil.upper.map((q) => q[1]);
  const loX = airfoil.lower.map((q) => q[0]);
  const loY = airfoil.lower.map((q) => q[1]);
  const P0: V2 = [upperCutX, interp(upperCutX, upX, upY)];
  const P4: V2 = [lowerCutX, interp(lowerCutX, loX, loY)];
  const P1: V2 = add(P0, scale(forwardTangent(upX, upY, upperCutX), noseTangentUpper));
  const P3: V2 = add(P4, scale(forwardTangent(loX, loY, lowerCutX), noseTangentLower));
  return [P0, P1, noseTip, P3, P4];
}

export function buildFlaperon(airfoil: AirfoilCoords, shape: FlaperonShape, bluntThickness: number): V2[] {
  const { upperCutX, lowerCutX } = shape;
  const upX = airfoil.upper.map((q) => q[0]);
  const upY = airfoil.upper.map((q) => q[1]);
  const loX = airfoil.lower.map((q) => q[0]);
  const loY = airfoil.lower.map((q) => q[1]);

  // Blunt-truncate the TE: gap = yUp - yLo decreases to 0 at the TE; find the x
  // (nearest the TE) where it equals bluntThickness. Reverse for ascending interp.
  const scan = linspace(0.6, 1.0, 2001);
  const gap = scan.map((x) => interp(x, upX, upY) - interp(x, loX, loY));
  const xTE = interp(bluntThickness, gap.slice().reverse(), scan.slice().reverse());

  const ctrl = flaperonControlPoints(airfoil, shape);
  const [P0, , , , P4] = ctrl;
  const cps = ctrl.map((p, i) => ({ x: p[0], y: p[1], id: String(i) }));
  const nose: V2[] = evaluateBSpline(cps, 3, 41).map((p) => [p.x, p.y]);

  // Airfoil surfaces aft of each cut, blunt-truncated at xTE.
  const upperAft: V2[] = [P0, ...airfoil.upper.filter((q) => q[0] > upperCutX && q[0] < xTE), [xTE, interp(xTE, upX, upY)]];
  const lowerAft: V2[] = [P4, ...airfoil.lower.filter((q) => q[0] > lowerCutX && q[0] < xTE), [xTE, interp(xTE, loX, loY)]];

  return [
    ...upperAft.slice().reverse(),   // upper TE corner → upper cut (P0)
    ...nose.slice(1, -1),            // nose interior (P0..P4 endpoints dropped — they're the cuts)
    ...lowerAft,                     // lower cut (P4) → lower TE corner
  ];
}

/** The vane's 5 logical control points: [te, ...points] (the TE + 4 others). */
export function vaneControlPoints(shape: VaneShape): V2[] {
  return [shape.te, ...shape.points];
}

/**
 * Vane contour: an open clamped cubic B-spline. The TE location is split into the
 * spline's start and end points by the blunt-TE thickness, added orthogonal to the
 * averaged direction of the TE's two adjacent control points (`points[0]` and
 * `points[3]`); the 4 control points shape the loop in between. The spline passes
 * through the two TE corners (gap = `bluntThickness`) and approximates the 4 points.
 * Selig-ordered: TE-corner → P1 → … → P4 → other-TE-corner.
 */
export function buildVane(shape: VaneShape, bluntThickness: number): V2[] {
  const { te, points } = shape;
  const [p1, p2, p3, p4] = points;
  // averaged direction toward the two control points adjacent to the TE
  const dir = unit(add(unit(sub(p1, te)), unit(sub(p4, te))));
  let n: V2 = [-dir[1], dir[0]];                       // orthogonal to that direction
  if (dot(sub(p1, te), n) < 0) n = [-n[0], -n[1]];     // orient toward p1's side
  const h = bluntThickness / 2;
  const teStart: V2 = [te[0] + h * n[0], te[1] + h * n[1]];   // start (adjacent to p1)
  const teEnd: V2 = [te[0] - h * n[0], te[1] - h * n[1]];     // end (adjacent to p4)
  const cps = [teStart, p1, p2, p3, p4, teEnd].map((p, i) => ({ x: p[0], y: p[1], id: String(i) }));
  return evaluateBSpline(cps, 3, 81).map((q) => [q.x, q.y]);
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
  /** The 5 vane B-spline control points, carried to the current deployment. */
  vaneControls: V2[];
  /** Flaperon contour at the current deployment. */
  flaperon: V2[];
  /** Flaperon hinge pivot, carried to the current deployment. */
  flaperonPivot: V2;
  /** The 5 flaperon nose control points, carried to the current deployment. */
  flaperonControls: V2[];
}

/**
 * Build the element contours for the configuration's operating point. The lead-axel
 * arc-length is `sLead = design.axelSpacing + operation.deploy` (deploy = 0 is
 * stowed). The main element is fixed; the vane and flaperon deploy together as one
 * rigid assembly riding the track. The flaperon is additionally rotated about its
 * hinge (relative to the assembly), then carried by the track motion along with its
 * pivot. All elements share the same blunt TE thickness (`main.bluntThickness`).
 */
export function buildConfiguration(cfg: HighLiftAirfoil): Configuration {
  const { main: m, design: d, operation: op } = cfg;
  const sLead = d.axelSpacing + op.deploy;
  const main = buildMainWing(m.coords, d.cove, m.bluntThickness);
  const vaneFlat = buildVane(d.vane, m.bluntThickness);
  const flaperonFlat = buildFlaperon(m.coords, d.flaperon, m.bluntThickness);
  const move = trackMotion(d.track, d.axelSpacing, sLead);
  const pivot = d.flaperonHinge.pivot;
  const angleRad = (op.flaperonAngleDeg * Math.PI) / 180.0;
  const deploy = (pts: V2[]): V2[] => rotate(pts, angleRad, pivot).map(move);
  return {
    main,
    vane: vaneFlat.map(move),
    vaneControls: vaneControlPoints(d.vane).map(move),
    flaperon: deploy(flaperonFlat),
    flaperonPivot: move(pivot),
    flaperonControls: deploy(flaperonControlPoints(m.coords, d.flaperon)),
  };
}
