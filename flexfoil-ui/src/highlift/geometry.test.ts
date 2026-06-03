/**
 * Behavioral tests for the high-lift geometry engine. The elements are now
 * parameterized (contour-following cove, airfoil-derived flaperon, free-form
 * B-spline vane), so they no longer match the old reference .dat fixtures; the
 * NACA primitives (naca4 / anchorAirfoil) are still exercised directly.
 */

import { describe, expect, it } from 'vitest';

import {
  anchorAirfoil,
  buildConfiguration,
  buildFlaperon,
  flaperonControlPoints,
  buildMainWing,
  buildVane,
  vaneControlPoints,
  naca4,
  trackPoint,
  trailingAxel,
} from './geometry';
import { DEFAULT_HIGH_LIFT_AIRFOIL as CFG } from './estolConfig';
import type { HighLiftAirfoil, V2 } from './estolConfig';

const norm = (a: V2, b: V2) => Math.hypot(a[0] - b[0], a[1] - b[1]);

const AXEL = CFG.design.axelSpacing;

/**
 * A config clone at a given lead-axel arc-length `sLead` (so `sLead = AXEL` is
 * stowed), with optional flaperon angle / hinge pivot overrides. Lets the
 * kinematic tests drive deployment the way the old `buildConfiguration(cfg, sLead)`
 * signature did, now that deploy lives in `operation`.
 */
function at(sLead: number, over: { flaperonAngleDeg?: number; pivot?: V2 } = {}): HighLiftAirfoil {
  return {
    ...CFG,
    design: { ...CFG.design, flaperonHinge: { pivot: over.pivot ?? CFG.design.flaperonHinge.pivot } },
    operation: {
      deploy: sLead - AXEL,
      flaperonAngleDeg: over.flaperonAngleDeg ?? CFG.operation.flaperonAngleDeg,
    },
  };
}

describe('geometry primitives', () => {
  it('naca4 produces 2*(n/2+1)-1 points with a blunt TE gap', () => {
    const pts = naca4('0012', 0.01, 160);
    expect(pts.length).toBe(161);
    // First and last points are the two TE corners; their gap is the blunt TE.
    const teGap = Math.hypot(pts[0][0] - pts[pts.length - 1][0], pts[0][1] - pts[pts.length - 1][1]);
    expect(teGap).toBeCloseTo(0.01, 3);
  });

  it('anchorAirfoil places LE and TE at the requested anchors', () => {
    const raw = naca4('2412', 0.005, 80);
    const le: V2 = [0.55, -0.048];
    const te: V2 = [0.73, 0.035];
    const a = anchorAirfoil(raw, le, te);
    // LE = min-x point of the raw foil maps to the target LE.
    let leIdx = 0;
    for (let i = 1; i < raw.length; i++) if (raw[i][0] < raw[leIdx][0]) leIdx = i;
    expect(a[leIdx][0]).toBeCloseTo(le[0], 9);
    expect(a[leIdx][1]).toBeCloseTo(le[1], 9);
    const teMid: V2 = [0.5 * (a[0][0] + a[a.length - 1][0]), 0.5 * (a[0][1] + a[a.length - 1][1])];
    expect(teMid[0]).toBeCloseTo(te[0], 9);
    expect(teMid[1]).toBeCloseTo(te[1], 9);
  });
});

describe('main cove cutout (contour-following)', () => {
  const cove = CFG.design.cove;
  const bt = CFG.main.bluntThickness;
  const up = CFG.main.coords.upper;
  const lo = CFG.main.coords.lower;
  const c = buildMainWing(CFG.main.coords, cove, bt);
  const near = (a: number, b: number, tol = 1e-9) => Math.abs(a - b) < tol;
  const has = (p: V2, tol = 1e-9) => c.some((q) => near(q[0], p[0], tol) && near(q[1], p[1], tol));

  it('the ceiling parallels the upper surface, offset down by bt (coveX..upperCutX)', () => {
    const samples = up.filter((q) => q[0] > cove.coveX && q[0] < cove.upperCutX);
    expect(samples.length).toBeGreaterThan(0);
    for (const q of samples) expect(has([q[0], q[1] - bt])).toBe(true);
  });

  it('the floor parallels the lower surface, offset up by bt (coveX..lowerCutX)', () => {
    const samples = lo.filter((q) => q[0] > cove.coveX && q[0] < cove.lowerCutX);
    expect(samples.length).toBeGreaterThan(0);
    for (const q of samples) expect(has([q[0], q[1] + bt])).toBe(true);
  });

  it('keeps the airfoil outer surfaces forward of the cuts', () => {
    const u = up.find((q) => q[0] > 0.1 && q[0] < cove.upperCutX)!;
    const l = lo.find((q) => q[0] > 0.1 && q[0] < cove.lowerCutX)!;
    expect(has(u)).toBe(true);
    expect(has(l)).toBe(true);
  });

  it('has a blunt TE of bt at both cuts', () => {
    for (const xc of [cove.upperCutX, cove.lowerCutX]) {
      const ys = c.filter((p) => near(p[0], xc)).map((p) => p[1]).sort((a, b) => a - b);
      expect(ys.some((y, i) => i > 0 && Math.abs(y - ys[i - 1] - bt) < 1e-9)).toBe(true);
    }
  });

  it('places the back wall forward of both cuts', () => {
    expect(cove.coveX).toBeLessThan(cove.lowerCutX);
    expect(cove.coveX).toBeLessThan(cove.upperCutX);
  });
});

describe('flaperon (airfoil-derived, B-spline nose)', () => {
  const shape = CFG.design.flaperon;
  const bt = CFG.main.bluntThickness;
  const up = CFG.main.coords.upper;
  const lo = CFG.main.coords.lower;
  const c = buildFlaperon(CFG.main.coords, shape, bt);
  const near = (a: number, b: number, tol = 1e-9) => Math.abs(a - b) < tol;
  const has = (p: V2, tol = 1e-9) => c.some((q) => near(q[0], p[0], tol) && near(q[1], p[1], tol));
  // local linear interp matching geometry.ts's interp (ascending xp)
  const lerp = (x: number, pts: V2[]) => {
    let i = 1;
    while (i < pts.length && pts[i][0] < x) i++;
    const [x0, y0] = pts[i - 1];
    const [x1, y1] = pts[i];
    return y0 + ((y1 - y0) * (x - x0)) / (x1 - x0);
  };

  it('the surfaces aft of the cuts are exactly the airfoil', () => {
    const u = up.filter((q) => q[0] > shape.upperCutX && q[0] < 0.96);
    const l = lo.filter((q) => q[0] > shape.lowerCutX && q[0] < 0.96);
    expect(u.length + l.length).toBeGreaterThan(2);
    for (const q of [...u, ...l]) expect(has(q)).toBe(true);
  });

  it('the nose passes through both cut points (on the airfoil)', () => {
    expect(has([shape.upperCutX, lerp(shape.upperCutX, up)], 1e-6)).toBe(true);
    expect(has([shape.lowerCutX, lerp(shape.lowerCutX, lo)], 1e-6)).toBe(true);
  });

  it('exposes 5 nose control points (P0/P4 on the cuts, P2 = noseTip)', () => {
    const cp = flaperonControlPoints(CFG.main.coords, shape);
    expect(cp.length).toBe(5);
    expect(cp[0][0]).toBeCloseTo(shape.upperCutX, 9);
    expect(cp[0][1]).toBeCloseTo(lerp(shape.upperCutX, up), 9);
    expect(cp[4][0]).toBeCloseTo(shape.lowerCutX, 9);
    expect(cp[4][1]).toBeCloseTo(lerp(shape.lowerCutX, lo), 9);
    expect(cp[2]).toEqual(shape.noseTip);
  });

  it('the nose extends forward of both cuts', () => {
    const minX = Math.min(...c.map((p) => p[0]));
    expect(minX).toBeLessThan(Math.min(shape.upperCutX, shape.lowerCutX));
  });

  it('has a blunt TE of bt at the trailing edge', () => {
    const maxX = Math.max(...c.map((p) => p[0]));
    const ys = c.filter((p) => near(p[0], maxX)).map((p) => p[1]).sort((a, b) => a - b);
    expect(ys[ys.length - 1] - ys[0]).toBeCloseTo(bt, 6);
  });

  it('the nose leaves the upper cut along the airfoil tangent (C1)', () => {
    const i0 = c.findIndex((p) => near(p[0], shape.upperCutX, 1e-6) && near(p[1], lerp(shape.upperCutX, up), 1e-6));
    expect(i0).toBeGreaterThanOrEqual(0);
    const u2 = (v: V2): V2 => { const m = Math.hypot(v[0], v[1]); return [v[0] / m, v[1] / m]; };
    const cos = (a: V2, b: V2) => a[0] * b[0] + a[1] * b[1];
    const dNose = u2([c[i0 + 1][0] - c[i0][0], c[i0 + 1][1] - c[i0][1]]);    // into the nose
    const h = 1e-3;                                                          // airfoil forward tangent at the cut
    const dAir = u2([-2 * h, lerp(shape.upperCutX - h, up) - lerp(shape.upperCutX + h, up)]);
    const dTip = u2([shape.noseTip[0] - c[i0][0], shape.noseTip[1] - c[i0][1]]); // straight-to-tip (non-tangent case)
    // The nose direction at the join aligns with the airfoil tangent (cos≈1; the
    // small residual is finite-sample secant error), and far more than with the
    // straight-to-tip direction it would take without the tangency constraint.
    expect(cos(dNose, dAir)).toBeGreaterThan(0.99);
    expect(cos(dNose, dAir)).toBeGreaterThan(cos(dNose, dTip));
  });
});

describe('vane (free-form closed B-spline + blunt TE)', () => {
  const bt = CFG.main.bluntThickness;
  const shape = CFG.design.vane;
  const c = buildVane(shape, bt);

  it('exposes 5 control points with [0] = te', () => {
    const cp = vaneControlPoints(shape);
    expect(cp.length).toBe(5);
    expect(cp[0]).toEqual(shape.te);
  });

  it('is an open contour with a blunt TE gap of bt', () => {
    expect(c.length).toBeGreaterThan(20);
    expect(norm(c[0], c[c.length - 1])).toBeCloseTo(bt, 6);   // the two TE corners, exactly bt apart
  });

  it('stays within the control-polygon bounding box (B-spline approximates inward)', () => {
    const cp = vaneControlPoints(shape);
    const pad = bt;   // the TE corners sit ±bt/2 off the te control point
    const xmin = Math.min(...cp.map((p) => p[0])) - pad, xmax = Math.max(...cp.map((p) => p[0])) + pad;
    const ymin = Math.min(...cp.map((p) => p[1])) - pad, ymax = Math.max(...cp.map((p) => p[1])) + pad;
    for (const p of c) {
      expect(p[0]).toBeGreaterThanOrEqual(xmin);
      expect(p[0]).toBeLessThanOrEqual(xmax);
      expect(p[1]).toBeGreaterThanOrEqual(ymin);
      expect(p[1]).toBeLessThanOrEqual(ymax);
    }
  });
});

describe('flap-track kinematics (1-DOF)', () => {
  // A well-formed local track (linear segment longer than the axel bar) so the
  // mechanism invariants are exactly testable. The production default deliberately
  // keeps axelSpacing > linearLength (the assembly is always partly on the arc, so
  // it never purely translates) — a valid design, just not what these tests probe.
  const KAXEL = 0.1;
  const KTRACK = { anchor: [0.6, -0.05] as V2, angleDeg: -6, linearLength: 0.35, arcRadius: -0.135, length: 0.5 };
  const KCFG: HighLiftAirfoil = { ...CFG, design: { ...CFG.design, axelSpacing: KAXEL, track: KTRACK } };
  const kat = (sLead: number): HighLiftAirfoil => ({
    ...KCFG,
    operation: { deploy: sLead - KAXEL, flaperonAngleDeg: 0 },
  });
  const vaneFlat = buildVane(CFG.design.vane, CFG.main.bluntThickness);
  const leIdx = (() => {
    let k = 0;
    for (let i = 1; i < vaneFlat.length; i++) if (vaneFlat[i][0] < vaneFlat[k][0]) k = i;
    return k;
  })();
  const angle = (pts: V2[]) => Math.atan2(pts[0][1] - pts[leIdx][1], pts[0][0] - pts[leIdx][0]);

  it('deployment is the identity when stowed (deploy = 0)', () => {
    const { vane } = buildConfiguration(kat(KAXEL));
    for (let i = 0; i < vaneFlat.length; i++) {
      expect(vane[i][0]).toBeCloseTo(vaneFlat[i][0], 9);
      expect(vane[i][1]).toBeCloseTo(vaneFlat[i][1], 9);
    }
  });

  it('the two axels stay exactly axelSpacing apart through deployment', () => {
    for (const sLead of [KAXEL, 0.2, 0.4, 0.6]) {
      const lead = trackPoint(KTRACK, sLead);
      const trail = trackPoint(KTRACK, trailingAxel(KTRACK, KAXEL, sLead));
      expect(norm(lead, trail)).toBeCloseTo(KAXEL, 9);
    }
  });

  it('pure translation while both axels are on the linear segment', () => {
    const sLead = KTRACK.linearLength - 1e-3; // both axels still on the line ⇒ no rotation
    const { vane } = buildConfiguration(kat(sLead));
    const t0: V2 = [vane[0][0] - vaneFlat[0][0], vane[0][1] - vaneFlat[0][1]];
    for (let i = 1; i < vaneFlat.length; i++) {
      expect(vane[i][0] - vaneFlat[i][0]).toBeCloseTo(t0[0], 9);
      expect(vane[i][1] - vaneFlat[i][1]).toBeCloseTo(t0[1], 9);
    }
  });

  it('develops TE-down rotation once on the arc (negative radius)', () => {
    const { vane } = buildConfiguration(kat(0.6));
    // Chord vector (LE→TE) rotates clockwise (TE down) relative to stowed.
    expect(angle(vane) - angle(vaneFlat)).toBeLessThan(-0.3);
  });
});

describe('flaperon hinge (relative to the assembly)', () => {
  // Orientation from two fixed material points (TE corner → LE), so the measure
  // tracks the same vertices regardless of rotation. LE is the contour midpoint.
  const mid = Math.floor(buildConfiguration(at(AXEL)).flaperon.length / 2);
  const flaperonAngle = (pts: V2[]) => Math.atan2(pts[0][1] - pts[mid][1], pts[0][0] - pts[mid][0]);

  it('rotates the flaperon by exactly the relative angle, at any deployment', () => {
    for (const sLead of [AXEL, 0.2, 0.35]) {
      const f0 = buildConfiguration(at(sLead, { flaperonAngleDeg: 0 })).flaperon;
      const fA = buildConfiguration(at(sLead, { flaperonAngleDeg: -15 })).flaperon;
      expect(((flaperonAngle(fA) - flaperonAngle(f0)) * 180) / Math.PI).toBeCloseTo(-15, 6);
    }
  });

  it('leaves the vane unaffected by the flaperon angle', () => {
    const v0 = buildConfiguration(at(0.25, { flaperonAngleDeg: 0 })).vane;
    const vA = buildConfiguration(at(0.25, { flaperonAngleDeg: -20 })).vane;
    for (let i = 0; i < v0.length; i++) {
      expect(vA[i][0]).toBeCloseTo(v0[i][0], 12);
      expect(vA[i][1]).toBeCloseTo(v0[i][1], 12);
    }
  });

  it('reports the pivot at its stowed location when stowed, carried otherwise', () => {
    const pivot: V2 = [0.72, -0.02];
    const stowed = buildConfiguration(at(AXEL, { pivot, flaperonAngleDeg: -25 })).flaperonPivot;
    expect(stowed[0]).toBeCloseTo(0.72, 9);
    expect(stowed[1]).toBeCloseTo(-0.02, 9);
    // deployed: the pivot is carried off its stowed spot by the track motion.
    const deployed = buildConfiguration(at(0.3, { pivot, flaperonAngleDeg: -25 })).flaperonPivot;
    expect(norm(deployed, stowed)).toBeGreaterThan(0.05);
  });
});
