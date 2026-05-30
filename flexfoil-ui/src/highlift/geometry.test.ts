/**
 * Validates the TypeScript geometry engine against the reference Selig .dat
 * files produced by ~/himalaya/geometry/airfoils/build_estol_geometry.py.
 *
 * The fixtures in ./__fixtures__/ are verbatim copies of main_wing.dat,
 * vane_stowed.dat and flap_stowed.dat. They are written to 6 decimal places,
 * so a max rounding error of 5e-7 is expected; we assert agreement to 1e-6.
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

import {
  anchorAirfoil,
  buildConfiguration,
  buildMainWing,
  buildNacaElement,
  naca4,
  trackPoint,
  trailingAxel,
} from './geometry';
import { DEFAULT_ESTOL_CONFIG as CFG } from './estolConfig';
import type { V2 } from './estolConfig';

const norm = (a: V2, b: V2) => Math.hypot(a[0] - b[0], a[1] - b[1]);

const TOL = 1e-6;

function loadFixture(name: string): V2[] {
  const path = fileURLToPath(new URL(`./__fixtures__/${name}`, import.meta.url));
  const lines = readFileSync(path, 'utf8').split('\n');
  const pts: V2[] = [];
  for (const line of lines.slice(1)) {
    const trimmed = line.trim();
    if (!trimmed) continue;
    const [x, y] = trimmed.split(/\s+/).map(Number);
    pts.push([x, y]);
  }
  return pts;
}

function expectContourMatches(got: V2[], expected: V2[], tol = TOL): void {
  expect(got.length).toBe(expected.length);
  let maxErr = 0;
  for (let i = 0; i < expected.length; i++) {
    maxErr = Math.max(maxErr, Math.abs(got[i][0] - expected[i][0]), Math.abs(got[i][1] - expected[i][1]));
  }
  // Surfaced in the failure message to make a near-miss easy to diagnose.
  expect(maxErr, `max coordinate error ${maxErr.toExponential(3)}`).toBeLessThan(tol);
}

const chord = (le: V2, te: V2) => Math.hypot(te[0] - le[0], te[1] - le[1]);

describe('high-lift geometry engine vs reference build_estol_geometry.py', () => {
  it('coved main wing matches main_wing.dat', () => {
    const got = buildMainWing(CFG.mainAirfoil, CFG.mainCutouts, CFG.bluntThickness);
    expectContourMatches(got, loadFixture('main_wing.dat'));
  });

  it('NACA 9621 vane matches vane_stowed.dat', () => {
    const te = CFG.bluntThickness / chord(CFG.vane.stowedLe, CFG.vane.stowedTe);
    const got = buildNacaElement(CFG.vane.naca, te, CFG.vane.stowedLe, CFG.vane.stowedTe);
    expectContourMatches(got, loadFixture('vane_stowed.dat'));
  });

  it('NACA 6311 aft flap matches flap_stowed.dat', () => {
    const te = CFG.bluntThickness / chord(CFG.aftFlap.stowedLe, CFG.aftFlap.stowedTe);
    const got = buildNacaElement(CFG.aftFlap.naca, te, CFG.aftFlap.stowedLe, CFG.aftFlap.stowedTe);
    expectContourMatches(got, loadFixture('flap_stowed.dat'));
  });
});

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

describe('flap-track kinematics (1-DOF)', () => {
  const { track, axelSpacing: d } = CFG;
  // Same TE thickness buildConfiguration derives, so the reference foil matches.
  const vaneTe = CFG.bluntThickness / norm(CFG.vane.stowedTe, CFG.vane.stowedLe);
  const vaneFlat = buildNacaElement(CFG.vane.naca, vaneTe, CFG.vane.stowedLe, CFG.vane.stowedTe);
  const leIdx = (() => {
    let k = 0;
    for (let i = 1; i < vaneFlat.length; i++) if (vaneFlat[i][0] < vaneFlat[k][0]) k = i;
    return k;
  })();
  const angle = (pts: V2[]) => Math.atan2(pts[0][1] - pts[leIdx][1], pts[0][0] - pts[leIdx][0]);

  it('deployment is the identity when stowed (sLead = axelSpacing)', () => {
    const { vane } = buildConfiguration(CFG, d);
    for (let i = 0; i < vaneFlat.length; i++) {
      expect(vane[i][0]).toBeCloseTo(vaneFlat[i][0], 9);
      expect(vane[i][1]).toBeCloseTo(vaneFlat[i][1], 9);
    }
  });

  it('the two axels stay exactly axelSpacing apart through deployment', () => {
    for (const sLead of [d, 0.2, 0.35, 0.5]) {
      const lead = trackPoint(track, sLead);
      const trail = trackPoint(track, trailingAxel(track, d, sLead));
      expect(norm(lead, trail)).toBeCloseTo(d, 9);
    }
  });

  it('pure translation while both axels are on the linear segment', () => {
    const sLead = track.linearLength - 1e-3; // both axels still on the line ⇒ no rotation
    const { vane } = buildConfiguration(CFG, sLead);
    const t0: V2 = [vane[0][0] - vaneFlat[0][0], vane[0][1] - vaneFlat[0][1]];
    for (let i = 1; i < vaneFlat.length; i++) {
      expect(vane[i][0] - vaneFlat[i][0]).toBeCloseTo(t0[0], 9);
      expect(vane[i][1] - vaneFlat[i][1]).toBeCloseTo(t0[1], 9);
    }
  });

  it('develops TE-down rotation once on the arc (negative radius)', () => {
    const { vane } = buildConfiguration(CFG, 0.5);
    // Chord vector (LE→TE) rotates clockwise (TE down) relative to stowed.
    expect(angle(vane) - angle(vaneFlat)).toBeLessThan(-0.3);
  });
});
