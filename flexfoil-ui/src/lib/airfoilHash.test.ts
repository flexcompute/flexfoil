import { describe, it, expect } from 'vitest';
import {
  canonicalGeometryString,
  computeAirfoilHash,
  computeGeometryHash,
  isIdentityPlacement,
  type GeometryElement,
  type Placement,
} from './airfoilHash';

// Shared cross-language fixture. The identical fixture and the identical
// expected canonical strings appear in
// packages/flexfoil-python/tests/test_geometry_hash.py — if either side's
// canonicalisation drifts, one of the two pinned strings stops matching.
const MAIN = [
  { x: 0.0, y: 0.0 },
  { x: 1.0, y: 0.0 },
  { x: 0.5, y: 0.06 },
  { x: 0.25, y: -0.03 },
];

const FLAP = [
  { x: 0.8, y: 0.0 },
  { x: 1.2, y: 0.0 },
  { x: 1.0, y: 0.02 },
];

const MAIN_COORDS =
  '0.00000000,0.00000000;1.00000000,0.00000000;0.50000000,0.06000000;0.25000000,-0.03000000';
const FLAP_COORDS = '0.80000000,0.00000000;1.20000000,0.00000000;1.00000000,0.02000000';

// Digest of MAIN alone, as produced by the coordinate-only implementation that
// shipped before placement was folded in. Existing cached runs are keyed on
// this; it must never change.
const MAIN_LEGACY_DIGEST =
  'ef89b0ef74a644048ea6c5572e54ef740084c204088a0e3807200db3751990e2';

const PLACEMENT: Placement = {
  rotation: 25,
  pivot: { x: 0.8, y: 0.01 },
  translation: { x: 0.02, y: -0.05 },
  scale: 1,
};

// The placement/element encoding is the part that MUST agree byte-for-byte
// with Python. This exact literal is asserted on both sides.
const PLACEMENT_SUFFIX = '@25.00000000,0.80000000,0.01000000,0.02000000,-0.05000000,1.00000000';
const ELEMENT_SEP = '#';

describe('computeAirfoilHash — backward compatibility', () => {
  it('pins the legacy digest for a single element with no placement', async () => {
    expect(await computeAirfoilHash(MAIN)).toBe(MAIN_LEGACY_DIGEST);
  });

  it('canonicalises a single unplaced element to coordinates alone', () => {
    expect(canonicalGeometryString([{ panels: MAIN }])).toBe(MAIN_COORDS);
  });

  it('treats an identity placement as absent', async () => {
    const identity: Placement = {
      rotation: 0,
      pivot: { x: 0, y: 0 },
      translation: { x: 0, y: 0 },
      scale: 1,
    };
    expect(await computeAirfoilHash(MAIN, identity)).toBe(MAIN_LEGACY_DIGEST);
    expect(await computeAirfoilHash(MAIN, {})).toBe(MAIN_LEGACY_DIGEST);
    expect(await computeAirfoilHash(MAIN, null)).toBe(MAIN_LEGACY_DIGEST);
  });

  it('treats a pivot-only placement as identity (a pivot alone moves nothing)', async () => {
    expect(await computeAirfoilHash(MAIN, { pivot: { x: 0.75, y: 0.02 } })).toBe(
      MAIN_LEGACY_DIGEST
    );
    expect(isIdentityPlacement({ pivot: { x: 0.75, y: 0.02 } })).toBe(true);
  });

  it('treats a single-element list as equivalent to a bare panel list', async () => {
    expect(await computeGeometryHash([{ panels: MAIN }])).toBe(MAIN_LEGACY_DIGEST);
  });

  it('collapses -0 so that the digest matches Python formatting', async () => {
    const negZero = await computeAirfoilHash(MAIN, { translation: { x: -0, y: 0.1 } });
    const posZero = await computeAirfoilHash(MAIN, { translation: { x: 0, y: 0.1 } });
    expect(negZero).toBe(posZero);
  });
});

describe('computeAirfoilHash — placement sensitivity', () => {
  it('distinguishes configurations differing only in translation', async () => {
    const a = await computeAirfoilHash(MAIN, { translation: { x: 0.01, y: 0 } });
    const b = await computeAirfoilHash(MAIN, { translation: { x: 0.02, y: 0 } });
    expect(a).not.toBe(b);
    expect(a).not.toBe(MAIN_LEGACY_DIGEST);
  });

  it('distinguishes configurations differing only in rotation', async () => {
    const a = await computeAirfoilHash(MAIN, { rotation: 10, pivot: { x: 0.75, y: 0 } });
    const b = await computeAirfoilHash(MAIN, { rotation: 20, pivot: { x: 0.75, y: 0 } });
    expect(a).not.toBe(b);
  });

  it('distinguishes configurations differing only in pivot when rotated', async () => {
    const a = await computeAirfoilHash(MAIN, { rotation: 10, pivot: { x: 0.7, y: 0 } });
    const b = await computeAirfoilHash(MAIN, { rotation: 10, pivot: { x: 0.75, y: 0 } });
    expect(a).not.toBe(b);
  });

  it('distinguishes configurations differing only in scale', async () => {
    const a = await computeAirfoilHash(MAIN, { scale: 1.0000001 });
    const b = await computeAirfoilHash(MAIN, { scale: 0.3 });
    expect(a).not.toBe(b);
  });

  it('resolves placement differences below the 8-decimal precision to the same key', async () => {
    // Documented precision limit: differences finer than 1e-8 are the same key.
    const a = await computeAirfoilHash(MAIN, { translation: { x: 0.1, y: 0 } });
    const b = await computeAirfoilHash(MAIN, { translation: { x: 0.1 + 1e-12, y: 0 } });
    expect(a).toBe(b);
  });
});

describe('computeGeometryHash — element count and order', () => {
  const main: GeometryElement = { panels: MAIN };
  const flap: GeometryElement = { panels: FLAP, placement: PLACEMENT };

  it('encodes the placement suffix exactly as Python does', () => {
    expect(canonicalGeometryString([main, flap])).toBe(
      MAIN_COORDS + ELEMENT_SEP + FLAP_COORDS + PLACEMENT_SUFFIX
    );
  });

  it('distinguishes element count', async () => {
    const one = await computeGeometryHash([main]);
    const two = await computeGeometryHash([main, main]);
    expect(one).not.toBe(two);
  });

  it('distinguishes element order', async () => {
    const ab = await computeGeometryHash([main, { panels: FLAP }]);
    const ba = await computeGeometryHash([{ panels: FLAP }, main]);
    expect(ab).not.toBe(ba);
  });

  it('distinguishes assemblies that share shapes but differ in placement only', async () => {
    const closed = await computeGeometryHash([
      main,
      { panels: FLAP, placement: { translation: { x: 0, y: -0.01 } } },
    ]);
    const open = await computeGeometryHash([
      main,
      { panels: FLAP, placement: { translation: { x: 0, y: -0.03 } } },
    ]);
    expect(closed).not.toBe(open);
  });

  it('is stable across repeated calls', async () => {
    expect(await computeGeometryHash([main, flap])).toBe(
      await computeGeometryHash([main, flap])
    );
  });

  it('hashes an empty assembly without throwing', async () => {
    expect(await computeGeometryHash([])).toHaveLength(64);
  });
});
