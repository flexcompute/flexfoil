import { describe, it, expect } from 'vitest';
import {
  canonicalGeometryString,
  computeAirfoilHash,
  computeGeometryHash,
  isIdentityPlacement,
  isLegacySingleElement,
  sharedCanonicalGeometryString,
  type GeometryElement,
  type Placement,
} from './airfoilHash';

// ---------------------------------------------------------------------------
// Shared cross-language fixtures.
//
// The identical fixtures, canonical strings and digest literals appear in
// packages/flexfoil-python/tests/test_geometry_hash.py. The shared digests are
// the deliverable: they are the *same hex strings* on both sides, so if either
// implementation's canonicalisation drifts, one of the two suites fails.
// ---------------------------------------------------------------------------

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

// On this side the frozen legacy node separator and the shared node separator
// are both ";", so one string serves both encodings. Python's differ ("|" vs
// ";") and its test file therefore pins two.
const MAIN_COORDS =
  '0.00000000,0.00000000;1.00000000,0.00000000;0.50000000,0.06000000;0.25000000,-0.03000000';
const FLAP_COORDS = '0.80000000,0.00000000;1.20000000,0.00000000;1.00000000,0.02000000';

// Digest of MAIN alone, as produced by the coordinate-only implementation that
// shipped before placement was folded in. Existing cached runs are keyed on
// this; it must never change. Not comparable with Python's legacy digest — see
// the module notes in airfoilHash.ts.
const MAIN_LEGACY_DIGEST =
  'ef89b0ef74a644048ea6c5572e54ef740084c204088a0e3807200db3751990e2';

const PLACEMENT: Placement = {
  rotation: 25,
  pivot: { x: 0.8, y: 0.01 },
  translation: { x: 0.02, y: -0.05 },
  scale: 1,
};

const SHARED_SCHEMA = 'ffgeom1:';
const PLACEMENT_SUFFIX = '@25.00000000,0.80000000,0.01000000,0.02000000,-0.05000000,1.00000000';
const ELEMENT_SEP = '#';

// Vector 1 — a two-element assembly with a non-identity placement on the
// second element. Pinned byte-for-byte against Python.
const SHARED_ASSEMBLY_CANONICAL =
  SHARED_SCHEMA + MAIN_COORDS + ELEMENT_SEP + FLAP_COORDS + PLACEMENT_SUFFIX;
const SHARED_ASSEMBLY_DIGEST =
  'cbec85f2a62ddd0613c392aec6b3ecc751665d6096a2206e67bed2e2b1c5e53f';

// Vector 2 — float-formatting edge cases, which is where the two languages
// diverge if the formatter is not shared: negative zero in a coordinate and in
// a placement field, and 0.001953125 (an odd multiple of 1/512, the smallest
// family of doubles that lands exactly on a half at the eighth decimal, where
// Python's round-half-to-even and JavaScript's round-half-away-from-zero
// disagree). Both must render 0.00195313.
const EDGE = [
  { x: -0, y: 0 },
  { x: 0.001953125, y: -0.001953125 },
  { x: 0.5, y: -0 },
];
const EDGE_PLACEMENT: Placement = {
  rotation: -0,
  pivot: { x: 0, y: 0 },
  translation: { x: 0.001953125, y: -0 },
  scale: 1.001953125,
};
const SHARED_EDGE_CANONICAL =
  SHARED_SCHEMA +
  '0.00000000,0.00000000;0.00195313,-0.00195313;0.50000000,0.00000000' +
  '@0.00000000,0.00000000,0.00000000,0.00195313,0.00000000,1.00195313';
const SHARED_EDGE_DIGEST =
  '3b3c8bf81354ec30ee2fe3f151771c7ed2791eff9a70d76cb2080b5807ed9c98';

describe('computeAirfoilHash — backward compatibility', () => {
  it('pins the legacy digest for a single element with no placement', async () => {
    expect(await computeAirfoilHash(MAIN)).toBe(MAIN_LEGACY_DIGEST);
  });

  it('canonicalises a single unplaced element to coordinates alone', () => {
    expect(canonicalGeometryString([{ panels: MAIN, placement: null }])).toBe(MAIN_COORDS);
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
    expect(await computeGeometryHash([{ panels: MAIN, placement: null }])).toBe(
      MAIN_LEGACY_DIGEST
    );
  });

  it('routes only the lone unplaced element to the legacy encoding', () => {
    expect(isLegacySingleElement([{ panels: MAIN, placement: null }])).toBe(true);
    expect(isLegacySingleElement([{ panels: MAIN, placement: {} }])).toBe(true);
    expect(isLegacySingleElement([{ panels: MAIN, placement: PLACEMENT }])).toBe(false);
    expect(
      isLegacySingleElement([
        { panels: MAIN, placement: null },
        { panels: FLAP, placement: null },
      ])
    ).toBe(false);
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
  const main: GeometryElement = { panels: MAIN, placement: null };
  const flap: GeometryElement = { panels: FLAP, placement: PLACEMENT };

  it('distinguishes element count', async () => {
    const one = await computeGeometryHash([main]);
    const two = await computeGeometryHash([main, main]);
    expect(one).not.toBe(two);
  });

  it('distinguishes element order', async () => {
    const bareFlap: GeometryElement = { panels: FLAP, placement: null };
    const ab = await computeGeometryHash([main, bareFlap]);
    const ba = await computeGeometryHash([bareFlap, main]);
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

// ---------------------------------------------------------------------------
// Shared cross-language vectors. These digests are pinned as the same literal
// hex strings in packages/flexfoil-python/tests/test_geometry_hash.py.
// ---------------------------------------------------------------------------

describe('shared assembly digest — cross-language vectors', () => {
  const main: GeometryElement = { panels: MAIN, placement: null };
  const flap: GeometryElement = { panels: FLAP, placement: PLACEMENT };

  it('canonicalises the assembly vector exactly as Python does', () => {
    expect(sharedCanonicalGeometryString([main, flap])).toBe(SHARED_ASSEMBLY_CANONICAL);
    expect(canonicalGeometryString([main, flap])).toBe(SHARED_ASSEMBLY_CANONICAL);
  });

  it('pins the assembly vector digest shared with Python', async () => {
    expect(await computeGeometryHash([main, flap])).toBe(SHARED_ASSEMBLY_DIGEST);
  });

  it('canonicalises the float-formatting vector exactly as Python does', () => {
    expect(
      canonicalGeometryString([{ panels: EDGE, placement: EDGE_PLACEMENT }])
    ).toBe(SHARED_EDGE_CANONICAL);
  });

  it('pins the float-formatting vector digest shared with Python', async () => {
    expect(await computeGeometryHash([{ panels: EDGE, placement: EDGE_PLACEMENT }])).toBe(
      SHARED_EDGE_DIGEST
    );
  });

  it('normalises signed zero in coordinates as well as placement fields', async () => {
    const negZeroCoord = await computeGeometryHash([
      { panels: [{ x: -0, y: -0 }], placement: PLACEMENT },
    ]);
    const posZeroCoord = await computeGeometryHash([
      { panels: [{ x: 0, y: 0 }], placement: PLACEMENT },
    ]);
    expect(negZeroCoord).toBe(posZeroCoord);

    const negZeroField = await computeAirfoilHash(MAIN, { translation: { x: -0, y: 0.1 } });
    const posZeroField = await computeAirfoilHash(MAIN, { translation: { x: 0, y: 0.1 } });
    expect(negZeroField).toBe(posZeroField);
  });

  it('keeps the sign of a non-zero negative that rounds to zero', () => {
    // Both languages emit "-0.00000000" here; only exact -0 is normalised.
    expect(
      sharedCanonicalGeometryString([
        { panels: [{ x: -1e-12, y: 0 }], placement: PLACEMENT },
      ])
    ).toContain('-0.00000000,0.00000000');
  });
});
