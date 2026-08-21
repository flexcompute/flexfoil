import { describe, it, expect } from 'vitest';
import {
  findLeadingEdgeIndex,
  sampleMidpointY,
  splitContourSurfaces,
  SurfaceSplitError,
  trailingEdgePoint,
  type SurfaceSplitPoint,
} from './surfaceSplit';

/** NACA 4-digit half-thickness. */
function halfThickness(x: number, t = 0.12): number {
  return (
    5 *
    t *
    (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1015 * x ** 4)
  );
}

/** NACA 4-digit mean camber line (max camber `m` at 40% chord). */
function camberLine(x: number, m: number): number {
  const p = 0.4;
  return x < p
    ? (m / (p * p)) * (2 * p * x - x * x)
    : (m / ((1 - p) * (1 - p))) * (1 - 2 * p + 2 * p * x - x * x);
}

/**
 * Build a contour in the conventional TE -> upper -> LE -> lower -> TE order.
 *
 * The upper surface owns `nUpper` nodes (index 0 is the TE, index `nUpper - 1`
 * is the LE) and the lower surface owns the remaining `nLower` nodes, so the
 * true leading-edge index is always `nUpper - 1` and the paneling density of
 * the two surfaces can be varied independently.
 */
function buildContour(nUpper: number, nLower: number, camber = 0): SurfaceSplitPoint[] {
  const pts: SurfaceSplitPoint[] = [];
  for (let i = 0; i < nUpper; i++) {
    const x = 0.5 * (1 + Math.cos((Math.PI * i) / (nUpper - 1)));
    pts.push({ x, y: camberLine(x, camber) + halfThickness(x) });
  }
  for (let j = 1; j <= nLower; j++) {
    const x = 0.5 * (1 - Math.cos((Math.PI * j) / nLower));
    pts.push({ x, y: camberLine(x, camber) - halfThickness(x) });
  }
  return pts;
}

/** Rotate/offset a contour rigidly (nose-down rotation about `pivot`). */
function transform(
  pts: SurfaceSplitPoint[],
  degrees: number,
  dx: number,
  dy: number,
  pivot = 0.5,
): SurfaceSplitPoint[] {
  const a = (degrees * Math.PI) / 180;
  const c = Math.cos(a);
  const s = Math.sin(a);
  return pts.map(({ x, y }) => ({
    x: pivot + (x - pivot) * c - y * s + dx,
    y: (x - pivot) * s + y * c + dy,
  }));
}

/**
 * Place a copy of a contour as a second element: scaled to `chord`, rotated by
 * `degrees` about its leading edge and moved so that leading edge sits at
 * (`x0`, `y0`) - i.e. what a deflected flap looks like behind a main element.
 */
function placeElement(
  pts: SurfaceSplitPoint[],
  chord: number,
  degrees: number,
  x0: number,
  y0: number,
): SurfaceSplitPoint[] {
  const a = (degrees * Math.PI) / 180;
  const c = Math.cos(a);
  const s = Math.sin(a);
  return pts.map(({ x, y }) => ({
    x: x0 + chord * (x * c - y * s),
    y: y0 + chord * (x * s + y * c),
  }));
}

/** Run `fn` and return whatever it threw, or `null` when it did not throw. */
function captureError(fn: () => unknown): unknown {
  try {
    fn();
    return null;
  } catch (e) {
    return e;
  }
}

/** The split this module replaces: bare index arithmetic on the node count. */
function midpointSplit(nodeCount: number, panelCount: number) {
  const mid = Math.floor(nodeCount / 2);
  const upper: number[] = [];
  const lower: number[] = [];
  for (let i = 0; i < panelCount; i++) {
    if (i < mid) upper.push(i);
    else lower.push(i);
  }
  return { upper, lower };
}

function indexRange(start: number, end: number): number[] {
  const out: number[] = [];
  for (let i = start; i < end; i++) out.push(i);
  return out;
}

/** Naive "leading edge is the leftmost node" rule, for contrast. */
function minXIndex(pts: SurfaceSplitPoint[]): number {
  return pts.reduce((best, p, i) => (p.x < pts[best].x ? i : best), 0);
}

describe('findLeadingEdgeIndex', () => {
  it('finds the LE of a conventional section', () => {
    const pts = buildContour(21, 20);
    expect(findLeadingEdgeIndex(pts)).toBe(20);
  });

  it('is invariant to rotation and translation, where minimum x is not', () => {
    const pts = buildContour(21, 20);
    // A flap-sized deflection plus an arbitrary offset, i.e. what an element of
    // a deployed multi-element configuration actually looks like.
    const moved = transform(pts, 25, -3, 7.5);
    expect(findLeadingEdgeIndex(moved)).toBe(20);
    // The naive "leftmost node" rule disagrees at this deflection - that is why
    // the LE is taken as the furthest node from the TE instead.
    expect(minXIndex(moved)).not.toBe(20);
    // Small rotations are benign for minimum x; the failure is not immediate,
    // which is exactly what makes it easy to miss.
    expect(minXIndex(transform(pts, 5, 0, 0))).toBe(20);
  });

  it('takes the TE as the midpoint of the first and last node (blunt TE)', () => {
    const pts: SurfaceSplitPoint[] = [
      { x: 1, y: 0.01 },
      { x: 0.5, y: 0.06 },
      { x: 0, y: 0 },
      { x: 0.5, y: -0.06 },
      { x: 1, y: -0.01 },
    ];
    expect(trailingEdgePoint(pts)).toEqual({ x: 1, y: 0 });
    expect(findLeadingEdgeIndex(pts)).toBe(2);
  });

  it('does not throw on an empty contour', () => {
    expect(findLeadingEdgeIndex([])).toBe(0);
  });
});

describe('splitContourSurfaces - conventional single-element section', () => {
  it('reproduces the midpoint split when the LE really is at the midpoint', () => {
    // 41 nodes, LE at node 20 === Math.floor(41 / 2): the one case the old
    // index arithmetic got right, so the visual result must be unchanged.
    const pts = buildContour(21, 20);
    expect(pts.length).toBe(41);
    expect(Math.floor(pts.length / 2)).toBe(20);

    const split = splitContourSurfaces(pts);
    const old = midpointSplit(pts.length, pts.length - 1);

    expect(split.leadingEdgeIndex).toBe(20);
    expect(split.reversed).toBe(false);
    expect(split.upperIndices).toEqual(old.upper);
    expect(split.lowerIndices).toEqual(old.lower);
    expect(split.upperNodeSlice).toEqual([0, 21]);
    expect(split.lowerNodeSlice).toEqual([20, 41]);
  });

  it('partitions every panel exactly once', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts);
    expect([...split.upperIndices, ...split.lowerIndices]).toEqual(
      indexRange(0, pts.length - 1),
    );
  });

  it('puts positive-y nodes on the upper surface and negative-y on the lower', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts);
    const upperNodes = pts.slice(...split.upperNodeSlice);
    const lowerNodes = pts.slice(...split.lowerNodeSlice);
    expect(upperNodes.every((p) => p.y >= 0)).toBe(true);
    expect(lowerNodes.every((p) => p.y <= 0)).toBe(true);
  });

  it('corrects the off-by-one on an even node count', () => {
    // 40 nodes: the LE sits at 19, but Math.floor(40 / 2) is 20. Real XFOIL
    // paneling hits this (testdata/naca0012_xfoil_paneled.dat has 160 nodes
    // with its LE at index 79), so the old split mislabelled one panel.
    const pts = buildContour(20, 20);
    expect(pts.length).toBe(40);
    const split = splitContourSurfaces(pts);
    expect(split.leadingEdgeIndex).toBe(19);
    expect(split.upperIndices).toEqual(indexRange(0, 19));
    expect(midpointSplit(40, 39).upper).toEqual(indexRange(0, 20));
  });
});

describe('splitContourSurfaces - non-uniform paneling', () => {
  it('splits at the real LE when the LE is nowhere near the midpoint', () => {
    // 20 upper nodes against 6 lower nodes: LE at 19, midpoint says 13.
    const pts = buildContour(20, 6);
    expect(pts.length).toBe(26);

    const split = splitContourSurfaces(pts);
    expect(split.leadingEdgeIndex).toBe(19);
    expect(split.upperIndices).toEqual(indexRange(0, 19));
    expect(split.lowerIndices).toEqual(indexRange(19, 25));

    // The old behaviour, for the record: six panels of the upper surface were
    // being drawn as the lower surface, with no error reported.
    const old = midpointSplit(pts.length, pts.length - 1);
    expect(old.upper).toEqual(indexRange(0, 13));
    const misassigned = old.lower.filter((i) => split.upperIndices.includes(i));
    expect(misassigned).toEqual(indexRange(13, 19));
  });

  it('keeps upper and lower surfaces on the correct side of the chord', () => {
    const pts = buildContour(20, 6);
    const split = splitContourSurfaces(pts);
    expect(pts.slice(...split.upperNodeSlice).every((p) => p.y >= 0)).toBe(true);
    expect(pts.slice(...split.lowerNodeSlice).every((p) => p.y <= 0)).toBe(true);

    // Under the old split, "upper" leaked onto the lower surface.
    const old = midpointSplit(pts.length, pts.length - 1);
    expect(pts.slice(0, Math.floor(pts.length / 2) + 1).every((p) => p.y >= 0)).toBe(true);
    expect(old.lower.some((i) => pts[i].y > 0)).toBe(true);
  });
});

describe('splitContourSurfaces - rotated and offset section', () => {
  it('splits a cambered section rotated nose-down and translated', () => {
    const pts = buildContour(18, 11, 0.04);
    const moved = transform(pts, 25, 12, -4);
    const split = splitContourSurfaces(moved);

    expect(split.leadingEdgeIndex).toBe(17);
    expect(split.reversed).toBe(false);
    expect(split.upperIndices).toEqual(indexRange(0, 17));
    expect(split.lowerIndices).toEqual(indexRange(17, 28));

    // Same partition as the un-transformed geometry: the split follows the
    // section, not the axes.
    expect(splitContourSurfaces(pts).upperIndices).toEqual(split.upperIndices);
  });
});

describe('splitContourSurfaces - reversed contour ordering', () => {
  it('detects TE -> lower -> LE -> upper -> TE and swaps the surfaces', () => {
    const forward = buildContour(21, 20, 0.03);
    const reversed = [...forward].reverse();

    const split = splitContourSurfaces(reversed);
    expect(split.reversed).toBe(true);
    expect(split.leadingEdgeIndex).toBe(20);
    // First segment is now the lower surface.
    expect(split.lowerIndices).toEqual(indexRange(0, 20));
    expect(split.upperIndices).toEqual(indexRange(20, 40));
    expect(split.upperNodeSlice).toEqual([20, 41]);
    expect(split.lowerNodeSlice).toEqual([0, 21]);

    const upperNodes = reversed.slice(...split.upperNodeSlice);
    const lowerNodes = reversed.slice(...split.lowerNodeSlice);
    // The camber line is positive everywhere, so compare against it rather
    // than against y = 0.
    expect(upperNodes.every((p) => p.y >= camberLine(p.x, 0.03) - 1e-12)).toBe(true);
    expect(lowerNodes.every((p) => p.y <= camberLine(p.x, 0.03) + 1e-12)).toBe(true);
  });

  it('still detects reversed ordering after rotation', () => {
    const reversed = transform([...buildContour(21, 20)].reverse(), 20, 0, 0);
    expect(splitContourSurfaces(reversed).reversed).toBe(true);
  });

  it('keeps the last node-based sample readable when the upper surface comes second', () => {
    // The viscous path reports Cp at nodes, so the last sample is node n - 1 and
    // has no node after it. On a reversed contour that sample belongs to the
    // upper surface, so a panel-midpoint lookup written as
    // `(nodes[i].y + nodes[i + 1].y) / 2` reads past the end of the array -
    // which is what `sampleMidpointY` exists to prevent, on both surfaces.
    const reversed = [...buildContour(21, 20)].reverse();
    const split = splitContourSurfaces(reversed, { sampleCount: reversed.length });
    const last = reversed.length - 1;

    expect(split.reversed).toBe(true);
    expect(split.upperIndices).toContain(last);
    expect(reversed[last + 1]).toBeUndefined();

    for (const i of [...split.upperIndices, ...split.lowerIndices]) {
      expect(Number.isFinite(sampleMidpointY(reversed, i))).toBe(true);
    }
    expect(sampleMidpointY(reversed, last)).toBe(reversed[last].y);
  });
});

describe('splitContourSurfaces - more than one contour in the array', () => {
  /** A main element plus a deflected flap, concatenated into one array. */
  function mainPlusFlap(): SurfaceSplitPoint[] {
    const main = buildContour(21, 20);
    const flap = placeElement(buildContour(11, 10), 0.3, -25, 0.85, -0.03);
    expect(main.length).toBe(41);
    expect(flap.length).toBe(21);
    return [...main, ...flap];
  }

  it('reports a concatenated main-plus-flap array rather than splitting it', () => {
    const pts = mainPlusFlap();
    expect(pts.length).toBe(62);

    // Left to the geometric derivation, the array yields the main element's LE,
    // which puts every flap node on one surface: a single curve that jumps
    // between two elements. Two contours have no one upper/lower pair.
    expect(findLeadingEdgeIndex(pts)).toBe(20);

    const error = captureError(() => splitContourSurfaces(pts));
    expect(error).toBeInstanceOf(SurfaceSplitError);
    expect((error as SurfaceSplitError).code).toBe('multiple-contours');
    expect((error as SurfaceSplitError).message).toContain('node 41');
    expect((error as SurfaceSplitError).message).toContain(
      'per-element splitting is not yet supported',
    );
  });

  it('reports a slat-main-flap array, where the first closure found is not the split', () => {
    // Three elements in file order, as testdata/mda_30p_30n_trimmed.dat lists
    // them. The first element's apparent closure can fall a node early, so
    // every closure of the first loop has to be tried before concluding that
    // the array holds one contour.
    const slat = placeElement(buildContour(9, 8), 0.15, -30, -0.1, 0.01);
    const main = buildContour(21, 20);
    const flap = placeElement(buildContour(11, 10), 0.3, -25, 0.85, -0.03);
    const pts = [...slat, ...main, ...flap];

    const error = captureError(() => splitContourSurfaces(pts));
    expect((error as SurfaceSplitError).code).toBe('multiple-contours');
    expect((error as SurfaceSplitError).message).toContain(`of ${pts.length}`);
  });

  it('reports it whatever the caller says the leading edge is', () => {
    const pts = mainPlusFlap();
    for (const options of [{}, { leadingEdgeIndex: 20 }, { sampleCount: pts.length }]) {
      const error = captureError(() => splitContourSurfaces(pts, options));
      expect((error as SurfaceSplitError).code).toBe('multiple-contours');
    }
  });

  it('leaves single elements alone, including fine paneling and a blunt TE', () => {
    const blunt = buildContour(81, 80);
    blunt[0] = { x: 1, y: 0.01 };
    blunt[blunt.length - 1] = { x: 1, y: -0.01 };

    const singles: SurfaceSplitPoint[][] = [
      buildContour(21, 20),
      buildContour(20, 6),
      buildContour(41, 40, 0.04),
      // 161 nodes, i.e. XFOIL-scale paneling: the nodes either side of the TE
      // are a few thousandths of a chord apart, which is the case most likely
      // to read as a closure.
      buildContour(81, 80),
      blunt,
      transform(buildContour(21, 20), 25, -3, 7.5),
      [...buildContour(21, 20)].reverse(),
    ];

    for (const pts of singles) {
      expect(captureError(() => splitContourSurfaces(pts))).toBeNull();
    }
  });
});

describe('splitContourSurfaces - explicit leading-edge index', () => {
  it('prefers a caller-supplied LE index over the derived one', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts, { leadingEdgeIndex: 7 });
    expect(split.leadingEdgeIndex).toBe(7);
    expect(split.upperIndices).toEqual(indexRange(0, 7));
    expect(split.lowerIndices).toEqual(indexRange(7, 40));
    expect(split.upperNodeSlice).toEqual([0, 8]);
    expect(split.lowerNodeSlice).toEqual([7, 41]);
  });

  it('honours an explicit index on a contour where derivation would disagree', () => {
    const pts = buildContour(20, 6);
    expect(splitContourSurfaces(pts).leadingEdgeIndex).toBe(19);
    expect(splitContourSurfaces(pts, { leadingEdgeIndex: 19 }).leadingEdgeIndex).toBe(19);
    expect(splitContourSurfaces(pts, { leadingEdgeIndex: 4 }).leadingEdgeIndex).toBe(4);
  });

  it('derives the LE only when no index is supplied', () => {
    const pts = buildContour(21, 20);
    for (const options of [{}, { leadingEdgeIndex: null }, { leadingEdgeIndex: undefined }]) {
      expect(splitContourSurfaces(pts, options).leadingEdgeIndex).toBe(20);
    }
  });

  it('rejects an index that does not split the contour, rather than clamping it', () => {
    // 41 nodes and 40 panel samples, so only [1, 39] divides the contour in two.
    // Clamping such an index - or quietly deriving one instead - would answer a
    // caller's wrong index with a wrong split and no error.
    const pts = buildContour(21, 20);
    for (const bad of [-1, 0, 40, 99, 3.5, NaN]) {
      const error = captureError(() => splitContourSurfaces(pts, { leadingEdgeIndex: bad }));
      expect(error).toBeInstanceOf(SurfaceSplitError);
      expect((error as SurfaceSplitError).code).toBe('invalid-leading-edge-index');
      expect((error as SurfaceSplitError).message).toContain('[1, 39]');
    }
  });

  it('rejects an index that no sample reaches', () => {
    // Ten samples of a 41-node contour cover panels 0-9, so an LE at node 20
    // would leave the second surface empty.
    const pts = buildContour(21, 20);
    const error = captureError(() =>
      splitContourSurfaces(pts, { leadingEdgeIndex: 20, sampleCount: 10 }),
    );
    expect((error as SurfaceSplitError).code).toBe('invalid-leading-edge-index');
    expect((error as SurfaceSplitError).message).toContain('[1, 9]');
  });
});

describe('splitContourSurfaces - sample count', () => {
  it('defaults to one sample per panel', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts);
    expect(split.upperIndices.concat(split.lowerIndices)).toEqual(indexRange(0, 40));
  });

  it('handles one sample per node, as the viscous solver reports Cp', () => {
    // The viscous path returns Cp at nodes (n samples) while the inviscid path
    // returns it at panel midpoints (n - 1). The extra trailing sample stays on
    // the lower surface instead of being dropped; the node slice stops at the
    // end of the contour, so callers must guard their own node lookups.
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts, { sampleCount: pts.length });
    expect(split.leadingEdgeIndex).toBe(20);
    expect(split.upperIndices).toEqual(indexRange(0, 20));
    expect(split.lowerIndices).toEqual(indexRange(20, 41));
    expect(split.lowerNodeSlice).toEqual([20, 41]);
  });

  it('respects a shorter sampleCount than the node count implies', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts, { sampleCount: 30 });
    expect(split.leadingEdgeIndex).toBe(20);
    expect(split.upperIndices).toEqual(indexRange(0, 20));
    expect(split.lowerIndices).toEqual(indexRange(20, 30));
    expect(split.lowerNodeSlice).toEqual([20, 31]);
  });

  it('keeps both surfaces non-empty when there are fewer samples than the LE index', () => {
    const pts = buildContour(21, 20);
    const split = splitContourSurfaces(pts, { sampleCount: 5 });
    expect(split.leadingEdgeIndex).toBe(4);
    expect(split.upperIndices).toEqual(indexRange(0, 4));
    expect(split.lowerIndices).toEqual([4]);
  });
});

describe('splitContourSurfaces - degenerate input', () => {

  it('never emits an empty surface for a splittable contour', () => {
    const pts: SurfaceSplitPoint[] = [
      { x: 1, y: 0 },
      { x: 0, y: 0.05 },
      { x: 1, y: -0.01 },
    ];
    const split = splitContourSurfaces(pts);
    expect(split.upperIndices.length).toBeGreaterThan(0);
    expect(split.lowerIndices.length).toBeGreaterThan(0);
  });

  it('degrades to a single surface when there is not enough contour', () => {
    const split = splitContourSurfaces([
      { x: 1, y: 0 },
      { x: 0, y: 0 },
    ]);
    expect(split.upperIndices).toEqual([0]);
    expect(split.lowerIndices).toEqual([]);
    expect(split.reversed).toBe(false);

    const empty = splitContourSurfaces([]);
    expect(empty.upperIndices).toEqual([]);
    expect(empty.lowerIndices).toEqual([]);
  });

  it('handles a zero-area contour (flat plate) without flipping the surfaces', () => {
    const pts: SurfaceSplitPoint[] = [
      { x: 1, y: 0 },
      { x: 0.5, y: 0 },
      { x: 0, y: 0 },
      { x: 0.5, y: 0 },
      { x: 1, y: 0 },
    ];
    const split = splitContourSurfaces(pts);
    expect(split.leadingEdgeIndex).toBe(2);
    expect(split.reversed).toBe(false);
    expect(split.upperIndices).toEqual([0, 1]);
    expect(split.lowerIndices).toEqual([2, 3]);
  });
});
