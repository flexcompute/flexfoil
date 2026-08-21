/**
 * surfaceSplit - locate the leading edge of an airfoil contour and split its
 * surface samples into upper and lower.
 *
 * Why this module exists: splitting a contour at the array midpoint
 * (`Math.floor(n / 2)`) is only correct when the paneling happens to be
 * symmetric about the leading edge. It is off by one for any even node count,
 * arbitrarily wrong for non-uniform paneling, and catastrophically wrong for a
 * concatenated multi-element array - in every case producing plausible-looking
 * but incorrect upper/lower curves with no error reported.
 *
 * The split is therefore derived from the geometry:
 *
 * - The leading edge is the node at maximum distance from the trailing edge
 *   (the midpoint of the first and last node). This is equivalent to the naive
 *   "minimum x" rule for a conventional section, but stays correct for
 *   cambered, offset and rotated sections - a flap element deflected 25 deg or
 *   more already breaks the minimum-x rule.
 * - Which side of the LE is the upper surface is decided from the contour's
 *   orientation (signed area), which is invariant under rotation and
 *   translation, rather than from array position.
 *
 * A caller that already knows the leading edge - e.g. a solver reporting a
 * per-element LE index for a multi-element configuration - can pass it in via
 * `leadingEdgeIndex` and skip the derivation entirely.
 */

/** Minimal point shape needed to split a contour. */
export interface SurfaceSplitPoint {
  x: number;
  y: number;
}

export interface SurfaceSplitOptions {
  /**
   * Leading-edge node index reported by the producer of the contour (for
   * example a solver that knows the per-element LE). Preferred over the
   * geometric derivation whenever it is a usable interior index
   * (`0 < i < nodes.length - 1`); ignored otherwise, because clamping a
   * nonsensical index would silently produce a wrong split - exactly the
   * failure mode this module exists to remove.
   */
  leadingEdgeIndex?: number | null;
  /**
   * How many per-station samples the caller holds (e.g. `cp.length`). Defaults
   * to `nodes.length - 1`, the panel count of a contour of `nodes`.
   *
   * The count is taken as given rather than clamped, because solvers do not
   * agree on where surface quantities live: this project's inviscid path
   * reports Cp at panel midpoints (`nodes.length - 1` samples) while the
   * viscous path reports it at nodes (`nodes.length` samples). A caller passing
   * more samples than there are panels therefore gets trailing indices with no
   * panel geometry behind them, and must guard its own node lookups.
   */
  sampleCount?: number;
}

export interface SurfaceSplit {
  /** Node index of the leading edge used for the split. */
  leadingEdgeIndex: number;
  /**
   * True when the contour runs TE -> lower -> LE -> upper -> TE, i.e. the
   * reverse of the usual TE -> upper -> LE -> lower -> TE ordering.
   */
  reversed: boolean;
  /** Sample indices on the upper surface, in contour order. */
  upperIndices: number[];
  /** Sample indices on the lower surface, in contour order. */
  lowerIndices: number[];
  /** `[start, end)` node slice spanning the upper surface, in contour order. */
  upperNodeSlice: [number, number];
  /** `[start, end)` node slice spanning the lower surface, in contour order. */
  lowerNodeSlice: [number, number];
}

/** Relative tolerance (against chord^2) below which a contour has no usable orientation. */
const AREA_EPS_REL = 1e-9;

function indexRange(start: number, end: number): number[] {
  const out: number[] = [];
  for (let i = start; i < end; i++) out.push(i);
  return out;
}

/**
 * Trailing-edge reference point: the midpoint of the first and last node.
 *
 * For a sharp (closed) contour the two coincide and this is the TE itself; for
 * a blunt TE it is the centre of the TE base.
 */
export function trailingEdgePoint(nodes: readonly SurfaceSplitPoint[]): SurfaceSplitPoint {
  if (nodes.length === 0) return { x: 0, y: 0 };
  const first = nodes[0];
  const last = nodes[nodes.length - 1];
  return { x: (first.x + last.x) / 2, y: (first.y + last.y) / 2 };
}

/**
 * Index of the leading-edge node: the node furthest from the trailing edge.
 *
 * Ties resolve to the first such node. Returns the raw argmax - callers that
 * need a splittable index should use `splitSurfaces`, which clamps it to the
 * interior of the contour.
 */
export function findLeadingEdgeIndex(nodes: readonly SurfaceSplitPoint[]): number {
  const n = nodes.length;
  if (n === 0) return 0;

  const te = trailingEdgePoint(nodes);
  let bestIdx = 0;
  let bestD2 = -1;
  for (let i = 0; i < n; i++) {
    const dx = nodes[i].x - te.x;
    const dy = nodes[i].y - te.y;
    const d2 = dx * dx + dy * dy;
    if (d2 > bestD2) {
      bestD2 = d2;
      bestIdx = i;
    }
  }
  return bestIdx;
}

/** Shoelace signed area of the contour, closing last -> first. Positive = counter-clockwise. */
function signedArea(nodes: readonly SurfaceSplitPoint[]): number {
  const n = nodes.length;
  let acc = 0;
  for (let i = 0; i < n; i++) {
    const p = nodes[i];
    const q = nodes[(i + 1) % n];
    acc += p.x * q.y - q.x * p.y;
  }
  return acc / 2;
}

/**
 * Mean offset of a node range from the chord line, measured along the chord
 * normal that points "up" for a contour whose chord runs TE -> LE.
 */
function meanChordOffset(
  nodes: readonly SurfaceSplitPoint[],
  start: number,
  end: number,
  te: SurfaceSplitPoint,
  le: SurfaceSplitPoint,
): number {
  // Chord vector TE -> LE; the "up" normal is that vector rotated by -90 deg,
  // which reduces to +y for a conventional nose-left/TE-right section.
  const nx = le.y - te.y;
  const ny = -(le.x - te.x);
  let acc = 0;
  let count = 0;
  for (let i = start; i < end; i++) {
    acc += (nodes[i].x - te.x) * nx + (nodes[i].y - te.y) * ny;
    count++;
  }
  return count > 0 ? acc / count : 0;
}

/**
 * True when the first segment of the contour (node 0 to the LE) is the lower
 * surface rather than the upper one.
 *
 * Decided from the contour orientation, which is rotation- and
 * translation-invariant: TE -> upper -> LE -> lower -> TE traverses an airfoil
 * counter-clockwise, so a clockwise contour is reverse-ordered. Degenerate
 * zero-area contours (e.g. a flat plate) fall back to comparing the two
 * segments across the chord line.
 */
function isReversedOrder(nodes: readonly SurfaceSplitPoint[], leIndex: number): boolean {
  const te = trailingEdgePoint(nodes);
  const le = nodes[leIndex];
  const chord2 = (le.x - te.x) ** 2 + (le.y - te.y) ** 2;
  const area = signedArea(nodes);

  if (Math.abs(area) > AREA_EPS_REL * Math.max(chord2, Number.MIN_VALUE)) {
    return area < 0;
  }

  const firstOffset = meanChordOffset(nodes, 0, leIndex + 1, te, le);
  const secondOffset = meanChordOffset(nodes, leIndex, nodes.length, te, le);
  return secondOffset > firstOffset;
}

/**
 * Split an airfoil contour into upper and lower surfaces at its real leading
 * edge.
 *
 * `nodes` is the contour in traversal order, conventionally
 * TE -> upper -> LE -> lower -> TE (the reverse ordering is detected and
 * handled). Samples `[0, leadingEdgeIndex)` go to the first surface and
 * `[leadingEdgeIndex, sampleCount)` to the second, so for panel-midpoint data
 * the LE panel boundary is the divide and for node data the LE node itself is
 * the first sample of the second surface.
 *
 * Contours too short to have two surfaces are returned as a single upper
 * surface, so callers degrade to one curve rather than to garbage.
 */
export function splitSurfaces(
  nodes: readonly SurfaceSplitPoint[],
  options: SurfaceSplitOptions = {},
): SurfaceSplit {
  const n = nodes.length;
  const requested = options.sampleCount ?? Math.max(0, n - 1);
  const sampleCount = Math.max(0, requested);

  if (n < 3 || sampleCount < 2) {
    return {
      leadingEdgeIndex: 0,
      reversed: false,
      upperIndices: indexRange(0, sampleCount),
      lowerIndices: [],
      upperNodeSlice: [0, Math.min(sampleCount + 1, n)],
      lowerNodeSlice: [n, n],
    };
  }

  const explicit = options.leadingEdgeIndex;
  const useExplicit =
    typeof explicit === 'number' &&
    Number.isInteger(explicit) &&
    explicit > 0 &&
    explicit < n - 1;

  const rawLe = useExplicit ? explicit : findLeadingEdgeIndex(nodes);
  // Both surfaces must own at least one sample, and the upper surface must stay
  // inside the node array, for the split to mean anything.
  const leIndex = Math.min(Math.max(rawLe, 1), n - 2, sampleCount - 1);

  const reversed = isReversedOrder(nodes, leIndex);

  const firstIndices = indexRange(0, leIndex);
  const secondIndices = indexRange(leIndex, sampleCount);
  const firstNodes: [number, number] = [0, leIndex + 1];
  const secondNodes: [number, number] = [leIndex, Math.min(sampleCount + 1, n)];

  return reversed
    ? {
        leadingEdgeIndex: leIndex,
        reversed,
        upperIndices: secondIndices,
        lowerIndices: firstIndices,
        upperNodeSlice: secondNodes,
        lowerNodeSlice: firstNodes,
      }
    : {
        leadingEdgeIndex: leIndex,
        reversed,
        upperIndices: firstIndices,
        lowerIndices: secondIndices,
        upperNodeSlice: firstNodes,
        lowerNodeSlice: secondNodes,
      };
}
