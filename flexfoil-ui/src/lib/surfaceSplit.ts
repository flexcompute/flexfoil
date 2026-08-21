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
 *
 * Scope: one contour per call. Several elements concatenated into a single
 * array have no single upper/lower pair to split into, so that input is
 * reported as a `SurfaceSplitError` with code `multiple-contours` instead of
 * being split as though it were one element. Splitting per element needs
 * element offsets from the producer of the geometry, which is a later
 * workstream.
 */

/** Minimal point shape needed to split a contour. */
export interface SurfaceSplitPoint {
  x: number;
  y: number;
}

/** Why a contour could not be split. */
export type SurfaceSplitErrorCode = 'invalid-leading-edge-index' | 'multiple-contours';

/**
 * Raised instead of returning a split that would not describe the input.
 *
 * Callers are expected to surface the message: the whole point of the module is
 * that an unsplittable contour is reported rather than drawn as two curves that
 * look like surfaces but are not.
 */
export class SurfaceSplitError extends Error {
  readonly code: SurfaceSplitErrorCode;

  constructor(code: SurfaceSplitErrorCode, message: string) {
    super(message);
    this.name = 'SurfaceSplitError';
    this.code = code;
  }
}

export interface SurfaceSplitOptions {
  /**
   * Leading-edge node index reported by the producer of the contour (for
   * example a solver that knows the per-element LE). Used exactly as given when
   * it is a splittable interior index, i.e. an integer in
   * `[1, min(nodes.length - 2, sampleCount - 1)]`.
   *
   * Any other value - including `0`, the last node, or an index past the last
   * sample - raises `SurfaceSplitError`. It is neither clamped nor quietly
   * replaced by the geometric derivation, because either would turn a caller's
   * wrong index into a wrong split with no error reported, which is the failure
   * mode this module exists to remove. `null` and `undefined` mean "not
   * supplied" and select the geometric derivation.
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

/**
 * Fraction of a loop's own size within which its two ends must meet for the
 * loop to count as closed.
 */
const CLOSURE_TOL_REL = 0.02;

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
 * need a splittable index should use `splitContourSurfaces`, which clamps the
 * derived index to the interior of the contour.
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
 * Every index at which the traverse starting at `start` comes back to its own
 * first node, i.e. every way the nodes from `start` on could close a loop.
 *
 * Closure is measured against the size of the loop so far (its bounding-box
 * diagonal), which makes the test independent of chord, position and units. It
 * also means a monotone run of nodes - whose ends are as far apart as the run is
 * long - never reads as closed, however finely it is paneled.
 */
function closureCandidates(nodes: readonly SurfaceSplitPoint[], start: number): number[] {
  const n = nodes.length;
  const first = nodes[start];
  const found: number[] = [];

  let minX = first.x;
  let maxX = first.x;
  let minY = first.y;
  let maxY = first.y;

  for (let j = start + 1; j < n; j++) {
    const p = nodes[j];
    minX = Math.min(minX, p.x);
    maxX = Math.max(maxX, p.x);
    minY = Math.min(minY, p.y);
    maxY = Math.max(maxY, p.y);

    // A loop needs three distinct nodes plus the one that closes it.
    if (j - start < 3) continue;
    const extent = Math.hypot(maxX - minX, maxY - minY);
    if (extent <= 0) continue;
    if (Math.hypot(p.x - first.x, p.y - first.y) <= CLOSURE_TOL_REL * extent) found.push(j);
  }

  return found;
}

/**
 * Index of the first node of a second contour, or `null` when the nodes look
 * like a single contour.
 *
 * A single contour closes once, at its own end. Several elements concatenated
 * into one array leave a signature instead: the traverse closes a loop part way
 * through and the nodes after it close a loop of their own. Every closure of the
 * first loop is tried, because a contour whose trailing edge is paneled finely
 * can appear to close a node or two before it really does.
 *
 * This is a cheap detection test rather than per-element splitting: elements
 * whose own trailing edges are blunt relative to their chord can still read as
 * one contour. It is here so that the common concatenation is reported instead
 * of being drawn as a single curve spanning two elements.
 */
function findSecondContourStart(nodes: readonly SurfaceSplitPoint[]): number | null {
  const n = nodes.length;
  // Two loops need four nodes each, counting the node that closes them.
  if (n < 8) return null;

  for (const firstClose of closureCandidates(nodes, 0)) {
    const secondStart = firstClose + 1;
    if (secondStart > n - 4) break;
    if (closureCandidates(nodes, secondStart).length > 0) return secondStart;
  }

  return null;
}

/**
 * `y` at the midpoint of the panel a sample sits on.
 *
 * Node-based surface data has one sample more than there are panels (see
 * `sampleCount`), so its last sample has no following node and falls back to
 * the node itself instead of reading past the end of the array. A sample with
 * no node at all returns `NaN`, which plots as a gap.
 */
export function sampleMidpointY(nodes: readonly SurfaceSplitPoint[], index: number): number {
  const a = nodes[index];
  if (!a) return Number.NaN;
  const b = nodes[index + 1] ?? a;
  return (a.y + b.y) / 2;
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
 *
 * Throws `SurfaceSplitError` when `leadingEdgeIndex` is supplied but not
 * splittable, and when the nodes hold more than one contour.
 */
export function splitContourSurfaces(
  nodes: readonly SurfaceSplitPoint[],
  options: SurfaceSplitOptions = {},
): SurfaceSplit {
  const n = nodes.length;
  const requested = options.sampleCount ?? Math.max(0, n - 1);
  const sampleCount = Math.max(0, requested);

  const explicit = options.leadingEdgeIndex;
  const maxLeIndex = Math.min(n - 2, sampleCount - 1);
  if (explicit != null && (!Number.isInteger(explicit) || explicit < 1 || explicit > maxLeIndex)) {
    throw new SurfaceSplitError(
      'invalid-leading-edge-index',
      `leadingEdgeIndex ${explicit} does not split a contour of ${n} node(s) ` +
        `with ${sampleCount} sample(s): expected an integer in [1, ${maxLeIndex}].`,
    );
  }

  const secondContourStart = findSecondContourStart(nodes);
  if (secondContourStart !== null) {
    throw new SurfaceSplitError(
      'multiple-contours',
      `multiple contours detected (a second one starts at node ${secondContourStart} ` +
        `of ${n}); per-element splitting is not yet supported.`,
    );
  }

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

  // An explicit index is validated above rather than clamped, so it arrives
  // here already inside [1, maxLeIndex]. A derived index is clamped
  // defensively: the argmax can land on a contour endpoint for a degenerate
  // contour, and both surfaces must own at least one sample - and the upper
  // surface must stay inside the node array - for the split to mean anything.
  const leIndex =
    explicit != null ? explicit : Math.min(Math.max(findLeadingEdgeIndex(nodes), 1), maxLeIndex);

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
