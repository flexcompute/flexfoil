import type { AirfoilPoint } from '../types';

export interface ParsedAirfoilFile {
  name: string;
  coordinates: AirfoilPoint[];
}

/**
 * Result of parsing a `.dat` file, which may describe more than one element
 * (slat / main / flap). `elements` holds every element contour in file order;
 * `coordinates` is `elements[0]`, chosen explicitly for the single-element
 * consumers rather than a concatenation of every block into one multi-loop
 * contour.
 */
export interface ParsedAirfoilDat extends ParsedAirfoilFile {
  elements: AirfoilPoint[][];
}

export interface ImportedAirfoilState extends ParsedAirfoilFile {
  panels: AirfoilPoint[];
}

/**
 * Fewest points that can plausibly describe one element contour. Blocks shorter
 * than this are treated as stray fragments, not elements (see `foldBlocks`).
 */
const MIN_ELEMENT_POINTS = 3;

/**
 * XFOIL / MSES element separator sentinel. A `999.0  999.0` line marks an
 * element boundary; per XFOIL convention any coordinate pair whose both values
 * reach the sentinel is a separator, never a real point.
 */
const ELEMENT_SEPARATOR_SENTINEL = 999;

/** Chord-normalised tolerance for calling a block's first and last point coincident. */
const CLOSURE_TOLERANCE = 1e-3;

/** Tolerance for a backwards step in x while testing surface monotonicity. */
const MONOTONE_TOLERANCE = 1e-4;

/** Minimum x-range overlap (as a fraction of the wider block's span) for two Lednicer surfaces. */
const MIN_LEDNICER_X_OVERLAP = 0.5;

/**
 * Fewest points for a block to be accepted as a separate element on the strength
 * of an annotation boundary alone. A coarse element contour still has to resolve
 * both of its surfaces.
 */
const MIN_ANNOTATION_SPLIT_POINTS = 8;

/**
 * Chordwise distance travelled along a block divided by the block's chordwise
 * span, above which the block is taken to be a contour rather than a single
 * surface. A contour wraps round the leading edge, so x doubles back and the
 * ratio is about 2 (the three elements of `30p-30n.dat` score 1.99 to 2.32); a
 * surface runs one way only and scores about 1.
 */
const MIN_CONTOUR_TRAVERSAL_RATIO = 1.5;

type RepanelAirfoil = (coordinates: AirfoilPoint[], nPanels: number) => AirfoilPoint[];

function fallbackAirfoilName(fileName: string): string {
  return fileName
    .replace(/\.[^.]+$/, '')
    .replace(/[_-]+/g, ' ')
    .trim() || 'Imported Airfoil';
}

/** Matches Lednicer-style count lines like "61 61", "61.0  61.0", "17. 17." */
function isLikelyCountLine(line: string): boolean {
  const trimmed = line.trim();
  if (!/^\d+\.?\d*\s+\d+\.?\d*$/.test(trimmed)) return false;

  const [a, b] = trimmed.split(/\s+/).map(Number);
  return a > 2 && b > 2;
}

/**
 * Matches an XFOIL/MSES element separator line such as `999.0 999.0`, tolerating
 * whitespace and format variation (`999. 999.`, `1000.0  1000.0`, commas).
 */
function isElementSeparatorLine(line: string): boolean {
  const parts = line.trim().replace(/,/g, ' ').split(/\s+/);
  if (parts.length < 2) return false;

  const x = Number(parts[0]);
  const y = Number(parts[1]);
  if (!Number.isFinite(x) || !Number.isFinite(y)) return false;

  return x >= ELEMENT_SEPARATOR_SENTINEL && y >= ELEMENT_SEPARATOR_SENTINEL;
}

function parseCoordinateLine(line: string): AirfoilPoint | null {
  const trimmed = line.trim();
  if (!trimmed || isLikelyCountLine(trimmed)) {
    return null;
  }

  const parts = trimmed.replace(/,/g, ' ').split(/\s+/);
  if (parts.length < 2) {
    return null;
  }

  const x = Number(parts[0]);
  const y = Number(parts[1]);
  if (!Number.isFinite(x) || !Number.isFinite(y)) {
    return null;
  }

  if (x > 1.5 || x < -0.5 || y > 1.5 || y < -1.5) {
    return null;
  }

  return { x, y };
}

function findLeadingEdgeIndex(coordinates: AirfoilPoint[]): number {
  let leIndex = 0;
  for (let i = 1; i < coordinates.length; i += 1) {
    if (coordinates[i].x < coordinates[leIndex].x) {
      leIndex = i;
    }
  }
  return leIndex;
}

function annotateSurfaces(coordinates: AirfoilPoint[]): AirfoilPoint[] {
  const leIndex = findLeadingEdgeIndex(coordinates);

  return coordinates.map((point, index) => ({
    ...point,
    surface: index <= leIndex ? 'upper' : 'lower',
  }));
}

/**
 * What stated an element boundary in a `.dat` file.
 *  - `sentinel`: an XFOIL `999.0 999.0` line — the file declares a boundary.
 *  - `blank`: a blank line or a Lednicer count line.
 *  - `annotation`: a comment or other text line. Honoured as a boundary only
 *    when the resulting blocks look like separate parts (see `splitIntoElements`).
 */
type SeparatorKind = 'sentinel' | 'blank' | 'annotation';

/** One lexical item of a `.dat` file's coordinate section. */
type DatToken =
  | { kind: 'point'; point: AirfoilPoint }
  | { kind: 'separator'; separator: SeparatorKind };

/** Split the file into a header plus a stream of points and element separators. */
function tokenizeDat(text: string): { header: string | null; tokens: DatToken[] } {
  const tokens: DatToken[] = [];
  let header: string | null = null;
  let sawPoint = false;

  for (const line of text.split(/\r?\n/)) {
    const trimmed = line.trim();

    // Checked before coordinate parsing so the sentinel is recognised as a
    // declared element boundary. (`parseCoordinateLine` would reject it as an
    // out-of-range coordinate, leaving only an implicit boundary behind.)
    if (isElementSeparatorLine(trimmed)) {
      tokens.push({ kind: 'separator', separator: 'sentinel' });
      continue;
    }

    const coord = parseCoordinateLine(trimmed);
    if (coord) {
      tokens.push({ kind: 'point', point: coord });
      sawPoint = true;
      continue;
    }

    // Any other line ends the current block. A blank line is the block separator
    // in Lednicer and multi-element files. A comment or other text line may be a
    // boundary too (`30p-30n.dat` in the bundled library labels its elements
    // `# Slat` / `# Main Element` / `# Flap` and uses nothing else) or may
    // annotate one contour; which it is is decided in `splitIntoElements`.
    const structural = !trimmed || isLikelyCountLine(trimmed);
    tokens.push({ kind: 'separator', separator: structural ? 'blank' : 'annotation' });

    if (structural) continue;

    // First non-blank, non-count, non-coordinate line is the header
    if (header === null && !sawPoint) {
      header = trimmed;
    }
  }

  return { header, tokens };
}

/**
 * Fold a token stream into element blocks.
 *
 * A run of fewer than `MIN_ELEMENT_POINTS` points cannot be an element contour,
 * and real `.dat` files do contain stray blank lines mid-contour (`cap21c.dat`
 * in the bundled library separates its closing trailing-edge point with two
 * blank lines). Such a fragment is re-joined to its neighbouring block so no
 * coordinate is lost and no spurious one-point "element" appears. A fragment
 * after an explicit `999.0` separator is kept as its own block: the file
 * declared an element there, so a degenerate element should be reported, not
 * absorbed.
 *
 * `annotationsSplit` selects whether annotation lines act as boundaries; see
 * `splitIntoElements`, which folds the same tokens both ways.
 */
function foldBlocks(tokens: DatToken[], annotationsSplit: boolean): AirfoilPoint[][] {
  const raw: { points: AirfoilPoint[]; explicit: boolean }[] = [];
  let current: AirfoilPoint[] = [];
  // Whether the separator that opened `current` was an explicit 999.0 line.
  let currentExplicit = false;

  for (const token of tokens) {
    if (token.kind === 'point') {
      current.push(token.point);
      continue;
    }

    if (token.separator === 'annotation' && !annotationsSplit) {
      continue;
    }
    const explicit = token.separator === 'sentinel';

    if (current.length > 0) {
      raw.push({ points: current, explicit: currentExplicit });
      current = [];
      currentExplicit = explicit;
    } else {
      // Consecutive separators: an explicit one still marks the next block.
      currentExplicit = currentExplicit || explicit;
    }
  }
  if (current.length > 0) {
    raw.push({ points: current, explicit: currentExplicit });
  }

  const blocks: AirfoilPoint[][] = [];
  // Fragment held over because there is no preceding block to join it to.
  let carry: AirfoilPoint[] = [];

  for (const { points, explicit } of raw) {
    const block = carry.length > 0 ? [...carry, ...points] : points;
    carry = [];

    if (block.length < MIN_ELEMENT_POINTS && !explicit) {
      const previous = blocks[blocks.length - 1];
      if (previous) {
        previous.push(...block);
      } else {
        carry = block;
      }
      continue;
    }

    blocks.push(block);
  }
  if (carry.length > 0) {
    blocks.push(carry);
  }

  return blocks;
}

function xRange(block: AirfoilPoint[]): { min: number; max: number } {
  let min = block[0].x;
  let max = block[0].x;
  for (const p of block) {
    if (p.x < min) min = p.x;
    if (p.x > max) max = p.x;
  }
  return { min, max };
}

/** True when a block's first and last point coincide, i.e. it is a closed loop. */
function isClosedContour(block: AirfoilPoint[]): boolean {
  const first = block[0];
  const last = block[block.length - 1];
  return Math.hypot(last.x - first.x, last.y - first.y) <= CLOSURE_TOLERANCE;
}

/** True when x never steps backwards, i.e. the block is a single LE→TE surface. */
function isMonotoneSurface(block: AirfoilPoint[]): boolean {
  for (let i = 1; i < block.length; i += 1) {
    if (block[i].x - block[i - 1].x < -MONOTONE_TOLERANCE) return false;
  }
  return true;
}

/**
 * True when a block looks like an element contour in its own right, rather than
 * one annotated stretch of a longer contour.
 *
 * A contour that encloses an element wraps round its leading edge, so the
 * chordwise distance travelled along it is about twice its chordwise span. A
 * single surface (`# upper surface` up to `# lower surface`) runs one way only,
 * so the two are about equal. See `MIN_CONTOUR_TRAVERSAL_RATIO`.
 */
function looksLikeElementContour(block: AirfoilPoint[]): boolean {
  if (block.length < MIN_ANNOTATION_SPLIT_POINTS) return false;

  const { min, max } = xRange(block);
  const span = max - min;
  if (span <= 0) return false;

  let travelled = 0;
  for (let i = 1; i < block.length; i += 1) {
    travelled += Math.abs(block[i].x - block[i - 1].x);
  }
  return travelled / span >= MIN_CONTOUR_TRAVERSAL_RATIO;
}

/** Overlap of two blocks' x-ranges, as a fraction of the wider block's span. */
function xOverlapFraction(a: AirfoilPoint[], b: AirfoilPoint[]): number {
  const ra = xRange(a);
  const rb = xRange(b);
  const overlap = Math.min(ra.max, rb.max) - Math.max(ra.min, rb.min);
  const span = Math.max(ra.max - ra.min, rb.max - rb.min);
  if (span <= 0) return 0;
  return Math.max(0, overlap) / span;
}

/**
 * Distinguish a genuine Lednicer file (two blocks that are the two SURFACES of
 * ONE airfoil) from a genuine two-ELEMENT configuration (e.g. main + flap).
 * Both can start near x = 0, so "starts at the LE" alone is not enough: on that
 * test alone a two-element file converts as though it were one airfoil, its first
 * block reversed and the second appended.
 *
 * A Lednicer surface runs LE→TE and stops there, so:
 *  - it is NOT a closed contour (first point ≠ last point), and
 *  - x never steps backwards along it, and
 *  - both surfaces cover the same chord, so their x-ranges overlap almost fully.
 *
 * An element of a multi-element configuration is a closed (or near-closed) loop:
 * x doubles back on itself and the first/last points coincide. Elements are also
 * usually offset in x (a flap sits aft of the main), so their x-ranges overlap
 * little. Failing any Lednicer test means "treat the blocks as elements", which
 * is the safe direction: elements are kept separate rather than glued together.
 */
function looksLikeLednicer(groups: AirfoilPoint[][]): boolean {
  if (groups.length !== 2) return false;

  const [upper, lower] = groups;
  if (upper.length < 2 || lower.length < 2) return false;

  // Both surfaces are listed from the leading edge outwards.
  if (upper[0].x >= 0.1 || lower[0].x >= 0.1) return false;

  if (isClosedContour(upper) || isClosedContour(lower)) return false;
  if (!isMonotoneSurface(upper) || !isMonotoneSurface(lower)) return false;

  return xOverlapFraction(upper, lower) >= MIN_LEDNICER_X_OVERLAP;
}

/**
 * Convert Lednicer groups (upper LE→TE then lower LE→TE) to Selig order,
 * or return null when the groups are not Lednicer surfaces.
 */
function lednicerToSelig(groups: AirfoilPoint[][]): AirfoilPoint[] | null {
  if (!looksLikeLednicer(groups)) return null;

  const [upper, lower] = groups;

  // Reverse upper (TE→LE) and append lower (LE→TE), skipping duplicate LE point
  const reversed = [...upper].reverse();
  return [...reversed, ...lower.slice(1)];
}

/**
 * Whether honouring annotation lines as element boundaries produced blocks that
 * really do look like separate parts of the file, rather than one contour cut up
 * by its own annotations.
 */
function annotationSplitIsCredible(blocks: AirfoilPoint[][]): boolean {
  // Two Lednicer surfaces are collapsed into one element downstream, so honouring
  // an annotation between them changes nothing but keeps Lednicer detection
  // working for a file whose only surface boundary is a comment line.
  if (looksLikeLednicer(blocks)) return true;

  return blocks.every(looksLikeElementContour);
}

/**
 * Group the token stream into element blocks, deciding whether annotation lines
 * are element boundaries or decoration.
 *
 * Blank lines and `999.0` sentinels always separate elements. A comment line is
 * ambiguous: `30p-30n.dat` in the bundled library labels its three elements
 * `# Slat` / `# Main Element` / `# Flap` and uses nothing else as a boundary,
 * while a single-element file may just as reasonably carry `# upper surface`
 * annotations mid-contour. Annotations are therefore honoured as boundaries only
 * when the blocks they produce stand up as separate parts; otherwise they are
 * decoration and the blocks stay joined.
 */
function splitIntoElements(tokens: DatToken[]): AirfoilPoint[][] {
  const joined = foldBlocks(tokens, false);
  const split = foldBlocks(tokens, true);

  if (split.length > joined.length && annotationSplitIsCredible(split)) {
    return split;
  }
  return joined;
}

export function parseAirfoilDat(text: string, fileName: string): ParsedAirfoilDat {
  const { header, tokens } = tokenizeDat(text);
  const groups = splitIntoElements(tokens);

  // Lednicer's two blocks are two surfaces of one element, so they collapse to a
  // single contour. Otherwise every block is its own element, kept separate
  // rather than concatenated into one multi-loop contour.
  const fromLednicer = lednicerToSelig(groups);
  const blocks = fromLednicer ? [fromLednicer] : groups;

  if (blocks.length === 0 || blocks[0].length < 3) {
    throw new Error('Expected at least 3 airfoil coordinates.');
  }

  // Later elements must be real contours too: an explicit `999.0` separator can
  // declare a degenerate element, and keeping it would corrupt the geometry.
  blocks.slice(1).forEach((block, index) => {
    if (block.length < 3) {
      throw new Error(
        `Element ${index + 2} has only ${block.length} coordinate(s); expected at least 3.`,
      );
    }
  });

  const primary = blocks[0];
  const leIndex = findLeadingEdgeIndex(primary);
  if (leIndex === 0 || leIndex === primary.length - 1) {
    throw new Error('Expected a full airfoil loop ordered TE -> upper -> LE -> lower -> TE.');
  }

  const elements = blocks.map(annotateSurfaces);

  return {
    name: header ?? fallbackAirfoilName(fileName),
    // Single-element consumers get element 0 explicitly, not a flattened blend.
    coordinates: elements[0],
    elements,
  };
}

export function prepareImportedAirfoil(
  parsed: ParsedAirfoilFile,
  nPanels: number,
  repanelAirfoil?: RepanelAirfoil,
): ImportedAirfoilState {
  let panels = parsed.coordinates;

  if (repanelAirfoil) {
    const repaneled = repanelAirfoil(parsed.coordinates, nPanels);
    if (repaneled.length > 0) {
      panels = annotateSurfaces(repaneled);
    }
  }

  return {
    ...parsed,
    panels,
  };
}
