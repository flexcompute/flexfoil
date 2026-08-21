/**
 * Deterministic SHA-256 hash of assembly geometry.
 * Used as the geometry component of the solver cache key.
 *
 * The hash must distinguish every geometry that can produce a different flow
 * solution. For a multi-element assembly that means *placement* as well as
 * element shape: two configurations built from the same element shapes but
 * with a different gap/overlap/deflection are entirely different aerodynamic
 * problems, and hashing coordinates alone would serve one's cached results for
 * the other.
 *
 * Two distinct digests live here, and which one you get matters:
 *
 * 1. LEGACY SINGLE-ELEMENT DIGEST — one element, no (or identity) placement.
 *    This is what this module has always returned, so every run already cached
 *    in a user's browser database is keyed on it and its encoding is frozen:
 *    nodes joined with ";", all 64 hex characters of the digest. It is
 *    deliberately *not* comparable with the Python digest for the same
 *    coordinates — packages/flexfoil-python/src/flexfoil/airfoil.py froze its
 *    own legacy encoding ("|" between nodes, truncated to 16 hex characters)
 *    against the runs already in users' ~/.flexfoil/runs.db. Unifying the two
 *    invalidates one of the caches, so it belongs with the database-migration
 *    workstream.
 *
 * 2. SHARED ASSEMBLY DIGEST — everything else: more than one element, or an
 *    element carrying a non-identity placement. Nothing is cached against this
 *    yet, so it is defined to be byte-for-byte identical to the Python
 *    implementation: same schema tag, separators, field order, float
 *    formatting (see `fmtShared`) and no truncation. The shared test vectors
 *    in airfoilHash.test.ts pin the exact digests, and the identical literals
 *    are pinned in packages/flexfoil-python/tests/test_geometry_hash.py, so a
 *    change to one side alone fails CI.
 *
 * Canonical forms:
 *
 *     legacy = "x,y" per node, 8 decimals, joined with ";"
 *     shared = "ffgeom1:" + element_0 ["#" element_1 ["#" element_2 ...]]
 *              element_i = <coords_i>["@" <placement_i>]
 *              coords_i  = "x,y" per node, 8 decimals, joined with ";"
 *              placement = rot,pivot_x,pivot_y,trans_x,trans_y,scale
 *
 * `canonicalGeometryString` and `computeGeometryHash` select between the two on
 * the shape of their input (see `isLegacySingleElement`). The shared form
 * carries a schema tag, so the two canonical strings can never coincide and one
 * geometry always has exactly one digest.
 */

// Legacy encoding. Changing this invalidates every run already cached against
// the old digest; the pinned test in airfoilHash.test.ts exists to catch that.
// (It happens to be the same character as SHARED_NODE_SEP below; Python's two
// separators differ, which is why both are named here.)
const LEGACY_NODE_SEP = ';';

// Shared encoding. Nothing is cached against these yet, but changing one
// without the matching change in
// packages/flexfoil-python/src/flexfoil/airfoil.py splits the two
// implementations' cache keys; bump the schema tag if the form has to change.
// Both sides' pinned shared vectors exist to catch that.
const SHARED_SCHEMA = 'ffgeom1:';
const SHARED_NODE_SEP = ';';
const ELEMENT_SEP = '#';
const PLACEMENT_PREFIX = '@';

// Used by both encodings, so a change here invalidates the cache *and* has to
// be mirrored on the Python side.
const PRECISION = 8;

export interface HashPoint {
  x: number;
  y: number;
}

/**
 * Rigid placement of one element inside a multi-element assembly.
 * `rotation` is in degrees about `pivot` (in the element's own coordinates),
 * applied before `translation`. `scale` is uniform. All fields default to the
 * identity, which hashes as if the placement were absent.
 */
export interface Placement {
  pivot?: HashPoint;
  rotation?: number;
  translation?: HashPoint;
  scale?: number;
}

/**
 * One element of an assembly: its own panel nodes plus where it sits.
 *
 * `placement` is a required property on purpose. Inside an assembly the
 * placement is part of the aerodynamic problem, so an element built without one
 * would key its results under the wrong geometry; requiring the property means
 * a caller has to write `placement: null` deliberately for an element that
 * really does sit in its own coordinates.
 */
export interface GeometryElement {
  panels: HashPoint[];
  placement: Placement | null;
}

/** True if the placement moves nothing (a pivot on its own is a no-op). */
export function isIdentityPlacement(placement?: Placement | null): boolean {
  if (!placement) return true;
  return (
    (placement.rotation ?? 0) === 0 &&
    (placement.translation?.x ?? 0) === 0 &&
    (placement.translation?.y ?? 0) === 0 &&
    (placement.scale ?? 1) === 1
  );
}

/**
 * Format one scalar for the shared canonical form. `toFixed` is the reference
 * behaviour for both languages here; Python's `_shared_fmt` reproduces it.
 *
 * `-0` collapses to `0` because `toFixed` strips the sign before rounding, so
 * JavaScript never emits "-0.00000000" for a true negative zero while Python's
 * `%.8f` does. Non-zero negatives keep their sign on both sides
 * (`(-1e-12).toFixed(8)` is "-0.00000000" too), so only exact zero is
 * normalised.
 */
function fmtShared(value: number): string {
  return (value === 0 ? 0 : value).toFixed(PRECISION);
}

/**
 * Frozen legacy coordinate encoding: raw `toFixed`, ";" between nodes.
 *
 * Deliberately not routed through `fmtShared` — the raw format is what the
 * already-cached digests were produced with. One consequence is that a
 * coordinate of exactly -0 renders "0.00000000" here (`toFixed` drops the sign)
 * while Python's frozen legacy encoding renders "-0.00000000" for the same
 * input. The two legacy encodings were frozen independently and are not
 * comparable in any case (different node separator, different digest length);
 * the shared form below normalises signed zero on both sides.
 */
function legacyCanonicalCoords(panels: HashPoint[]): string {
  return panels
    .map(p => `${p.x.toFixed(PRECISION)},${p.y.toFixed(PRECISION)}`)
    .join(LEGACY_NODE_SEP);
}

function sharedCanonicalCoords(panels: HashPoint[]): string {
  return panels.map(p => `${fmtShared(p.x)},${fmtShared(p.y)}`).join(SHARED_NODE_SEP);
}

function sharedCanonicalPlacement(placement: Placement | null): string {
  if (isIdentityPlacement(placement)) return '';
  const p = placement as Placement;
  return (
    PLACEMENT_PREFIX +
    [
      p.rotation ?? 0,
      p.pivot?.x ?? 0,
      p.pivot?.y ?? 0,
      p.translation?.x ?? 0,
      p.translation?.y ?? 0,
      p.scale ?? 1,
    ]
      .map(v => fmtShared(v))
      .join(',')
  );
}

/**
 * True for the one input shape the frozen legacy digest covers: exactly one
 * element, with no placement or an identity placement.
 */
export function isLegacySingleElement(elements: GeometryElement[]): boolean {
  return elements.length === 1 && isIdentityPlacement(elements[0].placement);
}

/**
 * Cross-language canonical string for an assembly. Byte-for-byte identical to
 * `shared_canonical_geometry` in
 * packages/flexfoil-python/src/flexfoil/airfoil.py.
 */
export function sharedCanonicalGeometryString(elements: GeometryElement[]): string {
  return (
    SHARED_SCHEMA +
    elements
      .map(e => sharedCanonicalCoords(e.panels) + sharedCanonicalPlacement(e.placement))
      .join(ELEMENT_SEP)
  );
}

/**
 * Canonical string for an assembly: the frozen legacy coordinate string for a
 * lone unplaced element, the shared form for everything else.
 */
export function canonicalGeometryString(elements: GeometryElement[]): string {
  if (isLegacySingleElement(elements)) return legacyCanonicalCoords(elements[0].panels);
  return sharedCanonicalGeometryString(elements);
}

async function sha256Hex(canonical: string): Promise<string> {
  const data = new TextEncoder().encode(canonical);
  const hashBuffer = await crypto.subtle.digest('SHA-256', data);
  const hashArray = Array.from(new Uint8Array(hashBuffer));
  return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
}

/**
 * Geometry half of the run-cache key. Returns the frozen legacy digest for a
 * single element with no (or identity) placement — the shape every
 * already-cached run has — and the cross-language shared digest for everything
 * else. Both are 64 hex characters on this side; Python truncates only its
 * legacy digest.
 */
export async function computeGeometryHash(elements: GeometryElement[]): Promise<string> {
  return sha256Hex(canonicalGeometryString(elements));
}

/**
 * Single-element convenience wrapper.
 *
 * Omitting `placement` hashes the element as sitting in its own coordinates
 * and yields the frozen legacy digest, which is what a lone airfoil wants. The
 * current callers — `hashPanels` in src/stores/runStore.ts and the three call
 * sites in src/lib/sweepEngine.ts — all hash lone airfoils and pass nothing.
 * Code that hashes an element *of an assembly* must go through
 * `computeGeometryHash`, whose `GeometryElement.placement` is a required
 * property, so the placement has to be stated rather than defaulted.
 */
export async function computeAirfoilHash(
  panels: HashPoint[],
  placement?: Placement | null
): Promise<string> {
  return computeGeometryHash([{ panels, placement: placement ?? null }]);
}
