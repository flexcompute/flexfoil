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
 * Canonical form:
 *
 *     element_0 [# element_1 [# element_2 ...]]
 *     element_i = <coords_i>[@<placement_i>]
 *     coords_i  = "x,y" per node, 8 decimal places, joined with ";"
 *     placement = rot,pivot_x,pivot_y,trans_x,trans_y,scale  (8 decimals each)
 *
 * The placement suffix is emitted only when the placement is present *and*
 * non-identity, and the element separator only appears between elements, so a
 * single element with no (or identity) placement canonicalises to exactly the
 * coordinate string this module has always produced — every run already cached
 * in a user's browser database stays valid.
 *
 * NOTE (known divergence, pre-existing): the Python implementation in
 * packages/flexfoil-python/src/flexfoil/airfoil.py joins nodes with "|" and
 * truncates the digest to 16 hex chars, so Python-side and browser-side hashes
 * for the *same* coordinates have never matched. The element/placement
 * encoding below is byte-for-byte identical to Python's, and both sides pin
 * the exact shared suffix in their tests; unifying the coordinate half is a
 * cache-invalidating change and belongs with the database migration.
 */

const COORDINATE_PRECISION = 8;
const ELEMENT_SEP = '#';
const PLACEMENT_PREFIX = '@';

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

/** One element of an assembly: its own panel nodes plus where it sits. */
export interface GeometryElement {
  panels: HashPoint[];
  placement?: Placement | null;
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
 * Format a placement scalar. `-0` collapses to `0` so that JavaScript and
 * Python (whose `%.8f` would otherwise emit "-0.00000000") agree.
 */
function fmt(value: number): string {
  return (value === 0 ? 0 : value).toFixed(COORDINATE_PRECISION);
}

function canonicalCoords(panels: HashPoint[]): string {
  return panels
    .map(p => `${p.x.toFixed(COORDINATE_PRECISION)},${p.y.toFixed(COORDINATE_PRECISION)}`)
    .join(';');
}

function canonicalPlacement(placement?: Placement | null): string {
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
      .map(fmt)
      .join(',')
  );
}

/** Canonical string for an assembly (see the module notes above). */
export function canonicalGeometryString(elements: GeometryElement[]): string {
  return elements
    .map(e => canonicalCoords(e.panels) + canonicalPlacement(e.placement))
    .join(ELEMENT_SEP);
}

async function sha256Hex(canonical: string): Promise<string> {
  const data = new TextEncoder().encode(canonical);
  const hashBuffer = await crypto.subtle.digest('SHA-256', data);
  const hashArray = Array.from(new Uint8Array(hashBuffer));
  return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
}

/**
 * Placement-aware geometry hash. Reduces to the historical coordinate-only
 * digest for a single element with no (or identity) placement.
 */
export async function computeGeometryHash(elements: GeometryElement[]): Promise<string> {
  return sha256Hex(canonicalGeometryString(elements));
}

export async function computeAirfoilHash(
  panels: HashPoint[],
  placement?: Placement | null
): Promise<string> {
  return computeGeometryHash([{ panels, placement }]);
}
