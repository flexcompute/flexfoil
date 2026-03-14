/**
 * Deterministic SHA-256 hash of airfoil panel coordinates.
 * Used as the geometry component of the solver cache key.
 */

const COORDINATE_PRECISION = 8;

export async function computeAirfoilHash(
  panels: { x: number; y: number }[]
): Promise<string> {
  const canonical = panels
    .map(p => `${p.x.toFixed(COORDINATE_PRECISION)},${p.y.toFixed(COORDINATE_PRECISION)}`)
    .join(';');

  const data = new TextEncoder().encode(canonical);
  const hashBuffer = await crypto.subtle.digest('SHA-256', data);
  const hashArray = Array.from(new Uint8Array(hashBuffer));
  return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
}
