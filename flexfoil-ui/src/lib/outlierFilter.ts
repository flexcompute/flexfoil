/**
 * IQR-based outlier filtering for polar data and run rows.
 *
 * Uses the standard 1.5×IQR fence:
 *   lower = Q1 − 1.5·IQR
 *   upper = Q3 + 1.5·IQR
 *
 * Points outside the fence on *any* active axis are excluded.
 */

export interface Fence {
  lower: number;
  upper: number;
}

function quartiles(sorted: number[]): { q1: number; q3: number } {
  const n = sorted.length;
  if (n < 4) return { q1: sorted[0], q3: sorted[n - 1] };
  const q1Idx = (n - 1) * 0.25;
  const q3Idx = (n - 1) * 0.75;
  const q1 =
    sorted[Math.floor(q1Idx)] +
    (q1Idx % 1) * (sorted[Math.ceil(q1Idx)] - sorted[Math.floor(q1Idx)]);
  const q3 =
    sorted[Math.floor(q3Idx)] +
    (q3Idx % 1) * (sorted[Math.ceil(q3Idx)] - sorted[Math.floor(q3Idx)]);
  return { q1, q3 };
}

/** Build a single IQR fence from raw numeric values. Returns null if
 *  there aren't enough values or the IQR is zero (constant data). */
export function buildFence(values: number[]): Fence | null {
  const finite = values.filter(Number.isFinite);
  if (finite.length < 4) return null;
  finite.sort((a, b) => a - b);
  const { q1, q3 } = quartiles(finite);
  const iqr = q3 - q1;
  if (iqr === 0) return null;
  return { lower: q1 - 1.5 * iqr, upper: q3 + 1.5 * iqr };
}

/** Check whether a single value sits inside a fence. Null/undefined/NaN
 *  values and null fences are treated as "pass" (not excluded). */
export function isInlier(value: number | null | undefined, fence: Fence | null): boolean {
  if (fence === null) return true;
  if (value == null || !Number.isFinite(value)) return true;
  return value >= fence.lower && value <= fence.upper;
}

/**
 * Generic outlier filter: keep items where every extracted value
 * sits within the IQR fence for its corresponding axis.
 */
export function filterOutliers<T>(
  items: T[],
  extractors: Array<(item: T) => number | null | undefined>,
): T[] {
  if (items.length < 4) return items;

  const fences = extractors.map((extract) => {
    const values = items.map(extract).filter((v): v is number => v != null && Number.isFinite(v));
    return buildFence(values);
  });

  if (fences.every((f) => f === null)) return items;

  return items.filter((item) =>
    extractors.every((extract, i) => isInlier(extract(item), fences[i])),
  );
}
