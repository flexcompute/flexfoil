/**
 * Intelligent polar-curve detection for the PlotBuilder.
 *
 * Groups RunRows into coherent polar curves by compound key
 * (airfoil_hash + Re + Mach + Ncrit), splits on alpha gaps,
 * and generates labels that only show distinguishing parameters.
 */

import type { RunRow } from '../types';

export interface PolarGroup {
  key: string;
  label: string;
  rows: RunRow[];
  isSinglePoint: boolean;
}

// ─── Reynolds formatting ────────────────────────────────────────────

export function formatReynolds(re: number): string {
  if (re >= 1e6 && re % 1e6 === 0) return `${re / 1e6}M`;
  if (re >= 1e6) return `${+(re / 1e6).toFixed(2)}M`;
  if (re >= 1e3 && re % 1e3 === 0) return `${re / 1e3}k`;
  if (re >= 1e3) return `${+(re / 1e3).toFixed(1)}k`;
  return String(re);
}

// ─── Alpha-gap splitting ────────────────────────────────────────────

/**
 * Split a sorted array of rows into sub-arrays wherever the gap between
 * consecutive alpha values exceeds `gapMultiplier` times the median step.
 * Returns the original array wrapped in a single-element array if no gaps
 * are detected or the array has fewer than 3 points.
 */
export function splitByAlphaGap(
  rows: RunRow[],
  gapMultiplier = 3,
): RunRow[][] {
  if (rows.length < 3) return [rows];

  const sorted = [...rows].sort((a, b) => a.alpha - b.alpha);
  const steps: number[] = [];
  for (let i = 1; i < sorted.length; i++) {
    steps.push(sorted[i].alpha - sorted[i - 1].alpha);
  }

  const medianStep = median(steps);
  if (medianStep <= 0) return [sorted];

  const threshold = medianStep * gapMultiplier;
  const segments: RunRow[][] = [[sorted[0]]];

  for (let i = 1; i < sorted.length; i++) {
    const gap = sorted[i].alpha - sorted[i - 1].alpha;
    if (gap > threshold) {
      segments.push([sorted[i]]);
    } else {
      segments[segments.length - 1].push(sorted[i]);
    }
  }

  return segments.length > 1 ? segments : [sorted];
}

function median(values: number[]): number {
  if (values.length === 0) return 0;
  const s = [...values].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 === 0 ? (s[mid - 1] + s[mid]) / 2 : s[mid];
}

// ─── Smart label generation ─────────────────────────────────────────

interface GroupMeta {
  airfoilName: string;
  reynolds: number;
  mach: number;
  ncrit: number;
}

/**
 * Build labels that only include parameters that actually vary across
 * the provided groups. When every group shares the same airfoil, the
 * airfoil name is omitted, etc.
 */
function buildSmartLabels(groups: GroupMeta[]): string[] {
  if (groups.length === 0) return [];
  if (groups.length === 1) return [groups[0].airfoilName];

  const uniqueAirfoils = new Set(groups.map(g => g.airfoilName));
  const uniqueRe = new Set(groups.map(g => g.reynolds));
  const uniqueMach = new Set(groups.map(g => g.mach));
  const uniqueNcrit = new Set(groups.map(g => g.ncrit));

  const showAirfoil = uniqueAirfoils.size > 1;
  const showRe = uniqueRe.size > 1;
  const showMach = uniqueMach.size > 1;
  const showNcrit = uniqueNcrit.size > 1;

  // If nothing varies (identical groups), fall back to airfoil name
  if (!showAirfoil && !showRe && !showMach && !showNcrit) {
    return groups.map(g => g.airfoilName);
  }

  return groups.map(g => {
    const parts: string[] = [];
    if (showAirfoil) parts.push(g.airfoilName);
    if (showRe) parts.push(`Re=${formatReynolds(g.reynolds)}`);
    if (showMach) parts.push(`M=${g.mach}`);
    if (showNcrit) parts.push(`Nc=${g.ncrit}`);

    if (parts.length === 0) return g.airfoilName;
    if (showAirfoil && parts.length > 1) {
      const name = parts.shift()!;
      return `${name} @ ${parts.join(', ')}`;
    }
    return parts.join(', ');
  });
}

// ─── Main detection ─────────────────────────────────────────────────

function polarCompoundKey(row: RunRow): string {
  return `${row.airfoil_hash}|${row.reynolds}|${row.mach}|${row.ncrit}`;
}

/**
 * Detect polar groups from a flat array of RunRows.
 *
 * 1. Group by compound key (airfoil_hash, Re, Mach, Ncrit)
 * 2. Within each group, sort by `sortField` (default: alpha)
 * 3. Split on alpha discontinuities
 * 4. Generate smart labels
 */
export function detectPolarGroups(
  rows: RunRow[],
  sortField: keyof RunRow = 'alpha',
): PolarGroup[] {
  if (rows.length === 0) return [];

  // Step 1: group by compound key
  const compoundGroups = new Map<string, RunRow[]>();
  for (const row of rows) {
    const key = polarCompoundKey(row);
    if (!compoundGroups.has(key)) compoundGroups.set(key, []);
    compoundGroups.get(key)!.push(row);
  }

  // Step 2: sort each group, split on alpha gaps
  const allSegments: { key: string; meta: GroupMeta; rows: RunRow[]; subIndex: number; totalSubs: number }[] = [];

  for (const [key, groupRows] of compoundGroups) {
    const sorted = [...groupRows].sort((a, b) => {
      const va = a[sortField];
      const vb = b[sortField];
      if (typeof va === 'number' && typeof vb === 'number') return va - vb;
      return 0;
    });

    const segments = splitByAlphaGap(sorted);
    const meta: GroupMeta = {
      airfoilName: sorted[0].airfoil_name,
      reynolds: sorted[0].reynolds,
      mach: sorted[0].mach,
      ncrit: sorted[0].ncrit,
    };

    for (let i = 0; i < segments.length; i++) {
      allSegments.push({
        key: segments.length > 1 ? `${key}__seg${i}` : key,
        meta,
        rows: segments[i],
        subIndex: i,
        totalSubs: segments.length,
      });
    }
  }

  // Step 3: build smart labels
  const labels = buildSmartLabels(allSegments.map(s => s.meta));

  return allSegments.map((seg, i) => {
    let label = labels[i];
    if (seg.totalSubs > 1) {
      label += ` (${seg.subIndex + 1})`;
    }
    return {
      key: seg.key,
      label,
      rows: seg.rows,
      isSinglePoint: seg.rows.length === 1,
    };
  });
}
