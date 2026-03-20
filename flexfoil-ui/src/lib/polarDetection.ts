/**
 * Intelligent polar-curve detection for the PlotBuilder.
 *
 * Groups RunRows into coherent polar curves by compound key
 * (airfoil_hash + Re + Mach + Ncrit), splits on alpha gaps,
 * and generates labels that only show distinguishing parameters.
 */

import type { RunRow } from '../types';
import { ALL_INPUT_PARAMETER_FIELDS } from './plotFields';
import { valueKey } from './plotStyling';

export interface PolarGroup {
  key: string;
  label: string;
  rows: RunRow[];
  isSinglePoint: boolean;
}

export interface DetectGroupOptions {
  sortField?: keyof RunRow;
  plottedFields?: (keyof RunRow)[];
  encodingFields?: Array<keyof RunRow | ''>;
  invariantFields?: (keyof RunRow)[];
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
  return splitByNumericGap(rows, 'alpha', gapMultiplier);
}

function median(values: number[]): number {
  if (values.length === 0) return 0;
  const s = [...values].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 === 0 ? (s[mid - 1] + s[mid]) / 2 : s[mid];
}

// ─── Smart label generation ─────────────────────────────────────────

interface GroupMeta {
  airfoil_name: string;
  alpha: number;
  reynolds: number;
  mach: number;
  ncrit: number;
  n_panels: number;
  max_iter: number;
  solver_mode: RunRow['solver_mode'];
  flap_deflection: number | null;
  flap_hinge_x: number | null;
}

/**
 * Build labels that only include parameters that actually vary across
 * the provided groups. When every group shares the same airfoil, the
 * airfoil name is omitted, etc.
 */
function buildSmartLabels(
  groups: GroupMeta[],
  excludeFromLabels: Set<string> = new Set(),
): string[] {
  if (groups.length === 0) return [];
  if (groups.length === 1) return [groups[0].airfoil_name];

  const labelFields: Array<keyof GroupMeta> = [
    'airfoil_name',
    'alpha',
    'reynolds',
    'mach',
    'ncrit',
    'flap_deflection',
    'flap_hinge_x',
    'n_panels',
    'max_iter',
    'solver_mode',
  ];
  const varyingFields = labelFields.filter((field) => {
    if (excludeFromLabels.has(field)) return false;
    return new Set(groups.map((group) => valueKey(group[field]))).size > 1;
  });

  return groups.map((group) => {
    const parts = varyingFields.map((field) => formatLabelPart(field, group[field]));
    if (parts.length === 0) return group.airfoil_name;
    if (varyingFields[0] === 'airfoil_name' && parts.length > 1) {
      const [name, ...rest] = parts;
      return `${name} @ ${rest.join(', ')}`;
    }
    return parts.join(', ');
  });
}

// ─── Main detection ─────────────────────────────────────────────────

export function splitByNumericGap(
  rows: RunRow[],
  field: keyof RunRow,
  gapMultiplier = 5,
): RunRow[][] {
  if (rows.length < 3) return [sortRows(rows, field)];

  const sorted = sortRows(rows, field);
  const steps: number[] = [];
  for (let i = 1; i < sorted.length; i++) {
    const prev = sorted[i - 1][field];
    const current = sorted[i][field];
    if (typeof prev !== 'number' || typeof current !== 'number') {
      return [sorted];
    }
    steps.push(current - prev);
  }

  const medianStep = median(steps);
  if (medianStep <= 0) return [sorted];

  const threshold = medianStep * gapMultiplier;
  const segments: RunRow[][] = [[sorted[0]]];

  for (let i = 1; i < sorted.length; i++) {
    const prev = sorted[i - 1][field];
    const current = sorted[i][field];
    if (typeof prev !== 'number' || typeof current !== 'number') {
      segments[segments.length - 1].push(sorted[i]);
      continue;
    }
    const gap = current - prev;
    if (gap > threshold) {
      segments.push([sorted[i]]);
    } else {
      segments[segments.length - 1].push(sorted[i]);
    }
  }

  return segments.length > 1 ? segments : [sorted];
}

function sortRows(rows: RunRow[], field: keyof RunRow): RunRow[] {
  return [...rows].sort((a, b) => {
    const va = a[field];
    const vb = b[field];
    if (typeof va === 'number' && typeof vb === 'number') return va - vb;
    return String(va ?? '').localeCompare(String(vb ?? ''));
  });
}

function buildInvariantKey(row: RunRow, invariantFields: (keyof RunRow)[]): string {
  return invariantFields
    .map((field) => `${String(field)}=${valueKey(row[field])}`)
    .join('|');
}

function buildGroupMeta(row: RunRow): GroupMeta {
  return {
    airfoil_name: row.airfoil_name,
    alpha: row.alpha,
    reynolds: row.reynolds,
    mach: row.mach,
    ncrit: row.ncrit,
    n_panels: row.n_panels,
    max_iter: row.max_iter,
    solver_mode: row.solver_mode,
    flap_deflection: row.flap_deflection,
    flap_hinge_x: row.flap_hinge_x,
  };
}

function formatLabelPart(field: keyof GroupMeta, value: GroupMeta[keyof GroupMeta]): string {
  switch (field) {
    case 'airfoil_name':
      return String(value);
    case 'alpha':
      return `α=${typeof value === 'number' ? value : String(value)}`;
    case 'reynolds':
      return `Re=${typeof value === 'number' ? formatReynolds(value) : String(value)}`;
    case 'mach':
      return `M=${String(value)}`;
    case 'ncrit':
      return `Nc=${String(value)}`;
    case 'flap_deflection':
      return value != null ? `δ=${value}°` : '';
    case 'flap_hinge_x':
      return value != null ? `x/c=${value}` : '';
    case 'n_panels':
      return `Panels=${String(value)}`;
    case 'max_iter':
      return `Iter=${String(value)}`;
    case 'solver_mode':
      return value === 'inviscid' ? 'Inviscid' : 'Viscous';
    default:
      return String(value);
  }
}

/**
 * Detect polar groups from a flat array of RunRows.
 *
 * Core principle: a group is "runs where every input parameter is the same
 * except the one on the X-axis (the swept variable)."
 *
 * 1. Take all input parameters, remove the X-axis field (sortField) and
 *    any explicitly plotted / encoding fields.
 * 2. Group by the remaining invariant fields.
 * 3. Optionally split on large gaps in the sort field (only for larger groups).
 * 4. Generate labels from only the fields that vary across resulting groups.
 */
export function detectSmartRunGroups(
  rows: RunRow[],
  options: DetectGroupOptions = {},
): PolarGroup[] {
  if (rows.length === 0) return [];

  const {
    sortField = 'alpha',
    plottedFields = [],
    encodingFields = [],
    invariantFields,
  } = options;

  // Build grouping fields: all input parameters minus the swept variable
  const baseFields = invariantFields ?? ALL_INPUT_PARAMETER_FIELDS;
  const excludedFields = new Set<keyof RunRow>([
    sortField,
    ...plottedFields,
    ...encodingFields.filter((field): field is keyof RunRow => Boolean(field)),
  ]);
  const groupFields = baseFields.filter((field) => !excludedFields.has(field));

  const compoundGroups = new Map<string, RunRow[]>();
  for (const row of rows) {
    const key = buildInvariantKey(row, groupFields);
    if (!compoundGroups.has(key)) compoundGroups.set(key, []);
    compoundGroups.get(key)!.push(row);
  }

  const allSegments: { key: string; meta: GroupMeta; rows: RunRow[]; subIndex: number; totalSubs: number }[] = [];

  for (const [key, groupRows] of [...compoundGroups.entries()].sort(([a], [b]) => a.localeCompare(b))) {
    const sorted = sortRows(groupRows, sortField);
    const meta = buildGroupMeta(sorted[0]);

    // Only gap-split larger groups (small ones are likely a single intentional sweep)
    // and only when the sort field is numeric (gap detection on categorical is meaningless)
    const shouldSplit = sorted.length >= 6
      && typeof sorted[0]?.[sortField] === 'number';
    const segments = shouldSplit
      ? splitByNumericGap(sorted, sortField)
      : [sorted];

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

  const labels = buildSmartLabels(
    allSegments.map(s => s.meta),
    new Set([...(plottedFields as string[]), String(sortField)]),
  );

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

export function detectPolarGroups(
  rows: RunRow[],
  sortField: keyof RunRow = 'alpha',
): PolarGroup[] {
  return detectSmartRunGroups(rows, { sortField, plottedFields: [sortField] });
}
