/**
 * Custom AG Grid aggregation functions for aerodynamic sweep analysis.
 *
 * Provides argmax/argmin factories that return the value of one column
 * at the row where another column reaches its extremum within a group.
 *
 * Example: alpha_stall = argmax(cl) applied to the alpha column
 *          → returns the alpha value at the row where CL is maximized.
 *
 * All functions respect AG Grid column filters (Converged Only, Exclude
 * Outliers, etc.) by traversing only childrenAfterFilter, so manually
 * flagged outliers and filtered-out rows are excluded from aggregation.
 */

import type { IAggFunc, IAggFuncParams, IRowNode } from 'ag-grid-community';
import type { RunRow } from '../types';

type NumericKey = {
  [K in keyof RunRow]: RunRow[K] extends number | null ? K : never;
}[keyof RunRow];

/**
 * Collect all filtered leaf data rows beneath a group node.
 * Uses childrenAfterFilter so grid column filters (converged, is_outlier,
 * etc.) are respected — filtered-out rows never contribute to aggregation.
 */
function getFilteredLeaves(node: IRowNode<RunRow>): RunRow[] {
  if (!node.group && node.data) return [node.data];

  const result: RunRow[] = [];
  const children = node.childrenAfterFilter;
  if (!children) return result;
  for (const child of children) {
    if (child.group) {
      result.push(...getFilteredLeaves(child));
    } else if (child.data) {
      result.push(child.data);
    }
  }
  return result;
}

/**
 * Create an agg function that returns the value of the *current* column
 * at the filtered leaf row where `optimizeField` is maximized.
 */
export function makeArgmax(optimizeField: NumericKey): IAggFunc<RunRow> {
  return (params: IAggFuncParams<RunRow>) => {
    const leaves = getFilteredLeaves(params.rowNode);
    if (leaves.length === 0) return null;

    const colField = params.colDef.field as keyof RunRow | undefined;
    let bestVal = -Infinity;
    let bestResult: number | null = null;

    for (const data of leaves) {
      const v = data[optimizeField];
      if (typeof v === 'number' && Number.isFinite(v) && v > bestVal) {
        bestVal = v;
        bestResult = colField ? (data[colField] as number | null) : null;
      }
    }

    return bestResult;
  };
}

/**
 * Create an agg function that returns the value of the *current* column
 * at the filtered leaf row where `optimizeField` is minimized.
 */
export function makeArgmin(optimizeField: NumericKey): IAggFunc<RunRow> {
  return (params: IAggFuncParams<RunRow>) => {
    const leaves = getFilteredLeaves(params.rowNode);
    if (leaves.length === 0) return null;

    const colField = params.colDef.field as keyof RunRow | undefined;
    let bestVal = Infinity;
    let bestResult: number | null = null;

    for (const data of leaves) {
      const v = data[optimizeField];
      if (typeof v === 'number' && Number.isFinite(v) && v < bestVal) {
        bestVal = v;
        bestResult = colField ? (data[colField] as number | null) : null;
      }
    }

    return bestResult;
  };
}

/** Median aggregation (not built into AG Grid). */
export const medianAgg: IAggFunc<RunRow> = (params: IAggFuncParams<RunRow>) => {
  const nums = params.values.filter(
    (v): v is number => typeof v === 'number' && Number.isFinite(v),
  );
  if (nums.length === 0) return null;
  nums.sort((a, b) => a - b);
  const mid = Math.floor(nums.length / 2);
  return nums.length % 2 === 0 ? (nums[mid - 1] + nums[mid]) / 2 : nums[mid];
};

/**
 * All custom aggregation functions to register on the AG Grid instance.
 *
 * Names are chosen to be user-readable in the column sidebar dropdown.
 */
export const CUSTOM_AGG_FUNCS: Record<string, IAggFunc<RunRow>> = {
  'at max(CL)': makeArgmax('cl'),
  'at min(CD)': makeArgmin('cd'),
  'at max(L/D)': makeArgmax('ld'),
  'at min(CL)': makeArgmin('cl'),
  'at max(CD)': makeArgmax('cd'),
  'at min(L/D)': makeArgmin('ld'),
  median: medianAgg,
};

/**
 * The full set of agg function names available for numeric columns,
 * including both AG Grid built-ins and our custom ones.
 */
export const ALL_AGG_FUNC_NAMES: string[] = [
  'sum',
  'min',
  'max',
  'count',
  'avg',
  'first',
  'last',
  ...Object.keys(CUSTOM_AGG_FUNCS),
];
