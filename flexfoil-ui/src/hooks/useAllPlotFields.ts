import { useMemo } from 'react';
import {
  PLOT_FIELDS,
  NUMERIC_PLOT_FIELDS,
  ENCODING_PLOT_FIELDS,
  type PlotFieldMeta,
} from '../lib/plotFields';
import { useCustomColumnStore, compileCustomColumns } from '../stores/customColumnStore';
import type { RunRow } from '../types';

export interface AllPlotFields {
  allFields: PlotFieldMeta[];
  numericFields: PlotFieldMeta[];
  encodingFields: PlotFieldMeta[];
  customEvaluators: Map<string, (row: RunRow) => number | null>;
  getFieldLabel: (key: string) => string;
  getFieldValue: (row: RunRow, key: string) => unknown;
}

/**
 * Returns the merged list of built-in + custom plot fields,
 * plus evaluator functions for custom columns.
 */
export function useAllPlotFields(): AllPlotFields {
  const customColumns = useCustomColumnStore((s) => s.columns);

  const customEvaluators = useMemo(
    () => compileCustomColumns(customColumns),
    [customColumns],
  );

  const customMetas: PlotFieldMeta[] = useMemo(
    () =>
      customColumns.map((col) => ({
        key: `custom:${col.id}` as keyof RunRow,
        label: col.name,
        kind: 'numeric' as const,
      })),
    [customColumns],
  );

  const allFields = useMemo(
    () => [...PLOT_FIELDS, ...customMetas],
    [customMetas],
  );

  const numericFields = useMemo(
    () => [...NUMERIC_PLOT_FIELDS, ...customMetas],
    [customMetas],
  );

  const encodingFields = useMemo(
    () => [...ENCODING_PLOT_FIELDS, ...customMetas],
    [customMetas],
  );

  const getFieldLabel = useMemo(() => {
    const labelMap = new Map<string, string>();
    for (const f of allFields) labelMap.set(f.key as string, f.label);
    return (key: string) => labelMap.get(key) ?? key;
  }, [allFields]);

  const getFieldValue = useMemo(() => {
    return (row: RunRow, key: string): unknown => {
      if (key.startsWith('custom:')) {
        const colId = key.slice('custom:'.length);
        const evalFn = customEvaluators.get(colId);
        return evalFn ? evalFn(row) : null;
      }
      return row[key as keyof RunRow];
    };
  }, [customEvaluators]);

  return {
    allFields,
    numericFields,
    encodingFields,
    customEvaluators,
    getFieldLabel,
    getFieldValue,
  };
}
