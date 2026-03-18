/**
 * DataExplorerPanel — AG Grid Enterprise table + Plotly correlogram.
 *
 * Two views toggled by tabs at the top:
 *   - Table: AG Grid with enterprise features (Thread-style)
 *   - Correlogram: Plotly SPLOM scatter-plot matrix (Thread-style)
 *
 * Filter state flows from the grid to the run store so the
 * PlotBuilder can consume the filtered subset.
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { AgGridReact } from 'ag-grid-react';
import {
  themeQuartz,
  colorSchemeDark,
  type ColDef,
  type GridReadyEvent,
  type FilterChangedEvent,
  type CellValueChangedEvent,
  type GridApi,
  type CellSelectionOptions,
} from 'ag-grid-community';
import Plot from 'react-plotly.js';
import { useRunStore } from '../../stores/runStore';
import { AUTO_GROUP_KEY, useRouteUiStore } from '../../stores/routeUiStore';
import { useTheme } from '../../contexts/ThemeContext';
import { useOnboarding } from '../../onboarding/useOnboarding';
import type { DataSource, RunRow } from '../../types';
import {
  ENCODING_PLOT_FIELDS,
  NUMERIC_PLOT_FIELDS,
  getPlotFieldLabel,
  isNumericPlotField,
} from '../../lib/plotFields';
import { detectSmartRunGroups } from '../../lib/polarDetection';
import { buildMarkerEncoding, colorForKey } from '../../lib/plotStyling';
import { buildFence, isInlier, type Fence } from '../../lib/outlierFilter';
import { useCustomColumnStore, compileCustomColumns } from '../../stores/customColumnStore';
import { useAllPlotFields } from '../../hooks/useAllPlotFields';
import { CustomColumnEditor } from '../CustomColumnEditor';

// ────────────────────────────────────────────────────────────
// Helpers
// ────────────────────────────────────────────────────────────

function pearsonR(xs: number[], ys: number[]): number | null {
  const n = Math.min(xs.length, ys.length);
  if (n < 2) return null;
  let sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;
  for (let i = 0; i < n; i++) {
    sx += xs[i]; sy += ys[i];
    sxy += xs[i] * ys[i];
    sx2 += xs[i] * xs[i];
    sy2 += ys[i] * ys[i];
  }
  const denom = Math.sqrt((n * sx2 - sx * sx) * (n * sy2 - sy * sy));
  if (denom === 0) return null;
  return (n * sxy - sx * sy) / denom;
}

function getLabel(key: keyof RunRow): string {
  return getPlotFieldLabel(key);
}

function axisRef(prefix: 'x' | 'y', index: number): string {
  return index === 1 ? prefix : `${prefix}${index}`;
}

function axisLayoutKey(prefix: 'x' | 'y', index: number): string {
  return index === 1 ? `${prefix}axis` : `${prefix}axis${index}`;
}

function computeNumericRange(
  rows: RunRow[],
  key: keyof RunRow,
  getter?: (row: RunRow, key: string) => unknown,
): [number, number] {
  const values = rows
    .map((row) => getter ? getter(row, key as string) : row[key])
    .filter((value): value is number => typeof value === 'number' && Number.isFinite(value));
  if (values.length === 0) return [0, 1];
  const min = Math.min(...values);
  const max = Math.max(...values);
  const pad = (max - min) * 0.05 || 0.5;
  return [min - pad, max + pad];
}

function buildRunHoverDetails(run: RunRow): Array<string | number> {
  return [
    run.airfoil_name,
    run.alpha,
    run.reynolds,
    run.mach,
    run.ncrit,
    run.n_panels,
    run.max_iter,
    run.solver_mode,
    run.session_id ?? 'n/a',
    run.created_at,
  ];
}

const RUN_HOVER_TEMPLATE = [
  'Run #%{meta}',
  'Airfoil=%{customdata[0]}',
  'Alpha=%{customdata[1]}',
  'Re=%{customdata[2]}',
  'Mach=%{customdata[3]}',
  'Ncrit=%{customdata[4]}',
  'Panels=%{customdata[5]}',
  'MaxIter=%{customdata[6]}',
  'Solver=%{customdata[7]}',
  'Session=%{customdata[8]}',
  'Created=%{customdata[9]}',
].join('<br>');

// ────────────────────────────────────────────────────────────
// AG Grid column definitions
// ────────────────────────────────────────────────────────────

function buildColumnDefs(): ColDef<RunRow>[] {
  return [
    { field: 'id', headerName: 'ID', width: 70, hide: true },
    { field: 'airfoil_name', headerName: 'Airfoil', pinned: 'left', width: 140, editable: true, enableRowGroup: true, chartDataType: 'category' as const },
    { field: 'alpha', headerName: 'α (°)', width: 90, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(2) ?? '' },
    { field: 'reynolds', headerName: 'Re', width: 120, chartDataType: 'series' as const, valueFormatter: p => p.value != null ? p.value.toExponential(2) : '' },
    { field: 'mach', headerName: 'Mach', width: 80, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(3) ?? '' },
    { field: 'ncrit', headerName: 'Ncrit', width: 80, chartDataType: 'series' as const },
    { field: 'cl', headerName: 'CL', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(5) ?? '' },
    { field: 'cd', headerName: 'CD', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(6) ?? '' },
    { field: 'cm', headerName: 'CM', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(5) ?? '' },
    { field: 'ld', headerName: 'L/D', width: 100, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(3) ?? '' },
    { field: 'converged', headerName: 'Converged', width: 100, chartDataType: 'category' as const, enableRowGroup: true,
      valueFormatter: p => p.value ? 'Yes' : 'No',
      cellStyle: p => ({ color: p.value ? 'var(--accent-primary)' : 'var(--accent-danger)' }),
    },
    { field: 'x_tr_upper', headerName: 'Xtr Upper', width: 100, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(4) ?? '' },
    { field: 'x_tr_lower', headerName: 'Xtr Lower', width: 100, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(4) ?? '' },
    { field: 'iterations', headerName: 'Iter', width: 70, chartDataType: 'series' as const },
    { field: 'residual', headerName: 'Residual', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toExponential(2) ?? '' },
    { field: 'n_panels', headerName: 'Panels', width: 85, chartDataType: 'series' as const },
    { field: 'max_iter', headerName: 'MaxIter', width: 85, hide: true, chartDataType: 'series' as const },
    { field: 'solver_mode', headerName: 'Solver', width: 90, chartDataType: 'category' as const, enableRowGroup: true },
    { field: 'session_id', headerName: 'Session', width: 130, enableRowGroup: true, hide: true, chartDataType: 'category' as const },
    { field: 'created_at', headerName: 'Time', width: 170, sort: 'desc', chartDataType: 'category' as const },
    { field: 'airfoil_hash', headerName: 'Hash', width: 130, hide: true, chartDataType: 'category' as const },
  ];
}

type GridFilterModel = ReturnType<GridApi<RunRow>['getFilterModel']>;

// ────────────────────────────────────────────────────────────
// Component
// ────────────────────────────────────────────────────────────

export function DataExplorerPanel() {
  const { isDark } = useTheme();
  const { startTour } = useOnboarding();
  const {
    allRuns,
    filteredRuns,
    setFilteredRuns,
    clearAll,
    exportDb,
    importDb,
    restoreRunById,
    renameRun,
    selectedRunId,
  } = useRunStore();

  const view = useRouteUiStore((state) => state.dataExplorerView);
  const setView = useRouteUiStore((state) => state.setDataExplorerView);
  const splomKeys = useRouteUiStore((state) => state.dataExplorerSplomKeys);
  const setSplomKeys = useRouteUiStore((state) => state.setDataExplorerSplomKeys);
  const colorBy = useRouteUiStore((state) => state.dataExplorerColorBy);
  const setColorBy = useRouteUiStore((state) => state.setDataExplorerColorBy);
  const filterModel = useRouteUiStore((state) => state.dataExplorerFilterModel);
  const setFilterModel = useRouteUiStore((state) => state.setDataExplorerFilterModel);
  const [dataSource, setDataSource] = useState<DataSource>('full');
  const [showColumnEditor, setShowColumnEditor] = useState(false);
  const customColumns = useCustomColumnStore((s) => s.columns);
  const { numericFields: mergedNumericFields, encodingFields: mergedEncodingFields, getFieldLabel: mergedGetLabel, getFieldValue } = useAllPlotFields();
  const [groupBy, setGroupBy] = useState<keyof RunRow | '' | typeof AUTO_GROUP_KEY>(AUTO_GROUP_KEY);
  const [sizeBy, setSizeBy] = useState<keyof RunRow | ''>('');
  const [symbolBy, setSymbolBy] = useState<keyof RunRow | ''>('');

  // ── AG Grid refs & config ──
  const gridRef = useRef<AgGridReact<RunRow>>(null);
  const apiRef = useRef<GridApi<RunRow> | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [confirmClear, setConfirmClear] = useState(false);

  const gridTheme = useMemo(
    () => isDark ? themeQuartz.withPart(colorSchemeDark) : themeQuartz,
    [isDark],
  );
  const compiledCustomCols = useMemo(() => compileCustomColumns(customColumns), [customColumns]);
  const columnDefs = useMemo(() => {
    const base = buildColumnDefs();
    const custom: ColDef<RunRow>[] = customColumns.map((col) => {
      const evalFn = compiledCustomCols.get(col.id);
      return {
        headerName: col.name,
        colId: `custom_${col.id}`,
        width: 110,
        chartDataType: 'series' as const,
        valueGetter: (params: { data?: RunRow }) => {
          if (!params.data || !evalFn) return null;
          return evalFn(params.data);
        },
        valueFormatter: (p: { value?: number | null }) => p.value?.toFixed(4) ?? '',
      };
    });
    return [...base, ...custom];
  }, [customColumns, compiledCustomCols]);
  const defaultColDef = useMemo<ColDef<RunRow>>(() => ({
    sortable: true, filter: true, resizable: true,
    enableRowGroup: true, enableValue: true,
  }), []);
  const cellSelection = useMemo<CellSelectionOptions>(() => ({ handle: { mode: 'fill' } }), []);

  const onGridReady = useCallback((e: GridReadyEvent<RunRow>) => {
    apiRef.current = e.api;
    if (filterModel) {
      e.api.setFilterModel(filterModel as GridFilterModel);
    }
  }, [filterModel]);
  const onFilterChanged = useCallback((e: FilterChangedEvent<RunRow>) => {
    const filtered: RunRow[] = [];
    e.api.forEachNodeAfterFilterAndSort(node => { if (node.data) filtered.push(node.data); });
    setFilterModel(e.api.getFilterModel());
    setFilteredRuns(filtered);
  }, [setFilterModel, setFilteredRuns]);
  const onCellValueChanged = useCallback((e: CellValueChangedEvent<RunRow>) => {
    if (e.colDef.field === 'airfoil_name' && e.data && e.newValue !== e.oldValue) {
      const trimmed = String(e.newValue ?? '').trim();
      if (trimmed) {
        renameRun(e.data.id, trimmed);
      }
    }
  }, [renameRun]);

  const outlierFilter = useRouteUiStore((state) => state.outlierFilterEnabled);
  const setOutlierFilter = useRouteUiStore((state) => state.setOutlierFilterEnabled);

  // ── Correlogram ──
  const data = dataSource === 'full' ? allRuns : filteredRuns;
  const successOnly = useMemo(() => data.filter((r) => r.success), [data]);

  const toggleSplomKey = useCallback((key: keyof RunRow) => {
    if (splomKeys.includes(key)) {
      setSplomKeys(splomKeys.length <= 2 ? splomKeys : splomKeys.filter(k => k !== key));
      return;
    }
    setSplomKeys([...splomKeys, key]);
  }, [setSplomKeys, splomKeys]);

  useEffect(() => {
    if (!apiRef.current) return;
    apiRef.current.setFilterModel((filterModel as GridFilterModel) ?? null);
  }, [filterModel]);

  const groupedRows = useMemo(() => {
    if (successOnly.length === 0) return [];
    if (groupBy === AUTO_GROUP_KEY) {
      return detectSmartRunGroups(successOnly, {
        sortField: 'alpha',
        plottedFields: splomKeys,
        encodingFields: [colorBy, sizeBy, symbolBy],
      });
    }
    if (!groupBy) {
      return [{
        key: 'all',
        label: `All (${successOnly.length})`,
        rows: successOnly,
        isSinglePoint: successOnly.length === 1,
      }];
    }

    const manualGroups = new Map<string, RunRow[]>();
    for (const row of successOnly) {
      const key = String(row[groupBy] ?? 'null');
      if (!manualGroups.has(key)) manualGroups.set(key, []);
      manualGroups.get(key)!.push(row);
    }

    return [...manualGroups.entries()]
      .sort(([a], [b]) => a.localeCompare(b))
      .map(([key, rows]) => ({
        key: `${String(groupBy)}=${key}`,
        label: `${getLabel(groupBy)}: ${key}`,
        rows,
        isSinglePoint: rows.length === 1,
      }));
  }, [successOnly, groupBy, splomKeys, colorBy, sizeBy, symbolBy]);

  const splomResult = useMemo(() => {
    if (successOnly.length === 0 || splomKeys.length < 2) return null;
    const gv = getFieldValue;
    const gl = mergedGetLabel;

    const n = splomKeys.length;
    const gridGap = n > 4 ? 0.012 : 0.02;
    const cellSize = (1 - gridGap * (n - 1)) / n;
    const annotations: Plotly.Annotations[] = [];
    const gridColor = isDark ? '#333' : '#ddd';
    const traces: Plotly.Data[] = [];
    const axisOverrides: Record<string, unknown> = {};

    // Pre-compute one IQR fence per SPLOM column so each cell only
    // applies its own X/Y fences rather than filtering globally.
    const fenceMap: Record<string, Fence | null> = {};
    if (outlierFilter) {
      for (const key of splomKeys) {
        const vals = successOnly
          .map((r) => gv(r, key as string))
          .filter((v): v is number => typeof v === 'number' && Number.isFinite(v));
        fenceMap[key as string] = buildFence(vals);
      }
    }

    const passesOutlierFence = (run: RunRow, ...keys: string[]): boolean => {
      if (!outlierFilter) return true;
      return keys.every((k) => isInlier(gv(run, k) as number | null, fenceMap[k] ?? null));
    };

    // Axis ranges should reflect the fenced data when filtering is on
    const ranges = Object.fromEntries(
      splomKeys.map((key) => {
        if (outlierFilter) {
          const fence = fenceMap[key as string];
          const values = successOnly
            .map((r) => gv(r, key as string))
            .filter((v): v is number => typeof v === 'number' && Number.isFinite(v) && isInlier(v, fence));
          if (values.length === 0) return [key, [0, 1]];
          const min = Math.min(...values);
          const max = Math.max(...values);
          const pad = (max - min) * 0.05 || 0.5;
          return [key, [min - pad, max + pad]];
        }
        return [key, computeNumericRange(successOnly, key, gv)];
      }),
    ) as Record<keyof RunRow, [number, number]>;

    let colorScaleShown = false;

    for (let row = 0; row < n; row++) {
      for (let col = 0; col < n; col++) {
        const axisIndex = row * n + col + 1;
        const xRef = axisRef('x', axisIndex);
        const yRef = axisRef('y', axisIndex);
        const xLayoutKey = axisLayoutKey('x', axisIndex);
        const yLayoutKey = axisLayoutKey('y', axisIndex);
        const domainX: [number, number] = [
          col * (cellSize + gridGap),
          col * (cellSize + gridGap) + cellSize,
        ];
        const domainY: [number, number] = [
          1 - (row + 1) * cellSize - row * gridGap,
          1 - row * (cellSize + gridGap),
        ];
        const xKey = splomKeys[col];
        const yKey = splomKeys[row];
        const xStr = xKey as string;
        const yStr = yKey as string;

        axisOverrides[xLayoutKey] = {
          domain: domainX,
          anchor: yRef,
          range: ranges[xKey],
          gridcolor: gridColor,
          zerolinecolor: gridColor,
          showline: true,
          mirror: true,
          tickfont: { size: 9 },
          showticklabels: row === n - 1,
          title: row === n - 1 ? { text: gl(xStr), font: { size: 10 } } : undefined,
        };
        axisOverrides[yLayoutKey] = {
          domain: domainY,
          anchor: xRef,
          gridcolor: gridColor,
          zerolinecolor: gridColor,
          showline: true,
          mirror: true,
          tickfont: { size: 9 },
          showticklabels: col === 0,
          title: col === 0 ? { text: gl(yStr), font: { size: 10 } } : undefined,
        };

        if (row === col) {
          // Diagonal: histogram — filter on the single variable
          groupedRows.forEach((group) => {
            const rows = group.rows.filter(
              (run) => gv(run, xStr) != null && passesOutlierFence(run, xStr),
            );
            if (rows.length === 0) return;
            traces.push({
              type: 'histogram',
              x: rows.map((run) => gv(run, xStr)) as number[],
              xaxis: xRef,
              yaxis: yRef,
              name: group.label,
              legendgroup: group.key,
              marker: { color: colorForKey(group.key || xStr), opacity: groupedRows.length > 1 ? 0.55 : 0.85 },
              showlegend: row === 0 && col === 0 && groupedRows.length > 1,
              hovertemplate: `${gl(xStr)}=%{x}<br>Count=%{y}<extra>${group.label}</extra>`,
            });
          });
          continue;
        }

        // Off-diagonal: scatter — filter on both X and Y for this cell
        const cellRows = successOnly.filter(
          (run) => gv(run, xStr) != null && gv(run, yStr) != null && passesOutlierFence(run, xStr, yStr),
        );
        groupedRows.forEach((group, groupIndex) => {
          const rows = group.rows.filter(
            (run) => gv(run, xStr) != null && gv(run, yStr) != null && passesOutlierFence(run, xStr, yStr),
          );
          if (rows.length === 0) return;
          const { marker } = buildMarkerEncoding({
            rows,
            colorBy,
            sizeBy,
            symbolBy,
            defaultColor: colorForKey(group.key || `${xStr}|${yStr}`),
            defaultSize: 5,
            opacity: 0.65,
            minSize: 3,
            maxSize: 11,
            showColorScale: !colorScaleShown && groupIndex === 0,
          });
          if (colorBy && isNumericPlotField(colorBy) && !colorScaleShown && groupIndex === 0) {
            colorScaleShown = true;
          }
          traces.push({
            type: 'scatter',
            mode: 'markers',
            x: rows.map((run) => gv(run, xStr)) as number[],
            y: rows.map((run) => gv(run, yStr)) as number[],
            xaxis: xRef,
            yaxis: yRef,
            name: group.label,
            legendgroup: group.key,
            marker,
            customdata: rows.map(buildRunHoverDetails),
            showlegend: false,
            hovertemplate: `${gl(xStr)}=%{x}<br>${gl(yStr)}=%{y}<br>${RUN_HOVER_TEMPLATE}<extra>${group.label}</extra>`,
          });
        });

        const xValues = cellRows.map((run) => gv(run, xStr) as number);
        const yValues = cellRows.map((run) => gv(run, yStr) as number);
        const r = pearsonR(xValues, yValues);
        if (r != null) {
          annotations.push({
            text: `r=${r.toFixed(2)}`,
            font: { size: 9, color: Math.abs(r) > 0.7 ? colorForKey('corr-strong') : (isDark ? '#888' : '#666') },
            showarrow: false,
            xref: `${xLayoutKey} domain` as Plotly.XAxisName,
            yref: `${yLayoutKey} domain` as Plotly.YAxisName,
            x: 0.95,
            y: 0.95,
            xanchor: 'right',
            yanchor: 'top',
          } as Plotly.Annotations);
        }
      }
    }

    const layout: Partial<Plotly.Layout> = {
      autosize: true,
      margin: { l: 80, r: 30, t: 20, b: groupedRows.length > 1 ? 90 : 60 },
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: isDark ? '#ccc' : '#333', size: 10 },
      showlegend: groupedRows.length > 1,
      legend: { orientation: 'h', y: -0.12 },
      barmode: groupedRows.length > 1 ? 'overlay' : undefined,
      dragmode: 'zoom' as const,
      hovermode: 'closest' as const,
      annotations, ...axisOverrides,
    };

    return { traces, layout };
  }, [successOnly, splomKeys, colorBy, sizeBy, symbolBy, isDark, groupedRows, getFieldValue, mergedGetLabel, outlierFilter]);

  const handlePlotClick = useCallback((event: Readonly<Plotly.PlotMouseEvent>) => {
    const rawRunId = event.points?.[0]?.customdata;
    if (typeof rawRunId === 'number') {
      restoreRunById(rawRunId);
    }
  }, [restoreRunById]);

  // ── Context menu (right-click → Open in Plot Builder) ──
  const setPlotXField = useRouteUiStore((state) => state.setPlotXField);
  const setPlotYField = useRouteUiStore((state) => state.setPlotYField);
  const setPlotChartType = useRouteUiStore((state) => state.setPlotChartType);
  const applyRouteActivePanel = useRouteUiStore((state) => state.applyRouteActivePanel);

  const hoveredCellRef = useRef<{ xKey: keyof RunRow; yKey: keyof RunRow } | null>(null);
  const [contextMenu, setContextMenu] = useState<{ x: number; y: number; xKey: keyof RunRow; yKey: keyof RunRow } | null>(null);

  const handlePlotHover = useCallback((event: Readonly<Plotly.PlotHoverEvent>) => {
    const point = event.points?.[0];
    if (!point) return;
    const xAxisId = (point as any).xaxis?._id as string | undefined;
    if (!xAxisId) return;
    const axisNum = xAxisId === 'x' ? 1 : parseInt(xAxisId.replace('x', ''), 10);
    const n = splomKeys.length;
    const col = (axisNum - 1) % n;
    const row = Math.floor((axisNum - 1) / n);
    if (col < n && row < n) {
      hoveredCellRef.current = { xKey: splomKeys[col], yKey: splomKeys[row] };
    }
  }, [splomKeys]);

  const handleContextMenu = useCallback((e: React.MouseEvent) => {
    const cell = hoveredCellRef.current;
    if (!cell || cell.xKey === cell.yKey) return;
    e.preventDefault();
    const rect = (e.currentTarget as HTMLElement).getBoundingClientRect();
    setContextMenu({ x: e.clientX - rect.left, y: e.clientY - rect.top, xKey: cell.xKey, yKey: cell.yKey });
  }, []);

  const handleOpenInPlotBuilder = useCallback(() => {
    if (!contextMenu) return;
    setPlotXField(contextMenu.xKey);
    setPlotYField(contextMenu.yKey);
    setPlotChartType('scatter');
    applyRouteActivePanel('plot-builder');
    setContextMenu(null);
  }, [contextMenu, setPlotXField, setPlotYField, setPlotChartType, applyRouteActivePanel]);

  useEffect(() => {
    if (!contextMenu) return;
    const dismiss = () => setContextMenu(null);
    window.addEventListener('click', dismiss);
    window.addEventListener('scroll', dismiss, true);
    return () => {
      window.removeEventListener('click', dismiss);
      window.removeEventListener('scroll', dismiss, true);
    };
  }, [contextMenu]);

  // ── Styles ──
  const btnStyle: React.CSSProperties = {
    padding: '2px 8px', fontSize: '10px', background: 'transparent',
    border: '1px solid var(--border-color)', borderRadius: '3px',
    color: 'var(--text-secondary)', cursor: 'pointer',
  };
  const tabStyle = (active: boolean): React.CSSProperties => ({
    padding: '4px 12px', fontSize: '11px', fontWeight: active ? 600 : 400,
    background: active ? (isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.06)') : 'transparent',
    border: 'none', borderBottom: active ? '2px solid var(--accent-primary)' : '2px solid transparent',
    color: active ? 'var(--text-primary)' : 'var(--text-secondary)',
    cursor: 'pointer',
  });
  const chipStyle = (active: boolean): React.CSSProperties => ({
    padding: '2px 7px', fontSize: '10px', borderRadius: '10px',
    border: `1px solid ${active ? colorForKey('accent-active') : 'var(--border-color)'}`,
    background: active ? (isDark ? 'rgba(0,212,170,0.15)' : 'rgba(0,212,170,0.1)') : 'transparent',
    color: active ? colorForKey('accent-active') : 'var(--text-secondary)', cursor: 'pointer', whiteSpace: 'nowrap',
  });
  const selectStyle: React.CSSProperties = {
    padding: '3px 6px', fontSize: '11px',
    background: 'var(--bg-tertiary)', border: '1px solid var(--border-color)',
    borderRadius: '3px', color: 'var(--text-primary)', minWidth: '110px',
  };
  const labelStyle: React.CSSProperties = {
    fontSize: '10px',
    color: 'var(--text-muted)',
    marginBottom: '2px',
  };

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
      {/* Tab bar + toolbar */}
      <div style={{
        display: 'flex', alignItems: 'center', gap: '4px',
        padding: '2px 10px 0', borderBottom: '1px solid var(--border-color)',
        flexShrink: 0, flexWrap: 'wrap',
      }}>
        {/* View tabs */}
        <button data-tour="de-tab-table" onClick={() => setView('table')} style={tabStyle(view === 'table')}>Table</button>
        <button data-tour="de-tab-correlogram" onClick={() => setView('correlogram')} style={tabStyle(view === 'correlogram')}>Correlogram</button>
        <button
          onClick={() => startTour('dataExplorer', true)}
          title="Data Explorer guide"
          style={{
            padding: '1px 6px', fontSize: '11px', fontWeight: 600,
            background: 'transparent', border: '1px solid var(--border-color)',
            borderRadius: '50%', color: 'var(--text-muted)', cursor: 'pointer',
            lineHeight: '16px', width: '20px', height: '20px',
            display: 'flex', alignItems: 'center', justifyContent: 'center',
          }}
        >?</button>

        <span style={{ flex: 1 }} />

        {/* Common actions */}
        <span style={{ fontSize: '11px', color: 'var(--text-muted)', fontWeight: 600 }}>{allRuns.length} runs</span>
        <button
          onClick={() => setShowColumnEditor(true)}
          title="Custom computed columns"
          style={{
            ...btnStyle,
            fontWeight: 600,
            display: 'flex',
            alignItems: 'center',
            gap: '3px',
          }}
        >
          <span style={{ fontSize: '13px', lineHeight: 1 }}>+</span>
          Columns{customColumns.length > 0 ? ` (${customColumns.length})` : ''}
        </button>
        {view === 'table' && (
          <>
            <button
              onClick={() => {
                apiRef.current?.setFilterModel(null);
                setFilterModel(null);
              }}
              style={btnStyle}
            >
              Clear Filters
            </button>
            <button data-tour="de-export-csv" onClick={() => apiRef.current?.exportDataAsCsv({ fileName: 'flexfoil-runs.csv' })} style={btnStyle}>Export CSV</button>
          </>
        )}
        <button
          data-tour="de-export-db"
          onClick={() => {
            const data = exportDb();
            const blob = new Blob([new Uint8Array(data) as BlobPart], { type: 'application/octet-stream' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url; a.download = 'flexfoil-runs.sqlite'; a.click();
            URL.revokeObjectURL(url);
          }}
          style={btnStyle}
        >Export DB</button>
        <button data-tour="de-import-db" onClick={() => fileInputRef.current?.click()} style={btnStyle}>Import DB</button>
        <input ref={fileInputRef} type="file" accept=".sqlite,.db" style={{ display: 'none' }}
          onChange={async (e) => {
            const file = e.target.files?.[0];
            if (!file) return;
            await importDb(new Uint8Array(await file.arrayBuffer()));
            e.target.value = '';
          }}
        />
        <button
          onClick={() => { if (confirmClear) { clearAll(); setConfirmClear(false); } else setConfirmClear(true); }}
          onBlur={() => setConfirmClear(false)}
          style={{ ...btnStyle, background: confirmClear ? 'var(--accent-danger)' : undefined, borderColor: confirmClear ? 'var(--accent-danger)' : undefined, color: confirmClear ? '#fff' : undefined }}
        >{confirmClear ? 'Confirm Clear' : 'Clear All'}</button>
      </div>

      {/* ── Table view ── */}
      {view === 'table' && (
        <div data-tour="de-table-view" style={{ flex: 1, minHeight: 0 }}>
          <AgGridReact<RunRow>
            ref={gridRef}
            theme={gridTheme}
            rowData={allRuns}
            columnDefs={columnDefs}
            defaultColDef={defaultColDef}
            onGridReady={onGridReady}
            onFilterChanged={onFilterChanged}
            onCellValueChanged={onCellValueChanged}
            enableCharts
            enableAdvancedFilter
            cellSelection={cellSelection}
            sideBar
            rowGroupPanelShow="always"
            enableCellTextSelection
            getRowId={(params) => String(params.data.id)}
            autoSizeStrategy={{ type: 'fitCellContents' }}
            statusBar={{
              statusPanels: [
                { statusPanel: 'agTotalAndFilteredRowCountComponent' },
                { statusPanel: 'agAggregationComponent' },
              ],
            }}
          />
        </div>
      )}

      {/* ── Correlogram view ── */}
      {view === 'correlogram' && (
        <div style={{ flex: 1, minHeight: 0, display: 'flex', flexDirection: 'column' }}>
          {/* Column & plot selectors */}
          <div data-tour="de-correlogram-controls" style={{
            padding: '6px 10px', borderBottom: '1px solid var(--border-color)',
            display: 'flex', gap: '6px', flexWrap: 'wrap', alignItems: 'center', flexShrink: 0,
          }}>
            <div data-tour="de-column-chips" style={{ display: 'flex', flexDirection: 'column', minWidth: '280px', flex: 1 }}>
              <span style={labelStyle}>Columns</span>
              <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap' }}>
                {mergedNumericFields.map((f) => (
                  <button key={f.key as string} onClick={() => toggleSplomKey(f.key)} style={chipStyle(splomKeys.includes(f.key))}>
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Data</span>
              <select value={dataSource} onChange={(e) => setDataSource(e.target.value as DataSource)} style={selectStyle}>
                <option value="full">All ({allRuns.length})</option>
                <option value="filtered">Grid Filter ({filteredRuns.length})</option>
              </select>
            </div>

            <div data-tour="de-encoding-controls" style={{ display: 'flex', gap: '6px', flexWrap: 'wrap', alignItems: 'center' }}>
            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Group By</span>
              <select
                value={groupBy as string}
                onChange={(e) => setGroupBy((e.target.value || '') as keyof RunRow | '' | typeof AUTO_GROUP_KEY)}
                style={selectStyle}
              >
                <option value={AUTO_GROUP_KEY}>Auto (Smart)</option>
                <option value="">None</option>
                {mergedEncodingFields.map((f) => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Color</span>
              <select
                value={colorBy as string}
                onChange={(e) => setColorBy((e.target.value || '') as keyof RunRow | '')}
                style={selectStyle}
              >
                <option value="">Color Auto</option>
                {mergedEncodingFields.map((f) => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Marker Size</span>
              <select
                value={sizeBy as string}
                onChange={(e) => setSizeBy((e.target.value || '') as keyof RunRow | '')}
                style={selectStyle}
              >
                <option value="">None</option>
                {mergedEncodingFields.map((f) => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Marker Type</span>
              <select
                value={symbolBy as string}
                onChange={(e) => setSymbolBy((e.target.value || '') as keyof RunRow | '')}
                style={selectStyle}
              >
                <option value="">None</option>
                {mergedEncodingFields.map((f) => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <label style={{
              display: 'flex', alignItems: 'center', gap: '4px', fontSize: '11px',
              color: 'var(--text-secondary)', cursor: 'pointer', userSelect: 'none',
              alignSelf: 'flex-end', paddingBottom: '3px',
            }}>
              <input
                type="checkbox"
                checked={outlierFilter}
                onChange={(e) => setOutlierFilter(e.target.checked)}
                style={{ margin: 0, accentColor: 'var(--accent-primary)' }}
              />
              Remove Outliers
            </label>
            </div>
          </div>

          {/* SPLOM plot */}
          <div data-tour="de-correlogram-plot" style={{ flex: 1, minHeight: 0, position: 'relative' }} onContextMenu={handleContextMenu}>
            {successOnly.length === 0 ? (
              <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', color: 'var(--text-muted)', fontSize: '13px' }}>
                No data. Run some analyses in the Solve panel.
              </div>
            ) : splomKeys.length < 2 ? (
              <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', color: 'var(--text-muted)', fontSize: '13px' }}>
                Select at least 2 columns.
              </div>
            ) : (
              <Plot
                data={(splomResult?.traces ?? []) as Plotly.Data[]}
                layout={(splomResult?.layout ?? {}) as Partial<Plotly.Layout>}
                config={{ responsive: true, displayModeBar: true, displaylogo: false }}
                onClick={handlePlotClick}
                onHover={handlePlotHover}
                useResizeHandler
                style={{ width: '100%', height: '100%' }}
              />
            )}
            {contextMenu && (
              <div
                style={{
                  position: 'absolute',
                  left: contextMenu.x,
                  top: contextMenu.y,
                  zIndex: 1000,
                  background: isDark ? '#1c241f' : '#fff',
                  border: `1px solid ${isDark ? '#3a463d' : '#d2dbd0'}`,
                  borderRadius: '6px',
                  boxShadow: isDark ? '0 4px 16px rgba(0,0,0,0.4)' : '0 4px 16px rgba(0,0,0,0.12)',
                  padding: '4px 0',
                  minWidth: '180px',
                }}
              >
                <button
                  onClick={handleOpenInPlotBuilder}
                  style={{
                    display: 'block',
                    width: '100%',
                    padding: '7px 14px',
                    background: 'transparent',
                    border: 'none',
                    borderRadius: 0,
                    textAlign: 'left',
                    fontSize: '12px',
                    color: isDark ? '#f0f2f0' : '#111813',
                    cursor: 'pointer',
                  }}
                  onMouseEnter={(e) => { e.currentTarget.style.background = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.06)'; }}
                  onMouseLeave={(e) => { e.currentTarget.style.background = 'transparent'; }}
                >
                  Open in Plot Builder
                  <span style={{ fontSize: '10px', color: 'var(--text-muted)', marginLeft: '8px' }}>
                    {getLabel(contextMenu.xKey)} vs {getLabel(contextMenu.yKey)}
                  </span>
                </button>
              </div>
            )}
          </div>
          <div style={{
            padding: '4px 10px',
            borderTop: '1px solid var(--border-color)',
            fontSize: '10px',
            color: 'var(--text-muted)',
            flexShrink: 0,
          }}>
            Diagonal cells show distributions; click any scatter point to load its historical flowfield.
            {selectedRunId != null ? ` Current restored run: #${selectedRunId}` : ''}
          </div>
        </div>
      )}

      {showColumnEditor && (
        <CustomColumnEditor onClose={() => setShowColumnEditor(false)} />
      )}
    </div>
  );
}
