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

import { useCallback, useMemo, useRef, useState } from 'react';
import { AgGridReact } from 'ag-grid-react';
import {
  themeQuartz,
  colorSchemeDark,
  type ColDef,
  type GridReadyEvent,
  type FilterChangedEvent,
  type GridApi,
  type CellSelectionOptions,
} from 'ag-grid-community';
import Plot from 'react-plotly.js';
import { useRunStore } from '../../stores/runStore';
import { useTheme } from '../../contexts/ThemeContext';
import type { RunRow } from '../../types';

// ────────────────────────────────────────────────────────────
// Shared constants
// ────────────────────────────────────────────────────────────

type ViewMode = 'table' | 'correlogram';

interface NumericField { key: keyof RunRow; label: string }

const NUMERIC_FIELDS: NumericField[] = [
  { key: 'alpha', label: 'α (°)' },
  { key: 'cl', label: 'CL' },
  { key: 'cd', label: 'CD' },
  { key: 'cm', label: 'CM' },
  { key: 'reynolds', label: 'Re' },
  { key: 'mach', label: 'Mach' },
  { key: 'ncrit', label: 'Ncrit' },
  { key: 'x_tr_upper', label: 'Xtr Up' },
  { key: 'x_tr_lower', label: 'Xtr Lo' },
  { key: 'iterations', label: 'Iter' },
  { key: 'residual', label: 'Residual' },
  { key: 'n_panels', label: 'Panels' },
];

const COLOR_FIELDS: NumericField[] = [
  ...NUMERIC_FIELDS,
  { key: 'airfoil_name', label: 'Airfoil' },
  { key: 'converged', label: 'Converged' },
  { key: 'session_id', label: 'Session' },
];

const DEFAULT_SPLOM_KEYS: (keyof RunRow)[] = ['alpha', 'cl', 'cd', 'cm'];

const COLORS = [
  '#00d4aa', '#ff6b6b', '#4ecdc4', '#ffe66d', '#a29bfe',
  '#fd79a8', '#6c5ce7', '#00cec9', '#fab1a0', '#74b9ff',
];

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
  return NUMERIC_FIELDS.find(f => f.key === key)?.label
    ?? COLOR_FIELDS.find(f => f.key === key)?.label
    ?? String(key);
}

// ────────────────────────────────────────────────────────────
// AG Grid column definitions
// ────────────────────────────────────────────────────────────

function buildColumnDefs(): ColDef<RunRow>[] {
  return [
    { field: 'id', headerName: 'ID', width: 70, hide: true },
    { field: 'airfoil_name', headerName: 'Airfoil', pinned: 'left', width: 140, enableRowGroup: true, chartDataType: 'category' as const },
    { field: 'alpha', headerName: 'α (°)', width: 90, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(2) ?? '' },
    { field: 'reynolds', headerName: 'Re', width: 120, chartDataType: 'series' as const, valueFormatter: p => p.value != null ? p.value.toExponential(2) : '' },
    { field: 'mach', headerName: 'Mach', width: 80, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(3) ?? '' },
    { field: 'ncrit', headerName: 'Ncrit', width: 80, chartDataType: 'series' as const },
    { field: 'cl', headerName: 'CL', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(5) ?? '' },
    { field: 'cd', headerName: 'CD', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(6) ?? '' },
    { field: 'cm', headerName: 'CM', width: 110, chartDataType: 'series' as const, valueFormatter: p => p.value?.toFixed(5) ?? '' },
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
    { field: 'session_id', headerName: 'Session', width: 130, enableRowGroup: true, hide: true, chartDataType: 'category' as const },
    { field: 'created_at', headerName: 'Time', width: 170, sort: 'desc', chartDataType: 'category' as const },
    { field: 'airfoil_hash', headerName: 'Hash', width: 130, hide: true, chartDataType: 'category' as const },
  ];
}

// ────────────────────────────────────────────────────────────
// Component
// ────────────────────────────────────────────────────────────

export function DataExplorerPanel() {
  const { isDark } = useTheme();
  const { allRuns, setFilteredRuns, clearAll, exportDb, importDb } = useRunStore();

  const [view, setView] = useState<ViewMode>('table');
  const [splomKeys, setSplomKeys] = useState<(keyof RunRow)[]>(DEFAULT_SPLOM_KEYS);
  const [colorBy, setColorBy] = useState<keyof RunRow | ''>('');

  // ── AG Grid refs & config ──
  const gridRef = useRef<AgGridReact<RunRow>>(null);
  const apiRef = useRef<GridApi<RunRow> | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [confirmClear, setConfirmClear] = useState(false);

  const gridTheme = useMemo(
    () => isDark ? themeQuartz.withPart(colorSchemeDark) : themeQuartz,
    [isDark],
  );
  const columnDefs = useMemo(buildColumnDefs, []);
  const defaultColDef = useMemo<ColDef<RunRow>>(() => ({
    sortable: true, filter: true, resizable: true,
    enableRowGroup: true, enableValue: true,
  }), []);
  const cellSelection = useMemo<CellSelectionOptions>(() => ({ handle: { mode: 'fill' } }), []);

  const onGridReady = useCallback((e: GridReadyEvent<RunRow>) => { apiRef.current = e.api; }, []);
  const onFilterChanged = useCallback((e: FilterChangedEvent<RunRow>) => {
    const filtered: RunRow[] = [];
    e.api.forEachNodeAfterFilterAndSort(node => { if (node.data) filtered.push(node.data); });
    setFilteredRuns(filtered);
  }, [setFilteredRuns]);

  // ── Correlogram ──
  const successOnly = useMemo(() => allRuns.filter(r => r.success), [allRuns]);

  const toggleSplomKey = useCallback((key: keyof RunRow) => {
    setSplomKeys(prev => {
      if (prev.includes(key)) return prev.length <= 2 ? prev : prev.filter(k => k !== key);
      return [...prev, key];
    });
  }, []);

  const splomResult = useMemo(() => {
    if (successOnly.length === 0 || splomKeys.length < 2) return null;

    const dimensions = splomKeys.map(key => ({
      label: getLabel(key),
      values: successOnly.map(r => (r[key] as number) ?? 0),
    }));

    let markerConfig: Record<string, unknown> = { size: 3, opacity: 0.6, color: COLORS[0] };
    if (colorBy) {
      const vals = successOnly.map(r => r[colorBy]);
      const isNumeric = vals.every(v => v == null || typeof v === 'number');
      if (isNumeric) {
        markerConfig = { ...markerConfig, color: vals.map(v => (v as number) ?? 0), colorscale: 'Viridis', colorbar: { title: getLabel(colorBy), thickness: 15, len: 0.5 } };
      } else {
        const unique = [...new Set(vals.map(v => String(v ?? 'null')))];
        const colorMap = Object.fromEntries(unique.map((v, i) => [v, COLORS[i % COLORS.length]]));
        markerConfig = { ...markerConfig, color: vals.map(v => colorMap[String(v ?? 'null')]) };
      }
    }

    const trace = { type: 'splom' as const, dimensions, marker: markerConfig, diagonal: { visible: true }, showupperhalf: true, showlowerhalf: true };

    const n = splomKeys.length;
    const annotations: Plotly.Annotations[] = [];
    for (let row = 0; row < n; row++) {
      for (let col = 0; col < n; col++) {
        if (row === col) continue;
        const r = pearsonR(dimensions[col].values, dimensions[row].values);
        if (r == null) continue;
        const axIdx = row * n + col;
        const xKey = axIdx === 0 ? 'xaxis' : `xaxis${axIdx + 1}`;
        const yKey = axIdx === 0 ? 'yaxis' : `yaxis${axIdx + 1}`;
        annotations.push({
          text: `r=${r.toFixed(2)}`,
          font: { size: 9, color: Math.abs(r) > 0.7 ? COLORS[0] : (isDark ? '#888' : '#666') },
          showarrow: false,
          xref: `${xKey} domain` as Plotly.XAxisName,
          yref: `${yKey} domain` as Plotly.YAxisName,
          x: 0.95, y: 0.95, xanchor: 'right', yanchor: 'top',
        } as Plotly.Annotations);
      }
    }

    const gridColor = isDark ? '#333' : '#ddd';
    const axisOverrides: Record<string, unknown> = {};
    for (let i = 0; i < n * n; i++) {
      const suffix = i === 0 ? '' : `${i + 1}`;
      axisOverrides[`xaxis${suffix}`] = { gridcolor: gridColor, zerolinecolor: gridColor };
      axisOverrides[`yaxis${suffix}`] = { gridcolor: gridColor, zerolinecolor: gridColor };
    }

    const layout: Partial<Plotly.Layout> = {
      autosize: true,
      margin: { l: 60, r: 20, t: 20, b: 40 },
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: isDark ? '#ccc' : '#333', size: 10 },
      showlegend: false, dragmode: 'select' as const,
      annotations, ...axisOverrides,
    };

    return { traces: [trace], layout };
  }, [successOnly, splomKeys, colorBy, isDark]);

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
    border: `1px solid ${active ? COLORS[0] : 'var(--border-color)'}`,
    background: active ? (isDark ? 'rgba(0,212,170,0.15)' : 'rgba(0,212,170,0.1)') : 'transparent',
    color: active ? COLORS[0] : 'var(--text-secondary)', cursor: 'pointer', whiteSpace: 'nowrap',
  });

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
      {/* Tab bar + toolbar */}
      <div style={{
        display: 'flex', alignItems: 'center', gap: '4px',
        padding: '2px 10px 0', borderBottom: '1px solid var(--border-color)',
        flexShrink: 0, flexWrap: 'wrap',
      }}>
        {/* View tabs */}
        <button onClick={() => setView('table')} style={tabStyle(view === 'table')}>Table</button>
        <button onClick={() => setView('correlogram')} style={tabStyle(view === 'correlogram')}>Correlogram</button>

        <span style={{ flex: 1 }} />

        {/* Common actions */}
        <span style={{ fontSize: '11px', color: 'var(--text-muted)', fontWeight: 600 }}>{allRuns.length} runs</span>
        {view === 'table' && (
          <>
            <button onClick={() => apiRef.current?.setFilterModel(null)} style={btnStyle}>Clear Filters</button>
            <button onClick={() => apiRef.current?.exportDataAsCsv({ fileName: 'flexfoil-runs.csv' })} style={btnStyle}>Export CSV</button>
          </>
        )}
        <button
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
        <button onClick={() => fileInputRef.current?.click()} style={btnStyle}>Import DB</button>
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
        <div style={{ flex: 1, minHeight: 0 }}>
          <AgGridReact<RunRow>
            ref={gridRef}
            theme={gridTheme}
            rowData={allRuns}
            columnDefs={columnDefs}
            defaultColDef={defaultColDef}
            onGridReady={onGridReady}
            onFilterChanged={onFilterChanged}
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
          {/* Column & color-by selectors */}
          <div style={{
            padding: '6px 10px', borderBottom: '1px solid var(--border-color)',
            display: 'flex', gap: '6px', flexWrap: 'wrap', alignItems: 'center', flexShrink: 0,
          }}>
            <span style={{ fontSize: '10px', color: 'var(--text-muted)', marginRight: '2px' }}>Columns:</span>
            {NUMERIC_FIELDS.map(f => (
              <button key={f.key as string} onClick={() => toggleSplomKey(f.key)} style={chipStyle(splomKeys.includes(f.key))}>
                {f.label}
              </button>
            ))}
            <span style={{ borderLeft: '1px solid var(--border-color)', height: '16px', margin: '0 4px' }} />
            <span style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Color:</span>
            <select
              value={colorBy as string}
              onChange={e => setColorBy((e.target.value || '') as keyof RunRow | '')}
              style={{
                padding: '2px 6px', fontSize: '10px',
                background: 'var(--bg-tertiary)', border: '1px solid var(--border-color)',
                borderRadius: '3px', color: 'var(--text-primary)',
              }}
            >
              <option value="">None</option>
              {COLOR_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
            </select>
          </div>

          {/* SPLOM plot */}
          <div style={{ flex: 1, minHeight: 0, position: 'relative' }}>
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
                useResizeHandler
                style={{ width: '100%', height: '100%' }}
              />
            )}
          </div>
        </div>
      )}
    </div>
  );
}
