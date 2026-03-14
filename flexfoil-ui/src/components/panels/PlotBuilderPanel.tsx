/**
 * PlotBuilderPanel — "plot anything vs anything" with Plotly.
 *
 * Inspired by Flexcompute-Thread's PlotBuilder. The user picks X/Y
 * axes from all numeric columns in the run database, optionally
 * groups traces by a categorical field, and selects a chart type.
 *
 * Data source can be switched between the full dataset and the subset
 * currently visible in the AG Grid (filtered rows).
 */

import { useState, useMemo, useCallback } from 'react';
import Plot from 'react-plotly.js';
import { useRunStore } from '../../stores/runStore';
import { useTheme } from '../../contexts/ThemeContext';
import { detectPolarGroups } from '../../lib/polarDetection';
import type { RunRow } from '../../types';

const AUTO_GROUP = '__auto__';

type ChartType = 'scatter' | 'line' | 'bar' | 'histogram';
type DataSource = 'full' | 'filtered';
type AxisScale = 'linear' | 'log';

interface AxisField {
  key: keyof RunRow;
  label: string;
}

const NUMERIC_FIELDS: AxisField[] = [
  { key: 'alpha', label: 'α (°)' },
  { key: 'cl', label: 'CL' },
  { key: 'cd', label: 'CD' },
  { key: 'cm', label: 'CM' },
  { key: 'reynolds', label: 'Re' },
  { key: 'mach', label: 'Mach' },
  { key: 'ncrit', label: 'Ncrit' },
  { key: 'x_tr_upper', label: 'Xtr Upper' },
  { key: 'x_tr_lower', label: 'Xtr Lower' },
  { key: 'iterations', label: 'Iterations' },
  { key: 'residual', label: 'Residual' },
  { key: 'n_panels', label: 'N Panels' },
];

const GROUP_FIELDS: AxisField[] = [
  { key: 'airfoil_name', label: 'Airfoil' },
  { key: 'reynolds', label: 'Re' },
  { key: 'mach', label: 'Mach' },
  { key: 'ncrit', label: 'Ncrit' },
  { key: 'converged', label: 'Converged' },
  { key: 'session_id', label: 'Session' },
  { key: 'n_panels', label: 'N Panels' },
];

const CHART_TYPES: { value: ChartType; label: string }[] = [
  { value: 'scatter', label: 'Scatter' },
  { value: 'line', label: 'Line' },
  { value: 'bar', label: 'Bar' },
  { value: 'histogram', label: 'Histogram' },
];

const COLORS = [
  '#00d4aa', '#ff6b6b', '#4ecdc4', '#ffe66d', '#a29bfe',
  '#fd79a8', '#6c5ce7', '#00cec9', '#fab1a0', '#74b9ff',
  '#55efc4', '#fdcb6e', '#e17055', '#0984e3', '#d63031',
];

export function PlotBuilderPanel() {
  const { isDark } = useTheme();
  const { allRuns, filteredRuns } = useRunStore();

  const [dataSource, setDataSource] = useState<DataSource>('full');
  const [chartType, setChartType] = useState<ChartType>('scatter');
  const [xField, setXField] = useState<keyof RunRow>('alpha');
  const [yField, setYField] = useState<keyof RunRow>('cl');
  const [groupBy, setGroupBy] = useState<keyof RunRow | '' | typeof AUTO_GROUP>(AUTO_GROUP);
  const [xScale, setXScale] = useState<AxisScale>('linear');
  const [yScale, setYScale] = useState<AxisScale>('linear');

  const data = dataSource === 'full' ? allRuns : filteredRuns;
  const successOnly = useMemo(() => data.filter(r => r.success), [data]);

  const getLabel = useCallback((key: keyof RunRow) => {
    return NUMERIC_FIELDS.find(f => f.key === key)?.label
      ?? GROUP_FIELDS.find(f => f.key === key)?.label
      ?? key;
  }, []);

  const traces = useMemo(() => {
    if (successOnly.length === 0) return [];

    if (chartType === 'histogram') {
      return [{
        x: successOnly.map(r => r[xField]).filter(v => v != null) as number[],
        type: 'histogram' as const,
        name: getLabel(xField),
        marker: { color: COLORS[0] },
      }];
    }

    // --- Auto (Polars) mode: compound-key grouping with gap detection ---
    if (groupBy === AUTO_GROUP) {
      const polarGroups = detectPolarGroups(successOnly, xField);

      return polarGroups.map((pg, i) => {
        const x = pg.rows.map(r => r[xField]).filter(v => v != null) as number[];
        const y = pg.rows.map(r => r[yField]).filter(v => v != null) as number[];
        const mode = pg.isSinglePoint ? 'markers' as const
          : chartType === 'scatter' ? 'markers' as const
          : chartType === 'line' ? 'lines+markers' as const
          : undefined;
        return {
          x, y,
          type: (chartType === 'bar' ? 'bar' : 'scatter') as 'bar' | 'scatter',
          mode,
          name: pg.label,
          marker: { color: COLORS[i % COLORS.length], size: 6 },
          line: { width: 2 },
        };
      });
    }

    // --- No grouping ---
    if (!groupBy) {
      const x = successOnly.map(r => r[xField]).filter(v => v != null) as number[];
      const y = successOnly.map(r => r[yField]).filter(v => v != null) as number[];
      const mode = chartType === 'scatter' ? 'markers' as const
        : chartType === 'line' ? 'lines+markers' as const
        : undefined;
      return [{
        x, y,
        type: (chartType === 'bar' ? 'bar' : 'scatter') as 'bar' | 'scatter',
        mode,
        name: `${getLabel(yField)} vs ${getLabel(xField)}`,
        marker: { color: COLORS[0], size: 6 },
        line: { width: 2 },
      }];
    }

    // --- Manual single-field grouping ---
    const groups = new Map<string, RunRow[]>();
    for (const row of successOnly) {
      const key = String(row[groupBy] ?? 'null');
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key)!.push(row);
    }

    return Array.from(groups.entries()).map(([key, rows], i) => {
      const x = rows.map(r => r[xField]).filter(v => v != null) as number[];
      const y = rows.map(r => r[yField]).filter(v => v != null) as number[];
      const mode = chartType === 'scatter' ? 'markers' as const
        : chartType === 'line' ? 'lines+markers' as const
        : undefined;
      return {
        x, y,
        type: (chartType === 'bar' ? 'bar' : 'scatter') as 'bar' | 'scatter',
        mode,
        name: `${getLabel(groupBy as keyof RunRow)}: ${key}`,
        marker: { color: COLORS[i % COLORS.length], size: 6 },
        line: { width: 2 },
      };
    });
  }, [successOnly, chartType, xField, yField, groupBy, getLabel]);

  const layout = useMemo(() => ({
    autosize: true,
    margin: { l: 60, r: 20, t: 30, b: 50 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: isDark ? '#ccc' : '#333', size: 11 },
    xaxis: {
      title: getLabel(xField),
      type: xScale as 'linear' | 'log',
      gridcolor: isDark ? '#333' : '#ddd',
      zerolinecolor: isDark ? '#555' : '#999',
    },
    yaxis: {
      title: chartType !== 'histogram' ? getLabel(yField) : 'Count',
      type: yScale as 'linear' | 'log',
      gridcolor: isDark ? '#333' : '#ddd',
      zerolinecolor: isDark ? '#555' : '#999',
    },
    legend: { orientation: 'h' as const, y: -0.2 },
    showlegend: groupBy === AUTO_GROUP ? traces.length > 1 : !!groupBy,
  }), [xField, yField, xScale, yScale, chartType, groupBy, traces.length, getLabel, isDark]);

  const config = useMemo(() => ({
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'] as Plotly.ModeBarDefaultButtons[],
  }), []);

  const selectStyle: React.CSSProperties = {
    padding: '3px 6px', fontSize: '11px',
    background: 'var(--bg-tertiary)', border: '1px solid var(--border-color)',
    borderRadius: '3px', color: 'var(--text-primary)', minWidth: 0, flex: 1,
  };
  const labelStyle: React.CSSProperties = { fontSize: '10px', color: 'var(--text-muted)', marginBottom: '2px' };

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column', overflow: 'hidden' }}>
      {/* Controls */}
      <div style={{
        padding: '8px 10px', borderBottom: '1px solid var(--border-color)',
        display: 'flex', flexWrap: 'wrap', gap: '8px', alignItems: 'flex-end', flexShrink: 0,
      }}>
        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Data</span>
          <select value={dataSource} onChange={e => setDataSource(e.target.value as DataSource)} style={selectStyle}>
            <option value="full">All ({allRuns.length})</option>
            <option value="filtered">Grid Filter ({filteredRuns.length})</option>
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Chart</span>
          <select value={chartType} onChange={e => setChartType(e.target.value as ChartType)} style={selectStyle}>
            {CHART_TYPES.map(ct => <option key={ct.value} value={ct.value}>{ct.label}</option>)}
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>X Axis</span>
          <div style={{ display: 'flex', gap: '4px' }}>
            <select value={xField as string} onChange={e => setXField(e.target.value as keyof RunRow)} style={selectStyle}>
              {NUMERIC_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
            </select>
            <select value={xScale} onChange={e => setXScale(e.target.value as AxisScale)} style={{ ...selectStyle, flex: 'none', width: '55px' }}>
              <option value="linear">Lin</option>
              <option value="log">Log</option>
            </select>
          </div>
        </div>

        {chartType !== 'histogram' && (
          <div style={{ display: 'flex', flexDirection: 'column' }}>
            <span style={labelStyle}>Y Axis</span>
            <div style={{ display: 'flex', gap: '4px' }}>
              <select value={yField as string} onChange={e => setYField(e.target.value as keyof RunRow)} style={selectStyle}>
                {NUMERIC_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
              <select value={yScale} onChange={e => setYScale(e.target.value as AxisScale)} style={{ ...selectStyle, flex: 'none', width: '55px' }}>
                <option value="linear">Lin</option>
                <option value="log">Log</option>
              </select>
            </div>
          </div>
        )}

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Group By</span>
          <select value={groupBy as string} onChange={e => setGroupBy((e.target.value || '') as keyof RunRow | '' | typeof AUTO_GROUP)} style={selectStyle}>
            <option value={AUTO_GROUP}>Auto (Polars)</option>
            <option value="">None</option>
            {GROUP_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
          </select>
        </div>
      </div>

      {/* Plot area */}
      <div style={{ flex: 1, minHeight: 0, position: 'relative' }}>
        {successOnly.length === 0 ? (
          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', color: 'var(--text-muted)', fontSize: '13px' }}>
            No data. Run some analyses in the Solve panel.
          </div>
        ) : (
          <Plot
            data={traces as Plotly.Data[]}
            layout={layout as Partial<Plotly.Layout>}
            config={config}
            useResizeHandler
            style={{ width: '100%', height: '100%' }}
          />
        )}
      </div>
    </div>
  );
}
