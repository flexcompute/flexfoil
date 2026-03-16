import { useState, useMemo, useCallback } from 'react';
import Plot from 'react-plotly.js';
import { useRunStore } from '../../stores/runStore';
import { AUTO_GROUP_KEY, useRouteUiStore } from '../../stores/routeUiStore';
import { useTheme } from '../../contexts/ThemeContext';
import { detectSmartRunGroups } from '../../lib/polarDetection';
import {
  ENCODING_PLOT_FIELDS,
  NUMERIC_PLOT_FIELDS,
  getPlotFieldLabel,
  isNumericPlotField,
} from '../../lib/plotFields';
import { buildMarkerEncoding, colorForKey, colorForValue } from '../../lib/plotStyling';
import type { AxisScale, ChartType, DataSource, RunRow } from '../../types';

const CHART_TYPES: { value: ChartType; label: string }[] = [
  { value: 'scatter', label: 'Scatter' },
  { value: 'line', label: 'Line' },
  { value: 'bar', label: 'Bar' },
  { value: 'histogram', label: 'Histogram' },
];

export function PlotBuilderPanel() {
  const { isDark } = useTheme();
  const { allRuns, filteredRuns, restoreRunById, selectedRunId } = useRunStore();

  const dataSource = useRouteUiStore((state) => state.plotDataSource);
  const setDataSource = useRouteUiStore((state) => state.setPlotDataSource);
  const chartType = useRouteUiStore((state) => state.plotChartType);
  const setChartType = useRouteUiStore((state) => state.setPlotChartType);
  const xField = useRouteUiStore((state) => state.plotXField);
  const setXField = useRouteUiStore((state) => state.setPlotXField);
  const yField = useRouteUiStore((state) => state.plotYField);
  const setYField = useRouteUiStore((state) => state.setPlotYField);
  const groupBy = useRouteUiStore((state) => state.plotGroupBy);
  const setGroupBy = useRouteUiStore((state) => state.setPlotGroupBy);
  const xScale = useRouteUiStore((state) => state.plotXScale);
  const setXScale = useRouteUiStore((state) => state.setPlotXScale);
  const yScale = useRouteUiStore((state) => state.plotYScale);
  const setYScale = useRouteUiStore((state) => state.setPlotYScale);
  const [colorBy, setColorBy] = useState<keyof RunRow | ''>('');
  const [sizeBy, setSizeBy] = useState<keyof RunRow | ''>('');
  const [symbolBy, setSymbolBy] = useState<keyof RunRow | ''>('');

  const data = dataSource === 'full' ? allRuns : filteredRuns;
  const successOnly = useMemo(() => data.filter(r => r.success), [data]);
  const getLabel = useCallback((key: keyof RunRow) => getPlotFieldLabel(key), []);

  const groups = useMemo(() => {
    if (successOnly.length === 0) return [];
    if (groupBy === AUTO_GROUP_KEY) {
      return detectSmartRunGroups(successOnly, {
        sortField: xField,
        plottedFields: chartType === 'histogram' ? [xField] : [xField, yField],
        encodingFields: [colorBy, sizeBy, symbolBy],
      });
    }
    if (!groupBy) {
      return [{
        key: 'all',
        label: chartType === 'histogram' ? getLabel(xField) : `${getLabel(yField)} vs ${getLabel(xField)}`,
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
  }, [successOnly, chartType, xField, yField, colorBy, sizeBy, symbolBy, groupBy, getLabel]);

  const traces = useMemo(() => {
    return groups.reduce<Plotly.Data[]>((allTraces, group, groupIndex) => {
      const defaultColor = colorForKey(group.key);
      if (chartType === 'histogram') {
        const rows = group.rows.filter((row) => row[xField] != null);
        if (rows.length === 0) return allTraces;
        allTraces.push({
          x: rows.map((row) => row[xField]) as number[],
          type: 'histogram' as const,
          name: group.label,
          marker: { color: defaultColor, opacity: 0.8 },
          customdata: rows.map((row) => row.id),
          hovertemplate: `%{x}<extra>${group.label}</extra>`,
        });
        return allTraces;
      }

      const rows = group.rows.filter((row) => row[xField] != null && row[yField] != null);
      if (rows.length === 0) return allTraces;
      const x = rows.map((row) => row[xField]) as number[];
      const y = rows.map((row) => row[yField]) as number[];
      const mode = group.isSinglePoint ? 'markers' as const
        : chartType === 'scatter' ? 'markers' as const
        : chartType === 'line' ? 'lines+markers' as const
        : undefined;
      const { marker, lineColor } = buildMarkerEncoding({
        rows,
        colorBy,
        sizeBy,
        symbolBy,
        defaultColor,
        defaultSize: 7,
        opacity: 0.85,
        minSize: 4,
        maxSize: 14,
        showColorScale: groupIndex === 0,
      });
      const barMarker = colorBy
        ? {
            color: rows.map((row) => (
              colorBy && isNumericPlotField(colorBy)
                ? row[colorBy] as number
                : colorForValue(row[colorBy])
            )),
            ...(colorBy && isNumericPlotField(colorBy)
              ? { colorscale: 'Viridis', colorbar: { title: getLabel(colorBy), thickness: 12, len: 0.45 } }
              : {}),
          }
        : { color: defaultColor };
      allTraces.push({
        x,
        y,
        type: (chartType === 'bar' ? 'bar' : 'scatter') as 'bar' | 'scatter',
        mode,
        name: group.label,
        marker: chartType === 'bar' ? barMarker : marker,
        line: { width: 2, color: lineColor },
        customdata: rows.map((row) => row.id),
      });
      return allTraces;
    }, []);
  }, [groups, chartType, xField, yField, colorBy, sizeBy, symbolBy, getLabel]);

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
    showlegend: traces.length > 1,
  }), [xField, yField, xScale, yScale, chartType, traces.length, getLabel, isDark]);

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
  const handlePlotClick = useCallback((event: Readonly<Plotly.PlotMouseEvent>) => {
    const rawRunId = event.points?.[0]?.customdata;
    if (typeof rawRunId === 'number') {
      restoreRunById(rawRunId);
    }
  }, [restoreRunById]);

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
              {NUMERIC_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
                {NUMERIC_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
          <select value={groupBy as string} onChange={e => setGroupBy((e.target.value || '') as keyof RunRow | '' | typeof AUTO_GROUP_KEY)} style={selectStyle}>
            <option value={AUTO_GROUP_KEY}>Auto (Smart)</option>
            <option value="">None</option>
            {ENCODING_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Color</span>
          <select value={colorBy as string} onChange={e => setColorBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
            <option value="">Color Auto</option>
            {ENCODING_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Marker Size</span>
          <select value={sizeBy as string} onChange={e => setSizeBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
            <option value="">None</option>
            {ENCODING_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Marker Type</span>
          <select value={symbolBy as string} onChange={e => setSymbolBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
            <option value="">None</option>
            {ENCODING_PLOT_FIELDS.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
            onClick={handlePlotClick}
            useResizeHandler
            style={{ width: '100%', height: '100%' }}
          />
        )}
      </div>
      <div style={{
        padding: '4px 10px',
        borderTop: '1px solid var(--border-color)',
        fontSize: '10px',
        color: 'var(--text-muted)',
        flexShrink: 0,
      }}>
        Click a plotted point to restore its historical flowfield.
        {selectedRunId != null ? ` Current restored run: #${selectedRunId}` : ''}
      </div>
    </div>
  );
}
