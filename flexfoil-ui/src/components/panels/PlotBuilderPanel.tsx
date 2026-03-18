import { useState, useMemo, useCallback, useEffect, useRef } from 'react';
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
import { useSavedPlotStore, type SavedPlot } from '../../stores/savedPlotStore';
import { useAllPlotFields } from '../../hooks/useAllPlotFields';
import { filterOutliers } from '../../lib/outlierFilter';
import type { AxisScale, ChartType, DataSource, RunRow } from '../../types';

const CHART_TYPES: { value: ChartType; label: string }[] = [
  { value: 'scatter', label: 'Scatter' },
  { value: 'line', label: 'Line' },
  { value: 'bar', label: 'Bar' },
  { value: 'histogram', label: 'Histogram' },
];

const CHART_TYPE_COLORS: Record<string, string> = {
  scatter: '#22c55e',
  line: '#3b82f6',
  bar: '#f59e0b',
  histogram: '#8b5cf6',
};

type PlotBuilderView = 'build' | 'saved';

export function PlotBuilderPanel() {
  const plotAreaRef = useRef<HTMLDivElement | null>(null);
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

  const [view, setView] = useState<PlotBuilderView>('build');
  const [saveTitle, setSaveTitle] = useState('');
  const [showSaveInput, setShowSaveInput] = useState(false);
  const [saveStatus, setSaveStatus] = useState<'idle' | 'saved'>('idle');

  const { savedPlots, addPlot, deletePlot } = useSavedPlotStore();
  const [confirmDeleteId, setConfirmDeleteId] = useState<string | null>(null);

  const {
    numericFields: mergedNumericFields,
    encodingFields: mergedEncodingFields,
    getFieldLabel: mergedGetLabel,
    getFieldValue,
  } = useAllPlotFields();

  const outlierFilter = useRouteUiStore((state) => state.outlierFilterEnabled);
  const setOutlierFilter = useRouteUiStore((state) => state.setOutlierFilterEnabled);

  const data = dataSource === 'full' ? allRuns : filteredRuns;
  const successOnlyRaw = useMemo(() => data.filter(r => r.success), [data]);
  const getLabel = useCallback((key: keyof RunRow) => mergedGetLabel(key as string), [mergedGetLabel]);

  const successOnly = useMemo(() => {
    if (!outlierFilter) return successOnlyRaw;
    const extractors: Array<(row: RunRow) => number | null | undefined> = [
      (row) => getFieldValue(row, xField as string) as number | null,
    ];
    if (chartType !== 'histogram') {
      extractors.push((row) => getFieldValue(row, yField as string) as number | null);
    }
    return filterOutliers(successOnlyRaw, extractors);
  }, [successOnlyRaw, outlierFilter, xField, yField, chartType, getFieldValue]);
  const outlierCount = successOnlyRaw.length - successOnly.length;

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

  const gfv = getFieldValue;
  const traces = useMemo(() => {
    const xKey = xField as string;
    const yKey = yField as string;
    return groups.reduce<Plotly.Data[]>((allTraces, group, groupIndex) => {
      const defaultColor = colorForKey(group.key);
      if (chartType === 'histogram') {
        const rows = group.rows.filter((row) => gfv(row, xKey) != null);
        if (rows.length === 0) return allTraces;
        allTraces.push({
          x: rows.map((row) => gfv(row, xKey)) as number[],
          type: 'histogram' as const,
          name: group.label,
          marker: { color: defaultColor, opacity: 0.8 },
          customdata: rows.map((row) => row.id),
          hovertemplate: `%{x}<extra>${group.label}</extra>`,
        });
        return allTraces;
      }

      const rows = group.rows.filter((row) => gfv(row, xKey) != null && gfv(row, yKey) != null);
      if (rows.length === 0) return allTraces;
      const x = rows.map((row) => gfv(row, xKey)) as number[];
      const y = rows.map((row) => gfv(row, yKey)) as number[];
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
  }, [groups, chartType, xField, yField, colorBy, sizeBy, symbolBy, getLabel, gfv]);

  const layout = useMemo(() => ({
    autosize: true,
    margin: { l: 60, r: 20, t: 30, b: 50 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: isDark ? '#ccc' : '#333', size: 11 },
    xaxis: {
      title: { text: getLabel(xField), font: { size: 12 } },
      type: xScale as 'linear' | 'log',
      gridcolor: isDark ? '#333' : '#ddd',
      zerolinecolor: isDark ? '#555' : '#999',
    },
    yaxis: {
      title: { text: chartType !== 'histogram' ? getLabel(yField) : 'Count', font: { size: 12 } },
      type: yScale as 'linear' | 'log',
      gridcolor: isDark ? '#333' : '#ddd',
      zerolinecolor: isDark ? '#555' : '#999',
    },
    legend: { orientation: 'h' as const, y: -0.2 },
    showlegend: traces.length > 1,
    hovermode: 'closest' as const,
  }), [xField, yField, xScale, yScale, chartType, traces.length, getLabel, isDark]);

  const plotConfig = useMemo(() => ({
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

  const handleSavePlot = useCallback(() => {
    const title = saveTitle.trim() || `${getLabel(yField)} vs ${getLabel(xField)}`;
    addPlot({
      title,
      chartType,
      xField,
      yField,
      groupBy,
      colorBy,
      sizeBy,
      symbolBy,
      xScale,
      yScale,
      dataSource,
    });
    setSaveTitle('');
    setShowSaveInput(false);
    setSaveStatus('saved');
    setTimeout(() => setSaveStatus('idle'), 2000);
  }, [saveTitle, getLabel, yField, xField, addPlot, chartType, groupBy, colorBy, sizeBy, symbolBy, xScale, yScale, dataSource]);

  const handleLoadPlot = useCallback((plot: SavedPlot) => {
    setChartType(plot.chartType);
    setXField(plot.xField);
    setYField(plot.yField);
    setGroupBy(plot.groupBy);
    setColorBy(plot.colorBy);
    setSizeBy(plot.sizeBy);
    setSymbolBy(plot.symbolBy);
    setXScale(plot.xScale);
    setYScale(plot.yScale);
    setDataSource(plot.dataSource);
    setView('build');
  }, [setChartType, setXField, setYField, setGroupBy, setXScale, setYScale, setDataSource]);

  useEffect(() => {
    const container = plotAreaRef.current;
    if (!container) return;

    let frameId: number | null = null;
    const requestPlotResize = () => {
      if (frameId != null) {
        cancelAnimationFrame(frameId);
      }
      frameId = requestAnimationFrame(() => {
        window.dispatchEvent(new Event('resize'));
      });
    };

    const observer = new ResizeObserver(() => {
      requestPlotResize();
    });

    observer.observe(container);
    requestPlotResize();

    return () => {
      observer.disconnect();
      if (frameId != null) {
        cancelAnimationFrame(frameId);
      }
    };
  }, []);

  const tabStyle = (active: boolean): React.CSSProperties => ({
    padding: '4px 12px', fontSize: '11px', fontWeight: active ? 600 : 400,
    background: active ? (isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.06)') : 'transparent',
    border: 'none', borderBottom: active ? '2px solid var(--accent-primary)' : '2px solid transparent',
    color: active ? 'var(--text-primary)' : 'var(--text-secondary)',
    cursor: 'pointer', borderRadius: 0,
  });

  const saveBtnStyle: React.CSSProperties = {
    padding: '3px 10px', fontSize: '11px', fontWeight: 600,
    background: saveStatus === 'saved' ? 'var(--accent-primary)' : 'transparent',
    border: `1px solid ${saveStatus === 'saved' ? 'var(--accent-primary)' : 'var(--border-color)'}`,
    borderRadius: '4px',
    color: saveStatus === 'saved' ? 'var(--bg-primary)' : 'var(--text-secondary)',
    cursor: 'pointer', whiteSpace: 'nowrap',
  };

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column', overflow: 'hidden' }}>
      {/* Tab bar */}
      <div style={{
        display: 'flex', alignItems: 'center', gap: '4px',
        padding: '2px 10px 0', borderBottom: '1px solid var(--border-color)',
        flexShrink: 0,
      }}>
        <button data-tour="pb-tab-build" onClick={() => setView('build')} style={tabStyle(view === 'build')}>Build</button>
        <button data-tour="pb-tab-saved" onClick={() => setView('saved')} style={tabStyle(view === 'saved')}>
          Saved{savedPlots.length > 0 ? ` (${savedPlots.length})` : ''}
        </button>
        <span style={{ flex: 1 }} />
        {view === 'build' && (
          showSaveInput ? (
            <div style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
              <input
                autoFocus
                placeholder={`${getLabel(yField)} vs ${getLabel(xField)}`}
                value={saveTitle}
                onChange={(e) => setSaveTitle(e.target.value)}
                onKeyDown={(e) => { if (e.key === 'Enter') handleSavePlot(); if (e.key === 'Escape') setShowSaveInput(false); }}
                style={{
                  padding: '2px 6px', fontSize: '11px', width: '160px',
                  background: 'var(--bg-tertiary)', border: '1px solid var(--border-color)',
                  borderRadius: '3px', color: 'var(--text-primary)',
                }}
              />
              <button onClick={handleSavePlot} style={{ ...saveBtnStyle, background: 'var(--accent-primary)', color: 'var(--bg-primary)', borderColor: 'var(--accent-primary)' }}>
                Save
              </button>
              <button onClick={() => setShowSaveInput(false)} style={{ ...saveBtnStyle, padding: '3px 6px' }}>
                Cancel
              </button>
            </div>
          ) : (
            <button data-tour="pb-save-plot" onClick={() => setShowSaveInput(true)} style={saveBtnStyle}>
              {saveStatus === 'saved' ? 'Saved!' : 'Save Plot'}
            </button>
          )
        )}
      </div>

      {/* ── Build view ── */}
      {view === 'build' && (
        <>
          {/* Controls */}
          <div data-tour="pb-controls" style={{
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
                  {mergedNumericFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
                    {mergedNumericFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
                {mergedEncodingFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Color</span>
              <select value={colorBy as string} onChange={e => setColorBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
                <option value="">Color Auto</option>
                {mergedEncodingFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Marker Size</span>
              <select value={sizeBy as string} onChange={e => setSizeBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
                <option value="">None</option>
                {mergedEncodingFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
              </select>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column' }}>
              <span style={labelStyle}>Marker Type</span>
              <select value={symbolBy as string} onChange={e => setSymbolBy((e.target.value || '') as keyof RunRow | '')} style={selectStyle}>
                <option value="">None</option>
                {mergedEncodingFields.map(f => <option key={f.key as string} value={f.key as string}>{f.label}</option>)}
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
              {outlierFilter && outlierCount > 0 && (
                <span style={{ color: 'var(--accent-warning, #f59e0b)', fontSize: '10px' }}>
                  ({outlierCount})
                </span>
              )}
            </label>
          </div>

          {/* Plot area */}
          <div data-tour="pb-plot-area" ref={plotAreaRef} style={{ flex: 1, minHeight: 0, position: 'relative' }}>
            {successOnly.length === 0 ? (
              <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', color: 'var(--text-muted)', fontSize: '13px' }}>
                No data. Run some analyses in the Solve panel.
              </div>
            ) : (
              <Plot
                data={traces as Plotly.Data[]}
                layout={layout as Partial<Plotly.Layout>}
                config={plotConfig}
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
        </>
      )}

      {/* ── Saved plots view ── */}
      {view === 'saved' && (
        <div style={{ flex: 1, minHeight: 0, overflow: 'auto', padding: '12px' }}>
          {savedPlots.length === 0 ? (
            <div style={{
              display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center',
              height: '100%', gap: '8px', color: 'var(--text-muted)',
            }}>
              <span style={{ fontSize: '28px', opacity: 0.3 }}>📊</span>
              <span style={{ fontSize: '13px' }}>No saved plots yet.</span>
              <span style={{ fontSize: '11px' }}>Use the Build tab to create and save a plot.</span>
            </div>
          ) : (
            <div style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fill, minmax(200px, 1fr))',
              gap: '10px',
            }}>
              {savedPlots.map((plot) => (
                <SavedPlotCard
                  key={plot.id}
                  plot={plot}
                  isDark={isDark}
                  onLoad={() => handleLoadPlot(plot)}
                  onDelete={() => {
                    if (confirmDeleteId === plot.id) {
                      deletePlot(plot.id);
                      setConfirmDeleteId(null);
                    } else {
                      setConfirmDeleteId(plot.id);
                    }
                  }}
                  isConfirmingDelete={confirmDeleteId === plot.id}
                  onCancelDelete={() => setConfirmDeleteId(null)}
                  getLabel={getLabel}
                />
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

function SavedPlotCard({
  plot,
  isDark,
  onLoad,
  onDelete,
  isConfirmingDelete,
  onCancelDelete,
  getLabel,
}: {
  plot: SavedPlot;
  isDark: boolean;
  onLoad: () => void;
  onDelete: () => void;
  isConfirmingDelete: boolean;
  onCancelDelete: () => void;
  getLabel: (key: keyof RunRow) => string;
}) {
  const accentColor = CHART_TYPE_COLORS[plot.chartType] ?? '#00643c';
  const date = new Date(plot.createdAt);
  const timeStr = date.toLocaleDateString(undefined, {
    month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit',
  });

  return (
    <div
      onClick={onLoad}
      style={{
        borderRadius: '8px',
        border: `1px solid ${isDark ? '#2a302c' : '#d0d5d2'}`,
        background: isDark ? '#1c211e' : '#fff',
        cursor: 'pointer',
        transition: 'border-color 0.15s, box-shadow 0.15s, transform 0.15s',
        overflow: 'hidden',
      }}
      onMouseEnter={(e) => {
        e.currentTarget.style.borderColor = accentColor;
        e.currentTarget.style.boxShadow = isDark ? '0 4px 16px rgba(0,0,0,0.3)' : '0 4px 16px rgba(0,0,0,0.1)';
        e.currentTarget.style.transform = 'translateY(-1px)';
      }}
      onMouseLeave={(e) => {
        e.currentTarget.style.borderColor = isDark ? '#2a302c' : '#d0d5d2';
        e.currentTarget.style.boxShadow = 'none';
        e.currentTarget.style.transform = 'none';
      }}
    >
      <div style={{ height: 3, background: accentColor }} />

      <div style={{
        height: 64,
        background: isDark ? '#161a18' : '#f5f5f4',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '8px 12px',
        gap: '6px',
      }}>
        <span style={{ fontSize: '10px', color: 'var(--text-muted)', fontFamily: 'var(--font-mono)', textAlign: 'center', lineHeight: 1.4 }}>
          {plot.chartType === 'histogram' ? getLabel(plot.xField) : `${getLabel(plot.yField)} vs ${getLabel(plot.xField)}`}
          {plot.groupBy && plot.groupBy !== '__auto__' ? (
            <><br />grouped by {getLabel(plot.groupBy as keyof RunRow)}</>
          ) : plot.groupBy === '__auto__' ? (
            <><br />auto-grouped</>
          ) : null}
        </span>
      </div>

      <div style={{ padding: '8px 10px' }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <div style={{
            fontSize: '12px', fontWeight: 600,
            color: isDark ? '#f0f2f0' : '#111813',
            overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap',
            flex: 1, minWidth: 0,
          }}>
            {plot.title}
          </div>
          {isConfirmingDelete ? (
            <div style={{ display: 'flex', gap: '2px', flexShrink: 0 }}>
              <button
                onClick={(e) => { e.stopPropagation(); onDelete(); }}
                style={{
                  background: 'var(--accent-danger)', border: 'none', borderRadius: '3px',
                  padding: '2px 6px', fontSize: '10px', color: '#fff', cursor: 'pointer',
                }}
              >
                Delete
              </button>
              <button
                onClick={(e) => { e.stopPropagation(); onCancelDelete(); }}
                style={{
                  background: 'transparent', border: '1px solid var(--border-color)', borderRadius: '3px',
                  padding: '2px 5px', fontSize: '10px', color: 'var(--text-muted)', cursor: 'pointer',
                }}
              >
                No
              </button>
            </div>
          ) : (
            <button
              onClick={(e) => { e.stopPropagation(); onDelete(); }}
              title="Delete plot"
              style={{
                background: 'none', border: 'none', cursor: 'pointer',
                padding: '2px 4px', fontSize: '11px', opacity: 0.3,
                color: isDark ? '#f0f2f0' : '#111813', flexShrink: 0,
              }}
              onMouseEnter={(e) => { e.currentTarget.style.opacity = '1'; }}
              onMouseLeave={(e) => { e.currentTarget.style.opacity = '0.3'; }}
            >
              ✕
            </button>
          )}
        </div>
        <div style={{
          display: 'flex', justifyContent: 'space-between', alignItems: 'center',
          marginTop: '4px', fontSize: '10px', color: 'var(--text-muted)',
        }}>
          <span style={{
            padding: '1px 5px', borderRadius: '3px', fontSize: '9px', fontWeight: 600,
            background: isDark ? 'rgba(255,255,255,0.06)' : 'rgba(0,0,0,0.05)',
            color: accentColor,
          }}>
            {plot.chartType}
          </span>
          <span style={{ opacity: 0.6 }}>{timeStr}</span>
        </div>
      </div>
    </div>
  );
}
