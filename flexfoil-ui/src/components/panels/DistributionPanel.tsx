/**
 * DistributionPanel - Surface distribution plot (Cp, Cf, BL quantities)
 *
 * Plots surface quantities against x/c, y, or arc-length s for pinned runs.
 * Runs are pinned via checkbox selection in the Data Explorer table.
 */

import { useMemo, useRef, useEffect } from 'react';
import Plot from 'react-plotly.js';
import { useDistributionStore } from '../../stores/distributionStore';
import { useRunStore } from '../../stores/runStore';
import { useTheme } from '../../contexts/ThemeContext';
import { colorForKey } from '../../lib/plotStyling';
import {
  analyzeAirfoil,
  analyzeAirfoilInviscid,
  getBLDistribution,
  getBLVisualizationData,
  isWasmReady,
  type AnalysisResult,
  type BLDistribution,
} from '../../lib/wasm';
import type { DistributionQuantity, RunRow, SurfaceCoordinate } from '../../types';

const X_LABELS: Record<SurfaceCoordinate, string> = { x: 'x/c', y: 'y/c', s: 's (arc-length)' };
const Y_LABELS: Record<DistributionQuantity, string> = {
  cp: 'Cp',
  cf: 'Cf',
  delta_star: 'δ* (displacement thickness)',
  theta: 'θ (momentum thickness)',
  h: 'H (shape factor)',
  ue: 'ue / U∞',
};

const QUANTITY_OPTIONS: { value: DistributionQuantity; label: string }[] = [
  { value: 'cp', label: 'Cp' },
  { value: 'cf', label: 'Cf' },
  { value: 'delta_star', label: 'δ*' },
  { value: 'theta', label: 'θ' },
  { value: 'h', label: 'H' },
  { value: 'ue', label: 'ue' },
];

const COORD_OPTIONS: { value: SurfaceCoordinate; label: string }[] = [
  { value: 'x', label: 'x/c' },
  { value: 'y', label: 'y' },
  { value: 's', label: 's' },
];

function computeArcLength(points: { x: number; y: number }[]): number[] {
  if (points.length === 0) return [];
  const s = [0];
  for (let i = 1; i < points.length; i++) {
    const dx = points[i].x - points[i - 1].x;
    const dy = points[i].y - points[i - 1].y;
    s.push(s[i - 1] + Math.sqrt(dx * dx + dy * dy));
  }
  const total = s[s.length - 1] || 1;
  return s.map((v) => v / total);
}

interface SurfaceData {
  x_upper: number[];
  x_lower: number[];
  y_upper: number[];
  y_lower: number[];
  s_upper: number[];
  s_lower: number[];
  val_upper: number[];
  val_lower: number[];
}

function extractSurfaceData(
  run: RunRow,
  quantity: DistributionQuantity,
): SurfaceData | null {
  if (!run.geometry_snapshot?.panels || !isWasmReady()) return null;

  const panels = run.geometry_snapshot.panels;
  const isViscous = run.solver_mode === 'viscous';

  try {
    if (quantity === 'cp') {
      const result: AnalysisResult = isViscous
        ? analyzeAirfoil(panels, run.alpha, run.reynolds, run.mach, run.ncrit, run.max_iter)
        : analyzeAirfoilInviscid(panels, run.alpha);

      if (!result.success || !result.cp || !result.cp_x) return null;

      // cp_x and cp are at panel midpoints; split into upper/lower by index
      const n = panels.length;
      const mid = Math.floor(n / 2);
      const cpX = result.cp_x;
      const cp = result.cp;

      // Upper: indices 0..mid-1 (TE→LE on upper), Lower: mid..end (LE→TE on lower)
      const upperIdx: number[] = [];
      const lowerIdx: number[] = [];
      for (let i = 0; i < cpX.length; i++) {
        if (i < mid) upperIdx.push(i);
        else lowerIdx.push(i);
      }

      const upperPanels = panels.slice(0, mid + 1);
      const lowerPanels = panels.slice(mid);
      const sUpper = computeArcLength(upperPanels);
      const sLower = computeArcLength(lowerPanels);
      // midpoints of arc-length
      const sMidUpper = upperIdx.map((_, j) => (sUpper[j] + sUpper[j + 1]) / 2);
      const sMidLower = lowerIdx.map((_, j) => (sLower[j] + sLower[j + 1]) / 2);

      return {
        x_upper: upperIdx.map((i) => cpX[i]),
        x_lower: lowerIdx.map((i) => cpX[i]),
        y_upper: upperIdx.map((i) => (panels[i].y + panels[i + 1].y) / 2),
        y_lower: lowerIdx.map((i) => (panels[i].y + (panels[i + 1]?.y ?? panels[i].y)) / 2),
        s_upper: sMidUpper,
        s_lower: sMidLower,
        val_upper: upperIdx.map((i) => cp[i]),
        val_lower: lowerIdx.map((i) => cp[i]),
      };
    }

    // BL quantities require viscous solve
    if (!isViscous) return null;

    const bl: BLDistribution = getBLDistribution(
      panels,
      run.alpha,
      run.reynolds,
      run.mach,
      run.ncrit,
      run.max_iter,
    );

    if (!bl.success) return null;

    if (quantity === 'ue') {
      const blVis = getBLVisualizationData(
        panels, run.alpha, run.reynolds, run.mach, run.ncrit, run.max_iter,
      );
      if (!blVis.success) return null;

      return {
        x_upper: [...blVis.upper.x],
        x_lower: [...blVis.lower.x],
        y_upper: [...blVis.upper.y],
        y_lower: [...blVis.lower.y],
        s_upper: computeArcLength(blVis.upper.x.map((x, i) => ({ x, y: blVis.upper.y[i] }))),
        s_lower: computeArcLength(blVis.lower.x.map((x, i) => ({ x, y: blVis.lower.y[i] }))),
        val_upper: [...blVis.upper.ue],
        val_lower: [...blVis.lower.ue],
      };
    }

    const valKey = {
      cf: { upper: bl.cf_upper, lower: bl.cf_lower },
      delta_star: { upper: bl.delta_star_upper, lower: bl.delta_star_lower },
      theta: { upper: bl.theta_upper, lower: bl.theta_lower },
      h: { upper: bl.h_upper, lower: bl.h_lower },
    }[quantity];

    if (!valKey) return null;

    const upperPanelsForS = bl.x_upper.map((x) => ({ x, y: 0 }));
    const lowerPanelsForS = bl.x_lower.map((x) => ({ x, y: 0 }));

    return {
      x_upper: [...bl.x_upper],
      x_lower: [...bl.x_lower],
      y_upper: bl.x_upper.map(() => 0),
      y_lower: bl.x_lower.map(() => 0),
      s_upper: computeArcLength(upperPanelsForS),
      s_lower: computeArcLength(lowerPanelsForS),
      val_upper: [...valKey.upper],
      val_lower: [...valKey.lower],
    };
  } catch (e) {
    console.warn('[DistributionPanel] solver error for run', run.id, e);
    return null;
  }
}

function getXValues(data: SurfaceData, coord: SurfaceCoordinate, surface: 'upper' | 'lower'): number[] {
  return surface === 'upper'
    ? (coord === 'x' ? data.x_upper : coord === 'y' ? data.y_upper : data.s_upper)
    : (coord === 'x' ? data.x_lower : coord === 'y' ? data.y_lower : data.s_lower);
}

function buildRunLabel(run: RunRow): string {
  const parts = [run.airfoil_name, `α=${run.alpha}°`];
  if (run.solver_mode === 'viscous') {
    parts.push(`Re=${run.reynolds.toExponential(1)}`);
  }
  return parts.join(' ');
}

export function DistributionPanel() {
  const plotAreaRef = useRef<HTMLDivElement | null>(null);
  const { isDark } = useTheme();
  const allRuns = useRunStore((s) => s.allRuns);

  const pinnedRunIds = useDistributionStore((s) => s.pinnedRunIds);
  const xAxis = useDistributionStore((s) => s.distributionXAxis);
  const setXAxis = useDistributionStore((s) => s.setDistributionXAxis);
  const yAxis = useDistributionStore((s) => s.distributionYAxis);
  const setYAxis = useDistributionStore((s) => s.setDistributionYAxis);
  const showUpper = useDistributionStore((s) => s.showUpper);
  const setShowUpper = useDistributionStore((s) => s.setShowUpper);
  const showLower = useDistributionStore((s) => s.showLower);
  const setShowLower = useDistributionStore((s) => s.setShowLower);
  const clearPinned = useDistributionStore((s) => s.clearPinned);

  const pinnedRuns = useMemo(
    () => allRuns.filter((r) => pinnedRunIds.has(r.id) && r.success),
    [allRuns, pinnedRunIds],
  );

  const traces = useMemo(() => {
    const result: Plotly.Data[] = [];
    for (const run of pinnedRuns) {
      const data = extractSurfaceData(run, yAxis);
      if (!data) continue;

      const color = colorForKey(`dist_${run.id}`);
      const label = buildRunLabel(run);

      if (showUpper) {
        result.push({
          x: getXValues(data, xAxis, 'upper'),
          y: data.val_upper,
          type: 'scatter',
          mode: 'lines+markers',
          name: `${label} (upper)`,
          marker: { color, size: 3 },
          line: { color, width: 2 },
        });
      }

      if (showLower) {
        result.push({
          x: getXValues(data, xAxis, 'lower'),
          y: data.val_lower,
          type: 'scatter',
          mode: 'lines+markers',
          name: `${label} (lower)`,
          marker: { color, size: 3, symbol: 'diamond' },
          line: { color, width: 2, dash: 'dash' },
        });
      }
    }
    return result;
  }, [pinnedRuns, xAxis, yAxis, showUpper, showLower]);

  const isCp = yAxis === 'cp';

  const layout = useMemo(
    () => ({
      autosize: true,
      margin: { l: 60, r: 20, t: 30, b: 50 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: isDark ? '#ccc' : '#333', size: 11 },
      xaxis: {
        title: { text: X_LABELS[xAxis], font: { size: 12 } },
        gridcolor: isDark ? '#333' : '#ddd',
        zerolinecolor: isDark ? '#555' : '#999',
      },
      yaxis: {
        title: { text: Y_LABELS[yAxis], font: { size: 12 } },
        gridcolor: isDark ? '#333' : '#ddd',
        zerolinecolor: isDark ? '#555' : '#999',
        autorange: isCp ? ('reversed' as const) : (true as const),
      },
      legend: { orientation: 'h' as const, y: -0.2 },
      showlegend: traces.length > 1,
      hovermode: 'closest' as const,
    }),
    [xAxis, yAxis, isCp, traces.length, isDark],
  );

  const plotConfig = useMemo(
    () => ({
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      modeBarButtonsToRemove: ['lasso2d', 'select2d'] as Plotly.ModeBarDefaultButtons[],
    }),
    [],
  );

  useEffect(() => {
    const container = plotAreaRef.current;
    if (!container) return;

    let frameId: number | null = null;
    const requestResize = () => {
      if (frameId != null) cancelAnimationFrame(frameId);
      frameId = requestAnimationFrame(() => window.dispatchEvent(new Event('resize')));
    };
    const observer = new ResizeObserver(() => requestResize());
    observer.observe(container);
    requestResize();
    return () => {
      observer.disconnect();
      if (frameId != null) cancelAnimationFrame(frameId);
    };
  }, []);

  const selectStyle: React.CSSProperties = {
    padding: '3px 6px',
    fontSize: '11px',
    background: 'var(--bg-tertiary)',
    border: '1px solid var(--border-color)',
    borderRadius: '3px',
    color: 'var(--text-primary)',
    minWidth: 0,
  };
  const labelStyle: React.CSSProperties = { fontSize: '10px', color: 'var(--text-muted)', marginBottom: '2px' };

  return (
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column', overflow: 'hidden' }}>
      {/* Controls */}
      <div
        style={{
          padding: '8px 10px',
          borderBottom: '1px solid var(--border-color)',
          display: 'flex',
          flexWrap: 'wrap',
          gap: '8px',
          alignItems: 'flex-end',
          flexShrink: 0,
        }}
      >
        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>X Axis</span>
          <select
            value={xAxis}
            onChange={(e) => setXAxis(e.target.value as SurfaceCoordinate)}
            style={selectStyle}
          >
            {COORD_OPTIONS.map((o) => (
              <option key={o.value} value={o.value}>
                {o.label}
              </option>
            ))}
          </select>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <span style={labelStyle}>Quantity</span>
          <select
            value={yAxis}
            onChange={(e) => setYAxis(e.target.value as DistributionQuantity)}
            style={selectStyle}
          >
            {QUANTITY_OPTIONS.map((o) => (
              <option key={o.value} value={o.value}>
                {o.label}
              </option>
            ))}
          </select>
        </div>

        <label
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: '4px',
            fontSize: '11px',
            color: 'var(--text-secondary)',
            cursor: 'pointer',
            userSelect: 'none',
            alignSelf: 'flex-end',
            paddingBottom: '3px',
          }}
        >
          <input
            type="checkbox"
            checked={showUpper}
            onChange={(e) => setShowUpper(e.target.checked)}
            style={{ margin: 0, accentColor: 'var(--accent-primary)' }}
          />
          Upper
        </label>

        <label
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: '4px',
            fontSize: '11px',
            color: 'var(--text-secondary)',
            cursor: 'pointer',
            userSelect: 'none',
            alignSelf: 'flex-end',
            paddingBottom: '3px',
          }}
        >
          <input
            type="checkbox"
            checked={showLower}
            onChange={(e) => setShowLower(e.target.checked)}
            style={{ margin: 0, accentColor: 'var(--accent-primary)' }}
          />
          Lower
        </label>

        <span style={{ flex: 1 }} />

        <span
          style={{
            fontSize: '10px',
            color: 'var(--text-muted)',
            alignSelf: 'flex-end',
            paddingBottom: '3px',
          }}
        >
          {pinnedRuns.length} run{pinnedRuns.length !== 1 ? 's' : ''} pinned
        </span>

        {pinnedRuns.length > 0 && (
          <button
            onClick={clearPinned}
            style={{
              fontSize: '10px',
              padding: '2px 8px',
              background: 'transparent',
              border: '1px solid var(--border-color)',
              borderRadius: '3px',
              color: 'var(--text-muted)',
              cursor: 'pointer',
              alignSelf: 'flex-end',
              marginBottom: '1px',
            }}
          >
            Clear
          </button>
        )}
      </div>

      {/* Plot area */}
      <div ref={plotAreaRef} style={{ flex: 1, minHeight: 0, position: 'relative' }}>
        {pinnedRuns.length === 0 ? (
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              height: '100%',
              color: 'var(--text-muted)',
              fontSize: '13px',
              flexDirection: 'column',
              gap: '8px',
            }}
          >
            <span>No runs pinned.</span>
            <span style={{ fontSize: '11px' }}>
              Select runs with checkboxes in the Data Explorer table.
            </span>
          </div>
        ) : (
          <Plot
            data={traces as Plotly.Data[]}
            layout={layout as Partial<Plotly.Layout>}
            config={plotConfig}
            useResizeHandler
            style={{ width: '100%', height: '100%' }}
          />
        )}
      </div>

      {/* Footer */}
      <div
        style={{
          padding: '4px 10px',
          borderTop: '1px solid var(--border-color)',
          fontSize: '10px',
          color: 'var(--text-muted)',
          flexShrink: 0,
        }}
      >
        {isCp ? 'Cp axis is inverted (negative up) per convention.' : 'Select runs in Data Explorer to overlay distributions.'}
      </div>
    </div>
  );
}
