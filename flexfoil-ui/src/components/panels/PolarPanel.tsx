/**
 * PolarPanel - Aerodynamic polar plot visualization
 * 
 * Displays multiple polar series overlaid (Cl vs alpha, Cl vs Cd, etc.)
 * with configurable axes and a color-coded legend.
 */

import { useRef, useEffect, useCallback, useMemo, useState } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useRouteUiStore } from '../../stores/routeUiStore';
import { colorForKey } from '../../lib/plotStyling';
import { filterOutliers } from '../../lib/outlierFilter';
import type { AxisVariable, PolarPoint } from '../../types';

const PLOT_MARGIN = { top: 20, right: 20, bottom: 40, left: 50 };

const AXIS_LABELS: Record<AxisVariable, string> = {
  alpha: 'α (deg)',
  cl: 'Cl',
  cd: 'Cd',
  cm: 'Cm',
  ld: 'L/D',
  reynolds: 'Re',
  mach: 'Mach',
  ncrit: 'Ncrit',
  flapDeflection: 'Flap δ (deg)',
  flapHingeX: 'Flap x/c',
};

function getValue(point: PolarPoint, variable: AxisVariable): number {
  if (variable === 'ld') {
    const cd = point.cd ?? 0;
    return Math.abs(cd) > 1e-10 ? point.cl / cd : 0;
  }
  return (point as Record<string, number | undefined>)[variable] ?? 0;
}

function getAxisBounds(values: number[]): [number, number] {
  if (values.length === 0) return [0, 1];
  const min = Math.min(...values);
  const max = Math.max(...values);
  const padding = (max - min) * 0.1 || 0.5;
  return [min - padding, max + padding];
}

function generateTicks(min: number, max: number, count: number = 5): number[] {
  const range = max - min;
  const step = range / (count - 1);
  const magnitude = Math.pow(10, Math.floor(Math.log10(step)));
  const niceStep = Math.ceil(step / magnitude) * magnitude;
  
  const ticks: number[] = [];
  let tick = Math.floor(min / niceStep) * niceStep;
  while (tick <= max + niceStep * 0.5) {
    if (tick >= min - niceStep * 0.1) {
      ticks.push(tick);
    }
    tick += niceStep;
  }
  return ticks;
}

export function PolarPanel() {
  const { polarData, removePolar, clearAllPolars } = useAirfoilStore();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  const xAxis = useRouteUiStore((state) => state.polarXAxis);
  const setXAxis = useRouteUiStore((state) => state.setPolarXAxis);
  const yAxis = useRouteUiStore((state) => state.polarYAxis);
  const setYAxis = useRouteUiStore((state) => state.setPolarYAxis);
  const outlierFilter = useRouteUiStore((state) => state.outlierFilterEnabled);
  const setOutlierFilter = useRouteUiStore((state) => state.setOutlierFilterEnabled);
  const [canvasSize, setCanvasSize] = useState({ width: 400, height: 300 });

  const filteredPolarData = useMemo(() => {
    if (!outlierFilter) return polarData;
    return polarData.map((series) => ({
      ...series,
      points: filterOutliers(series.points, [
        (p) => getValue(p, xAxis),
        (p) => getValue(p, yAxis),
      ]),
    }));
  }, [polarData, outlierFilter, xAxis, yAxis]);

  const totalPoints = useMemo(
    () => filteredPolarData.reduce((n, s) => n + s.points.length, 0),
    [filteredPolarData],
  );
  const rawTotalPoints = useMemo(
    () => polarData.reduce((n, s) => n + s.points.length, 0),
    [polarData],
  );
  const outlierCount = rawTotalPoints - totalPoints;

  const { xBounds, yBounds } = useMemo(() => {
    const allX: number[] = [];
    const allY: number[] = [];
    for (const series of filteredPolarData) {
      for (const p of series.points) {
        allX.push(getValue(p, xAxis));
        allY.push(getValue(p, yAxis));
      }
    }
    return {
      xBounds: getAxisBounds(allX),
      yBounds: getAxisBounds(allY),
    };
  }, [filteredPolarData, xAxis, yAxis]);
  
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    
    const { width, height } = canvasSize;
    const plotWidth = width - PLOT_MARGIN.left - PLOT_MARGIN.right;
    const plotHeight = height - PLOT_MARGIN.top - PLOT_MARGIN.bottom;
    
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--bg-primary').trim() || '#0f0f0f';
    ctx.fillRect(0, 0, width, height);
    
    const toCanvasX = (x: number) =>
      PLOT_MARGIN.left + ((x - xBounds[0]) / (xBounds[1] - xBounds[0])) * plotWidth;
    const toCanvasY = (y: number) =>
      PLOT_MARGIN.top + (1 - (y - yBounds[0]) / (yBounds[1] - yBounds[0])) * plotHeight;
    
    // Grid
    ctx.strokeStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--foil-grid').trim() || '#333333';
    ctx.lineWidth = 1;
    
    const xTicks = generateTicks(xBounds[0], xBounds[1]);
    const yTicks = generateTicks(yBounds[0], yBounds[1]);
    
    for (const x of xTicks) {
      const cx = toCanvasX(x);
      ctx.beginPath();
      ctx.moveTo(cx, PLOT_MARGIN.top);
      ctx.lineTo(cx, height - PLOT_MARGIN.bottom);
      ctx.stroke();
    }
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      ctx.beginPath();
      ctx.moveTo(PLOT_MARGIN.left, cy);
      ctx.lineTo(width - PLOT_MARGIN.right, cy);
      ctx.stroke();
    }
    
    // Axes
    ctx.strokeStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--text-secondary').trim() || '#888888';
    ctx.globalAlpha = 0.3;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(PLOT_MARGIN.left, height - PLOT_MARGIN.bottom);
    ctx.lineTo(width - PLOT_MARGIN.right, height - PLOT_MARGIN.bottom);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(PLOT_MARGIN.left, PLOT_MARGIN.top);
    ctx.lineTo(PLOT_MARGIN.left, height - PLOT_MARGIN.bottom);
    ctx.stroke();
    ctx.globalAlpha = 1;
    
    // Tick labels
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--text-secondary').trim() || '#888888';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    for (const x of xTicks) {
      ctx.fillText(x.toFixed(1), toCanvasX(x), height - PLOT_MARGIN.bottom + 15);
    }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const y of yTicks) {
      ctx.fillText(y.toFixed(2), PLOT_MARGIN.left - 5, toCanvasY(y));
    }
    
    // Axis labels
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--text-primary').trim() || '#000000';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText(AXIS_LABELS[xAxis], width / 2, height - 15);
    ctx.save();
    ctx.translate(15, height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'bottom';
    ctx.fillText(AXIS_LABELS[yAxis], 0, 0);
    ctx.restore();
    
    // Draw each series
    if (totalPoints > 0) {
      filteredPolarData.forEach((series) => {
        const color = colorForKey(series.key);
        const pts = series.points;
        if (pts.length === 0) return;

        const xVals = pts.map(p => getValue(p, xAxis));
        const yVals = pts.map(p => getValue(p, yAxis));

        // Line
        ctx.strokeStyle = color;
        ctx.globalAlpha = 0.8;
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < pts.length; i++) {
          const x = toCanvasX(xVals[i]);
          const y = toCanvasY(yVals[i]);
          if (i === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.globalAlpha = 1;

        // Points
        ctx.fillStyle = color;
        for (let i = 0; i < pts.length; i++) {
          ctx.beginPath();
          ctx.arc(toCanvasX(xVals[i]), toCanvasY(yVals[i]), 3, 0, Math.PI * 2);
          ctx.fill();
        }
      });
    } else {
      ctx.fillStyle = getComputedStyle(document.documentElement)
        .getPropertyValue('--text-secondary').trim() || '#888888';
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No polar data. Run a polar sweep in the Solve panel.', width / 2, height / 2);
    }
  }, [canvasSize, filteredPolarData, totalPoints, xAxis, yAxis, xBounds, yBounds]);
  
  // Resize observer
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;
    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        setCanvasSize({ width: Math.max(200, width), height: Math.max(150, height) });
      }
    });
    observer.observe(container);
    const rect = container.getBoundingClientRect();
    setCanvasSize({ width: Math.max(200, rect.width), height: Math.max(150, rect.height) });
    return () => observer.disconnect();
  }, []);
  
  useEffect(() => { draw(); }, [draw]);
  
  return (
    <div className="panel-content" style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      {/* Controls */}
      <div style={{ 
        display: 'flex', 
        gap: '12px', 
        padding: '8px',
        borderBottom: '1px solid var(--border-color)',
        alignItems: 'center',
      }}>
        <label style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '12px' }}>
          X:
          <select 
            value={xAxis}
            onChange={(e) => setXAxis(e.target.value as AxisVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            {Object.entries(AXIS_LABELS).map(([k, v]) => (
              <option key={k} value={k}>{v}</option>
            ))}
          </select>
        </label>
        
        <label style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '12px' }}>
          Y:
          <select 
            value={yAxis}
            onChange={(e) => setYAxis(e.target.value as AxisVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            {Object.entries(AXIS_LABELS).map(([k, v]) => (
              <option key={k} value={k}>{v}</option>
            ))}
          </select>
        </label>
        
        <label style={{
          display: 'flex', alignItems: 'center', gap: '4px', fontSize: '11px',
          color: 'var(--text-secondary)', cursor: 'pointer', userSelect: 'none',
        }}>
          <input
            type="checkbox"
            checked={outlierFilter}
            onChange={(e) => setOutlierFilter(e.target.checked)}
            style={{ margin: 0, accentColor: 'var(--accent-primary)' }}
          />
          Remove Outliers
        </label>
        
        <span style={{ 
          marginLeft: 'auto', 
          fontSize: '11px', 
          color: 'var(--text-secondary)' 
        }}>
          {totalPoints} pts
          {outlierFilter && outlierCount > 0 && (
            <span style={{ color: 'var(--accent-warning, #f59e0b)' }}> ({outlierCount} removed)</span>
          )}
          {filteredPolarData.length > 1 && ` · ${filteredPolarData.length} series`}
        </span>
      </div>
      
      {/* Canvas */}
      <div 
        ref={containerRef}
        style={{ flex: 1, minHeight: 0, overflow: 'hidden' }}
      >
        <canvas
          ref={canvasRef}
          width={canvasSize.width}
          height={canvasSize.height}
          style={{ display: 'block', width: '100%', height: '100%' }}
        />
      </div>

      {/* Legend */}
      {filteredPolarData.length > 0 && (
        <div style={{
          padding: '6px 8px',
          borderTop: '1px solid var(--border-color)',
          display: 'flex',
          flexDirection: 'column',
          gap: '3px',
          maxHeight: '100px',
          overflowY: 'auto',
        }}>
          {filteredPolarData.map((series) => (
            <div
              key={series.key}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: '6px',
                fontSize: '10px',
              }}
            >
              <span style={{
                width: 10,
                height: 10,
                borderRadius: 2,
                background: colorForKey(series.key),
                flexShrink: 0,
              }} />
              <span style={{
                flex: 1,
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                whiteSpace: 'nowrap',
                color: 'var(--text-secondary)',
              }}>
                {series.label} ({series.points.length} pts)
              </span>
              <button
                onClick={() => removePolar(series.key)}
                style={{
                  background: 'transparent',
                  border: 'none',
                  color: 'var(--text-muted)',
                  cursor: 'pointer',
                  padding: '0 2px',
                  fontSize: '11px',
                  lineHeight: 1,
                }}
                title="Remove this polar"
              >
                ×
              </button>
            </div>
          ))}
          {filteredPolarData.length > 1 && (
            <button
              onClick={clearAllPolars}
              style={{
                fontSize: '10px',
                padding: '2px 6px',
                background: 'transparent',
                border: '1px solid var(--border-color)',
                borderRadius: '3px',
                color: 'var(--text-muted)',
                cursor: 'pointer',
                alignSelf: 'flex-end',
                marginTop: '2px',
              }}
            >
              Clear all
            </button>
          )}
        </div>
      )}
    </div>
  );
}
