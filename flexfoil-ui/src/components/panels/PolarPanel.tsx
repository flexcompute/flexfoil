/**
 * PolarPanel - Aerodynamic polar plot visualization
 * 
 * Displays polar data (Cl vs alpha, Cl vs Cd, Cl vs Cm, etc.) with configurable axes.
 * Supports both inviscid and viscous polars.
 */

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import type { PolarPoint } from '../../types';

type AxisVariable = 'alpha' | 'cl' | 'cd' | 'cm' | 'ld';

const AXIS_LABELS: Record<AxisVariable, string> = {
  alpha: 'α (deg)',
  cl: 'Cl',
  cd: 'Cd',
  cm: 'Cm',
  ld: 'L/D',
};

// Get value from polar point by variable name
function getValue(point: PolarPoint, variable: AxisVariable): number {
  switch (variable) {
    case 'alpha':
      return point.alpha;
    case 'cl':
      return point.cl;
    case 'cd':
      return point.cd ?? 0;
    case 'cm':
      return point.cm;
    case 'ld':
      return point.cd && point.cd > 0 ? point.cl / point.cd : 0;
    default:
      return 0;
  }
}

// Auto-scale axis bounds with some padding
function getAxisBounds(values: number[]): [number, number] {
  if (values.length === 0) return [0, 1];
  const min = Math.min(...values);
  const max = Math.max(...values);
  const padding = (max - min) * 0.1 || 0.5;
  return [min - padding, max + padding];
}

// Generate nice tick values
function generateTicks(min: number, max: number, count: number = 5): number[] {
  const range = max - min;
  if (range <= 0) return [min];
  
  const step = range / (count - 1);
  
  // Round step to nice value
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
  const { polarData, solverMode } = useAirfoilStore();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  // Default to Cl vs alpha, but Cl vs Cd for viscous
  const [xAxis, setXAxis] = useState<AxisVariable>('alpha');
  const [yAxis, setYAxis] = useState<AxisVariable>('cl');
  const [canvasSize, setCanvasSize] = useState({ width: 400, height: 300 });
  
  // Check if we have Cd data
  const hasCdData = useMemo(() => {
    return polarData.some(p => p.cd !== undefined && p.cd !== null);
  }, [polarData]);
  
  // Margins for axis labels
  const margin = { top: 20, right: 20, bottom: 40, left: 50 };
  
  // Extract data for current axis selection
  const { xValues, yValues } = useMemo(() => {
    const xValues = polarData.map(p => getValue(p, xAxis));
    const yValues = polarData.map(p => getValue(p, yAxis));
    return { xValues, yValues };
  }, [polarData, xAxis, yAxis]);
  
  // Compute axis bounds
  const { xBounds, yBounds } = useMemo(() => {
    return {
      xBounds: getAxisBounds(xValues),
      yBounds: getAxisBounds(yValues),
    };
  }, [xValues, yValues]);
  
  // Draw the plot
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    
    const { width, height } = canvasSize;
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    
    // Clear
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--bg-primary').trim() || '#0f0f0f';
    ctx.fillRect(0, 0, width, height);
    
    // Transform functions
    const toCanvasX = (x: number) => {
      return margin.left + ((x - xBounds[0]) / (xBounds[1] - xBounds[0])) * plotWidth;
    };
    const toCanvasY = (y: number) => {
      return margin.top + (1 - (y - yBounds[0]) / (yBounds[1] - yBounds[0])) * plotHeight;
    };
    
    // Draw grid
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.1)';
    ctx.lineWidth = 1;
    
    const xTicks = generateTicks(xBounds[0], xBounds[1]);
    const yTicks = generateTicks(yBounds[0], yBounds[1]);
    
    for (const x of xTicks) {
      const cx = toCanvasX(x);
      ctx.beginPath();
      ctx.moveTo(cx, margin.top);
      ctx.lineTo(cx, height - margin.bottom);
      ctx.stroke();
    }
    
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      ctx.beginPath();
      ctx.moveTo(margin.left, cy);
      ctx.lineTo(width - margin.right, cy);
      ctx.stroke();
    }
    
    // Draw axes
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
    ctx.lineWidth = 1;
    
    // X axis
    ctx.beginPath();
    ctx.moveTo(margin.left, height - margin.bottom);
    ctx.lineTo(width - margin.right, height - margin.bottom);
    ctx.stroke();
    
    // Y axis
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, height - margin.bottom);
    ctx.stroke();
    
    // Draw tick labels
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    
    for (const x of xTicks) {
      const cx = toCanvasX(x);
      const label = xAxis === 'cd' ? x.toFixed(4) : x.toFixed(1);
      ctx.fillText(label, cx, height - margin.bottom + 15);
    }
    
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      const label = yAxis === 'cd' ? y.toFixed(4) : y.toFixed(2);
      ctx.fillText(label, margin.left - 5, cy);
    }
    
    // Draw axis labels
    ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
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
    
    // Draw data points and line
    if (polarData.length > 0) {
      // Line
      ctx.strokeStyle = hasCdData ? 'rgba(0, 200, 150, 0.8)' : 'rgba(100, 150, 255, 0.8)';
      ctx.lineWidth = 2;
      ctx.beginPath();
      
      for (let i = 0; i < polarData.length; i++) {
        const x = toCanvasX(xValues[i]);
        const y = toCanvasY(yValues[i]);
        if (i === 0) {
          ctx.moveTo(x, y);
        } else {
          ctx.lineTo(x, y);
        }
      }
      ctx.stroke();
      
      // Points
      ctx.fillStyle = hasCdData ? '#00d4aa' : '#6b9fff';
      for (let i = 0; i < polarData.length; i++) {
        const x = toCanvasX(xValues[i]);
        const y = toCanvasY(yValues[i]);
        ctx.beginPath();
        ctx.arc(x, y, 4, 0, Math.PI * 2);
        ctx.fill();
      }
      
      // Highlight max L/D point if showing drag polar
      if (xAxis === 'cd' && yAxis === 'cl' && hasCdData) {
        let maxLD = 0;
        let maxLDIdx = 0;
        for (let i = 0; i < polarData.length; i++) {
          const ld = getValue(polarData[i], 'ld');
          if (ld > maxLD) {
            maxLD = ld;
            maxLDIdx = i;
          }
        }
        
        if (maxLD > 0) {
          const x = toCanvasX(xValues[maxLDIdx]);
          const y = toCanvasY(yValues[maxLDIdx]);
          
          ctx.strokeStyle = '#ffcc00';
          ctx.lineWidth = 2;
          ctx.beginPath();
          ctx.arc(x, y, 8, 0, Math.PI * 2);
          ctx.stroke();
          
          ctx.fillStyle = '#ffcc00';
          ctx.font = '10px sans-serif';
          ctx.textAlign = 'left';
          ctx.fillText(`Max L/D = ${maxLD.toFixed(1)}`, x + 12, y);
        }
      }
    } else {
      // No data message
      ctx.fillStyle = 'rgba(255, 255, 255, 0.5)';
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No polar data. Run a polar sweep in the Solve panel.', width / 2, height / 2);
    }
    
    // Mode indicator
    ctx.fillStyle = hasCdData ? '#00d4aa' : '#6b9fff';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(hasCdData ? 'Viscous' : 'Inviscid', width - margin.right, margin.top - 5);
  }, [canvasSize, polarData, xAxis, yAxis, xBounds, yBounds, xValues, yValues, hasCdData]);
  
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
    
    // Initial size
    const rect = container.getBoundingClientRect();
    setCanvasSize({ width: Math.max(200, rect.width), height: Math.max(150, rect.height) });
    
    return () => observer.disconnect();
  }, []);
  
  // Redraw when data/size changes
  useEffect(() => {
    draw();
  }, [draw]);
  
  return (
    <div className="panel-content" style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      {/* Controls */}
      <div style={{ 
        display: 'flex', 
        gap: '12px', 
        padding: '8px',
        borderBottom: '1px solid var(--border-color)',
        alignItems: 'center',
        flexWrap: 'wrap',
      }}>
        <label style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '12px' }}>
          X:
          <select 
            value={xAxis}
            onChange={(e) => setXAxis(e.target.value as AxisVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            <option value="alpha">Alpha</option>
            <option value="cl">Cl</option>
            <option value="cd" disabled={!hasCdData}>Cd {!hasCdData ? '(viscous only)' : ''}</option>
            <option value="cm">Cm</option>
            <option value="ld" disabled={!hasCdData}>L/D {!hasCdData ? '(viscous only)' : ''}</option>
          </select>
        </label>
        
        <label style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '12px' }}>
          Y:
          <select 
            value={yAxis}
            onChange={(e) => setYAxis(e.target.value as AxisVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            <option value="alpha">Alpha</option>
            <option value="cl">Cl</option>
            <option value="cd" disabled={!hasCdData}>Cd {!hasCdData ? '(viscous only)' : ''}</option>
            <option value="cm">Cm</option>
            <option value="ld" disabled={!hasCdData}>L/D {!hasCdData ? '(viscous only)' : ''}</option>
          </select>
        </label>
        
        {/* Quick presets */}
        <div style={{ display: 'flex', gap: '4px' }}>
          <button
            onClick={() => { setXAxis('alpha'); setYAxis('cl'); }}
            style={{ padding: '2px 6px', fontSize: '10px' }}
          >
            Cl-α
          </button>
          {hasCdData && (
            <>
              <button
                onClick={() => { setXAxis('cd'); setYAxis('cl'); }}
                style={{ padding: '2px 6px', fontSize: '10px' }}
              >
                Drag Polar
              </button>
              <button
                onClick={() => { setXAxis('alpha'); setYAxis('ld'); }}
                style={{ padding: '2px 6px', fontSize: '10px' }}
              >
                L/D-α
              </button>
            </>
          )}
        </div>
        
        <span style={{ 
          marginLeft: 'auto', 
          fontSize: '11px', 
          color: 'var(--text-secondary)' 
        }}>
          {polarData.length} points
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
    </div>
  );
}
