/**
 * PolarPanel - Aerodynamic polar plot visualization
 * 
 * Displays polar data (Cl vs alpha, Cl vs Cm, etc.) with configurable axes.
 */

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import type { PolarPoint } from '../../types';

type AxisVariable = 'alpha' | 'cl' | 'cd' | 'cm';

const AXIS_LABELS: Record<AxisVariable, string> = {
  alpha: 'α (deg)',
  cl: 'Cl',
  cd: 'Cd',
  cm: 'Cm',
};

// Get value from polar point by variable name
function getValue(point: PolarPoint, variable: AxisVariable): number {
  return point[variable] ?? 0;
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
  const { polarData } = useAirfoilStore();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  const [xAxis, setXAxis] = useState<AxisVariable>('alpha');
  const [yAxis, setYAxis] = useState<AxisVariable>('cl');
  const [canvasSize, setCanvasSize] = useState({ width: 400, height: 300 });
  
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
    ctx.strokeStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--foil-grid').trim() || '#333333';
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
    ctx.strokeStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--text-secondary').trim() || '#888888';
    ctx.globalAlpha = 0.3;
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
    
    ctx.globalAlpha = 1; // Reset alpha
    
    // Draw tick labels
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--text-secondary').trim() || '#888888';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    
    for (const x of xTicks) {
      const cx = toCanvasX(x);
      ctx.fillText(x.toFixed(1), cx, height - margin.bottom + 15);
    }
    
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      ctx.fillText(y.toFixed(2), margin.left - 5, cy);
    }
    
    // Draw axis labels
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
    
    // Draw data points and line
    if (polarData.length > 0) {
      // Line
      ctx.strokeStyle = getComputedStyle(document.documentElement)
        .getPropertyValue('--foil-line').trim() || '#00d4aa';
      ctx.globalAlpha = 0.8;
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
      
      ctx.globalAlpha = 1; // Reset alpha
      
      // Points
      ctx.fillStyle = getComputedStyle(document.documentElement)
        .getPropertyValue('--foil-line').trim() || '#00d4aa';
      for (let i = 0; i < polarData.length; i++) {
        const x = toCanvasX(xValues[i]);
        const y = toCanvasY(yValues[i]);
        ctx.beginPath();
        ctx.arc(x, y, 4, 0, Math.PI * 2);
        ctx.fill();
      }
    } else {
      // No data message
      ctx.fillStyle = getComputedStyle(document.documentElement)
        .getPropertyValue('--text-secondary').trim() || '#888888';
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('No polar data. Run a polar sweep in the Solve panel.', width / 2, height / 2);
    }
  }, [canvasSize, polarData, xAxis, yAxis, xBounds, yBounds, xValues, yValues]);
  
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
            <option value="cd">Cd</option>
            <option value="cm">Cm</option>
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
            <option value="cd">Cd</option>
            <option value="cm">Cm</option>
          </select>
        </label>
        
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
