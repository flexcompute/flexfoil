/**
 * BoundaryLayerPanel - Boundary layer visualization
 * 
 * Displays:
 * - Momentum thickness (θ)
 * - Displacement thickness (δ*)
 * - Shape factor (H / Hk)
 * - Skin friction (Cf)
 * - Edge velocity (Ue)
 * - Amplification factor (N)
 * - Transition location markers
 * - Wake region (extends x-axis beyond x/c = 1)
 */

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { getBLVisualizationData, isWasmReady, type BLVisualizationData } from '../../lib/wasm';

type PlotVariable = 'theta' | 'delta_star' | 'h' | 'hk' | 'cf' | 'ue' | 'ampl';

const VARIABLE_LABELS: Record<PlotVariable, string> = {
  theta: 'θ (momentum thickness)',
  delta_star: 'δ* (displacement thickness)',
  h: 'H (shape factor)',
  hk: 'Hk (kinematic shape factor)',
  cf: 'Cf (skin friction)',
  ue: 'Ue (edge velocity)',
  ampl: 'N (amplification factor)',
};

const VARIABLE_COLORS = {
  upper: '#00d4aa',
  lower: '#ff6b9d',
  wake: '#ffc832',
};

function getVariableData(
  data: BLVisualizationData,
  variable: PlotVariable,
  surface: 'upper' | 'lower',
): { x: number[]; y: number[] } {
  const side = data[surface];
  const x = side.x;
  let y: number[];

  switch (variable) {
    case 'theta': y = side.theta; break;
    case 'delta_star': y = side.delta_star; break;
    case 'h':
      y = side.delta_star.map((ds, i) => side.theta[i] > 0 ? ds / side.theta[i] : 0);
      break;
    case 'hk': y = side.hk; break;
    case 'cf': y = side.cf; break;
    case 'ue': y = side.ue; break;
    case 'ampl': y = side.ampl; break;
  }

  return { x, y };
}

function getWakeData(
  data: BLVisualizationData,
  variable: PlotVariable,
): { x: number[]; y: number[] } {
  const wake = data.wake;
  if (wake.x.length === 0) return { x: [], y: [] };

  switch (variable) {
    case 'theta': return { x: wake.x, y: wake.theta };
    case 'delta_star': return { x: wake.x, y: wake.delta_star };
    default: return { x: [], y: [] };
  }
}

export function BoundaryLayerPanel() {
  const { panels, displayAlpha, reynolds, mach, ncrit, maxIterations } = useAirfoilStore();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  const [plotVariable, setPlotVariable] = useState<PlotVariable>('h');
  const [showWakeRegion, setShowWakeRegion] = useState(false);
  const [canvasSize, setCanvasSize] = useState({ width: 400, height: 300 });
  const [blData, setBlData] = useState<BLVisualizationData | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  
  const margin = { top: 20, right: 20, bottom: 40, left: 60 };

  useEffect(() => {
    if (!isWasmReady() || panels.length < 3) {
      setBlData(null);
      return;
    }

    setIsLoading(true);
    
    requestAnimationFrame(() => {
      try {
        const data = getBLVisualizationData(panels, displayAlpha, reynolds, mach, ncrit, maxIterations);
        setBlData(data.success ? data : null);
      } catch (e) {
        console.error('BL computation failed:', e);
        setBlData(null);
      }
      setIsLoading(false);
    });
  }, [panels, displayAlpha, reynolds, mach, ncrit, maxIterations]);

  const hasWakeData = useMemo(() => {
    return blData !== null && blData.wake.x.length > 0;
  }, [blData]);

  const bounds = useMemo(() => {
    if (!blData) return { xMin: 0, xMax: 1, yMin: 0, yMax: 1 };

    const upperData = getVariableData(blData, plotVariable, 'upper');
    const lowerData = getVariableData(blData, plotVariable, 'lower');
    const wakeData = showWakeRegion ? getWakeData(blData, plotVariable) : { x: [], y: [] };
    
    const allX = [...upperData.x, ...lowerData.x, ...wakeData.x];
    const allY = [...upperData.y, ...lowerData.y, ...wakeData.y];
    
    if (allX.length === 0) {
      return { xMin: 0, xMax: 1, yMin: 0, yMax: 1 };
    }

    const yMin = Math.min(...allY, 0);
    const yMax = Math.max(...allY);
    const xMax = showWakeRegion && wakeData.x.length > 0
      ? Math.max(1, ...wakeData.x)
      : 1;
    
    const yPadding = (yMax - yMin) * 0.1 || 0.1;
    
    return {
      xMin: 0,
      xMax,
      yMin: yMin - yPadding,
      yMax: yMax + yPadding,
    };
  }, [blData, plotVariable, showWakeRegion]);

  const generateTicks = (min: number, max: number, count: number = 5): number[] => {
    const range = max - min;
    if (range <= 0) return [min];
    
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
  };

  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    
    const { width, height } = canvasSize;
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--bg-primary').trim() || '#0f0f0f';
    ctx.fillRect(0, 0, width, height);
    
    const toCanvasX = (x: number) => {
      return margin.left + ((x - bounds.xMin) / (bounds.xMax - bounds.xMin)) * plotWidth;
    };
    const toCanvasY = (y: number) => {
      return margin.top + (1 - (y - bounds.yMin) / (bounds.yMax - bounds.yMin)) * plotHeight;
    };
    
    // Grid
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.1)';
    ctx.lineWidth = 1;
    
    const xTicks = generateTicks(bounds.xMin, bounds.xMax);
    const yTicks = generateTicks(bounds.yMin, bounds.yMax);
    
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

    // Wake region shading (x > 1)
    if (showWakeRegion && bounds.xMax > 1) {
      const wakeStart = toCanvasX(1);
      const wakeEnd = toCanvasX(bounds.xMax);
      ctx.fillStyle = 'rgba(255, 200, 50, 0.05)';
      ctx.fillRect(wakeStart, margin.top, wakeEnd - wakeStart, plotHeight);
      
      // TE marker line
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(wakeStart, margin.top);
      ctx.lineTo(wakeStart, height - margin.bottom);
      ctx.stroke();
      ctx.setLineDash([]);
      
      ctx.fillStyle = 'rgba(255, 200, 50, 0.5)';
      ctx.font = '9px sans-serif';
      ctx.textAlign = 'left';
      ctx.fillText('wake', wakeStart + 3, margin.top + 10);
    }
    
    // Axes
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
    ctx.lineWidth = 1;
    
    ctx.beginPath();
    ctx.moveTo(margin.left, height - margin.bottom);
    ctx.lineTo(width - margin.right, height - margin.bottom);
    ctx.stroke();
    
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, height - margin.bottom);
    ctx.stroke();
    
    // Tick labels
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    
    for (const x of xTicks) {
      const cx = toCanvasX(x);
      ctx.fillText(x.toFixed(2), cx, height - margin.bottom + 15);
    }
    
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      const useExponential = plotVariable !== 'h' && plotVariable !== 'hk' && plotVariable !== 'ue' && plotVariable !== 'ampl';
      const yStr = useExponential ? y.toExponential(1) : y.toFixed(1);
      ctx.fillText(yStr, margin.left - 5, cy);
    }
    
    // Axis labels
    ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('x/c', width / 2, height - 15);
    
    ctx.save();
    ctx.translate(15, height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'bottom';
    ctx.fillText(VARIABLE_LABELS[plotVariable], 0, 0);
    ctx.restore();
    
    // Draw data
    if (blData) {
      const drawLine = (data: { x: number[]; y: number[] }, color: string, dashed = false) => {
        if (data.x.length === 0) return;
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        if (dashed) ctx.setLineDash([6, 3]);
        ctx.beginPath();
        for (let i = 0; i < data.x.length; i++) {
          const cx = toCanvasX(data.x[i]);
          const cy = toCanvasY(data.y[i]);
          if (i === 0) ctx.moveTo(cx, cy);
          else ctx.lineTo(cx, cy);
        }
        ctx.stroke();
        if (dashed) ctx.setLineDash([]);
      };

      drawLine(getVariableData(blData, plotVariable, 'upper'), VARIABLE_COLORS.upper);
      drawLine(getVariableData(blData, plotVariable, 'lower'), VARIABLE_COLORS.lower);

      if (showWakeRegion) {
        const wakeData = getWakeData(blData, plotVariable);
        drawLine(wakeData, VARIABLE_COLORS.wake, true);
      }
      
      // Transition markers
      if (blData.x_tr_upper < 1) {
        const xTr = toCanvasX(blData.x_tr_upper);
        ctx.strokeStyle = VARIABLE_COLORS.upper;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(xTr, margin.top);
        ctx.lineTo(xTr, height - margin.bottom);
        ctx.stroke();
        ctx.setLineDash([]);
        
        ctx.fillStyle = VARIABLE_COLORS.upper;
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(`Tr U: ${(blData.x_tr_upper * 100).toFixed(0)}%`, xTr, margin.top - 5);
      }
      
      if (blData.x_tr_lower < 1) {
        const xTr = toCanvasX(blData.x_tr_lower);
        ctx.strokeStyle = VARIABLE_COLORS.lower;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(xTr, margin.top);
        ctx.lineTo(xTr, height - margin.bottom);
        ctx.stroke();
        ctx.setLineDash([]);
        
        ctx.fillStyle = VARIABLE_COLORS.lower;
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(`Tr L: ${(blData.x_tr_lower * 100).toFixed(0)}%`, xTr, margin.top + 10);
      }
    } else if (isLoading) {
      ctx.fillStyle = 'rgba(255, 255, 255, 0.5)';
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('Computing...', width / 2, height / 2);
    }
    
    // Legend
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    
    ctx.fillStyle = VARIABLE_COLORS.upper;
    ctx.fillRect(width - 90, 10, 12, 12);
    ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
    ctx.fillText('Upper', width - 74, 20);
    
    ctx.fillStyle = VARIABLE_COLORS.lower;
    ctx.fillRect(width - 90, 26, 12, 12);
    ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
    ctx.fillText('Lower', width - 74, 36);

    if (showWakeRegion && hasWakeData) {
      ctx.fillStyle = VARIABLE_COLORS.wake;
      ctx.fillRect(width - 90, 42, 12, 12);
      ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
      ctx.fillText('Wake', width - 74, 52);
    }
  }, [canvasSize, blData, plotVariable, bounds, isLoading, showWakeRegion, hasWakeData]);
  
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
          Variable:
          <select 
            value={plotVariable}
            onChange={(e) => setPlotVariable(e.target.value as PlotVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            <option value="h">Shape Factor H</option>
            <option value="hk">Kinematic Hk</option>
            <option value="theta">Momentum θ</option>
            <option value="delta_star">Displacement δ*</option>
            <option value="cf">Skin Friction Cf</option>
            <option value="ue">Edge Velocity Ue</option>
            <option value="ampl">Amplification N</option>
          </select>
        </label>

        {hasWakeData && (
          <label style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '11px', color: 'var(--text-secondary)' }}>
            <input
              type="checkbox"
              checked={showWakeRegion}
              onChange={(e) => setShowWakeRegion(e.target.checked)}
              style={{ width: '12px', height: '12px' }}
            />
            Wake
          </label>
        )}
        
        {blData && (
          <span style={{ 
            marginLeft: 'auto', 
            fontSize: '11px', 
            color: 'var(--text-secondary)' 
          }}>
            α={displayAlpha.toFixed(1)}° | Re={(reynolds / 1e6).toFixed(2)}M
            {mach > 0 ? ` | M=${mach.toFixed(2)}` : ''}
          </span>
        )}
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
