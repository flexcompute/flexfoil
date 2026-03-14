/**
 * BoundaryLayerPanel - Boundary layer visualization
 * 
 * Displays:
 * - Momentum thickness (θ)
 * - Displacement thickness (δ*)
 * - Shape factor (H)
 * - Skin friction (Cf)
 * - Transition location markers
 */

import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { getBLDistribution, isWasmReady, type BLDistribution } from '../../lib/wasm';

type PlotVariable = 'theta' | 'delta_star' | 'h' | 'cf';

const VARIABLE_LABELS: Record<PlotVariable, string> = {
  theta: 'θ (momentum thickness)',
  delta_star: 'δ* (displacement thickness)',
  h: 'H (shape factor)',
  cf: 'Cf (skin friction)',
};

const VARIABLE_COLORS = {
  upper: '#00d4aa',
  lower: '#ff6b9d',
};

export function BoundaryLayerPanel() {
  const { panels, displayAlpha, reynolds } = useAirfoilStore();
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  const [plotVariable, setPlotVariable] = useState<PlotVariable>('h');
  const [canvasSize, setCanvasSize] = useState({ width: 400, height: 300 });
  const [blData, setBLData] = useState<BLDistribution | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  
  // Margins for axis labels
  const margin = { top: 20, right: 20, bottom: 40, left: 60 };

  // Fetch BL data when inputs change
  useEffect(() => {
    if (!isWasmReady() || panels.length < 3) {
      setBLData(null);
      return;
    }

    setIsLoading(true);
    
    // Use requestAnimationFrame to avoid blocking
    requestAnimationFrame(() => {
      try {
        const data = getBLDistribution(panels, displayAlpha, reynolds);
        setBLData(data);
      } catch (e) {
        console.error('BL computation failed:', e);
        setBLData(null);
      }
      setIsLoading(false);
    });
  }, [panels, displayAlpha, reynolds]);

  // Get data for current variable
  const getData = useCallback((variable: PlotVariable, surface: 'upper' | 'lower'): { x: number[]; y: number[] } => {
    if (!blData) return { x: [], y: [] };

    const x = surface === 'upper' ? blData.x_upper : blData.x_lower;
    let y: number[];

    switch (variable) {
      case 'theta':
        y = surface === 'upper' ? blData.theta_upper : blData.theta_lower;
        break;
      case 'delta_star':
        y = surface === 'upper' ? blData.delta_star_upper : blData.delta_star_lower;
        break;
      case 'h':
        y = surface === 'upper' ? blData.h_upper : blData.h_lower;
        break;
      case 'cf':
        y = surface === 'upper' ? blData.cf_upper : blData.cf_lower;
        break;
    }

    return { x, y };
  }, [blData]);

  // Compute axis bounds
  const bounds = useMemo(() => {
    const upperData = getData(plotVariable, 'upper');
    const lowerData = getData(plotVariable, 'lower');
    
    const allX = [...upperData.x, ...lowerData.x];
    const allY = [...upperData.y, ...lowerData.y];
    
    if (allX.length === 0) {
      return { xMin: 0, xMax: 1, yMin: 0, yMax: 1 };
    }

    const yMin = Math.min(...allY, 0);
    const yMax = Math.max(...allY);
    
    const yPadding = (yMax - yMin) * 0.1 || 0.1;
    
    return {
      xMin: 0,
      xMax: 1,
      yMin: yMin - yPadding,
      yMax: yMax + yPadding,
    };
  }, [getData, plotVariable]);

  // Generate nice tick values
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
      return margin.left + ((x - bounds.xMin) / (bounds.xMax - bounds.xMin)) * plotWidth;
    };
    const toCanvasY = (y: number) => {
      return margin.top + (1 - (y - bounds.yMin) / (bounds.yMax - bounds.yMin)) * plotHeight;
    };
    
    // Draw grid
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
      ctx.fillText(x.toFixed(2), cx, height - margin.bottom + 15);
    }
    
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const y of yTicks) {
      const cy = toCanvasY(y);
      const yStr = plotVariable === 'h' ? y.toFixed(1) : y.toExponential(1);
      ctx.fillText(yStr, margin.left - 5, cy);
    }
    
    // Draw axis labels
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
      // Draw upper surface
      const upperData = getData(plotVariable, 'upper');
      if (upperData.x.length > 0) {
        ctx.strokeStyle = VARIABLE_COLORS.upper;
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < upperData.x.length; i++) {
          const cx = toCanvasX(upperData.x[i]);
          const cy = toCanvasY(upperData.y[i]);
          if (i === 0) ctx.moveTo(cx, cy);
          else ctx.lineTo(cx, cy);
        }
        ctx.stroke();
      }
      
      // Draw lower surface
      const lowerData = getData(plotVariable, 'lower');
      if (lowerData.x.length > 0) {
        ctx.strokeStyle = VARIABLE_COLORS.lower;
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < lowerData.x.length; i++) {
          const cx = toCanvasX(lowerData.x[i]);
          const cy = toCanvasY(lowerData.y[i]);
          if (i === 0) ctx.moveTo(cx, cy);
          else ctx.lineTo(cx, cy);
        }
        ctx.stroke();
      }
      
      // Draw transition markers
      if (blData.x_tr_upper < 1) {
        const xTr = toCanvasX(blData.x_tr_upper);
        ctx.strokeStyle = VARIABLE_COLORS.upper;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(xTr, margin.top);
        ctx.lineTo(xTr, height - margin.bottom);
        ctx.stroke();
        ctx.setLineDash([]);
        
        // Label
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
        
        // Label
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
  }, [canvasSize, blData, plotVariable, bounds, getData, isLoading]);
  
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
          Variable:
          <select 
            value={plotVariable}
            onChange={(e) => setPlotVariable(e.target.value as PlotVariable)}
            style={{ padding: '4px', fontSize: '12px' }}
          >
            <option value="h">Shape Factor H</option>
            <option value="theta">Momentum θ</option>
            <option value="delta_star">Displacement δ*</option>
            <option value="cf">Skin Friction Cf</option>
          </select>
        </label>
        
        {blData && (
          <span style={{ 
            marginLeft: 'auto', 
            fontSize: '11px', 
            color: 'var(--text-secondary)' 
          }}>
            α={displayAlpha.toFixed(1)}° | Re={(reynolds / 1e6).toFixed(2)}M
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
