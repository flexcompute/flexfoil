/**
 * SSPPlot - Interactive SVG plot for SSP spacing control
 */

import React, { useRef, useCallback } from 'react';
import type { SpacingKnot } from '../../types';
import { DraggableKnot } from './DraggableKnot';
import { OutputPoints } from './OutputPoints';

interface SSPPlotProps {
  knots: SpacingKnot[];
  si: number[];
  onKnotChange: (index: number, knot: SpacingKnot) => void;
  onAddKnot: (S: number, F: number) => void;
  onRemoveKnot: (index: number) => void;
}

// Plot dimensions and margins (viewBox coordinates)
const WIDTH = 500;
const HEIGHT = 320;
const MARGIN = { top: 25, right: 25, bottom: 45, left: 45 };
const PLOT_WIDTH = WIDTH - MARGIN.left - MARGIN.right;
const PLOT_HEIGHT = HEIGHT - MARGIN.top - MARGIN.bottom;

// Scale configuration
const S_MIN = 0;
const S_MAX = 1;
const F_MIN = 0;
const F_MAX = 2.5;

// Convert data coordinates to SVG coordinates
function toSvgX(s: number): number {
  return MARGIN.left + ((s - S_MIN) / (S_MAX - S_MIN)) * PLOT_WIDTH;
}

function toSvgY(f: number): number {
  return MARGIN.top + PLOT_HEIGHT - ((f - F_MIN) / (F_MAX - F_MIN)) * PLOT_HEIGHT;
}

// Convert SVG coordinates to data coordinates
function toDataS(svgX: number): number {
  return S_MIN + ((svgX - MARGIN.left) / PLOT_WIDTH) * (S_MAX - S_MIN);
}

function toDataF(svgY: number): number {
  return F_MIN + ((MARGIN.top + PLOT_HEIGHT - svgY) / PLOT_HEIGHT) * (F_MAX - F_MIN);
}

// Generate axis ticks
function generateTicks(min: number, max: number, count: number): number[] {
  const step = (max - min) / count;
  const ticks: number[] = [];
  for (let i = 0; i <= count; i++) {
    ticks.push(min + i * step);
  }
  return ticks;
}

// Helper to get SVG coordinates from client coordinates
function getSvgCoordsFromEvent(
  svg: SVGSVGElement,
  clientX: number,
  clientY: number
): { svgX: number; svgY: number } {
  const rect = svg.getBoundingClientRect();
  const viewBox = svg.viewBox.baseVal;
  const scaleX = viewBox.width / rect.width;
  const scaleY = viewBox.height / rect.height;
  return {
    svgX: (clientX - rect.left) * scaleX,
    svgY: (clientY - rect.top) * scaleY
  };
}

export const SSPPlot: React.FC<SSPPlotProps> = ({
  knots,
  si,
  onKnotChange,
  onAddKnot,
  onRemoveKnot
}) => {
  const svgRef = useRef<SVGSVGElement>(null);

  // Handle adding knot from coordinates
  const addKnotFromCoords = useCallback((clientX: number, clientY: number, isShift: boolean) => {
    if (isShift) return;
    
    const svg = svgRef.current;
    if (!svg) return;
    
    const { svgX, svgY } = getSvgCoordsFromEvent(svg, clientX, clientY);
    
    // Check if click is within plot area
    if (svgX < MARGIN.left || svgX > WIDTH - MARGIN.right ||
        svgY < MARGIN.top || svgY > HEIGHT - MARGIN.bottom) {
      return;
    }
    
    const S = Math.max(0.01, Math.min(0.99, toDataS(svgX)));
    const F = Math.max(0.01, Math.min(F_MAX - 0.1, toDataF(svgY)));
    
    onAddKnot(S, F);
  }, [onAddKnot]);

  // Handle click on plot area to add new knot
  const handlePlotClick = useCallback((e: React.MouseEvent<SVGRectElement>) => {
    addKnotFromCoords(e.clientX, e.clientY, e.shiftKey);
  }, [addKnotFromCoords]);

  // Handle touch on plot area to add new knot (double tap)
  const lastTapRef = useRef<number>(0);
  const handlePlotTouch = useCallback((e: React.TouchEvent<SVGRectElement>) => {
    const now = Date.now();
    const DOUBLE_TAP_DELAY = 300;
    
    if (now - lastTapRef.current < DOUBLE_TAP_DELAY) {
      e.preventDefault();
      const touch = e.touches[0] || e.changedTouches[0];
      addKnotFromCoords(touch.clientX, touch.clientY, false);
    }
    lastTapRef.current = now;
  }, [addKnotFromCoords]);

  // Generate grid lines
  const sTicks = generateTicks(S_MIN, S_MAX, 10);
  const fTicks = generateTicks(F_MIN, F_MAX, 5);

  // Build path for the piecewise-linear function line
  const linePath = knots.length >= 2
    ? `M ${knots.map(k => `${toSvgX(k.S)},${toSvgY(k.F)}`).join(' L ')}`
    : '';

  return (
    <svg
      ref={svgRef}
      viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
      style={{ 
        width: '100%',
        height: 'auto',
        minHeight: '200px',
        borderRadius: '4px',
        backgroundColor: 'var(--bg-primary, #1a1a2e)',
        touchAction: 'none'
      }}
    >
      {/* Grid lines */}
      <g className="grid" stroke="var(--border-color, #333)" strokeWidth="0.5" opacity="0.5">
        {/* Vertical grid lines (S axis) */}
        {sTicks.map(s => (
          <line
            key={`v-${s}`}
            x1={toSvgX(s)}
            y1={MARGIN.top}
            x2={toSvgX(s)}
            y2={HEIGHT - MARGIN.bottom}
          />
        ))}
        {/* Horizontal grid lines (F axis) */}
        {fTicks.map(f => (
          <line
            key={`h-${f}`}
            x1={MARGIN.left}
            y1={toSvgY(f)}
            x2={WIDTH - MARGIN.right}
            y2={toSvgY(f)}
          />
        ))}
      </g>

      {/* Axes */}
      <g className="axes" stroke="var(--text-secondary, #888)" strokeWidth="1.5">
        {/* X axis (S) */}
        <line
          x1={MARGIN.left}
          y1={HEIGHT - MARGIN.bottom}
          x2={WIDTH - MARGIN.right}
          y2={HEIGHT - MARGIN.bottom}
        />
        {/* Y axis (F) */}
        <line
          x1={MARGIN.left}
          y1={MARGIN.top}
          x2={MARGIN.left}
          y2={HEIGHT - MARGIN.bottom}
        />
      </g>

      {/* Axis labels */}
      <g className="axis-labels" fontSize="10" fill="var(--text-muted, #666)">
        {/* S axis ticks and labels */}
        {sTicks.filter((_, i) => i % 2 === 0).map(s => (
          <g key={`sl-${s}`}>
            <line
              x1={toSvgX(s)}
              y1={HEIGHT - MARGIN.bottom}
              x2={toSvgX(s)}
              y2={HEIGHT - MARGIN.bottom + 4}
              stroke="var(--text-muted, #666)"
            />
            <text
              x={toSvgX(s)}
              y={HEIGHT - MARGIN.bottom + 14}
              textAnchor="middle"
            >
              {s.toFixed(1)}
            </text>
          </g>
        ))}
        {/* F axis ticks and labels */}
        {fTicks.map(f => (
          <g key={`fl-${f}`}>
            <line
              x1={MARGIN.left - 4}
              y1={toSvgY(f)}
              x2={MARGIN.left}
              y2={toSvgY(f)}
              stroke="var(--text-muted, #666)"
            />
            <text
              x={MARGIN.left - 8}
              y={toSvgY(f) + 3}
              textAnchor="end"
            >
              {f.toFixed(1)}
            </text>
          </g>
        ))}
      </g>

      {/* Axis titles */}
      <text
        x={WIDTH / 2}
        y={HEIGHT - 6}
        textAnchor="middle"
        fontSize="11"
        fontWeight="600"
        fill="var(--text-secondary, #aaa)"
      >
        s (arc length)
      </text>
      <text
        x={12}
        y={HEIGHT / 2}
        textAnchor="middle"
        fontSize="11"
        fontWeight="600"
        fill="var(--text-secondary, #aaa)"
        transform={`rotate(-90, 12, ${HEIGHT / 2})`}
      >
        F (spacing)
      </text>

      {/* Clickable/tappable plot area for adding knots */}
      <rect
        x={MARGIN.left}
        y={MARGIN.top}
        width={PLOT_WIDTH}
        height={PLOT_HEIGHT}
        fill="transparent"
        onClick={handlePlotClick}
        onTouchStart={handlePlotTouch}
        style={{ cursor: 'crosshair' }}
      />

      {/* Piecewise-linear function line */}
      <path
        d={linePath}
        fill="none"
        stroke="var(--accent-primary, #0df)"
        strokeWidth="2"
        strokeLinejoin="round"
      />

      {/* Output points along s-axis */}
      <OutputPoints si={si} toSvgX={toSvgX} yPosition={HEIGHT - MARGIN.bottom + 22} />

      {/* Draggable knot points */}
      {knots.map((knot, index) => (
        <DraggableKnot
          key={index}
          index={index}
          knot={knot}
          isFirst={index === 0}
          isLast={index === knots.length - 1}
          prevS={index > 0 ? knots[index - 1].S : 0}
          nextS={index < knots.length - 1 ? knots[index + 1].S : 1}
          svgRef={svgRef}
          toSvgX={toSvgX}
          toSvgY={toSvgY}
          toDataS={toDataS}
          toDataF={toDataF}
          fMax={F_MAX}
          onKnotChange={onKnotChange}
          onRemoveKnot={onRemoveKnot}
          canRemove={knots.length > 2}
        />
      ))}
    </svg>
  );
};

export default SSPPlot;
