/**
 * FoilSpacingPlot - Interactive SVG plot showing SSP knots around the airfoil
 * 
 * Instead of the traditional F vs S plot, this shows:
 * - The airfoil shape
 * - Knots positioned at their arc-length position (S) on the foil
 * - F value represented as normal displacement from the surface
 */

import React, { useRef, useCallback, useMemo } from 'react';
import type { SpacingKnot, AirfoilPoint } from '../../types';
import { pauseHistory, resumeHistory } from '../../stores/airfoilStore';

interface FoilSpacingPlotProps {
  knots: SpacingKnot[];
  coordinates: AirfoilPoint[];
  onKnotChange: (index: number, knot: SpacingKnot) => void;
  onAddKnot: (S: number, F: number) => void;
  onRemoveKnot: (index: number) => void;
}

// Plot dimensions
const WIDTH = 500;
const HEIGHT = 250;
const MARGIN = { top: 20, right: 20, bottom: 20, left: 20 };

// Scale for F (normal displacement) - how much the F value extends from surface
const F_SCALE = 0.15; // F=1 extends 0.15 chord lengths from surface
const F_MIN = 0.01;
const F_MAX = 2.5;

/**
 * Compute cumulative arc lengths for coordinates
 */
function computeArcLengths(coords: AirfoilPoint[]): number[] {
  if (coords.length === 0) return [];
  const s: number[] = [0];
  for (let i = 1; i < coords.length; i++) {
    const dx = coords[i].x - coords[i - 1].x;
    const dy = coords[i].y - coords[i - 1].y;
    s.push(s[i - 1] + Math.sqrt(dx * dx + dy * dy));
  }
  return s;
}

/**
 * Compute surface normal at each point (outward pointing)
 */
function computeNormals(coords: AirfoilPoint[]): { nx: number; ny: number }[] {
  if (coords.length < 2) return [];
  
  const normals: { nx: number; ny: number }[] = [];
  
  for (let i = 0; i < coords.length; i++) {
    let tx: number, ty: number;
    
    if (i === 0) {
      // Forward difference
      tx = coords[1].x - coords[0].x;
      ty = coords[1].y - coords[0].y;
    } else if (i === coords.length - 1) {
      // Backward difference
      tx = coords[i].x - coords[i - 1].x;
      ty = coords[i].y - coords[i - 1].y;
    } else {
      // Central difference
      tx = coords[i + 1].x - coords[i - 1].x;
      ty = coords[i + 1].y - coords[i - 1].y;
    }
    
    // Normalize tangent
    const len = Math.sqrt(tx * tx + ty * ty);
    if (len > 1e-10) {
      tx /= len;
      ty /= len;
    }
    
    // Normal is perpendicular to tangent (rotate 90 degrees CCW for outward)
    // For airfoil going TE->upper->LE->lower->TE, outward is -90 rotation
    normals.push({ nx: ty, ny: -tx });
  }
  
  return normals;
}

/**
 * Interpolate position and normal at arc-length fraction S
 */
function interpolateAtS(
  coords: AirfoilPoint[],
  arcLengths: number[],
  normals: { nx: number; ny: number }[],
  S: number
): { x: number; y: number; nx: number; ny: number } {
  if (coords.length === 0) {
    return { x: 0.5, y: 0, nx: 0, ny: 1 };
  }
  
  const totalLength = arcLengths[arcLengths.length - 1];
  const targetS = S * totalLength;
  
  // Find segment
  let i = 0;
  for (let j = 0; j < arcLengths.length - 1; j++) {
    if (targetS >= arcLengths[j] && targetS <= arcLengths[j + 1]) {
      i = j;
      break;
    }
    i = j;
  }
  
  // Interpolate within segment
  const segLen = arcLengths[i + 1] - arcLengths[i];
  const t = segLen > 1e-10 ? (targetS - arcLengths[i]) / segLen : 0;
  
  return {
    x: coords[i].x + t * (coords[i + 1].x - coords[i].x),
    y: coords[i].y + t * (coords[i + 1].y - coords[i].y),
    nx: normals[i].nx + t * (normals[i + 1].nx - normals[i].nx),
    ny: normals[i].ny + t * (normals[i + 1].ny - normals[i].ny),
  };
}

export const FoilSpacingPlot: React.FC<FoilSpacingPlotProps> = ({
  knots,
  coordinates,
  onKnotChange,
  onAddKnot,
  onRemoveKnot,
}) => {
  const svgRef = useRef<SVGSVGElement>(null);
  const draggingRef = useRef<{ index: number; isEndpoint: boolean } | null>(null);

  // Precompute arc lengths and normals
  const { arcLengths, normals } = useMemo(() => {
    const al = computeArcLengths(coordinates);
    const n = computeNormals(coordinates);
    return { 
      arcLengths: al, 
      normals: n, 
    };
  }, [coordinates]);

  // Compute view bounds to fit airfoil with some padding
  const viewBounds = useMemo(() => {
    if (coordinates.length === 0) {
      return { minX: -0.1, maxX: 1.1, minY: -0.3, maxY: 0.3 };
    }
    
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    
    for (const p of coordinates) {
      minX = Math.min(minX, p.x);
      maxX = Math.max(maxX, p.x);
      minY = Math.min(minY, p.y);
      maxY = Math.max(maxY, p.y);
    }
    
    // Add padding for knot visualization
    const padX = (maxX - minX) * 0.15;
    const padY = Math.max((maxY - minY) * 0.3, 0.2);
    
    return {
      minX: minX - padX,
      maxX: maxX + padX,
      minY: minY - padY,
      maxY: maxY + padY,
    };
  }, [coordinates]);

  // Convert airfoil coordinates to SVG coordinates
  const toSvgX = useCallback((x: number): number => {
    const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
    return MARGIN.left + ((x - viewBounds.minX) / (viewBounds.maxX - viewBounds.minX)) * plotWidth;
  }, [viewBounds]);

  const toSvgY = useCallback((y: number): number => {
    const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
    return MARGIN.top + plotHeight - ((y - viewBounds.minY) / (viewBounds.maxY - viewBounds.minY)) * plotHeight;
  }, [viewBounds]);

  // Convert SVG to airfoil coordinates
  const toAirfoilX = useCallback((svgX: number): number => {
    const plotWidth = WIDTH - MARGIN.left - MARGIN.right;
    return viewBounds.minX + ((svgX - MARGIN.left) / plotWidth) * (viewBounds.maxX - viewBounds.minX);
  }, [viewBounds]);

  const toAirfoilY = useCallback((svgY: number): number => {
    const plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
    return viewBounds.minY + ((MARGIN.top + plotHeight - svgY) / plotHeight) * (viewBounds.maxY - viewBounds.minY);
  }, [viewBounds]);

  // Get SVG coordinates from mouse/touch event
  const getSvgCoords = useCallback((clientX: number, clientY: number): { x: number; y: number } => {
    const svg = svgRef.current;
    if (!svg) return { x: 0, y: 0 };
    
    const rect = svg.getBoundingClientRect();
    const viewBox = svg.viewBox.baseVal;
    const scaleX = viewBox.width / rect.width;
    const scaleY = viewBox.height / rect.height;
    
    return {
      x: (clientX - rect.left) * scaleX,
      y: (clientY - rect.top) * scaleY,
    };
  }, []);

  // Find S value from airfoil position (closest point on curve)
  const findSFromPosition = useCallback((px: number, py: number): number => {
    if (coordinates.length < 2) return 0.5;
    
    let minDist = Infinity;
    let bestS = 0.5;
    
    // Search along the curve
    for (let i = 0; i <= 100; i++) {
      const s = i / 100;
      const pos = interpolateAtS(coordinates, arcLengths, normals, s);
      const dist = Math.hypot(px - pos.x, py - pos.y);
      if (dist < minDist) {
        minDist = dist;
        bestS = s;
      }
    }
    
    return bestS;
  }, [coordinates, arcLengths, normals]);

  // Handle mouse down on knot
  const handleKnotMouseDown = useCallback((e: React.MouseEvent, index: number, isEndpoint: boolean) => {
    e.stopPropagation();
    e.preventDefault();
    
    if (e.shiftKey || e.button === 2) {
      // Remove knot
      onRemoveKnot(index);
      return;
    }
    
    pauseHistory();
    draggingRef.current = { index, isEndpoint };
    
    const handleMouseMove = (moveEvent: MouseEvent) => {
      if (!draggingRef.current) return;
      
      const svgCoords = getSvgCoords(moveEvent.clientX, moveEvent.clientY);
      const airfoilX = toAirfoilX(svgCoords.x);
      const airfoilY = toAirfoilY(svgCoords.y);
      
      const { index: idx, isEndpoint: isEnd } = draggingRef.current;
      const currentKnot = knots[idx];
      
      // Get position on foil at current S
      const pos = interpolateAtS(coordinates, arcLengths, normals, currentKnot.S);
      
      // Calculate new F based on distance from foil along normal
      const dx = airfoilX - pos.x;
      const dy = airfoilY - pos.y;
      const normalDist = dx * pos.nx + dy * pos.ny;
      const newF = Math.max(F_MIN, Math.min(F_MAX, normalDist / F_SCALE));
      
      if (isEnd) {
        // Endpoints can only change F (stay at S=0 or S=1)
        onKnotChange(idx, { S: currentKnot.S, F: newF });
      } else {
        // Interior knots can move along surface
        let newS = findSFromPosition(airfoilX, airfoilY);
        
        // Constrain to stay between neighbors
        const prevS = idx > 0 ? knots[idx - 1].S + 0.01 : 0.01;
        const nextS = idx < knots.length - 1 ? knots[idx + 1].S - 0.01 : 0.99;
        newS = Math.max(prevS, Math.min(nextS, newS));
        
        onKnotChange(idx, { S: newS, F: newF });
      }
    };
    
    const handleMouseUp = () => {
      draggingRef.current = null;
      resumeHistory();
      window.removeEventListener('mousemove', handleMouseMove);
      window.removeEventListener('mouseup', handleMouseUp);
    };
    
    window.addEventListener('mousemove', handleMouseMove);
    window.addEventListener('mouseup', handleMouseUp);
  }, [knots, coordinates, arcLengths, normals, onKnotChange, onRemoveKnot, getSvgCoords, toAirfoilX, toAirfoilY, findSFromPosition]);

  // Handle click on plot area to add knot
  const handlePlotClick = useCallback((e: React.MouseEvent<SVGPathElement>) => {
    if (e.shiftKey) return;
    
    const svgCoords = getSvgCoords(e.clientX, e.clientY);
    const airfoilX = toAirfoilX(svgCoords.x);
    const airfoilY = toAirfoilY(svgCoords.y);
    
    // Find closest S position
    const S = findSFromPosition(airfoilX, airfoilY);
    
    // Check not too close to existing knots
    const MIN_GAP = 0.02;
    if (knots.some(k => Math.abs(k.S - S) < MIN_GAP)) return;
    
    // Calculate F from distance to surface
    const pos = interpolateAtS(coordinates, arcLengths, normals, S);
    const dx = airfoilX - pos.x;
    const dy = airfoilY - pos.y;
    const normalDist = dx * pos.nx + dy * pos.ny;
    const F = Math.max(F_MIN, Math.min(F_MAX, Math.abs(normalDist) / F_SCALE + 0.5));
    
    onAddKnot(S, F);
  }, [knots, coordinates, arcLengths, normals, onAddKnot, getSvgCoords, toAirfoilX, toAirfoilY, findSFromPosition]);

  // Build airfoil path
  const foilPath = useMemo(() => {
    if (coordinates.length < 2) return '';
    return `M ${coordinates.map(p => `${toSvgX(p.x)},${toSvgY(p.y)}`).join(' L ')}`;
  }, [coordinates, toSvgX, toSvgY]);

  // Compute knot positions with normal offsets
  const knotPositions = useMemo(() => {
    return knots.map(knot => {
      const pos = interpolateAtS(coordinates, arcLengths, normals, knot.S);
      const offset = knot.F * F_SCALE;
      return {
        // Base position on surface
        baseX: toSvgX(pos.x),
        baseY: toSvgY(pos.y),
        // Offset position showing F value
        knotX: toSvgX(pos.x + pos.nx * offset),
        knotY: toSvgY(pos.y + pos.ny * offset),
        // Normal direction in SVG space
        nx: pos.nx,
        ny: pos.ny,
      };
    });
  }, [knots, coordinates, arcLengths, normals, toSvgX, toSvgY]);

  // Build the piecewise linear "ribbon" showing spacing function
  const ribbonPath = useMemo(() => {
    if (knots.length < 2) return '';
    
    // Sample along the curve and interpolate F
    const samples: { x: number; y: number }[] = [];
    const numSamples = 100;
    
    for (let i = 0; i <= numSamples; i++) {
      const s = i / numSamples;
      
      // Linear interpolation of F between knots
      let F = 1;
      for (let j = 0; j < knots.length - 1; j++) {
        if (s >= knots[j].S && s <= knots[j + 1].S) {
          const t = (s - knots[j].S) / (knots[j + 1].S - knots[j].S);
          F = knots[j].F + t * (knots[j + 1].F - knots[j].F);
          break;
        }
      }
      
      const pos = interpolateAtS(coordinates, arcLengths, normals, s);
      const offset = F * F_SCALE;
      
      samples.push({
        x: toSvgX(pos.x + pos.nx * offset),
        y: toSvgY(pos.y + pos.ny * offset),
      });
    }
    
    return `M ${samples.map(p => `${p.x},${p.y}`).join(' L ')}`;
  }, [knots, coordinates, arcLengths, normals, toSvgX, toSvgY]);

  return (
    <svg
      ref={svgRef}
      viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
      style={{ 
        width: '100%',
        height: 'auto',
        minHeight: '180px',
        borderRadius: '4px',
        backgroundColor: 'var(--bg-primary, #1a1a2e)',
        touchAction: 'none',
      }}
      onContextMenu={(e) => e.preventDefault()}
    >
      {/* Click area for adding knots */}
      <rect
        x={MARGIN.left}
        y={MARGIN.top}
        width={WIDTH - MARGIN.left - MARGIN.right}
        height={HEIGHT - MARGIN.top - MARGIN.bottom}
        fill="transparent"
        style={{ cursor: 'crosshair' }}
      />

      {/* Airfoil shape */}
      <path
        d={foilPath}
        fill="none"
        stroke="var(--text-muted, #666)"
        strokeWidth="2"
        onClick={handlePlotClick}
        style={{ cursor: 'crosshair' }}
      />

      {/* Spacing function ribbon */}
      <path
        d={ribbonPath}
        fill="none"
        stroke="var(--accent-primary, #0df)"
        strokeWidth="2"
        strokeLinejoin="round"
        opacity="0.8"
      />

      {/* Knot stems (lines from surface to knot) */}
      {knotPositions.map((pos, i) => (
        <line
          key={`stem-${i}`}
          x1={pos.baseX}
          y1={pos.baseY}
          x2={pos.knotX}
          y2={pos.knotY}
          stroke="var(--accent-primary, #0df)"
          strokeWidth="1"
          strokeDasharray="3,2"
          opacity="0.6"
        />
      ))}

      {/* Base markers on surface */}
      {knotPositions.map((pos, i) => (
        <circle
          key={`base-${i}`}
          cx={pos.baseX}
          cy={pos.baseY}
          r={3}
          fill="var(--text-muted, #666)"
        />
      ))}

      {/* Draggable knot circles */}
      {knotPositions.map((pos, i) => {
        const isEndpoint = i === 0 || i === knots.length - 1;
        return (
          <g key={`knot-${i}`}>
            {/* Hit area */}
            <circle
              cx={pos.knotX}
              cy={pos.knotY}
              r={12}
              fill="transparent"
              style={{ cursor: 'grab' }}
              onMouseDown={(e) => handleKnotMouseDown(e, i, isEndpoint)}
            />
            {/* Visible knot */}
            <circle
              cx={pos.knotX}
              cy={pos.knotY}
              r={6}
              fill={isEndpoint ? 'var(--accent-success, #0f8)' : 'var(--accent-primary, #0df)'}
              stroke="var(--bg-primary, #1a1a2e)"
              strokeWidth="2"
              style={{ cursor: 'grab' }}
              onMouseDown={(e) => handleKnotMouseDown(e, i, isEndpoint)}
            />
            {/* F value label */}
            <text
              x={pos.knotX + 10}
              y={pos.knotY + 4}
              fill="var(--text-muted, #888)"
              fontSize="9"
              fontFamily="var(--font-mono)"
            >
              F={knots[i].F.toFixed(1)}
            </text>
          </g>
        );
      })}

      {/* Legend */}
      <text
        x={WIDTH - MARGIN.right - 5}
        y={MARGIN.top + 12}
        fill="var(--text-muted, #666)"
        fontSize="9"
        textAnchor="end"
      >
        F = normal displacement
      </text>
      <text
        x={WIDTH - MARGIN.right - 5}
        y={MARGIN.top + 24}
        fill="var(--text-muted, #666)"
        fontSize="9"
        textAnchor="end"
      >
        Higher F = denser mesh
      </text>
    </svg>
  );
};

export default FoilSpacingPlot;
