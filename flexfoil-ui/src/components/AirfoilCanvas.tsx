/**
 * AirfoilCanvas - Canvas 2D renderer for airfoil visualization
 * 
 * Renders:
 * - Airfoil curve (spline)
 * - Panel points
 * - Control points (surface, bezier handles, b-spline)
 * - Grid and axes
 * 
 * Supports:
 * - Pan and zoom
 * - Point dragging
 * - Point selection
 * - Show/hide toggles for points and grid
 */

import { useRef, useEffect, useCallback, useState } from 'react';
import { useAirfoilStore } from '../stores/airfoilStore';
import type { Point, ViewportState } from '../types';

// Constants
const POINT_RADIUS = 5;
const PANEL_POINT_RADIUS = 2.5;
const HANDLE_RADIUS = 4;
const CONTROL_RADIUS = 6;
const HIT_RADIUS = 10;

export function AirfoilCanvas() {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  
  // State from store
  const { 
    coordinates, 
    panels, 
    controlMode,
    bezierHandles,
    bsplineControlPoints,
    updatePoint,
    updateBSplineControlPoint,
  } = useAirfoilStore();

  // Viewport state
  const [viewport, setViewport] = useState<ViewportState>({
    center: { x: 0.5, y: 0 },
    zoom: 400,
    width: 800,
    height: 600,
  });

  // Display toggles
  const [showGrid, setShowGrid] = useState(true);
  const [showCurve, setShowCurve] = useState(true);
  const [showPoints, setShowPoints] = useState(true);
  const [showControls, setShowControls] = useState(true);

  // Interaction state
  const [isDragging, setIsDragging] = useState(false);
  const [isPanning, setIsPanning] = useState(false);
  const [dragTarget, setDragTarget] = useState<{ type: string; index: number | string } | null>(null);
  const [hoveredPoint, setHoveredPoint] = useState<{ type: string; index: number | string } | null>(null);
  const lastMousePos = useRef<Point>({ x: 0, y: 0 });

  // Convert airfoil coordinates to canvas coordinates
  const toCanvas = useCallback((p: Point): Point => {
    return {
      x: viewport.width / 2 + (p.x - viewport.center.x) * viewport.zoom,
      y: viewport.height / 2 - (p.y - viewport.center.y) * viewport.zoom,
    };
  }, [viewport]);

  // Convert canvas coordinates to airfoil coordinates
  const toAirfoil = useCallback((p: Point): Point => {
    return {
      x: viewport.center.x + (p.x - viewport.width / 2) / viewport.zoom,
      y: viewport.center.y - (p.y - viewport.height / 2) / viewport.zoom,
    };
  }, [viewport]);

  // Find point at canvas position
  const findPointAt = useCallback((canvasPos: Point): { type: string; index: number | string } | null => {
    // Check B-spline control points first (on top)
    if (controlMode === 'bspline') {
      for (const cp of bsplineControlPoints) {
        const cpCanvas = toCanvas(cp);
        const dist = Math.hypot(canvasPos.x - cpCanvas.x, canvasPos.y - cpCanvas.y);
        if (dist < HIT_RADIUS) {
          return { type: 'bspline', index: cp.id };
        }
      }
    }

    // Check bezier handles
    if (controlMode === 'bezier') {
      for (let i = 0; i < bezierHandles.length; i++) {
        const handle = bezierHandles[i];
        const hCanvas = toCanvas(handle.position);
        const dist = Math.hypot(canvasPos.x - hCanvas.x, canvasPos.y - hCanvas.y);
        if (dist < HIT_RADIUS) {
          return { type: 'bezier', index: i };
        }
      }
    }

    // Check surface points
    if (controlMode === 'surface') {
      for (let i = 0; i < coordinates.length; i++) {
        const pCanvas = toCanvas(coordinates[i]);
        const dist = Math.hypot(canvasPos.x - pCanvas.x, canvasPos.y - pCanvas.y);
        if (dist < HIT_RADIUS) {
          return { type: 'surface', index: i };
        }
      }
    }

    return null;
  }, [controlMode, coordinates, bezierHandles, bsplineControlPoints, toCanvas]);

  // Draw the canvas
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const { width, height } = viewport;

    // Clear canvas
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--bg-primary').trim() || '#0f0f0f';
    ctx.fillRect(0, 0, width, height);

    // Draw grid (if enabled)
    if (showGrid) {
      drawGrid(ctx);
      drawAxes(ctx);
    }

    // Draw B-spline control polygon (if in bspline mode and controls visible)
    if (showControls && controlMode === 'bspline' && bsplineControlPoints.length > 1) {
      ctx.beginPath();
      ctx.strokeStyle = 'rgba(0, 153, 255, 0.3)';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      const firstCP = toCanvas(bsplineControlPoints[0]);
      ctx.moveTo(firstCP.x, firstCP.y);
      for (let i = 1; i < bsplineControlPoints.length; i++) {
        const cp = toCanvas(bsplineControlPoints[i]);
        ctx.lineTo(cp.x, cp.y);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Draw airfoil curve (if enabled)
    if (showCurve && panels.length > 1) {
      ctx.beginPath();
      ctx.strokeStyle = getComputedStyle(document.documentElement)
        .getPropertyValue('--foil-line').trim() || '#00d4aa';
      ctx.lineWidth = 2;
      
      const first = toCanvas(panels[0]);
      ctx.moveTo(first.x, first.y);
      
      for (let i = 1; i < panels.length; i++) {
        const p = toCanvas(panels[i]);
        ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
    }

    // Draw bezier handles (if in bezier mode and controls visible)
    if (showControls && controlMode === 'bezier') {
      for (let i = 0; i < bezierHandles.length; i++) {
        const handle = bezierHandles[i];
        const hCanvas = toCanvas(handle.position);
        const pCanvas = toCanvas(coordinates[handle.pointIndex]);
        
        // Draw handle line
        ctx.beginPath();
        ctx.strokeStyle = 'rgba(255, 102, 102, 0.5)';
        ctx.lineWidth = 1;
        ctx.moveTo(pCanvas.x, pCanvas.y);
        ctx.lineTo(hCanvas.x, hCanvas.y);
        ctx.stroke();
        
        // Draw handle point
        const isHovered = hoveredPoint?.type === 'bezier' && hoveredPoint.index === i;
        ctx.beginPath();
        ctx.fillStyle = isHovered ? '#ffaa00' : '#ff6666';
        ctx.arc(hCanvas.x, hCanvas.y, HANDLE_RADIUS, 0, Math.PI * 2);
        ctx.fill();
      }
    }

    // Draw B-spline control points (if in bspline mode and controls visible)
    if (showControls && controlMode === 'bspline') {
      for (const cp of bsplineControlPoints) {
        const cpCanvas = toCanvas(cp);
        const isHovered = hoveredPoint?.type === 'bspline' && hoveredPoint.index === cp.id;
        
        ctx.beginPath();
        ctx.fillStyle = isHovered ? '#ffaa00' : '#0099ff';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2;
        ctx.arc(cpCanvas.x, cpCanvas.y, CONTROL_RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

    // Draw surface points (if in surface mode and controls visible)
    if (showControls && controlMode === 'surface') {
      for (let i = 0; i < coordinates.length; i++) {
        const p = toCanvas(coordinates[i]);
        const isHovered = hoveredPoint?.type === 'surface' && hoveredPoint.index === i;
        
        ctx.beginPath();
        ctx.fillStyle = isHovered ? '#ffaa00' : '#ffffff';
        ctx.arc(p.x, p.y, POINT_RADIUS, 0, Math.PI * 2);
        ctx.fill();
      }
    }

    // Draw panel points (small white dots to show spacing distribution)
    if (showPoints) {
      ctx.fillStyle = '#ffffff';
      ctx.strokeStyle = 'rgba(0, 0, 0, 0.6)';
      ctx.lineWidth = 1;
      for (const p of panels) {
        const pCanvas = toCanvas(p);
        ctx.beginPath();
        ctx.arc(pCanvas.x, pCanvas.y, PANEL_POINT_RADIUS, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    }

  }, [viewport, panels, coordinates, controlMode, bezierHandles, bsplineControlPoints, hoveredPoint, showGrid, showCurve, showPoints, showControls, toCanvas]);

  // Draw grid
  const drawGrid = useCallback((ctx: CanvasRenderingContext2D) => {
    const { width, height, zoom, center } = viewport;
    
    ctx.strokeStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--foil-grid').trim() || '#333333';
    ctx.lineWidth = 0.5;

    // Determine grid spacing based on zoom
    let gridSpacing = 0.1;
    if (zoom < 100) gridSpacing = 0.5;
    if (zoom < 50) gridSpacing = 1.0;
    if (zoom > 500) gridSpacing = 0.05;
    if (zoom > 1000) gridSpacing = 0.01;

    // Calculate visible range
    const left = center.x - width / (2 * zoom);
    const right = center.x + width / (2 * zoom);
    const top = center.y + height / (2 * zoom);
    const bottom = center.y - height / (2 * zoom);

    // Draw vertical lines
    const startX = Math.floor(left / gridSpacing) * gridSpacing;
    for (let x = startX; x <= right; x += gridSpacing) {
      const canvasX = toCanvas({ x, y: 0 }).x;
      ctx.beginPath();
      ctx.moveTo(canvasX, 0);
      ctx.lineTo(canvasX, height);
      ctx.stroke();
    }

    // Draw horizontal lines
    const startY = Math.floor(bottom / gridSpacing) * gridSpacing;
    for (let y = startY; y <= top; y += gridSpacing) {
      const canvasY = toCanvas({ x: 0, y }).y;
      ctx.beginPath();
      ctx.moveTo(0, canvasY);
      ctx.lineTo(width, canvasY);
      ctx.stroke();
    }
  }, [viewport, toCanvas]);

  // Draw axes
  const drawAxes = useCallback((ctx: CanvasRenderingContext2D) => {
    const { width, height } = viewport;
    const origin = toCanvas({ x: 0, y: 0 });

    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
    ctx.lineWidth = 1;

    // X axis
    ctx.beginPath();
    ctx.moveTo(0, origin.y);
    ctx.lineTo(width, origin.y);
    ctx.stroke();

    // Y axis
    ctx.beginPath();
    ctx.moveTo(origin.x, 0);
    ctx.lineTo(origin.x, height);
    ctx.stroke();
  }, [viewport, toCanvas]);

  // Handle resize
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const resizeObserver = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        setViewport((v) => ({ ...v, width, height }));
      }
    });

    resizeObserver.observe(container);
    
    // Initial size
    const rect = container.getBoundingClientRect();
    setViewport((v) => ({ ...v, width: rect.width, height: rect.height }));

    return () => resizeObserver.disconnect();
  }, []);

  // Redraw when state changes
  useEffect(() => {
    draw();
  }, [draw]);

  // Mouse event handlers
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    lastMousePos.current = canvasPos;

    // Check for point hit
    const hit = findPointAt(canvasPos);
    if (hit) {
      setIsDragging(true);
      setDragTarget(hit);
    } else {
      // Start panning
      setIsPanning(true);
    }
  }, [findPointAt]);

  const handleMouseMove = useCallback((e: React.MouseEvent) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilPos = toAirfoil(canvasPos);

    if (isDragging && dragTarget) {
      // Update point position
      if (dragTarget.type === 'surface' && typeof dragTarget.index === 'number') {
        updatePoint(dragTarget.index, { x: airfoilPos.x, y: airfoilPos.y });
      } else if (dragTarget.type === 'bspline' && typeof dragTarget.index === 'string') {
        updateBSplineControlPoint(dragTarget.index, { x: airfoilPos.x, y: airfoilPos.y });
      }
    } else if (isPanning) {
      // Pan viewport
      const dx = (canvasPos.x - lastMousePos.current.x) / viewport.zoom;
      const dy = (canvasPos.y - lastMousePos.current.y) / viewport.zoom;
      setViewport((v) => ({
        ...v,
        center: { x: v.center.x - dx, y: v.center.y + dy },
      }));
    } else {
      // Update hover state
      const hit = findPointAt(canvasPos);
      setHoveredPoint(hit);
    }

    lastMousePos.current = canvasPos;
  }, [isDragging, isPanning, dragTarget, viewport.zoom, toAirfoil, findPointAt, updatePoint, updateBSplineControlPoint]);

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
    setIsPanning(false);
    setDragTarget(null);
  }, []);

  const handleWheel = useCallback((e: React.WheelEvent) => {
    e.preventDefault();
    
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const canvasPos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    const airfoilBefore = toAirfoil(canvasPos);

    // Zoom
    const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
    const newZoom = Math.max(50, Math.min(5000, viewport.zoom * zoomFactor));

    // Adjust center to zoom toward mouse position
    setViewport((v) => {
      const newViewport = { ...v, zoom: newZoom };
      const airfoilAfter = {
        x: v.center.x + (canvasPos.x - v.width / 2) / newZoom,
        y: v.center.y - (canvasPos.y - v.height / 2) / newZoom,
      };
      return {
        ...newViewport,
        center: {
          x: v.center.x + (airfoilBefore.x - airfoilAfter.x),
          y: v.center.y + (airfoilBefore.y - airfoilAfter.y),
        },
      };
    });
  }, [viewport, toAirfoil]);

  // Reset view
  const handleResetView = useCallback(() => {
    setViewport((v) => ({
      ...v,
      center: { x: 0.5, y: 0 },
      zoom: Math.min(v.width, v.height) * 0.8,
    }));
  }, []);

  return (
    <div 
      ref={containerRef} 
      className="canvas-container"
      style={{ position: 'relative' }}
    >
      <canvas
        ref={canvasRef}
        width={viewport.width}
        height={viewport.height}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
        style={{ cursor: isDragging ? 'grabbing' : isPanning ? 'grabbing' : hoveredPoint ? 'pointer' : 'crosshair' }}
      />
      
      {/* Overlay controls */}
      <div
        style={{
          position: 'absolute',
          bottom: '12px',
          right: '12px',
          display: 'flex',
          gap: '8px',
          alignItems: 'center',
          zIndex: 10,
          pointerEvents: 'auto',
          background: 'rgba(15, 15, 15, 0.8)',
          padding: '6px 10px',
          borderRadius: '6px',
          border: '1px solid var(--border-color)',
        }}
      >
        {/* Display toggles */}
        <label style={{ 
          display: 'flex', 
          alignItems: 'center', 
          gap: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
          cursor: 'pointer',
        }}>
          <input 
            type="checkbox" 
            checked={showGrid} 
            onChange={(e) => setShowGrid(e.target.checked)}
            style={{ width: '12px', height: '12px' }}
          />
          Grid
        </label>
        <label style={{ 
          display: 'flex', 
          alignItems: 'center', 
          gap: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
          cursor: 'pointer',
        }}>
          <input 
            type="checkbox" 
            checked={showCurve} 
            onChange={(e) => setShowCurve(e.target.checked)}
            style={{ width: '12px', height: '12px' }}
          />
          Curve
        </label>
        <label style={{ 
          display: 'flex', 
          alignItems: 'center', 
          gap: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
          cursor: 'pointer',
        }}>
          <input 
            type="checkbox" 
            checked={showPoints} 
            onChange={(e) => setShowPoints(e.target.checked)}
            style={{ width: '12px', height: '12px' }}
          />
          Points
        </label>
        <label style={{ 
          display: 'flex', 
          alignItems: 'center', 
          gap: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
          cursor: 'pointer',
        }}>
          <input 
            type="checkbox" 
            checked={showControls} 
            onChange={(e) => setShowControls(e.target.checked)}
            style={{ width: '12px', height: '12px' }}
          />
          Controls
        </label>
        
        <div style={{ width: '1px', height: '16px', background: 'var(--border-color)' }} />
        
        <button onClick={handleResetView} style={{ padding: '4px 8px', fontSize: '11px' }}>
          Reset View
        </button>
        <span style={{ 
          background: 'var(--bg-tertiary)', 
          padding: '4px 8px', 
          borderRadius: '4px',
          fontSize: '11px',
          color: 'var(--text-secondary)',
        }}>
          Zoom: {viewport.zoom.toFixed(0)}x
        </span>
      </div>
    </div>
  );
}
