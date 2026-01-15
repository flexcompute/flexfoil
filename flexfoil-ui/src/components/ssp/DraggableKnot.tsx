/**
 * DraggableKnot - An interactive knot point for the SSP plot
 */

import React, { useState, useCallback, useEffect, useRef } from 'react';
import type { SpacingKnot } from '../../types';

interface DraggableKnotProps {
  index: number;
  knot: SpacingKnot;
  isFirst: boolean;
  isLast: boolean;
  prevS: number;
  nextS: number;
  svgRef: React.RefObject<SVGSVGElement | null>;
  toSvgX: (s: number) => number;
  toSvgY: (f: number) => number;
  toDataS: (svgX: number) => number;
  toDataF: (svgY: number) => number;
  fMax: number;
  onKnotChange: (index: number, knot: SpacingKnot) => void;
  onRemoveKnot: (index: number) => void;
  canRemove: boolean;
}

const KNOT_RADIUS = 8;
const MIN_S_GAP = 0.02;
const MIN_F = 0.01;
const LONG_PRESS_DURATION = 500;

export const DraggableKnot: React.FC<DraggableKnotProps> = ({
  index,
  knot,
  isFirst,
  isLast,
  prevS,
  nextS,
  svgRef,
  toSvgX,
  toSvgY,
  toDataS,
  toDataF,
  fMax,
  onKnotChange,
  onRemoveKnot,
  canRemove
}) => {
  const [isDragging, setIsDragging] = useState(false);
  const [isHovered, setIsHovered] = useState(false);
  const longPressTimerRef = useRef<number | null>(null);
  const hasMoved = useRef(false);

  // Helper to get SVG coordinates from client coordinates
  const getSvgCoords = useCallback((clientX: number, clientY: number) => {
    if (!svgRef.current) return { svgX: 0, svgY: 0 };
    const rect = svgRef.current.getBoundingClientRect();
    const svg = svgRef.current;
    const viewBox = svg.viewBox.baseVal;
    const scaleX = viewBox.width / rect.width;
    const scaleY = viewBox.height / rect.height;
    return {
      svgX: (clientX - rect.left) * scaleX,
      svgY: (clientY - rect.top) * scaleY
    };
  }, [svgRef]);

  // Clear long press timer
  const clearLongPress = useCallback(() => {
    if (longPressTimerRef.current !== null) {
      clearTimeout(longPressTimerRef.current);
      longPressTimerRef.current = null;
    }
  }, []);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    // Shift+click to remove knot
    if (e.shiftKey && canRemove && !isFirst && !isLast) {
      e.preventDefault();
      e.stopPropagation();
      onRemoveKnot(index);
      return;
    }
    
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, [index, canRemove, isFirst, isLast, onRemoveKnot]);

  // Touch start - begin drag and start long press timer
  const handleTouchStart = useCallback((e: React.TouchEvent) => {
    e.preventDefault();
    e.stopPropagation();
    hasMoved.current = false;
    setIsDragging(true);
    
    // Start long press timer for removal
    if (canRemove && !isFirst && !isLast) {
      longPressTimerRef.current = window.setTimeout(() => {
        if (!hasMoved.current) {
          setIsDragging(false);
          onRemoveKnot(index);
        }
      }, LONG_PRESS_DURATION);
    }
  }, [index, canRemove, isFirst, isLast, onRemoveKnot]);

  const handleMouseMove = useCallback((e: MouseEvent) => {
    if (!isDragging || !svgRef.current) return;

    const { svgX, svgY } = getSvgCoords(e.clientX, e.clientY);

    let newS = toDataS(svgX);
    let newF = toDataF(svgY);

    // Apply constraints

    // First and last knots can only move vertically
    if (isFirst) {
      newS = 0;
    } else if (isLast) {
      newS = 1;
    } else {
      // Enforce monotonicity: S must be between neighbors with gap
      newS = Math.max(prevS + MIN_S_GAP, Math.min(nextS - MIN_S_GAP, newS));
    }

    // F must be positive
    newF = Math.max(MIN_F, Math.min(fMax - 0.1, newF));

    onKnotChange(index, { S: newS, F: newF });
  }, [isDragging, svgRef, getSvgCoords, toDataS, toDataF, isFirst, isLast, prevS, nextS, fMax, index, onKnotChange]);

  // Touch move handler
  const handleTouchMove = useCallback((e: TouchEvent) => {
    if (!isDragging || !svgRef.current) return;
    
    hasMoved.current = true;
    clearLongPress();
    
    const touch = e.touches[0];
    const { svgX, svgY } = getSvgCoords(touch.clientX, touch.clientY);

    let newS = toDataS(svgX);
    let newF = toDataF(svgY);

    // Apply constraints
    if (isFirst) {
      newS = 0;
    } else if (isLast) {
      newS = 1;
    } else {
      newS = Math.max(prevS + MIN_S_GAP, Math.min(nextS - MIN_S_GAP, newS));
    }

    newF = Math.max(MIN_F, Math.min(fMax - 0.1, newF));

    onKnotChange(index, { S: newS, F: newF });
  }, [isDragging, svgRef, getSvgCoords, toDataS, toDataF, isFirst, isLast, prevS, nextS, fMax, index, onKnotChange, clearLongPress]);

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
  }, []);

  const handleTouchEnd = useCallback(() => {
    setIsDragging(false);
    clearLongPress();
  }, [clearLongPress]);

  // Add global event listeners for drag
  useEffect(() => {
    if (isDragging) {
      window.addEventListener('mousemove', handleMouseMove);
      window.addEventListener('mouseup', handleMouseUp);
      window.addEventListener('touchmove', handleTouchMove, { passive: false });
      window.addEventListener('touchend', handleTouchEnd);
      window.addEventListener('touchcancel', handleTouchEnd);
      return () => {
        window.removeEventListener('mousemove', handleMouseMove);
        window.removeEventListener('mouseup', handleMouseUp);
        window.removeEventListener('touchmove', handleTouchMove);
        window.removeEventListener('touchend', handleTouchEnd);
        window.removeEventListener('touchcancel', handleTouchEnd);
      };
    }
  }, [isDragging, handleMouseMove, handleMouseUp, handleTouchMove, handleTouchEnd]);

  // Cleanup long press timer on unmount
  useEffect(() => {
    return () => clearLongPress();
  }, [clearLongPress]);

  // Handle right-click to remove
  const handleContextMenu = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    if (canRemove && !isFirst && !isLast) {
      onRemoveKnot(index);
    }
  }, [index, canRemove, isFirst, isLast, onRemoveKnot]);

  const cx = toSvgX(knot.S);
  const cy = toSvgY(knot.F);

  // Determine cursor based on constraints
  let cursor = 'grab';
  if (isDragging) {
    cursor = 'grabbing';
  } else if (isFirst || isLast) {
    cursor = 'ns-resize';
  }

  // Color coding
  let fillColor = 'var(--accent-primary, #2196F3)';
  if (isFirst || isLast) {
    fillColor = 'var(--accent-success, #4CAF50)';
  }
  if (isHovered) {
    fillColor = isDragging ? '#1565C0' : '#42A5F5';
    if (isFirst || isLast) {
      fillColor = isDragging ? '#2E7D32' : '#66BB6A';
    }
  }

  return (
    <g>
      {/* Larger invisible hit area */}
      <circle
        cx={cx}
        cy={cy}
        r={KNOT_RADIUS + 12}
        fill="transparent"
        style={{ cursor, touchAction: 'none' }}
        onMouseDown={handleMouseDown}
        onTouchStart={handleTouchStart}
        onMouseEnter={() => setIsHovered(true)}
        onMouseLeave={() => setIsHovered(false)}
        onContextMenu={handleContextMenu}
      />
      {/* Visible knot circle */}
      <circle
        cx={cx}
        cy={cy}
        r={isDragging ? KNOT_RADIUS + 2 : KNOT_RADIUS}
        fill={fillColor}
        stroke={isDragging ? '#0D47A1' : 'var(--accent-primary, #1976D2)'}
        strokeWidth={isDragging ? 3 : 2}
        style={{ 
          cursor,
          transition: isDragging ? 'none' : 'all 0.1s ease',
          pointerEvents: 'none'
        }}
      />
      {/* Knot label on hover */}
      {isHovered && (
        <text
          x={cx}
          y={cy - KNOT_RADIUS - 8}
          textAnchor="middle"
          fontSize="11"
          fill="var(--text-primary, #333)"
          fontWeight="500"
        >
          ({knot.S.toFixed(2)}, {knot.F.toFixed(2)})
        </text>
      )}
    </g>
  );
};

export default DraggableKnot;
