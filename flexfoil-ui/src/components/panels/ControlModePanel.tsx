/**
 * ControlModePanel - Switch between surface, bezier, and b-spline modes
 */

import { useCallback } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import type { ControlMode, BSplineControlPoint } from '../../types';

export function ControlModePanel() {
  const { 
    controlMode, 
    setControlMode,
    coordinates,
    bsplineControlPoints,
    setBSplineControlPoints,
    bsplineDegree,
    setBSplineDegree,
    addBSplineControlPoint,
  } = useAirfoilStore();

  const handleModeChange = useCallback((mode: ControlMode) => {
    setControlMode(mode);
    
    // Initialize B-spline control points if switching to bspline mode
    if (mode === 'bspline' && bsplineControlPoints.length === 0) {
      // Create initial control points from airfoil shape
      // Use a simplified set of control points
      const numCP = 8;
      const newControlPoints: BSplineControlPoint[] = [];
      
      for (let i = 0; i < numCP; i++) {
        const t = i / (numCP - 1);
        const idx = Math.floor(t * (coordinates.length - 1));
        const p = coordinates[idx];
        
        // Offset control points slightly to show they're off-surface
        const offset = 0.02;
        newControlPoints.push({
          id: `cp-${i}`,
          x: p.x,
          y: p.y + (p.y >= 0 ? offset : -offset),
          weight: 1,
        });
      }
      
      setBSplineControlPoints(newControlPoints);
    }
  }, [setControlMode, coordinates, bsplineControlPoints.length, setBSplineControlPoints]);

  const handleAddControlPoint = useCallback(() => {
    const id = `cp-${Date.now()}`;
    // Add at center
    addBSplineControlPoint({
      id,
      x: 0.5,
      y: 0.1,
      weight: 1,
    });
  }, [addBSplineControlPoint]);

  return (
    <div className="panel">
      <div className="panel-header">Control Mode</div>
      <div className="panel-content">
        {/* Mode selection */}
        <div className="form-group">
          <div className="form-label">Edit Mode</div>
          <div className="control-mode-group">
            <button
              className={`control-mode-btn ${controlMode === 'surface' ? 'active' : ''}`}
              onClick={() => handleModeChange('surface')}
            >
              Surface
            </button>
            <button
              className={`control-mode-btn ${controlMode === 'bezier' ? 'active' : ''}`}
              onClick={() => handleModeChange('bezier')}
            >
              Bezier
            </button>
            <button
              className={`control-mode-btn ${controlMode === 'bspline' ? 'active' : ''}`}
              onClick={() => handleModeChange('bspline')}
            >
              B-Spline
            </button>
          </div>
        </div>

        {/* Mode description */}
        <div style={{ 
          padding: '8px 12px', 
          background: 'var(--bg-tertiary)', 
          borderRadius: '4px',
          fontSize: '12px',
          color: 'var(--text-secondary)',
          marginBottom: '12px',
        }}>
          {controlMode === 'surface' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>Surface Points</strong>
              <p style={{ margin: '4px 0 0' }}>
                Drag points directly on the airfoil surface. 
                Best for fine adjustments and manual edits.
              </p>
            </>
          )}
          {controlMode === 'bezier' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>Bezier Handles</strong>
              <p style={{ margin: '4px 0 0' }}>
                Control tangent directions at each point.
                Useful for smooth curve shaping.
              </p>
            </>
          )}
          {controlMode === 'bspline' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>B-Spline Control</strong>
              <p style={{ margin: '4px 0 0' }}>
                Off-surface control polygon. The curve is attracted 
                to control points but doesn't pass through them.
              </p>
            </>
          )}
        </div>

        {/* B-spline specific controls */}
        {controlMode === 'bspline' && (
          <>
            <div className="form-group">
              <div className="form-label">Degree</div>
              <div className="form-row">
                <input
                  type="range"
                  min={1}
                  max={5}
                  value={bsplineDegree}
                  onChange={(e) => setBSplineDegree(parseInt(e.target.value))}
                  style={{ flex: 1 }}
                />
                <span style={{ 
                  width: '30px', 
                  textAlign: 'right',
                  fontFamily: 'var(--font-mono)',
                }}>
                  {bsplineDegree}
                </span>
              </div>
            </div>

            <div className="form-group">
              <div className="form-label">
                Control Points ({bsplineControlPoints.length})
              </div>
              <button onClick={handleAddControlPoint} style={{ width: '100%' }}>
                + Add Control Point
              </button>
            </div>

            <div style={{ 
              maxHeight: '200px', 
              overflow: 'auto',
              border: '1px solid var(--border-color)',
              borderRadius: '4px',
            }}>
              {bsplineControlPoints.map((cp, i) => (
                <div
                  key={cp.id}
                  style={{
                    padding: '6px 8px',
                    borderBottom: '1px solid var(--border-color)',
                    fontSize: '11px',
                    fontFamily: 'var(--font-mono)',
                    display: 'flex',
                    justifyContent: 'space-between',
                  }}
                >
                  <span style={{ color: 'var(--accent-secondary)' }}>CP{i}</span>
                  <span>x: {cp.x.toFixed(3)}</span>
                  <span>y: {cp.y.toFixed(3)}</span>
                </div>
              ))}
            </div>
          </>
        )}

        {/* Instructions */}
        <div style={{ 
          marginTop: '12px',
          padding: '8px',
          background: 'var(--bg-primary)',
          borderRadius: '4px',
          fontSize: '11px',
          color: 'var(--text-muted)',
        }}>
          <div><strong>Controls:</strong></div>
          <div>• Click + drag to move points</div>
          <div>• Scroll to zoom</div>
          <div>• Drag canvas to pan</div>
        </div>
      </div>
    </div>
  );
}
