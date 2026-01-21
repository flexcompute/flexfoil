/**
 * ControlModePanel - Switch between parameters, camber-spline, and thickness-spline modes
 */

import { useCallback, useEffect, useMemo, useRef } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useVisualizationStore } from '../../stores/visualizationStore';
import type { ControlMode } from '../../types';

// Debounce delay for repaneling after slider changes (ms)
const REPANEL_DEBOUNCE_MS = 150;

export function ControlModePanel() {
  const { 
    controlMode, 
    setControlMode,
    coordinates,
    baseCoordinates,
    // Parameters mode
    thicknessScale,
    camberScale,
    setThicknessScale,
    setCamberScale,
    applyScaling,
    repanelCoordinates,
    // Camber spline mode
    camberControlPoints,
    addCamberControlPoint,
    initializeCamberControlPoints,
    // Thickness spline mode
    thicknessControlPoints,
    addThicknessControlPoint,
    initializeThicknessControlPoints,
  } = useAirfoilStore();
  
  // Check if the base airfoil is symmetric (no camber to scale)
  // A symmetric airfoil has max camber < 0.1% of chord
  const isSymmetric = useMemo(() => {
    if (baseCoordinates.length < 10) return false;
    
    // Quick check: find max |y| among upper surface points (x > 0.1)
    // For symmetric airfoils, the mean y at each x should be ~0
    let maxCamber = 0;
    const midChord = baseCoordinates.filter(p => p.x > 0.2 && p.x < 0.8);
    if (midChord.length < 4) return false;
    
    // Group by approximate x and check if mean y is near zero
    for (const p of midChord) {
      // Find corresponding point on opposite surface (approximate)
      const opposite = midChord.find(q => Math.abs(q.x - p.x) < 0.05 && q !== p);
      if (opposite) {
        const meanY = (p.y + opposite.y) / 2;
        maxCamber = Math.max(maxCamber, Math.abs(meanY));
      }
    }
    
    return maxCamber < 0.001; // Less than 0.1% chord
  }, [baseCoordinates]);

  const { showControls, setShowControls } = useVisualizationStore();

  // Auto-enable control point visibility when in spline modes
  useEffect(() => {
    if ((controlMode === 'camber-spline' || controlMode === 'thickness-spline') && !showControls) {
      setShowControls(true);
    }
  }, [controlMode, showControls, setShowControls]);

  const handleModeChange = useCallback((mode: ControlMode) => {
    setControlMode(mode);
    
    // Auto-enable control point visibility for spline modes
    if (mode === 'camber-spline' || mode === 'thickness-spline') {
      setShowControls(true);
    }
    
    // Initialize camber control points if switching to camber-spline mode
    if (mode === 'camber-spline' && camberControlPoints.length === 0 && coordinates.length > 0) {
      initializeCamberControlPoints();
    }
    
    // Initialize thickness control points if switching to thickness-spline mode
    if (mode === 'thickness-spline' && thicknessControlPoints.length === 0 && coordinates.length > 0) {
      initializeThicknessControlPoints();
    }
  }, [
    setControlMode, 
    setShowControls,
    coordinates, 
    camberControlPoints.length, 
    thicknessControlPoints.length,
    initializeCamberControlPoints,
    initializeThicknessControlPoints,
  ]);

  // Debounced repaneling for smooth slider interaction
  const repanelTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  
  const debouncedRepanel = useCallback(() => {
    if (repanelTimeoutRef.current) {
      clearTimeout(repanelTimeoutRef.current);
    }
    repanelTimeoutRef.current = setTimeout(() => {
      repanelCoordinates();
    }, REPANEL_DEBOUNCE_MS);
  }, [repanelCoordinates]);
  
  // Cleanup timeout on unmount
  useEffect(() => {
    return () => {
      if (repanelTimeoutRef.current) {
        clearTimeout(repanelTimeoutRef.current);
      }
    };
  }, []);
  
  // Real-time slider handlers - update shape immediately, debounce repaneling
  const handleThicknessChange = useCallback((value: number) => {
    setThicknessScale(value, true); // skipRepanel=true for immediate visual update
    debouncedRepanel();
  }, [setThicknessScale, debouncedRepanel]);
  
  const handleCamberChange = useCallback((value: number) => {
    setCamberScale(value, true); // skipRepanel=true for immediate visual update
    debouncedRepanel();
  }, [setCamberScale, debouncedRepanel]);

  const handleAddCamberPoint = useCallback(() => {
    const id = `camber-${Date.now()}`;
    // Add at center with zero camber
    addCamberControlPoint({
      id,
      x: 0.5,
      y: 0,
    });
  }, [addCamberControlPoint]);

  const handleAddThicknessPoint = useCallback(() => {
    const id = `thickness-${Date.now()}`;
    // Add at center with average thickness
    addThicknessControlPoint({
      id,
      x: 0.5,
      t: 0.06, // Typical half-thickness
    });
  }, [addThicknessControlPoint]);

  return (
    <div className="panel">
      <div className="panel-header">Control Mode</div>
      <div className="panel-content">
        {/* Mode selection */}
        <div className="form-group">
          <div className="form-label">Edit Mode</div>
          <div className="control-mode-group">
            <button
              className={`control-mode-btn ${controlMode === 'parameters' ? 'active' : ''}`}
              onClick={() => handleModeChange('parameters')}
              data-tour="control-mode-parameters"
            >
              Parameters
            </button>
            <button
              className={`control-mode-btn ${controlMode === 'camber-spline' ? 'active' : ''}`}
              onClick={() => handleModeChange('camber-spline')}
              data-tour="control-mode-camber"
            >
              Camber
            </button>
            <button
              className={`control-mode-btn ${controlMode === 'thickness-spline' ? 'active' : ''}`}
              onClick={() => handleModeChange('thickness-spline')}
              data-tour="control-mode-thickness"
            >
              Thickness
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
          {controlMode === 'parameters' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>Parameter Scaling</strong>
              <p style={{ margin: '4px 0 0' }}>
                Adjust thickness and camber using sliders.
                Scales the current airfoil proportionally.
              </p>
            </>
          )}
          {controlMode === 'camber-spline' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>Camber Line Editor</strong>
              <p style={{ margin: '4px 0 0' }}>
                Drag control points to reshape the camber line.
                Click to add, shift+click to remove points.
              </p>
            </>
          )}
          {controlMode === 'thickness-spline' && (
            <>
              <strong style={{ color: 'var(--text-primary)' }}>Thickness Editor</strong>
              <p style={{ margin: '4px 0 0' }}>
                Drag control points to reshape the thickness distribution.
                Click to add, shift+click to remove points.
              </p>
            </>
          )}
        </div>

        {/* Parameters mode controls */}
        {controlMode === 'parameters' && (
          <>
            <div className="form-group" data-tour="thickness-slider">
              <div className="form-label">Thickness Scale</div>
              <div className="form-row">
                <input
                  type="range"
                  min={0.1}
                  max={3.0}
                  step={0.01}
                  value={thicknessScale}
                  onChange={(e) => handleThicknessChange(parseFloat(e.target.value))}
                  style={{ flex: 1 }}
                />
                <span style={{ 
                  width: '50px', 
                  textAlign: 'right',
                  fontFamily: 'var(--font-mono)',
                }}>
                  {(thicknessScale * 100).toFixed(0)}%
                </span>
              </div>
            </div>

            <div className="form-group" data-tour="camber-slider">
              <div className="form-label">Camber Scale</div>
              <div className="form-row">
                <input
                  type="range"
                  min={0}
                  max={3.0}
                  step={0.01}
                  value={camberScale}
                  onChange={(e) => handleCamberChange(parseFloat(e.target.value))}
                  style={{ flex: 1, opacity: isSymmetric ? 0.5 : 1 }}
                  disabled={isSymmetric}
                />
                <span style={{ 
                  width: '50px', 
                  textAlign: 'right',
                  fontFamily: 'var(--font-mono)',
                  opacity: isSymmetric ? 0.5 : 1,
                }}>
                  {(camberScale * 100).toFixed(0)}%
                </span>
              </div>
              {isSymmetric && (
                <p style={{ 
                  fontSize: '11px', 
                  color: 'var(--text-muted)', 
                  marginTop: '4px',
                  fontStyle: 'italic',
                }}>
                  Symmetric airfoil has no camber to scale. Try a cambered airfoil like 2412.
                </p>
              )}
            </div>

            <button 
              onClick={applyScaling} 
              style={{ width: '100%', marginTop: '8px' }}
              disabled={thicknessScale === 1.0 && camberScale === 1.0}
            >
              Apply & Reset Sliders
            </button>
            
            <p style={{ 
              fontSize: '11px', 
              color: 'var(--text-muted)', 
              marginTop: '8px',
              textAlign: 'center',
            }}>
              Apply saves the scaled shape as the new base.
            </p>
          </>
        )}

        {/* Camber spline mode controls */}
        {controlMode === 'camber-spline' && (
          <>
            <div className="form-group">
              <div className="form-label">
                Camber Points ({camberControlPoints.length})
              </div>
              <button onClick={handleAddCamberPoint} style={{ width: '100%' }}>
                + Add Control Point
              </button>
            </div>

            <div style={{ 
              maxHeight: '200px', 
              overflow: 'auto',
              border: '1px solid var(--border-color)',
              borderRadius: '4px',
            }}>
              {camberControlPoints.map((cp, i) => (
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
                  <span style={{ color: 'var(--accent-primary)' }}>C{i}</span>
                  <span>x: {cp.x.toFixed(3)}</span>
                  <span>y: {(cp.y * 100).toFixed(2)}%</span>
                </div>
              ))}
            </div>
          </>
        )}

        {/* Thickness spline mode controls */}
        {controlMode === 'thickness-spline' && (
          <>
            <div className="form-group">
              <div className="form-label">
                Thickness Points ({thicknessControlPoints.length})
              </div>
              <button onClick={handleAddThicknessPoint} style={{ width: '100%' }}>
                + Add Control Point
              </button>
            </div>

            <div style={{ 
              maxHeight: '200px', 
              overflow: 'auto',
              border: '1px solid var(--border-color)',
              borderRadius: '4px',
            }}>
              {thicknessControlPoints.map((cp, i) => (
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
                  <span style={{ color: 'var(--accent-warning)' }}>T{i}</span>
                  <span>x: {cp.x.toFixed(3)}</span>
                  <span>t: {(cp.t * 100).toFixed(2)}%</span>
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
          {controlMode === 'parameters' ? (
            <>
              <div>• Drag sliders to scale</div>
              <div>• Apply to save changes</div>
            </>
          ) : (
            <>
              <div>• Click + drag to move points</div>
              <div>• Shift+click to remove</div>
            </>
          )}
          <div>• Scroll to zoom canvas</div>
          <div>• Drag canvas to pan</div>
        </div>
      </div>
    </div>
  );
}
