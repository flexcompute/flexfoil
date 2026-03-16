/**
 * SpacingPanel - Panel spacing control with Simple and Advanced modes
 * 
 * Simple mode: Curvature-based automatic spacing
 * Advanced mode: Full SSP (Position-Based Mesh Spacing) control
 */

import { useCallback, useMemo, useEffect } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useRouteUiStore } from '../../stores/routeUiStore';
import { SSPPlot, FoilSpacingPlot } from '../ssp';
import { computeSpacing } from '../../lib/ssp';
import { isWasmReady } from '../../lib/wasm';
import type { SpacingKnot } from '../../types';

export function SpacingPanel() {
  const { 
    spacingKnots, 
    setSpacingKnots,
    updateSpacingKnot,
    addSpacingKnot,
    removeSpacingKnot,
    nPanels, 
    setNPanels,
    curvatureWeight,
    setCurvatureWeight,
    repanel,
    repanelWithXfoil,
    panels,
    coordinates,
    spacingPanelMode,
    setSpacingPanelMode,
    sspInterpolation,
    setSSPInterpolation,
    sspVisualization,
    setSSPVisualization,
  } = useAirfoilStore();

  const liveUpdate = useRouteUiStore((state) => state.spacingLiveUpdate);
  const setLiveUpdate = useRouteUiStore((state) => state.setSpacingLiveUpdate);

  // Compute output spacing for SSP mode
  const outputSpacing = useMemo(() => {
    return computeSpacing(spacingKnots, nPanels);
  }, [spacingKnots, nPanels]);

  // Real-time update: repanel based on current mode
  useEffect(() => {
    if (liveUpdate && isWasmReady() && coordinates.length > 0) {
      if (spacingPanelMode === 'simple') {
        // Simple mode uses curvature-based paneling
        repanelWithXfoil();
      } else {
        // Advanced mode uses SSP + curvature blend
        repanel();
      }
    }
  }, [spacingKnots, nPanels, curvatureWeight, liveUpdate, repanel, repanelWithXfoil, coordinates.length, spacingPanelMode]);

  // Handle manual apply
  const handleApply = useCallback(() => {
    if (isWasmReady()) {
      if (spacingPanelMode === 'simple') {
        repanelWithXfoil();
      } else {
        repanel();
      }
    }
  }, [repanel, repanelWithXfoil, spacingPanelMode]);

  // Handle knot change from drag
  const handleKnotChange = useCallback((index: number, knot: SpacingKnot) => {
    updateSpacingKnot(index, knot);
  }, [updateSpacingKnot]);

  // Handle add knot from click
  const handleAddKnot = useCallback((S: number, F: number) => {
    const MIN_GAP = 0.02;
    const tooClose = spacingKnots.some(k => Math.abs(k.S - S) < MIN_GAP);
    if (tooClose) return;
    addSpacingKnot({ S, F });
  }, [spacingKnots, addSpacingKnot]);

  // Handle remove knot
  const handleRemoveKnot = useCallback((index: number) => {
    if (spacingKnots.length <= 2) return;
    if (index === 0 || index === spacingKnots.length - 1) return;
    removeSpacingKnot(index);
  }, [spacingKnots.length, removeSpacingKnot]);

  // Presets for Advanced mode
  const applyPreset = useCallback((preset: 'uniform' | 'cosine' | 'le-te' | 'le-fine' | 'te-fine') => {
    switch (preset) {
      case 'uniform':
        setSpacingKnots([
          { S: 0, F: 1 },
          { S: 1, F: 1 },
        ]);
        break;
      case 'cosine':
        setSpacingKnots([
          { S: 0, F: 1.5 },
          { S: 0.25, F: 0.4 },
          { S: 0.5, F: 1.5 },
          { S: 0.75, F: 0.4 },
          { S: 1, F: 1.5 },
        ]);
        break;
      case 'le-te':
        setSpacingKnots([
          { S: 0, F: 1.5 },
          { S: 0.25, F: 0.3 },
          { S: 0.5, F: 1.5 },
          { S: 0.75, F: 0.3 },
          { S: 1, F: 1.5 },
        ]);
        break;
      case 'le-fine':
        setSpacingKnots([
          { S: 0, F: 0.5 },
          { S: 0.35, F: 0.3 },
          { S: 0.5, F: 2.0 },
          { S: 0.65, F: 0.3 },
          { S: 1, F: 0.5 },
        ]);
        break;
      case 'te-fine':
        setSpacingKnots([
          { S: 0, F: 2.0 },
          { S: 0.15, F: 0.3 },
          { S: 0.5, F: 0.5 },
          { S: 0.85, F: 0.3 },
          { S: 1, F: 2.0 },
        ]);
        break;
    }
  }, [setSpacingKnots]);

  // Handle n adjustment
  const decrementN = useCallback(() => setNPanels(Math.max(10, nPanels - 5)), [nPanels, setNPanels]);
  const incrementN = useCallback(() => setNPanels(Math.min(200, nPanels + 5)), [nPanels, setNPanels]);

  return (
    <div className="panel">
      <div className="panel-header">
        <span>Spacing Control</span>
      </div>
      <div className="panel-content" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
        
        {/* Mode Toggle */}
        <div className="form-group">
          <div style={{ 
            display: 'flex', 
            gap: '4px',
            background: 'var(--bg-tertiary)',
            padding: '4px',
            borderRadius: '6px',
          }}>
            <button 
              onClick={() => setSpacingPanelMode('simple')}
              style={{ 
                flex: 1,
                padding: '8px 12px',
                fontSize: '11px',
                fontWeight: 600,
                background: spacingPanelMode === 'simple' ? 'var(--accent-primary)' : 'transparent',
                color: spacingPanelMode === 'simple' ? 'var(--bg-primary)' : 'var(--text-secondary)',
                border: 'none',
                borderRadius: '4px',
                cursor: 'pointer',
                transition: 'all 0.15s ease',
              }}
            >
              Simple
            </button>
            <button 
              onClick={() => setSpacingPanelMode('advanced')}
              style={{ 
                flex: 1,
                padding: '8px 12px',
                fontSize: '11px',
                fontWeight: 600,
                background: spacingPanelMode === 'advanced' ? 'var(--accent-primary)' : 'transparent',
                color: spacingPanelMode === 'advanced' ? 'var(--bg-primary)' : 'var(--text-secondary)',
                border: 'none',
                borderRadius: '4px',
                cursor: 'pointer',
                transition: 'all 0.15s ease',
              }}
            >
              Advanced
            </button>
          </div>
        </div>

        {/* Panel count controls - always visible */}
        <div className="form-group">
          <div className="form-label" style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <span>Output Points (n)</span>
            <span style={{ 
              fontFamily: 'var(--font-mono)', 
              fontSize: '14px',
              fontWeight: 600,
              color: 'var(--accent-primary)',
            }}>
              {nPanels}
            </span>
          </div>
          <div className="form-row" style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <button 
              onClick={decrementN}
              style={{ 
                width: '32px', 
                height: '32px', 
                fontSize: '18px',
                fontWeight: 'bold',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
              }}
              disabled={nPanels <= 10}
            >
              −
            </button>
            <input
              type="range"
              min={10}
              max={200}
              step={5}
              value={nPanels}
              onChange={(e) => setNPanels(parseInt(e.target.value))}
              style={{ flex: 1 }}
            />
            <button 
              onClick={incrementN}
              style={{ 
                width: '32px', 
                height: '32px', 
                fontSize: '18px',
                fontWeight: 'bold',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
              }}
              disabled={nPanels >= 200}
            >
              +
            </button>
          </div>
        </div>

        {/* ========== SIMPLE MODE ========== */}
        {spacingPanelMode === 'simple' && (
          <>
            <div style={{ 
              background: 'var(--bg-tertiary)',
              borderRadius: '6px',
              padding: '12px',
              border: '1px solid var(--border-color)',
            }}>
              <div style={{ 
                fontSize: '12px', 
                fontWeight: 600,
                color: 'var(--accent-success)',
                marginBottom: '8px',
              }}>
                Curvature-Based Spacing
              </div>
              <div style={{ 
                fontSize: '11px', 
                color: 'var(--text-muted)',
                lineHeight: 1.5,
              }}>
                Automatically clusters panels at high-curvature regions (leading edge) 
                and the trailing edge. This is the recommended setting for most analyses.
              </div>
            </div>

            <div style={{ 
              fontSize: '10px', 
              color: 'var(--text-muted)',
              marginTop: '4px',
            }}>
              Points: {coordinates.length} → {panels.length}
            </div>
          </>
        )}

        {/* ========== ADVANCED MODE ========== */}
        {spacingPanelMode === 'advanced' && (
          <>
            {/* SSP Visualization Toggle */}
            <div className="form-group">
              <div className="form-label">Visualization</div>
              <div style={{ 
                display: 'flex', 
                gap: '4px',
                background: 'var(--bg-tertiary)',
                padding: '3px',
                borderRadius: '4px',
              }}>
                <button 
                  onClick={() => setSSPVisualization('plot')}
                  style={{ 
                    flex: 1,
                    padding: '6px 8px',
                    fontSize: '10px',
                    fontWeight: 500,
                    background: sspVisualization === 'plot' ? 'var(--bg-secondary)' : 'transparent',
                    color: sspVisualization === 'plot' ? 'var(--text-primary)' : 'var(--text-muted)',
                    border: sspVisualization === 'plot' ? '1px solid var(--border-color)' : '1px solid transparent',
                    borderRadius: '3px',
                    cursor: 'pointer',
                  }}
                >
                  F vs S Plot
                </button>
                <button 
                  onClick={() => setSSPVisualization('foil')}
                  style={{ 
                    flex: 1,
                    padding: '6px 8px',
                    fontSize: '10px',
                    fontWeight: 500,
                    background: sspVisualization === 'foil' ? 'var(--bg-secondary)' : 'transparent',
                    color: sspVisualization === 'foil' ? 'var(--text-primary)' : 'var(--text-muted)',
                    border: sspVisualization === 'foil' ? '1px solid var(--border-color)' : '1px solid transparent',
                    borderRadius: '3px',
                    cursor: 'pointer',
                  }}
                >
                  Foil View
                </button>
              </div>
            </div>

            {/* SSP Plot / Foil View */}
            <div style={{ 
              background: 'var(--bg-primary)', 
              borderRadius: '6px',
              padding: '8px',
              border: '1px solid var(--border-color)',
            }}>
              {sspVisualization === 'plot' ? (
                <SSPPlot
                  knots={spacingKnots}
                  si={outputSpacing}
                  onKnotChange={handleKnotChange}
                  onAddKnot={handleAddKnot}
                  onRemoveKnot={handleRemoveKnot}
                />
              ) : (
                <FoilSpacingPlot
                  knots={spacingKnots}
                  coordinates={coordinates}
                  onKnotChange={handleKnotChange}
                  onAddKnot={handleAddKnot}
                  onRemoveKnot={handleRemoveKnot}
                />
              )}
            </div>

            {/* Curvature blend control */}
            <div className="form-group">
              <div className="form-label" style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                <span>Curvature Blend</span>
                <span style={{ 
                  fontFamily: 'var(--font-mono)', 
                  fontSize: '14px',
                  fontWeight: 600,
                  color: curvatureWeight > 0 ? 'var(--accent-warning)' : 'var(--text-muted)',
                }}>
                  {Math.round(curvatureWeight * 100)}%
                </span>
              </div>
              <div className="form-row" style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                <span style={{ fontSize: '10px', color: 'var(--text-muted)', width: '28px' }}>SSP</span>
                <input
                  type="range"
                  min={0}
                  max={100}
                  step={5}
                  value={curvatureWeight * 100}
                  onChange={(e) => setCurvatureWeight(parseInt(e.target.value) / 100)}
                  style={{ flex: 1 }}
                />
                <span style={{ fontSize: '10px', color: 'var(--text-muted)', width: '28px', textAlign: 'right' }}>κ</span>
              </div>
              <div style={{ 
                fontSize: '10px', 
                color: 'var(--text-muted)',
                marginTop: '4px',
              }}>
                Blend between manual SSP knots and automatic curvature-based clustering
              </div>
            </div>

            {/* Interpolation Mode */}
            <div className="form-group">
              <div className="form-label">F(s) Interpolation</div>
              <div style={{ 
                display: 'flex', 
                gap: '4px',
                background: 'var(--bg-tertiary)',
                padding: '3px',
                borderRadius: '4px',
              }}>
                <button 
                  onClick={() => setSSPInterpolation('linear')}
                  style={{ 
                    flex: 1,
                    padding: '6px 8px',
                    fontSize: '10px',
                    fontWeight: 500,
                    background: sspInterpolation === 'linear' ? 'var(--bg-secondary)' : 'transparent',
                    color: sspInterpolation === 'linear' ? 'var(--text-primary)' : 'var(--text-muted)',
                    border: sspInterpolation === 'linear' ? '1px solid var(--border-color)' : '1px solid transparent',
                    borderRadius: '3px',
                    cursor: 'pointer',
                  }}
                >
                  Linear
                </button>
                <button 
                  onClick={() => setSSPInterpolation('spline')}
                  style={{ 
                    flex: 1,
                    padding: '6px 8px',
                    fontSize: '10px',
                    fontWeight: 500,
                    background: sspInterpolation === 'spline' ? 'var(--bg-secondary)' : 'transparent',
                    color: sspInterpolation === 'spline' ? 'var(--text-primary)' : 'var(--text-muted)',
                    border: sspInterpolation === 'spline' ? '1px solid var(--border-color)' : '1px solid transparent',
                    borderRadius: '3px',
                    cursor: 'pointer',
                  }}
                  disabled // Not yet implemented
                  title="Spline interpolation coming soon"
                >
                  Spline
                </button>
              </div>
              <div style={{ 
                fontSize: '9px', 
                color: 'var(--text-muted)',
                marginTop: '4px',
              }}>
                Linear is Drela's original method; spline gives smoother transitions
              </div>
            </div>

            {/* Presets */}
            <div className="form-group">
              <div className="form-label">Quick Presets</div>
              <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
                <button onClick={() => applyPreset('uniform')} style={{ flex: '1 1 auto', fontSize: '10px', padding: '6px 8px' }}>
                  Uniform
                </button>
                <button onClick={() => applyPreset('cosine')} style={{ flex: '1 1 auto', fontSize: '10px', padding: '6px 8px' }}>
                  Cosine
                </button>
                <button onClick={() => applyPreset('le-te')} style={{ flex: '1 1 auto', fontSize: '10px', padding: '6px 8px' }}>
                  LE+TE
                </button>
                <button onClick={() => applyPreset('le-fine')} style={{ flex: '1 1 auto', fontSize: '10px', padding: '6px 8px' }}>
                  LE Fine
                </button>
                <button onClick={() => applyPreset('te-fine')} style={{ flex: '1 1 auto', fontSize: '10px', padding: '6px 8px' }}>
                  TE Fine
                </button>
              </div>
            </div>

            {/* Knot info */}
            <div className="form-group">
              <div className="form-label">Knots ({spacingKnots.length})</div>
              <div style={{ 
                fontSize: '10px', 
                color: 'var(--text-muted)',
                background: 'var(--bg-tertiary)',
                padding: '8px',
                borderRadius: '4px',
                fontFamily: 'var(--font-mono)',
                maxHeight: '100px',
                overflowY: 'auto',
              }}>
                {spacingKnots.map((k, i) => (
                  <div key={i} style={{ display: 'flex', justifyContent: 'space-between' }}>
                    <span style={{ color: i === 0 || i === spacingKnots.length - 1 ? 'var(--accent-success)' : 'var(--accent-primary)' }}>
                      {i}:
                    </span>
                    <span>S={k.S.toFixed(3)}</span>
                    <span>F={k.F.toFixed(2)}</span>
                  </div>
                ))}
              </div>
            </div>

            {/* Instructions */}
            <div style={{ 
              fontSize: '10px', 
              color: 'var(--text-muted)',
              lineHeight: 1.5,
              borderTop: '1px solid var(--border-color)',
              paddingTop: '8px',
            }}>
              <strong>Controls:</strong>
              <ul style={{ margin: '4px 0 0 0', paddingLeft: '16px' }}>
                <li>Drag knots to adjust spacing</li>
                <li>Click plot to add a knot</li>
                <li>Right-click or shift+click to remove</li>
                <li>Green knots (endpoints) move vertically only</li>
              </ul>
            </div>
          </>
        )}

        {/* Apply button and live update toggle */}
        <div className="form-group">
          <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
            <button 
              onClick={handleApply}
              style={{ 
                flex: 1,
                padding: '10px',
                fontWeight: 600,
                background: 'var(--accent-primary)',
                color: 'var(--bg-primary)',
                border: 'none',
                borderRadius: '6px',
                cursor: 'pointer',
              }}
              disabled={!isWasmReady()}
            >
              Apply Spacing
            </button>
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
                checked={liveUpdate} 
                onChange={(e) => setLiveUpdate(e.target.checked)}
              />
              Live
            </label>
          </div>
          {spacingPanelMode === 'advanced' && (
            <div style={{ 
              fontSize: '10px', 
              color: 'var(--text-muted)',
              marginTop: '4px',
            }}>
              Points: {coordinates.length} → {panels.length}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
