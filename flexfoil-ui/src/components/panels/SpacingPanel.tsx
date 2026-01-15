/**
 * SpacingPanel - Full SSP (Position-Based Mesh Spacing) control panel
 * Ported from spacing-frontend with full interactive plot
 */

import { useCallback, useMemo, useEffect, useState } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { SSPPlot } from '../ssp';
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
    repanel,
    panels,
    coordinates,
  } = useAirfoilStore();

  const [liveUpdate, setLiveUpdate] = useState(true);

  // Compute output spacing
  const outputSpacing = useMemo(() => {
    return computeSpacing(spacingKnots, nPanels);
  }, [spacingKnots, nPanels]);

  // Real-time update: repanel immediately when knots or nPanels change
  useEffect(() => {
    if (liveUpdate && isWasmReady() && coordinates.length > 0) {
      // No debounce for real-time updates
      repanel();
    }
  }, [spacingKnots, nPanels, liveUpdate, repanel, coordinates.length]);

  // Handle manual apply
  const handleApply = useCallback(() => {
    if (isWasmReady()) {
      repanel();
    }
  }, [repanel]);

  // Handle knot change from drag
  const handleKnotChange = useCallback((index: number, knot: SpacingKnot) => {
    updateSpacingKnot(index, knot);
  }, [updateSpacingKnot]);

  // Handle add knot from click
  const handleAddKnot = useCallback((S: number, F: number) => {
    // Check if too close to existing knots
    const MIN_GAP = 0.02;
    const tooClose = spacingKnots.some(k => Math.abs(k.S - S) < MIN_GAP);
    if (tooClose) return;
    
    addSpacingKnot({ S, F });
  }, [spacingKnots, addSpacingKnot]);

  // Handle remove knot
  const handleRemoveKnot = useCallback((index: number) => {
    // Must keep at least 2 knots
    if (spacingKnots.length <= 2) return;
    // Cannot remove first or last knot
    if (index === 0 || index === spacingKnots.length - 1) return;
    
    removeSpacingKnot(index);
  }, [spacingKnots.length, removeSpacingKnot]);

  // Presets
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
          { S: 0, F: 0.3 },
          { S: 0.5, F: 1.5 },
          { S: 1, F: 0.3 },
        ]);
        break;
      case 'le-te':
        setSpacingKnots([
          { S: 0, F: 0.5 },
          { S: 0.25, F: 1.2 },
          { S: 0.5, F: 0.3 },
          { S: 0.75, F: 1.2 },
          { S: 1, F: 0.5 },
        ]);
        break;
      case 'le-fine':
        setSpacingKnots([
          { S: 0, F: 1.0 },
          { S: 0.4, F: 1.5 },
          { S: 0.5, F: 0.2 },
          { S: 0.6, F: 1.5 },
          { S: 1, F: 1.0 },
        ]);
        break;
      case 'te-fine':
        setSpacingKnots([
          { S: 0, F: 0.3 },
          { S: 0.2, F: 1.5 },
          { S: 0.5, F: 1.0 },
          { S: 0.8, F: 1.5 },
          { S: 1, F: 0.3 },
        ]);
        break;
    }
  }, [setSpacingKnots]);

  // Handle n adjustment with keyboard
  const decrementN = useCallback(() => setNPanels(Math.max(10, nPanels - 5)), [nPanels, setNPanels]);
  const incrementN = useCallback(() => setNPanels(Math.min(200, nPanels + 5)), [nPanels, setNPanels]);

  return (
    <div className="panel">
      <div className="panel-header">SSP Spacing Control</div>
      <div className="panel-content" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
        {/* Interactive SSP Plot */}
        <div style={{ 
          background: 'var(--bg-primary)', 
          borderRadius: '6px',
          padding: '8px',
          border: '1px solid var(--border-color)',
        }}>
          <SSPPlot
            knots={spacingKnots}
            si={outputSpacing}
            onKnotChange={handleKnotChange}
            onAddKnot={handleAddKnot}
            onRemoveKnot={handleRemoveKnot}
          />
        </div>

        {/* Panel count controls */}
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

        {/* Presets */}
        <div className="form-group">
          <div className="form-label">Spacing Presets</div>
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
          <div style={{ 
            fontSize: '10px', 
            color: 'var(--text-muted)',
            marginTop: '4px',
          }}>
            Points: {coordinates.length} → {panels.length}
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
      </div>
    </div>
  );
}
