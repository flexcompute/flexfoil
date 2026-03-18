/**
 * GeometryDesignPanel - GDES operations: flap deflection, TE gap, LE radius, transforms.
 */

import React, { useCallback, useRef } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import {
  deflectFlap,
  setTeGap,
  setLeRadius,
  rotateAirfoil,
  scaleAirfoil as scaleAirfoilWasm,
  translateAirfoil,
  isWasmReady,
  repanelXfoil,
} from '../../lib/wasm';
import type { AirfoilPoint } from '../../types';

function applyAndRepanel(
  coords: { x: number; y: number }[],
  nPanels: number,
): { coordinates: AirfoilPoint[]; panels: AirfoilPoint[] } {
  const coordinates: AirfoilPoint[] = coords.map(p => ({ x: p.x, y: p.y }));
  let panels = coordinates;
  if (isWasmReady()) {
    try {
      const repaneled = repanelXfoil(coords, nPanels);
      if (repaneled.length > 0) {
        panels = repaneled.map(pt => ({ x: pt.x, y: pt.y }));
      }
    } catch { /* keep raw */ }
  }
  return { coordinates, panels };
}

export function GeometryDesignPanel() {
  const {
    coordinates, nPanels, geometryDesign, setGeometryDesign,
    setCoordinates, setPanels,
  } = useAirfoilStore();
  
  const baseRef = useRef(coordinates);
  
  const applyFlap = useCallback(() => {
    if (!isWasmReady()) return;
    const result = deflectFlap(coordinates, geometryDesign.flapHingeX, geometryDesign.flapDeflection);
    if (result.length > 0) {
      const { coordinates: c, panels: p } = applyAndRepanel(result, nPanels);
      useAirfoilStore.setState({ coordinates: c, panels: p, baseCoordinates: c });
    }
  }, [coordinates, nPanels, geometryDesign.flapHingeX, geometryDesign.flapDeflection]);

  const applyTeGap = useCallback(() => {
    if (!isWasmReady()) return;
    const result = setTeGap(coordinates, geometryDesign.teGap, geometryDesign.teGapBlend);
    if (result.length > 0) {
      const { coordinates: c, panels: p } = applyAndRepanel(result, nPanels);
      useAirfoilStore.setState({ coordinates: c, panels: p, baseCoordinates: c });
    }
  }, [coordinates, nPanels, geometryDesign.teGap, geometryDesign.teGapBlend]);

  const applyLeRadius = useCallback(() => {
    if (!isWasmReady()) return;
    const result = setLeRadius(coordinates, geometryDesign.leRadiusFactor);
    if (result.length > 0) {
      const { coordinates: c, panels: p } = applyAndRepanel(result, nPanels);
      useAirfoilStore.setState({ coordinates: c, panels: p, baseCoordinates: c });
    }
  }, [coordinates, nPanels, geometryDesign.leRadiusFactor]);

  return (
    <div className="geometry-design-panel" style={{ padding: '10px', overflow: 'auto', fontSize: '12px' }}>
      {/* Flap Deflection */}
      <fieldset style={{ border: '1px solid var(--border-color)', borderRadius: 4, padding: '8px', marginBottom: '10px' }}>
        <legend style={{ fontSize: '11px', fontWeight: 600, color: 'var(--text-secondary)', padding: '0 4px' }}>Flap Deflection</legend>
        <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' }}>
          <label style={{ color: 'var(--text-muted)', minWidth: '70px' }}>Hinge x/c:</label>
          <input
            type="range" min={0.5} max={0.95} step={0.01}
            value={geometryDesign.flapHingeX}
            onChange={e => setGeometryDesign({ flapHingeX: parseFloat(e.target.value) })}
            style={{ flex: 1 }}
          />
          <span style={{ minWidth: '35px', textAlign: 'right' }}>{geometryDesign.flapHingeX.toFixed(2)}</span>
        </div>
        <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' }}>
          <label style={{ color: 'var(--text-muted)', minWidth: '70px' }}>Deflection:</label>
          <input
            type="range" min={-30} max={30} step={0.5}
            value={geometryDesign.flapDeflection}
            onChange={e => setGeometryDesign({ flapDeflection: parseFloat(e.target.value) })}
            style={{ flex: 1 }}
          />
          <span style={{ minWidth: '35px', textAlign: 'right' }}>{geometryDesign.flapDeflection.toFixed(1)}°</span>
        </div>
        <button className="btn btn-xs btn-primary" onClick={applyFlap}>Apply Flap</button>
      </fieldset>

      {/* TE Gap */}
      <fieldset style={{ border: '1px solid var(--border-color)', borderRadius: 4, padding: '8px', marginBottom: '10px' }}>
        <legend style={{ fontSize: '11px', fontWeight: 600, color: 'var(--text-secondary)', padding: '0 4px' }}>Trailing-Edge Gap</legend>
        <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' }}>
          <label style={{ color: 'var(--text-muted)', minWidth: '70px' }}>Gap (% c):</label>
          <input
            type="range" min={0} max={0.05} step={0.001}
            value={geometryDesign.teGap}
            onChange={e => setGeometryDesign({ teGap: parseFloat(e.target.value) })}
            style={{ flex: 1 }}
          />
          <span style={{ minWidth: '35px', textAlign: 'right' }}>{(geometryDesign.teGap * 100).toFixed(1)}%</span>
        </div>
        <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' }}>
          <label style={{ color: 'var(--text-muted)', minWidth: '70px' }}>Blend:</label>
          <input
            type="range" min={0.2} max={1.0} step={0.05}
            value={geometryDesign.teGapBlend}
            onChange={e => setGeometryDesign({ teGapBlend: parseFloat(e.target.value) })}
            style={{ flex: 1 }}
          />
          <span style={{ minWidth: '35px', textAlign: 'right' }}>{geometryDesign.teGapBlend.toFixed(2)}</span>
        </div>
        <button className="btn btn-xs btn-primary" onClick={applyTeGap}>Apply TE Gap</button>
      </fieldset>

      {/* LE Radius */}
      <fieldset style={{ border: '1px solid var(--border-color)', borderRadius: 4, padding: '8px', marginBottom: '10px' }}>
        <legend style={{ fontSize: '11px', fontWeight: 600, color: 'var(--text-secondary)', padding: '0 4px' }}>Leading-Edge Radius</legend>
        <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px' }}>
          <label style={{ color: 'var(--text-muted)', minWidth: '70px' }}>Scale:</label>
          <input
            type="range" min={0.2} max={3.0} step={0.05}
            value={geometryDesign.leRadiusFactor}
            onChange={e => setGeometryDesign({ leRadiusFactor: parseFloat(e.target.value) })}
            style={{ flex: 1 }}
          />
          <span style={{ minWidth: '35px', textAlign: 'right' }}>{geometryDesign.leRadiusFactor.toFixed(2)}×</span>
        </div>
        <button className="btn btn-xs btn-primary" onClick={applyLeRadius}>Apply LE Radius</button>
      </fieldset>
    </div>
  );
}
