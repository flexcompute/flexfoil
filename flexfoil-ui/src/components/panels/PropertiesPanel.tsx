/**
 * PropertiesPanel - Coordinate display and airfoil properties
 */

import { useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';

export function PropertiesPanel() {
  const { name, coordinates, nPanels, controlMode, bsplineControlPoints } = useAirfoilStore();

  // Compute basic airfoil properties
  const properties = useMemo(() => {
    if (coordinates.length < 3) {
      return { chord: 0, maxThickness: 0, maxCamber: 0, leRadius: 0 };
    }

    // Find chord (max x - min x)
    const xVals = coordinates.map(p => p.x);
    const minX = Math.min(...xVals);
    const maxX = Math.max(...xVals);
    const chord = maxX - minX;

    // Find LE index (min x)
    const leIndex = coordinates.findIndex(p => p.x === minX);

    // Split into upper and lower surfaces
    const upper = coordinates.slice(0, leIndex + 1);
    const lower = coordinates.slice(leIndex);

    // Compute thickness and camber at various x positions
    let maxThickness = 0;
    let maxCamber = 0;

    for (let x = 0.1; x <= 0.9; x += 0.1) {
      // Find upper y at this x (interpolate)
      let yu = 0;
      for (let i = 0; i < upper.length - 1; i++) {
        if ((upper[i].x >= x && upper[i + 1].x <= x) || 
            (upper[i].x <= x && upper[i + 1].x >= x)) {
          const t = (x - upper[i].x) / (upper[i + 1].x - upper[i].x);
          yu = upper[i].y + t * (upper[i + 1].y - upper[i].y);
          break;
        }
      }

      // Find lower y at this x
      let yl = 0;
      for (let i = 0; i < lower.length - 1; i++) {
        if ((lower[i].x >= x && lower[i + 1].x <= x) ||
            (lower[i].x <= x && lower[i + 1].x >= x)) {
          const t = (x - lower[i].x) / (lower[i + 1].x - lower[i].x);
          yl = lower[i].y + t * (lower[i + 1].y - lower[i].y);
          break;
        }
      }

      const thickness = Math.abs(yu - yl);
      const camber = (yu + yl) / 2;

      if (thickness > maxThickness) maxThickness = thickness;
      if (Math.abs(camber) > Math.abs(maxCamber)) maxCamber = camber;
    }

    // Estimate LE radius (simplified)
    const leRadius = coordinates[leIndex] ? 
      Math.min(
        Math.abs(coordinates[Math.max(0, leIndex - 1)]?.y || 0),
        Math.abs(coordinates[Math.min(coordinates.length - 1, leIndex + 1)]?.y || 0)
      ) : 0;

    return {
      chord,
      maxThickness: maxThickness / chord * 100,
      maxCamber: maxCamber / chord * 100,
      leRadius: leRadius / chord * 100,
    };
  }, [coordinates]);

  return (
    <div className="panel">
      <div className="panel-header">Properties</div>
      <div className="panel-content">
        {/* Airfoil name */}
        <div className="form-group">
          <div className="form-label">Airfoil</div>
          <div style={{
            padding: '6px 10px',
            background: 'var(--bg-tertiary)',
            borderRadius: '4px',
            fontWeight: 600,
            color: 'var(--accent-primary)',
          }}>
            {name}
          </div>
        </div>

        {/* Geometric properties */}
        <div className="form-group">
          <div className="form-label">Geometry</div>
          <div style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1fr',
            gap: '8px',
          }}>
            <PropertyCard label="Chord" value="1.000" unit="c" />
            <PropertyCard 
              label="Max t/c" 
              value={properties.maxThickness.toFixed(1)} 
              unit="%" 
            />
            <PropertyCard 
              label="Max Camber" 
              value={properties.maxCamber.toFixed(2)} 
              unit="%" 
            />
            <PropertyCard 
              label="LE Radius" 
              value={properties.leRadius.toFixed(2)} 
              unit="%c" 
            />
          </div>
        </div>

        {/* Point counts */}
        <div className="form-group">
          <div className="form-label">Points</div>
          <div style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1fr',
            gap: '8px',
          }}>
            <PropertyCard label="Input" value={coordinates.length.toString()} />
            <PropertyCard label="Panels" value={nPanels.toString()} />
            {controlMode === 'bspline' && (
              <PropertyCard label="Control" value={bsplineControlPoints.length.toString()} />
            )}
          </div>
        </div>

        {/* Coordinate table */}
        <div className="form-group">
          <div className="form-label">Coordinates (first 10)</div>
          <div style={{
            maxHeight: '200px',
            overflow: 'auto',
            border: '1px solid var(--border-color)',
            borderRadius: '4px',
            fontFamily: 'var(--font-mono)',
            fontSize: '10px',
          }}>
            {/* Header */}
            <div style={{
              display: 'grid',
              gridTemplateColumns: '30px 1fr 1fr',
              gap: '4px',
              padding: '4px 8px',
              background: 'var(--bg-tertiary)',
              borderBottom: '1px solid var(--border-color)',
              position: 'sticky',
              top: 0,
            }}>
              <span>#</span>
              <span>x/c</span>
              <span>y/c</span>
            </div>
            
            {/* Data rows */}
            {coordinates.slice(0, 10).map((p, i) => (
              <div
                key={i}
                style={{
                  display: 'grid',
                  gridTemplateColumns: '30px 1fr 1fr',
                  gap: '4px',
                  padding: '2px 8px',
                  borderBottom: '1px solid var(--border-color)',
                }}
              >
                <span style={{ color: 'var(--text-muted)' }}>{i}</span>
                <span>{p.x.toFixed(5)}</span>
                <span style={{ color: p.y >= 0 ? 'var(--accent-primary)' : 'var(--accent-secondary)' }}>
                  {p.y.toFixed(5)}
                </span>
              </div>
            ))}
            {coordinates.length > 10 && (
              <div style={{
                padding: '4px 8px',
                textAlign: 'center',
                color: 'var(--text-muted)',
              }}>
                ... {coordinates.length - 10} more points
              </div>
            )}
          </div>
        </div>

        {/* Export */}
        <button
          onClick={() => {
            const data = coordinates.map(p => `${p.x.toFixed(6)} ${p.y.toFixed(6)}`).join('\n');
            const blob = new Blob([`${name}\n${data}`], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `${name.replace(/\s+/g, '_')}.dat`;
            a.click();
            URL.revokeObjectURL(url);
          }}
          style={{ width: '100%' }}
        >
          Export Coordinates (.dat)
        </button>
      </div>
    </div>
  );
}

function PropertyCard({ label, value, unit }: { label: string; value: string; unit?: string }) {
  return (
    <div style={{
      padding: '8px',
      background: 'var(--bg-tertiary)',
      borderRadius: '4px',
    }}>
      <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginBottom: '2px' }}>
        {label}
      </div>
      <div style={{ fontFamily: 'var(--font-mono)', fontWeight: 600 }}>
        {value}
        {unit && <span style={{ fontSize: '10px', color: 'var(--text-secondary)', marginLeft: '2px' }}>{unit}</span>}
      </div>
    </div>
  );
}
