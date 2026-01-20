import React, { useState, useMemo } from 'react';

/**
 * Interactive plot showing the Hk (kinematic shape parameter) correlation
 * used in XFOIL's boundary layer closure relations.
 * 
 * This demonstrates:
 * - The Hk = (H - 0.29·M²) / (1 + 0.113·M²) compressibility correction
 * - Typical Hk values for laminar (2.5-4.0) and turbulent (1.3-2.0) flows
 */

interface HkCorrelationPlotProps {
  /** Initial Mach number (0-0.7) */
  initialMach?: number;
  /** Plot height in pixels */
  height?: number;
}

export default function HkCorrelationPlot({ 
  initialMach = 0, 
  height = 300 
}: HkCorrelationPlotProps): JSX.Element {
  const [mach, setMach] = useState(initialMach);

  // Compute Hk from H for various H values
  const data = useMemo(() => {
    const hValues: number[] = [];
    const hkValues: number[] = [];
    const m2 = mach * mach;

    for (let h = 1.1; h <= 5.0; h += 0.1) {
      const hk = (h - 0.29 * m2) / (1 + 0.113 * m2);
      hValues.push(h);
      hkValues.push(hk);
    }

    return { hValues, hkValues };
  }, [mach]);

  // SVG dimensions
  const width = 400;
  const padding = { top: 20, right: 20, bottom: 40, left: 50 };
  const plotWidth = width - padding.left - padding.right;
  const plotHeight = height - padding.top - padding.bottom;

  // Scales
  const xMin = 1.0, xMax = 5.0;
  const yMin = 0.5, yMax = 5.0;
  
  const scaleX = (h: number) => padding.left + ((h - xMin) / (xMax - xMin)) * plotWidth;
  const scaleY = (hk: number) => padding.top + plotHeight - ((hk - yMin) / (yMax - yMin)) * plotHeight;

  // Generate path
  const pathD = data.hValues.map((h, i) => {
    const x = scaleX(h);
    const y = scaleY(data.hkValues[i]);
    return `${i === 0 ? 'M' : 'L'} ${x} ${y}`;
  }).join(' ');

  // Reference lines for typical flow regimes
  const turbulentHkRange = { min: 1.3, max: 2.0 };
  const laminarHkRange = { min: 2.5, max: 4.0 };

  return (
    <div style={{ 
      border: '1px solid var(--ifm-color-emphasis-300)', 
      borderRadius: '8px', 
      padding: '16px',
      backgroundColor: 'var(--ifm-background-surface-color)',
      marginBottom: '1rem'
    }}>
      <h4 style={{ margin: '0 0 12px 0' }}>Hk Compressibility Correction</h4>
      
      <div style={{ marginBottom: '12px' }}>
        <label style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
          <span>Mach Number: <strong>{mach.toFixed(2)}</strong></span>
          <input
            type="range"
            min="0"
            max="0.7"
            step="0.01"
            value={mach}
            onChange={(e) => setMach(parseFloat(e.target.value))}
            style={{ width: '150px' }}
          />
        </label>
      </div>

      <svg width={width} height={height} style={{ display: 'block', margin: '0 auto' }}>
        {/* Grid lines */}
        {[1, 2, 3, 4, 5].map(h => (
          <line
            key={`grid-x-${h}`}
            x1={scaleX(h)}
            y1={padding.top}
            x2={scaleX(h)}
            y2={padding.top + plotHeight}
            stroke="var(--ifm-color-emphasis-200)"
            strokeDasharray="2,2"
          />
        ))}
        {[1, 2, 3, 4].map(hk => (
          <line
            key={`grid-y-${hk}`}
            x1={padding.left}
            y1={scaleY(hk)}
            x2={padding.left + plotWidth}
            y2={scaleY(hk)}
            stroke="var(--ifm-color-emphasis-200)"
            strokeDasharray="2,2"
          />
        ))}

        {/* Turbulent regime shading */}
        <rect
          x={padding.left}
          y={scaleY(turbulentHkRange.max)}
          width={plotWidth}
          height={scaleY(turbulentHkRange.min) - scaleY(turbulentHkRange.max)}
          fill="rgba(59, 130, 246, 0.1)"
        />
        <text
          x={padding.left + 5}
          y={scaleY((turbulentHkRange.min + turbulentHkRange.max) / 2) + 4}
          fontSize="10"
          fill="var(--ifm-color-emphasis-600)"
        >
          Turbulent
        </text>

        {/* Laminar regime shading */}
        <rect
          x={padding.left}
          y={scaleY(laminarHkRange.max)}
          width={plotWidth}
          height={scaleY(laminarHkRange.min) - scaleY(laminarHkRange.max)}
          fill="rgba(34, 197, 94, 0.1)"
        />
        <text
          x={padding.left + 5}
          y={scaleY((laminarHkRange.min + laminarHkRange.max) / 2) + 4}
          fontSize="10"
          fill="var(--ifm-color-emphasis-600)"
        >
          Laminar
        </text>

        {/* Main curve */}
        <path
          d={pathD}
          fill="none"
          stroke="var(--ifm-color-primary)"
          strokeWidth="2"
        />

        {/* Identity line (M=0 reference) */}
        {mach > 0 && (
          <line
            x1={scaleX(xMin)}
            y1={scaleY(xMin)}
            x2={scaleX(yMax)}
            y2={scaleY(yMax)}
            stroke="var(--ifm-color-emphasis-400)"
            strokeWidth="1"
            strokeDasharray="4,4"
          />
        )}

        {/* Axes */}
        <line
          x1={padding.left}
          y1={padding.top + plotHeight}
          x2={padding.left + plotWidth}
          y2={padding.top + plotHeight}
          stroke="var(--ifm-color-emphasis-500)"
        />
        <line
          x1={padding.left}
          y1={padding.top}
          x2={padding.left}
          y2={padding.top + plotHeight}
          stroke="var(--ifm-color-emphasis-500)"
        />

        {/* Axis labels */}
        <text
          x={padding.left + plotWidth / 2}
          y={height - 5}
          textAnchor="middle"
          fontSize="12"
          fill="var(--ifm-color-emphasis-700)"
        >
          H (shape factor)
        </text>
        <text
          x={15}
          y={padding.top + plotHeight / 2}
          textAnchor="middle"
          fontSize="12"
          fill="var(--ifm-color-emphasis-700)"
          transform={`rotate(-90, 15, ${padding.top + plotHeight / 2})`}
        >
          Hk (kinematic)
        </text>

        {/* Tick labels */}
        {[1, 2, 3, 4, 5].map(h => (
          <text
            key={`tick-x-${h}`}
            x={scaleX(h)}
            y={padding.top + plotHeight + 15}
            textAnchor="middle"
            fontSize="10"
            fill="var(--ifm-color-emphasis-600)"
          >
            {h}
          </text>
        ))}
        {[1, 2, 3, 4].map(hk => (
          <text
            key={`tick-y-${hk}`}
            x={padding.left - 8}
            y={scaleY(hk) + 4}
            textAnchor="end"
            fontSize="10"
            fill="var(--ifm-color-emphasis-600)"
          >
            {hk}
          </text>
        ))}
      </svg>

      <div style={{ 
        fontSize: '12px', 
        color: 'var(--ifm-color-emphasis-600)',
        marginTop: '8px',
        textAlign: 'center'
      }}>
        Formula: Hk = (H - 0.29·M²) / (1 + 0.113·M²) &nbsp;|&nbsp; 
        {mach > 0 && <span>Dashed line shows M=0 reference</span>}
        {mach === 0 && <span>Adjust Mach to see compressibility effect</span>}
      </div>
    </div>
  );
}
