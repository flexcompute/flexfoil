/**
 * VisualizationPanel - Controls for flow visualization options
 * 
 * Contains:
 * - Display toggles (Grid, Curve, Panels, Points, Controls, Streamlines, Smoke)
 * - Streamline options (density, adaptive)
 * - Smoke options (density, particles per blob, spawn interval, max age)
 * - Flow speed control
 */

import { useVisualizationStore } from '../../stores/visualizationStore';

export function VisualizationPanel() {
  const {
    // Display toggles
    showGrid,
    showCurve,
    showPanels,
    showPoints,
    showControls,
    showStreamlines,
    showSmoke,
    setShowGrid,
    setShowCurve,
    setShowPanels,
    setShowPoints,
    setShowControls,
    setShowStreamlines,
    setShowSmoke,
    
    // Streamline options
    streamlineDensity,
    adaptiveStreamlines,
    setStreamlineDensity,
    setAdaptiveStreamlines,
    
    // Smoke options
    smokeDensity,
    smokeParticlesPerBlob,
    smokeSpawnInterval,
    smokeMaxAge,
    setSmokeDensity,
    setSmokeParticlesPerBlob,
    setSmokeSpawnInterval,
    setSmokeMaxAge,
    
    // Flow speed
    flowSpeed,
    setFlowSpeed,
    
    // Reset
    resetVisualization,
  } = useVisualizationStore();

  return (
    <div className="panel-content" style={{ 
      padding: '12px', 
      display: 'flex', 
      flexDirection: 'column', 
      gap: '16px',
      overflowY: 'auto',
      height: '100%',
    }}>
      {/* Display Toggles Section */}
      <section>
        <h4 style={{ 
          margin: '0 0 8px 0', 
          fontSize: '12px', 
          fontWeight: 600,
          color: 'var(--text-primary)',
          textTransform: 'uppercase',
          letterSpacing: '0.5px',
        }}>
          Display
        </h4>
        <div style={{ 
          display: 'grid', 
          gridTemplateColumns: 'repeat(2, 1fr)', 
          gap: '6px',
        }}>
          <ToggleItem label="Grid" checked={showGrid} onChange={setShowGrid} />
          <ToggleItem label="Curve" checked={showCurve} onChange={setShowCurve} />
          <ToggleItem label="Panels" checked={showPanels} onChange={setShowPanels} />
          <ToggleItem label="Points" checked={showPoints} onChange={setShowPoints} />
          <ToggleItem label="Controls" checked={showControls} onChange={setShowControls} />
        </div>
      </section>

      {/* Flow Visualization Section */}
      <section>
        <h4 style={{ 
          margin: '0 0 8px 0', 
          fontSize: '12px', 
          fontWeight: 600,
          color: 'var(--text-primary)',
          textTransform: 'uppercase',
          letterSpacing: '0.5px',
        }}>
          Flow Visualization
        </h4>
        <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
          <ToggleItem label="Streamlines" checked={showStreamlines} onChange={setShowStreamlines} />
          <ToggleItem label="Smoke" checked={showSmoke} onChange={setShowSmoke} />
        </div>
      </section>

      {/* Streamline Options Section */}
      {showStreamlines && (
        <section>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            fontSize: '12px', 
            fontWeight: 600,
            color: 'var(--text-primary)',
            textTransform: 'uppercase',
            letterSpacing: '0.5px',
          }}>
            Streamline Options
          </h4>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
            <SliderItem
              label="Density"
              value={streamlineDensity}
              min={10}
              max={150}
              step={5}
              onChange={setStreamlineDensity}
              formatValue={(v) => `${v} lines`}
            />
            <ToggleItem 
              label="Adaptive (zoom-based)" 
              checked={adaptiveStreamlines} 
              onChange={setAdaptiveStreamlines} 
            />
          </div>
        </section>
      )}

      {/* Smoke Options Section */}
      {showSmoke && (
        <section>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            fontSize: '12px', 
            fontWeight: 600,
            color: 'var(--text-primary)',
            textTransform: 'uppercase',
            letterSpacing: '0.5px',
          }}>
            Smoke Options
          </h4>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
            <SliderItem
              label="Spawn Points"
              value={smokeDensity}
              min={10}
              max={80}
              step={2}
              onChange={setSmokeDensity}
              formatValue={(v) => `${v}`}
            />
            <SliderItem
              label="Particles/Blob"
              value={smokeParticlesPerBlob}
              min={3}
              max={30}
              step={1}
              onChange={setSmokeParticlesPerBlob}
              formatValue={(v) => `${v}`}
            />
            <SliderItem
              label="Spawn Interval"
              value={smokeSpawnInterval}
              min={0.05}
              max={1.0}
              step={0.05}
              onChange={setSmokeSpawnInterval}
              formatValue={(v) => `${v.toFixed(2)}s`}
            />
            <SliderItem
              label="Max Age"
              value={smokeMaxAge}
              min={0.5}
              max={10.0}
              step={0.5}
              onChange={setSmokeMaxAge}
              formatValue={(v) => `${v.toFixed(1)}s`}
            />
          </div>
        </section>
      )}

      {/* Flow Speed Section */}
      {(showStreamlines || showSmoke) && (
        <section>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            fontSize: '12px', 
            fontWeight: 600,
            color: 'var(--text-primary)',
            textTransform: 'uppercase',
            letterSpacing: '0.5px',
          }}>
            Animation Speed
          </h4>
          <SliderItem
            label="Flow Speed"
            value={flowSpeed}
            min={0.1}
            max={5.0}
            step={0.1}
            onChange={setFlowSpeed}
            formatValue={(v) => `${v.toFixed(1)}x`}
          />
        </section>
      )}

      {/* Reset Button */}
      <div style={{ marginTop: 'auto', paddingTop: '8px' }}>
        <button
          onClick={resetVisualization}
          style={{
            width: '100%',
            padding: '8px 12px',
            fontSize: '12px',
            background: 'var(--bg-tertiary)',
            border: '1px solid var(--border-color)',
            borderRadius: '4px',
            color: 'var(--text-secondary)',
            cursor: 'pointer',
          }}
        >
          Reset to Defaults
        </button>
      </div>
    </div>
  );
}

/**
 * Toggle checkbox item component
 */
interface ToggleItemProps {
  label: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
}

function ToggleItem({ label, checked, onChange }: ToggleItemProps) {
  return (
    <label style={{ 
      display: 'flex', 
      alignItems: 'center', 
      gap: '6px',
      fontSize: '12px',
      color: 'var(--text-secondary)',
      cursor: 'pointer',
      padding: '4px 0',
    }}>
      <input 
        type="checkbox" 
        checked={checked} 
        onChange={(e) => onChange(e.target.checked)}
        style={{ 
          width: '14px', 
          height: '14px',
          accentColor: 'var(--accent-primary)',
        }}
      />
      {label}
    </label>
  );
}

/**
 * Slider item component with label and value display
 */
interface SliderItemProps {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (value: number) => void;
  formatValue?: (value: number) => string;
}

function SliderItem({ label, value, min, max, step, onChange, formatValue }: SliderItemProps) {
  const displayValue = formatValue ? formatValue(value) : value.toString();
  
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
      <div style={{ 
        display: 'flex', 
        justifyContent: 'space-between', 
        alignItems: 'center',
      }}>
        <span style={{ 
          fontSize: '11px', 
          color: 'var(--text-secondary)',
        }}>
          {label}
        </span>
        <span style={{ 
          fontSize: '11px', 
          color: 'var(--text-primary)',
          fontFamily: 'monospace',
          minWidth: '50px',
          textAlign: 'right',
        }}>
          {displayValue}
        </span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        style={{
          width: '100%',
          height: '4px',
          accentColor: 'var(--accent-primary)',
          cursor: 'pointer',
        }}
      />
    </div>
  );
}
