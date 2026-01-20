/**
 * VisualizationPanel - Controls for flow visualization options
 * 
 * Contains:
 * - Display toggles (Grid, Curve, Panels, Points, Controls, Streamlines, Smoke)
 * - Aerodynamic overlays (Cp pressure, Force vectors)
 * - Streamline options (density, adaptive)
 * - Smoke options (density, particles per blob, spawn interval, max age)
 * - Animation options (morphing enable, duration)
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
    showPsiContours,
    showCp,
    showForces,
    setShowGrid,
    setShowCurve,
    setShowPanels,
    setShowPoints,
    setShowControls,
    setShowStreamlines,
    setShowSmoke,
    setShowPsiContours,
    setShowCp,
    setShowForces,
    
    // Animation options
    enableMorphing,
    morphDuration,
    setEnableMorphing,
    setMorphDuration,
    
    // Streamline options
    streamlineDensity,
    adaptiveStreamlines,
    setStreamlineDensity,
    setAdaptiveStreamlines,
    
    // Smoke options (spawn interval and max age are now auto-calculated)
    smokeDensity,
    smokeParticlesPerBlob,
    setSmokeDensity,
    setSmokeParticlesPerBlob,
    
    // Flow speed
    flowSpeed,
    setFlowSpeed,
    
    // Cp options
    cpDisplayMode,
    cpBarScale,
    setCpDisplayMode,
    setCpBarScale,
    
    // Force options
    forceScale,
    setForceScale,
    
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
          <ToggleItem label="Stream Function (ψ)" checked={showPsiContours} onChange={setShowPsiContours} />
          <ToggleItem label="Smoke" checked={showSmoke} onChange={setShowSmoke} />
        </div>
      </section>

      {/* Aerodynamic Overlays Section */}
      <section>
        <h4 style={{ 
          margin: '0 0 8px 0', 
          fontSize: '12px', 
          fontWeight: 600,
          color: 'var(--text-primary)',
          textTransform: 'uppercase',
          letterSpacing: '0.5px',
        }}>
          Aerodynamic Overlays
        </h4>
        <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
          <ToggleItem label="Pressure (Cp)" checked={showCp} onChange={setShowCp} />
          <ToggleItem label="Force Vectors" checked={showForces} onChange={setShowForces} />
        </div>
      </section>

      {/* Cp Options Section */}
      {showCp && (
        <section>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            fontSize: '12px', 
            fontWeight: 600,
            color: 'var(--text-primary)',
            textTransform: 'uppercase',
            letterSpacing: '0.5px',
          }}>
            Cp Display Options
          </h4>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
              <span style={{ fontSize: '11px', color: 'var(--text-secondary)' }}>Mode</span>
              <select 
                value={cpDisplayMode}
                onChange={(e) => setCpDisplayMode(e.target.value as 'color' | 'bars' | 'both')}
                style={{
                  padding: '6px 8px',
                  fontSize: '12px',
                  background: 'var(--bg-tertiary)',
                  border: '1px solid var(--border-color)',
                  borderRadius: '4px',
                  color: 'var(--text-primary)',
                }}
              >
                <option value="color">Color Only</option>
                <option value="bars">Bars Only</option>
                <option value="both">Color + Bars</option>
              </select>
            </div>
            <SliderItem
              label="Bar Scale"
              value={cpBarScale}
              min={0.01}
              max={0.5}
              step={0.01}
              onChange={setCpBarScale}
              formatValue={(v) => `${(v * 100).toFixed(0)}%`}
            />
          </div>
        </section>
      )}

      {/* Force Vector Options Section */}
      {showForces && (
        <section>
          <h4 style={{ 
            margin: '0 0 8px 0', 
            fontSize: '12px', 
            fontWeight: 600,
            color: 'var(--text-primary)',
            textTransform: 'uppercase',
            letterSpacing: '0.5px',
          }}>
            Force Vector Options
          </h4>
          <SliderItem
            label="Scale"
            value={forceScale}
            min={0.05}
            max={0.5}
            step={0.01}
            onChange={setForceScale}
            formatValue={(v) => `${(v * 100).toFixed(0)}%`}
          />
        </section>
      )}

      {/* Animation Section */}
      <section>
        <h4 style={{ 
          margin: '0 0 8px 0', 
          fontSize: '12px', 
          fontWeight: 600,
          color: 'var(--text-primary)',
          textTransform: 'uppercase',
          letterSpacing: '0.5px',
        }}>
          Animation
        </h4>
        <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
          <ToggleItem label="Smooth Morphing" checked={enableMorphing} onChange={setEnableMorphing} />
          {enableMorphing && (
            <SliderItem
              label="Duration"
              value={morphDuration}
              min={50}
              max={1000}
              step={50}
              onChange={setMorphDuration}
              formatValue={(v) => `${v}ms`}
            />
          )}
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
            {/* Spawn Interval and Max Age are now auto-calculated based on visible domain */}
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
