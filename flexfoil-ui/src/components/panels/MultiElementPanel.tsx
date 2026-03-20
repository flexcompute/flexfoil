/**
 * MultiElementPanel - Configure and analyze multi-element airfoil configurations.
 *
 * Provides controls for:
 * - Adding/removing airfoil elements (bodies)
 * - Positioning each element (x, y, rotation)
 * - Running multi-element inviscid analysis
 * - Viewing per-body results
 */

import { useState, useCallback, useRef } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { isWasmReady } from '../../lib/wasm';
import type { MultiElementBody, ElementPosition } from '../../types';
import { parseAirfoilDat } from '../../lib/airfoilImport';

const fieldsetStyle: React.CSSProperties = {
  border: '1px solid var(--border-color)', borderRadius: 4, padding: '8px', marginBottom: '10px',
};
const legendStyle: React.CSSProperties = { fontSize: '11px', fontWeight: 600, color: 'var(--text-secondary)', padding: '0 4px' };

function PositionRow({ label, value, onChange, step = 0.01, min, max }: {
  label: string;
  value: number;
  onChange: (v: number) => void;
  step?: number;
  min?: number;
  max?: number;
}) {
  return (
    <div style={{ display: 'flex', gap: '4px', alignItems: 'center', marginBottom: '3px' }}>
      <span style={{ fontSize: '10px', color: 'var(--text-muted)', minWidth: '50px' }}>{label}</span>
      <input
        type="number"
        value={Number(value.toFixed(4))}
        onChange={(e) => onChange(parseFloat(e.target.value) || 0)}
        step={step}
        min={min}
        max={max}
        style={{ flex: 1, fontSize: '11px', fontFamily: 'var(--font-mono)' }}
      />
    </div>
  );
}

function BodyCard({ body, isSelected, onSelect, onUpdatePosition, onRename, onRemove }: {
  body: MultiElementBody;
  isSelected: boolean;
  onSelect: () => void;
  onUpdatePosition: (pos: Partial<ElementPosition>) => void;
  onRename: (name: string) => void;
  onRemove: () => void;
}) {
  return (
    <div
      onClick={onSelect}
      style={{
        border: `1px solid ${isSelected ? body.color : 'var(--border-color)'}`,
        borderRadius: 4,
        padding: '6px',
        marginBottom: '6px',
        background: isSelected ? 'var(--bg-tertiary)' : 'transparent',
        cursor: 'pointer',
      }}
    >
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '4px' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
          <span style={{
            width: 10, height: 10, borderRadius: '50%',
            background: body.color, display: 'inline-block', flexShrink: 0,
          }} />
          <input
            type="text"
            value={body.name}
            onChange={(e) => onRename(e.target.value)}
            onClick={(e) => e.stopPropagation()}
            style={{
              fontSize: '11px', fontWeight: 600, padding: '1px 3px',
              background: 'transparent', border: '1px solid transparent',
              borderRadius: 2, color: 'var(--text-secondary)', width: '90px',
            }}
            onFocus={(e) => { e.target.style.borderColor = 'var(--border-color)'; }}
            onBlur={(e) => { e.target.style.borderColor = 'transparent'; }}
          />
        </div>
        <button
          onClick={(e) => { e.stopPropagation(); onRemove(); }}
          title="Remove element"
          style={{
            background: 'transparent', border: 'none', color: 'var(--text-muted)',
            cursor: 'pointer', fontSize: '14px', padding: '0 2px',
          }}
        >
          ×
        </button>
      </div>

      {isSelected && (
        <div onClick={(e) => e.stopPropagation()} style={{ marginTop: '6px' }}>
          <PositionRow label="X offset" value={body.position.x} onChange={(v) => onUpdatePosition({ x: v })} />
          <PositionRow label="Y offset" value={body.position.y} onChange={(v) => onUpdatePosition({ y: v })} />
          <PositionRow label="Angle (°)" value={body.position.angle} onChange={(v) => onUpdatePosition({ angle: v })} step={1} min={-90} max={90} />
          <div style={{ fontSize: '9px', color: 'var(--text-muted)', marginTop: '2px' }}>
            {body.panels.length - 1} panels &middot; {body.coordinates.length} points
          </div>
        </div>
      )}
    </div>
  );
}

export function MultiElementPanel() {
  const {
    multiElement,
    addBody,
    removeBody,
    updateBodyPosition,
    renameBody,
    selectBody,
    promoteToMultiElement,
    runMultiElementAnalysis,
    displayAlpha,
    name,
    coordinates,
  } = useAirfoilStore();

  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleImportDat = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files;
    if (!files) return;

    for (const file of Array.from(files)) {
      try {
        const text = await file.text();
        const parsed = parseAirfoilDat(text, file.name);
        if (parsed && parsed.coordinates.length >= 3) {
          addBody(parsed.name, parsed.coordinates.map(p => ({ x: p.x, y: p.y })));
        }
      } catch (err) {
        console.warn(`Failed to import ${file.name}:`, err);
      }
    }
    e.target.value = '';
  }, [addBody]);

  const handleRunAnalysis = useCallback(() => {
    if (!isWasmReady()) return;
    setIsAnalyzing(true);
    setError(null);
    try {
      const result = runMultiElementAnalysis(displayAlpha);
      if (result && !result.success) {
        setError(result.error || 'Analysis failed');
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Unknown error');
    } finally {
      setIsAnalyzing(false);
    }
  }, [displayAlpha, runMultiElementAnalysis]);

  const { bodies, selectedBodyId, lastResult, enabled } = multiElement;

  return (
    <div className="panel">
      <div className="panel-header">Multi-Element</div>
      <div className="panel-content">

        {!enabled && (
          <div style={{ marginBottom: '10px' }}>
            <div style={{ fontSize: '11px', color: 'var(--text-muted)', marginBottom: '8px' }}>
              Analyze multi-element airfoil configurations (slat + main + flap) using an inviscid panel method.
            </div>
            <button
              onClick={() => {
                promoteToMultiElement();
              }}
              style={{ width: '100%', marginBottom: '4px' }}
              disabled={coordinates.length < 3}
            >
              Start from current airfoil
            </button>
            <div style={{ fontSize: '9px', color: 'var(--text-muted)', textAlign: 'center' }}>
              Current airfoil ({name}) becomes the first element
            </div>
          </div>
        )}

        {enabled && (
          <>
            {/* Element list */}
            <fieldset style={fieldsetStyle}>
              <legend style={legendStyle}>Elements ({bodies.length})</legend>
              {bodies.map((body) => (
                <BodyCard
                  key={body.id}
                  body={body}
                  isSelected={body.id === selectedBodyId}
                  onSelect={() => selectBody(body.id)}
                  onUpdatePosition={(pos) => updateBodyPosition(body.id, pos)}
                  onRename={(n) => renameBody(body.id, n)}
                  onRemove={() => removeBody(body.id)}
                />
              ))}

              <div style={{ display: 'flex', gap: '4px' }}>
                <button
                  onClick={() => fileInputRef.current?.click()}
                  style={{ flex: 1, fontSize: '11px' }}
                >
                  + Import .dat
                </button>
                <input
                  ref={fileInputRef}
                  type="file"
                  accept=".dat,.txt"
                  multiple
                  onChange={handleImportDat}
                  style={{ display: 'none' }}
                />
              </div>
            </fieldset>

            {/* Analysis */}
            <fieldset style={fieldsetStyle}>
              <legend style={legendStyle}>Analysis</legend>
              <div style={{ display: 'flex', gap: '4px', alignItems: 'center', marginBottom: '6px' }}>
                <span style={{ fontSize: '10px', color: 'var(--text-muted)', minWidth: '50px' }}>Alpha (°)</span>
                <span style={{ fontSize: '11px', fontFamily: 'var(--font-mono)' }}>{displayAlpha.toFixed(1)}</span>
              </div>
              <button
                onClick={handleRunAnalysis}
                disabled={isAnalyzing || !isWasmReady() || bodies.length < 1}
                style={{ width: '100%', marginBottom: '6px' }}
              >
                {isAnalyzing ? 'Analyzing...' : 'Run Multi-Element (Inviscid)'}
              </button>
              <div style={{ fontSize: '9px', color: 'var(--text-muted)' }}>
                Uses alpha from the Solve panel. Inviscid only.
              </div>
            </fieldset>

            {/* Results */}
            {lastResult && lastResult.success && (
              <fieldset style={fieldsetStyle}>
                <legend style={legendStyle}>Results</legend>
                <div style={{
                  display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '4px', marginBottom: '8px',
                }}>
                  <ResultCard label="CL total" value={lastResult.cl_total.toFixed(4)} />
                  <ResultCard label="CM total" value={lastResult.cm_total.toFixed(4)} />
                </div>
                <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginBottom: '4px' }}>Per-body:</div>
                {lastResult.per_body.map((pb, i) => (
                  <div key={i} style={{
                    display: 'flex', justifyContent: 'space-between', alignItems: 'center',
                    padding: '2px 4px', fontSize: '11px', fontFamily: 'var(--font-mono)',
                    borderBottom: '1px solid var(--border-color)',
                  }}>
                    <span style={{ color: bodies[i]?.color ?? 'var(--text-primary)' }}>{pb.name}</span>
                    <span>CL = {pb.cl.toFixed(4)}</span>
                  </div>
                ))}
              </fieldset>
            )}

            {error && (
              <div style={{
                padding: '6px', background: 'var(--bg-tertiary)',
                border: '1px solid var(--accent-danger)', borderRadius: 4,
                fontSize: '10px', color: 'var(--accent-danger)',
              }}>
                {error}
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
}

function ResultCard({ label, value }: { label: string; value: string }) {
  return (
    <div style={{
      padding: '6px', background: 'var(--bg-tertiary)', borderRadius: 4, textAlign: 'center',
    }}>
      <div style={{ fontSize: '9px', color: 'var(--text-muted)', marginBottom: '2px' }}>{label}</div>
      <div style={{ fontFamily: 'var(--font-mono)', fontWeight: 600, fontSize: '13px' }}>{value}</div>
    </div>
  );
}
