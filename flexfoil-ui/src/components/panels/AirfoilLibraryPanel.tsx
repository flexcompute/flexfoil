/**
 * AirfoilLibraryPanel - NACA generator and file import
 */

import { useCallback, useRef, useState } from 'react';
import { parseAirfoilDat } from '../../lib/airfoilImport';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useRouteUiStore } from '../../stores/routeUiStore';

export function AirfoilLibraryPanel() {
  const { name, setName, generateNaca4, importAirfoil, reset } = useAirfoilStore();
  const [editingName, setEditingName] = useState(false);
  const [draftName, setDraftName] = useState('');
  const nameInputRef = useRef<HTMLInputElement>(null);
  const nacaCode = useRouteUiStore((state) => state.libraryNacaCode);
  const nPoints = useRouteUiStore((state) => state.libraryNPoints);
  const setNacaCode = useRouteUiStore((state) => state.setLibraryNacaCode);
  const setNPoints = useRouteUiStore((state) => state.setLibraryNPoints);

  const handleGenerateNaca = useCallback(() => {
    // Parse NACA code (e.g., "2412" -> m=0.02, p=0.4, t=0.12)
    if (nacaCode.length !== 4) {
      alert('Please enter a valid 4-digit NACA code (e.g., 0012, 2412)');
      return;
    }

    const m = parseInt(nacaCode[0], 10) / 100; // Max camber
    const p = parseInt(nacaCode[1], 10) / 10;  // Position of max camber
    const t = parseInt(nacaCode.slice(2), 10) / 100; // Max thickness

    generateNaca4({ m, p, t, nPoints });
  }, [nacaCode, nPoints, generateNaca4]);

  const handleFileImport = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const text = String(event.target?.result ?? '');
        const parsed = parseAirfoilDat(text, file.name);
        importAirfoil(parsed.name, parsed.coordinates);
      } catch (error) {
        const message = error instanceof Error ? error.message : 'Unknown import error.';
        window.alert(`Could not import ${file.name}: ${message}`);
      } finally {
        e.target.value = '';
      }
    };
    reader.readAsText(file);
  }, [importAirfoil]);

  return (
    <div className="panel">
      <div className="panel-header">Airfoil Library</div>
      <div className="panel-content">
        {/* Current airfoil — click to rename */}
        <div className="form-group">
          <div className="form-label">Current Airfoil</div>
          {editingName ? (
            <input
              ref={nameInputRef}
              type="text"
              value={draftName}
              onChange={(e) => setDraftName(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === 'Enter') {
                  const trimmed = draftName.trim();
                  if (trimmed) setName(trimmed);
                  setEditingName(false);
                } else if (e.key === 'Escape') {
                  setEditingName(false);
                }
              }}
              onBlur={() => {
                const trimmed = draftName.trim();
                if (trimmed) setName(trimmed);
                setEditingName(false);
              }}
              autoFocus
              style={{
                padding: '7px 12px',
                background: 'var(--bg-tertiary)',
                borderRadius: '4px',
                fontWeight: 600,
                color: 'var(--accent-primary)',
                border: '1px solid var(--accent-primary)',
                width: '100%',
                boxSizing: 'border-box',
              }}
            />
          ) : (
            <div
              onClick={() => { setDraftName(name); setEditingName(true); }}
              title="Click to rename"
              style={{
                padding: '8px 12px',
                background: 'var(--bg-tertiary)',
                borderRadius: '4px',
                fontWeight: 600,
                color: 'var(--accent-primary)',
                cursor: 'text',
              }}
            >
              {name}
            </div>
          )}
        </div>

        {/* NACA Generator */}
        <div className="form-group">
          <div className="form-label">NACA 4-Series Generator</div>
          <div className="form-row" style={{ marginBottom: '8px' }}>
            <input
              type="text"
              value={nacaCode}
              onChange={(e) => setNacaCode(e.target.value.replace(/\D/g, '').slice(0, 4))}
              placeholder="0012"
              style={{ width: '80px', fontFamily: 'var(--font-mono)', textAlign: 'center' }}
            />
            <span style={{ color: 'var(--text-secondary)', fontSize: '12px' }}>
              Code
            </span>
          </div>
          <div className="form-row" style={{ marginBottom: '8px' }}>
            <input
              type="number"
              defaultValue={nPoints}
              key={nPoints}
              onBlur={(e) => {
                const val = parseInt(e.target.value);
                if (!isNaN(val) && val >= 10 && val <= 500) setNPoints(val);
                else e.target.value = String(nPoints);
              }}
              onKeyDown={(e) => { if (e.key === 'Enter') (e.target as HTMLInputElement).blur(); }}
              min={10}
              max={500}
              style={{ width: '80px' }}
            />
            <span style={{ color: 'var(--text-secondary)', fontSize: '12px' }}>
              Points
            </span>
          </div>
          <button onClick={handleGenerateNaca} className="primary" style={{ width: '100%' }}>
            Generate NACA {nacaCode || '____'}
          </button>
        </div>

        {/* NACA presets */}
        <div className="form-group">
          <div className="form-label">Quick Presets</div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
            {['0012', '2412', '4412', '0015', '2415', '6412'].map((code) => (
              <button
                key={code}
                onClick={() => {
                  setNacaCode(code);
                  const m = parseInt(code[0], 10) / 100;
                  const p = parseInt(code[1], 10) / 10;
                  const t = parseInt(code.slice(2), 10) / 100;
                  generateNaca4({ m, p, t, nPoints });
                }}
                style={{ 
                  padding: '4px 8px', 
                  fontSize: '11px',
                  fontFamily: 'var(--font-mono)',
                  background: code === nacaCode ? 'var(--accent-primary)' : undefined,
                  color: code === nacaCode ? 'var(--bg-primary)' : undefined,
                }}
              >
                {code}
              </button>
            ))}
          </div>
        </div>

        {/* File import */}
        <div className="form-group">
          <div className="form-label">Import from File</div>
          <label style={{ 
            display: 'block',
            padding: '8px',
            border: '1px dashed var(--border-color)',
            borderRadius: '4px',
            textAlign: 'center',
            cursor: 'pointer',
            fontSize: '12px',
            color: 'var(--text-secondary)',
          }}>
            <input
              type="file"
              accept=".dat,.txt,.csv"
              onChange={handleFileImport}
              style={{ display: 'none' }}
            />
            Click to import .dat file
          </label>
        </div>

        {/* Reset */}
        <button onClick={reset} style={{ width: '100%', marginTop: '12px' }} data-tour="library-reset">
          Reset to Default
        </button>
      </div>
    </div>
  );
}
