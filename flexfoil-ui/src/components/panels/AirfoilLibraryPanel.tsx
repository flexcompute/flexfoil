/**
 * AirfoilLibraryPanel - NACA generator and file import
 */

import { useState, useCallback } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';

export function AirfoilLibraryPanel() {
  const { name, generateNaca4, reset } = useAirfoilStore();
  
  // NACA 4-series parameters
  const [nacaCode, setNacaCode] = useState('0012');
  const [nPoints, setNPoints] = useState(50);

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
      const text = event.target?.result as string;
      // TODO: Parse airfoil coordinate file (Selig format)
      console.log('Imported file:', file.name, text.slice(0, 100));
      alert('File import not yet implemented. Use NACA generator for now.');
    };
    reader.readAsText(file);
  }, []);

  return (
    <div className="panel">
      <div className="panel-header">Airfoil Library</div>
      <div className="panel-content">
        {/* Current airfoil */}
        <div className="form-group">
          <div className="form-label">Current Airfoil</div>
          <div style={{ 
            padding: '8px 12px', 
            background: 'var(--bg-tertiary)', 
            borderRadius: '4px',
            fontWeight: 600,
            color: 'var(--accent-primary)',
          }}>
            {name}
          </div>
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
              value={nPoints}
              onChange={(e) => setNPoints(Math.max(10, Math.min(200, parseInt(e.target.value) || 50)))}
              min={10}
              max={200}
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
        <button onClick={reset} style={{ width: '100%', marginTop: '12px' }}>
          Reset to Default
        </button>
      </div>
    </div>
  );
}
