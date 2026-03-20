/**
 * AirfoilLibraryPanel - NACA generator, Selig database browser, and file import
 */

import { useCallback, useEffect, useRef, useState } from 'react';
import { parseAirfoilDat } from '../../lib/airfoilImport';
import {
  fetchSeligAirfoil,
  getRandomCatalogEntry,
  getSeligCatalog,
  type SeligCatalogEntry,
} from '../../lib/airfoilDatabase';
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

  // Selig database state
  const [catalog, setCatalog] = useState<SeligCatalogEntry[]>([]);
  const [seligSearch, setSeligSearch] = useState('');
  const [seligLoading, setSeligLoading] = useState<string | null>(null);
  const [seligError, setSeligError] = useState<string | null>(null);

  useEffect(() => {
    getSeligCatalog().then(setCatalog);
  }, []);

  const filteredCatalog = seligSearch.length >= 2
    ? catalog.filter((entry) => {
        const q = seligSearch.toLowerCase();
        return entry.file.toLowerCase().includes(q) || entry.name.toLowerCase().includes(q);
      })
    : [];

  const loadSeligAirfoil = useCallback(async (entry: SeligCatalogEntry) => {
    setSeligLoading(entry.file);
    setSeligError(null);
    try {
      const parsed = await fetchSeligAirfoil(entry.file);
      importAirfoil(parsed.name, parsed.coordinates);
    } catch (err) {
      setSeligError(err instanceof Error ? err.message : 'Failed to load airfoil');
    } finally {
      setSeligLoading(null);
    }
  }, [importAirfoil]);

  const handleRandomFoil = useCallback(async () => {
    if (catalog.length === 0) return;
    const entry = getRandomCatalogEntry(catalog);
    await loadSeligAirfoil(entry);
  }, [catalog, loadSeligAirfoil]);

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

        {/* Selig Database */}
        <div className="form-group" data-tour="library-selig">
          <div className="form-label" style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <span>Selig Database</span>
            <span style={{ fontSize: '10px', color: 'var(--text-tertiary)', fontWeight: 400 }}>
              {catalog.length > 0 ? `${catalog.length} airfoils` : ''}
            </span>
          </div>
          <input
            type="text"
            value={seligSearch}
            onChange={(e) => setSeligSearch(e.target.value)}
            placeholder="Search airfoils (e.g. e387, clark, naca)…"
            style={{ width: '100%', boxSizing: 'border-box', marginBottom: '6px' }}
          />
          {seligSearch.length >= 2 && (
            <div style={{
              maxHeight: '180px',
              overflowY: 'auto',
              border: '1px solid var(--border-color)',
              borderRadius: '4px',
              marginBottom: '6px',
            }}>
              {filteredCatalog.length === 0 ? (
                <div style={{ padding: '8px 10px', fontSize: '11px', color: 'var(--text-tertiary)' }}>
                  No matches
                </div>
              ) : (
                filteredCatalog.slice(0, 100).map((entry) => (
                  <button
                    key={entry.file}
                    onClick={() => loadSeligAirfoil(entry)}
                    disabled={seligLoading !== null}
                    style={{
                      display: 'block',
                      width: '100%',
                      textAlign: 'left',
                      padding: '5px 10px',
                      fontSize: '11px',
                      lineHeight: '1.4',
                      background: seligLoading === entry.file ? 'var(--bg-tertiary)' : 'transparent',
                      border: 'none',
                      borderBottom: '1px solid var(--border-color)',
                      cursor: seligLoading ? 'wait' : 'pointer',
                      color: 'var(--text-primary)',
                    }}
                  >
                    <span style={{ fontFamily: 'var(--font-mono)', color: 'var(--accent-primary)' }}>
                      {entry.file}
                    </span>
                    {' '}
                    <span style={{ color: 'var(--text-secondary)' }}>{entry.name}</span>
                    {seligLoading === entry.file && (
                      <span style={{ marginLeft: '6px', color: 'var(--text-tertiary)' }}>loading…</span>
                    )}
                  </button>
                ))
              )}
              {filteredCatalog.length > 100 && (
                <div style={{ padding: '6px 10px', fontSize: '10px', color: 'var(--text-tertiary)' }}>
                  {filteredCatalog.length - 100} more — narrow your search
                </div>
              )}
            </div>
          )}
          {seligError && (
            <div style={{ fontSize: '11px', color: 'var(--error)', marginBottom: '6px' }}>
              {seligError}
            </div>
          )}
          <button
            data-tour="library-random-foil"
            onClick={handleRandomFoil}
            disabled={catalog.length === 0 || seligLoading !== null}
            style={{
              width: '100%',
              fontStyle: catalog.length === 0 ? 'italic' : undefined,
            }}
          >
            {seligLoading ? 'Loading…' : 'Random Foil'}
          </button>
          <div style={{ fontSize: '10px', color: 'var(--text-tertiary)', marginTop: '6px', textAlign: 'center' }}>
            Data from the{' '}
            <a
              href="https://m-selig.ae.illinois.edu/ads/coord_database.html"
              target="_blank"
              rel="noopener noreferrer"
              style={{ color: 'var(--text-secondary)' }}
            >
              UIUC Airfoil Coordinates Database
            </a>
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
              accept=".dat,.txt,.csv,text/plain,text/*"
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
