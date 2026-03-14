/**
 * SolvePanel - Aerodynamic analysis controls
 * 
 * Provides:
 * - Single-point faithful viscous analysis
 * - Alpha polar generation
 */

import { useState, useCallback, useMemo, useEffect, useRef } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { analyzeAirfoil, isWasmReady, type AnalysisResult } from '../../lib/wasm';
import type { PolarPoint } from '../../types';

type RunMode = 'alpha' | 'cl';

export function SolvePanel() {
  const {
    panels,
    name,
    polarData,
    displayAlpha,
    reynolds,
    setDisplayAlpha,
    setReynolds,
    setPolarData,
    clearPolar
  } = useAirfoilStore();
  
  const [runMode, setRunMode] = useState<RunMode>('alpha');
  
  // Single-point inputs - initialize from store's displayAlpha
  const [targetAlpha, setTargetAlphaLocal] = useState(displayAlpha);
  const [targetCl, setTargetCl] = useState(0.5);
  
  // Debounce timer ref
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  
  // Sync local targetAlpha with store on mount (respect store's initial value)
  const initializedRef = useRef(false);
  useEffect(() => {
    if (!initializedRef.current) {
      initializedRef.current = true;
      // On first load, update store to match initial targetAlpha if different
      if (displayAlpha !== targetAlpha) {
        setDisplayAlpha(targetAlpha);
      }
    }
  }, [displayAlpha, targetAlpha, setDisplayAlpha]);
  
  // Wrapper to update both local state and store displayAlpha (debounced)
  const setTargetAlpha = useCallback((alpha: number) => {
    setTargetAlphaLocal(alpha);
    
    // Debounce the store update to avoid excessive recomputation while typing
    if (debounceRef.current) {
      clearTimeout(debounceRef.current);
    }
    debounceRef.current = setTimeout(() => {
      setDisplayAlpha(alpha);
    }, 300); // 300ms debounce
  }, [setDisplayAlpha]);
  
  // Cleanup debounce timer on unmount
  useEffect(() => {
    return () => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
    };
  }, []);
  
  // Polar settings
  const [alphaStart, setAlphaStart] = useState(-5);
  const [alphaEnd, setAlphaEnd] = useState(15);
  const [alphaStep, setAlphaStep] = useState(1);
  
  // Results
  const [result, setResult] = useState<AnalysisResult | null>(null);
  const [resultAlpha, setResultAlpha] = useState<number | null>(null); // Alpha that was solved for
  const [isRunning, setIsRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  // Use store's polarData instead of local state
  const polar = polarData;

  // Run single-point analysis
  const runAnalysis = useCallback(() => {
    if (!isWasmReady() || panels.length < 3) {
      setError('WASM not ready or insufficient geometry');
      return;
    }
    
    setIsRunning(true);
    setError(null);
    setResultAlpha(null);
    
    try {
      if (runMode === 'alpha') {
        // Direct alpha analysis
        const res = analyzeAirfoil(panels, targetAlpha, reynolds);
        setResult(res);
        setResultAlpha(targetAlpha);
        if (!res.success) {
          setError(res.error || 'Analysis failed');
        } else if (!res.converged) {
          setError(`Faithful viscous solver did not converge (residual ${res.residual.toExponential(2)})`);
        }
      } else {
        // Iterate to find alpha for target CL
        let alpha = 0;
        let cl = 0;
        const maxIter = 20;
        const tol = 0.001;
        
        for (let i = 0; i < maxIter; i++) {
          const res = analyzeAirfoil(panels, alpha, reynolds);
          if (!res.success) {
            setError(res.error || 'Analysis failed during CL iteration');
            setIsRunning(false);
            return;
          }
          
          cl = res.cl;
          const clError = targetCl - cl;
          
          if (Math.abs(clError) < tol) {
            setResult(res);
            setResultAlpha(alpha);
            // Update display alpha to show the found alpha
            setDisplayAlpha(alpha);
            break;
          }
          
          // Estimate cl_alpha (thin airfoil: ~2π per radian ≈ 0.11 per degree)
          const clAlpha = 0.11;
          alpha += clError / clAlpha;
          
          // Clamp to reasonable range
          alpha = Math.max(-20, Math.min(25, alpha));
          
          if (i === maxIter - 1) {
            // Last iteration - use whatever we got
            const finalRes = analyzeAirfoil(panels, alpha, reynolds);
            setResult(finalRes);
            setResultAlpha(alpha);
            // Update display alpha even if not fully converged
            setDisplayAlpha(alpha);
            if (Math.abs(targetCl - finalRes.cl) > 0.01) {
              setError(`Could not converge to CL=${targetCl.toFixed(3)}. Got CL=${finalRes.cl.toFixed(3)} at α=${alpha.toFixed(2)}°`);
            }
          }
        }
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Unknown error');
    }
    
    setIsRunning(false);
  }, [panels, runMode, targetAlpha, targetCl, reynolds, setDisplayAlpha]);

  // Run polar sweep
  const runPolar = useCallback(() => {
    if (!isWasmReady() || panels.length < 3) {
      setError('WASM not ready or insufficient geometry');
      return;
    }
    
    setIsRunning(true);
    setError(null);
    clearPolar();
    
    try {
      const points: PolarPoint[] = [];
      
      for (let alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep) {
        const res = analyzeAirfoil(panels, alpha, reynolds);
        if (res.success) {
          points.push({
            alpha,
            cl: res.cl,
            cd: res.cd,
            cm: res.cm,
          });
        }
      }
      
      setPolarData(points);
      
      if (points.length === 0) {
        setError('No valid points in polar');
      }
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Unknown error');
    }
    
    setIsRunning(false);
  }, [panels, alphaStart, alphaEnd, alphaStep, reynolds, clearPolar, setPolarData]);

  // Compute cl_alpha from polar
  const clAlpha = useMemo(() => {
    if (polar.length < 2) return null;
    
    // Linear regression in the linear region (-5 to 10 deg typically)
    const linearPoints = polar.filter(p => p.alpha >= -5 && p.alpha <= 10);
    if (linearPoints.length < 2) return null;
    
    const n = linearPoints.length;
    const sumX = linearPoints.reduce((s, p) => s + p.alpha, 0);
    const sumY = linearPoints.reduce((s, p) => s + p.cl, 0);
    const sumXY = linearPoints.reduce((s, p) => s + p.alpha * p.cl, 0);
    const sumX2 = linearPoints.reduce((s, p) => s + p.alpha * p.alpha, 0);
    
    const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    return slope; // per degree
  }, [polar]);

  return (
    <div className="panel">
      <div className="panel-header">Solve</div>
      <div className="panel-content">
        <div className="form-group" data-tour="solve-mode">
          <div className="form-label">Solver</div>
          <div style={{ fontSize: '12px', color: 'var(--text-primary)' }}>
            Faithful XFOIL viscous solver
          </div>
          <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '4px' }}>
            Frontend flow analysis now uses the faithful Rust reproduction path.
          </div>
        </div>

        <div className="form-group">
          <div className="form-label">Reynolds Number</div>
          <input
            type="number"
            value={reynolds}
            onChange={(e) => setReynolds(Math.max(1, Number(e.target.value) || 1))}
            step={100000}
          />
        </div>

        {/* Single Point Analysis */}
        <div className="form-group" data-tour="solve-alpha">
          <div className="form-label">Single Point</div>
          
          {/* Run mode toggle */}
          <div style={{ display: 'flex', gap: '4px', marginBottom: '8px' }}>
            <button
              onClick={() => setRunMode('alpha')}
              className={runMode === 'alpha' ? 'active' : ''}
              style={{ flex: 1 }}
            >
              Run to α
            </button>
            <button
              onClick={() => setRunMode('cl')}
              className={runMode === 'cl' ? 'active' : ''}
              style={{ flex: 1 }}
            >
              Run to CL
            </button>
          </div>
          
          {/* Input field */}
          <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
            <label style={{ fontSize: '12px', minWidth: '60px' }}>
              {runMode === 'alpha' ? 'Alpha (°):' : 'Target CL:'}
            </label>
            <input
              type="number"
              value={runMode === 'alpha' ? targetAlpha : targetCl}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (runMode === 'alpha') {
                  setTargetAlpha(val);
                } else {
                  setTargetCl(val);
                }
              }}
              step={runMode === 'alpha' ? 0.5 : 0.1}
              style={{ flex: 1 }}
            />
            <button
              onClick={runAnalysis}
              disabled={isRunning || !isWasmReady()}
              style={{ minWidth: '60px' }}
              data-tour="solve-run"
            >
              {isRunning ? '...' : 'Run'}
            </button>
          </div>
        </div>

        {/* Single Point Results */}
        {result && result.success && (
          <div className="form-group">
            <div className="form-label">Results</div>
            <div style={{
              display: 'grid',
              gridTemplateColumns: runMode === 'cl' ? 'repeat(4, 1fr)' : 'repeat(3, 1fr)',
              gap: '8px',
            }}>
              {runMode === 'cl' && resultAlpha !== null && (
                <ResultCard label="α (°)" value={resultAlpha.toFixed(2)} highlight />
              )}
              <ResultCard label="CL" value={result.cl.toFixed(4)} />
              <ResultCard label="CD" value={result.cd.toFixed(5)} />
              <ResultCard label="CM" value={result.cm.toFixed(4)} />
            </div>
            <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '8px' }}>
              {result.converged
                ? `Converged in ${result.iterations} iterations`
                : `Not converged after ${result.iterations} iterations (residual ${result.residual.toExponential(2)})`}
            </div>
          </div>
        )}

        {/* Polar Sweep */}
        <div className="form-group" style={{ marginTop: '16px', paddingTop: '16px', borderTop: '1px solid var(--border-color)' }} data-tour="solve-polar">
          <div className="form-label">Alpha Polar</div>
          
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '8px', marginBottom: '8px' }}>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Start (°)</label>
              <input
                type="number"
                value={alphaStart}
                onChange={(e) => setAlphaStart(parseFloat(e.target.value))}
                step={1}
              />
            </div>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>End (°)</label>
              <input
                type="number"
                value={alphaEnd}
                onChange={(e) => setAlphaEnd(parseFloat(e.target.value))}
                step={1}
              />
            </div>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Step (°)</label>
              <input
                type="number"
                value={alphaStep}
                onChange={(e) => setAlphaStep(parseFloat(e.target.value))}
                step={0.5}
                min={0.1}
              />
            </div>
          </div>
          
          <button
            onClick={runPolar}
            disabled={isRunning || !isWasmReady()}
            style={{ width: '100%' }}
          >
            {isRunning ? 'Running...' : 'Generate Polar'}
          </button>
        </div>

        {/* Polar Results */}
        {polar.length > 0 && (
          <div className="form-group">
            <div className="form-label">
              Polar Results ({polar.length} points)
              {clAlpha !== null && (
                <span style={{ fontWeight: 'normal', marginLeft: '8px' }}>
                  CL_α = {(clAlpha * 180 / Math.PI).toFixed(3)}/rad
                </span>
              )}
            </div>
            
            {/* Mini polar table */}
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
                gridTemplateColumns: '1fr 1fr 1fr',
                gap: '4px',
                padding: '4px 8px',
                background: 'var(--bg-tertiary)',
                borderBottom: '1px solid var(--border-color)',
                position: 'sticky',
                top: 0,
              }}>
                <span>α (°)</span>
                <span>CL</span>
                <span>CM</span>
              </div>
              
              {/* Data rows */}
              {polar.map((p, i) => (
                <div
                  key={i}
                  style={{
                    display: 'grid',
                    gridTemplateColumns: '1fr 1fr 1fr',
                    gap: '4px',
                    padding: '2px 8px',
                    borderBottom: '1px solid var(--border-color)',
                  }}
                >
                  <span>{p.alpha.toFixed(1)}</span>
                  <span style={{ color: p.cl >= 0 ? 'var(--accent-primary)' : 'var(--accent-secondary)' }}>
                    {p.cl.toFixed(4)}
                  </span>
                  <span>{p.cm.toFixed(4)}</span>
                </div>
              ))}
            </div>
            
            {/* Export polar */}
            <button
              onClick={() => {
                const header = `# ${name} - Inviscid Polar\n# Alpha(deg)  CL        CM\n`;
                const data = polar.map(p => 
                  `${p.alpha.toFixed(2).padStart(8)} ${p.cl.toFixed(6).padStart(10)} ${p.cm.toFixed(6).padStart(10)}`
                ).join('\n');
                const blob = new Blob([header + data], { type: 'text/plain' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = `${name.replace(/\s+/g, '_')}_polar.txt`;
                a.click();
                URL.revokeObjectURL(url);
              }}
              style={{ width: '100%', marginTop: '8px' }}
              data-tour="solve-export"
            >
              Export Polar
            </button>
          </div>
        )}

        {/* Error display */}
        {error && (
          <div style={{
            padding: '8px',
            background: 'var(--bg-tertiary)',
            border: '1px solid var(--accent-danger)',
            borderRadius: '4px',
            fontSize: '11px',
            color: 'var(--accent-danger)',
            marginTop: '8px',
          }}>
            {error}
          </div>
        )}

        {/* Status */}
        <div style={{
          marginTop: '12px',
          padding: '8px',
          background: 'var(--bg-tertiary)',
          borderRadius: '4px',
          fontSize: '10px',
          color: 'var(--text-muted)',
        }}>
          <div>Panels: {panels.length - 1}</div>
          <div>Solver: Faithful XFOIL viscous</div>
          <div>WASM: {isWasmReady() ? '✓ Ready' : '⏳ Loading...'}</div>
        </div>
      </div>
    </div>
  );
}

function ResultCard({ label, value, highlight }: { label: string; value: string; highlight?: boolean }) {
  return (
    <div style={{
      padding: '8px',
      background: highlight ? 'var(--accent-primary)' : 'var(--bg-tertiary)',
      borderRadius: '4px',
      textAlign: 'center',
    }}>
      <div style={{ fontSize: '10px', color: highlight ? 'var(--bg-primary)' : 'var(--text-muted)', marginBottom: '2px' }}>
        {label}
      </div>
      <div style={{ 
        fontFamily: 'var(--font-mono)', 
        fontWeight: 600, 
        fontSize: '14px',
        color: highlight ? 'var(--bg-primary)' : 'inherit',
      }}>
        {value}
      </div>
    </div>
  );
}
