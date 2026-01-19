/**
 * SolvePanel - Aerodynamic analysis controls
 * 
 * Provides:
 * - Single-point analysis (run at alpha or CL)
 * - Alpha polar generation
 * - Inviscid/Viscous mode switch
 * - Reynolds number control
 */

import { useState, useCallback, useMemo } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { 
  analyzeAirfoil, 
  analyzeViscousExtended,
  isWasmReady, 
  type AnalysisResult,
  type ViscousAnalysisResult,
  TurbulentModelNames,
  TurbulentModelDescriptions,
  type TurbulentModel
} from '../../lib/wasm';
import type { PolarPoint, SolverMode } from '../../types';

type RunMode = 'alpha' | 'cl';

export function SolvePanel() {
  const { 
    panels, 
    name, 
    polarData, 
    reynolds,
    solverMode,
    turbulentModel,
    nCrit,
    setDisplayAlpha, 
    setPolarData, 
    clearPolar,
    setReynolds,
    setSolverMode,
    setTurbulentModel,
    setNCrit,
    setViscousSolution,
    setBLData
  } = useAirfoilStore();
  
  // Solver settings
  const [runMode, setRunMode] = useState<RunMode>('alpha');
  
  // Single-point inputs
  const [targetAlpha, setTargetAlphaLocal] = useState(5.0);
  const [targetCl, setTargetCl] = useState(0.5);
  
  // Wrapper to update both local state and store displayAlpha
  const setTargetAlpha = useCallback((alpha: number) => {
    setTargetAlphaLocal(alpha);
    setDisplayAlpha(alpha);
  }, [setDisplayAlpha]);
  
  // Polar settings
  const [alphaStart, setAlphaStart] = useState(-5);
  const [alphaEnd, setAlphaEnd] = useState(15);
  const [alphaStep, setAlphaStep] = useState(1);
  
  // Results
  const [result, setResult] = useState<AnalysisResult | ViscousAnalysisResult | null>(null);
  const [isRunning, setIsRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  // Use store's polarData instead of local state
  const polar = polarData;

  // Helper to check if result is viscous
  const isViscousResult = (r: AnalysisResult | ViscousAnalysisResult | null): r is ViscousAnalysisResult => {
    return r !== null && 'cd' in r;
  };

  // Run single-point analysis
  const runAnalysis = useCallback(() => {
    if (!isWasmReady() || panels.length < 3) {
      setError('WASM not ready or insufficient geometry');
      return;
    }
    
    setIsRunning(true);
    setError(null);
    
    // Helper to run viscous analysis with current settings
    const runViscousAnalysis = (alpha: number) => 
      analyzeViscousExtended(panels, alpha, reynolds, turbulentModel, nCrit);
    
    try {
      if (runMode === 'alpha') {
        // Direct alpha analysis
        if (solverMode === 'viscous') {
          const res = runViscousAnalysis(targetAlpha);
          setResult(res);
          if (res.success) {
            setViscousSolution({
              cl: res.cl,
              cd: res.cd,
              cd_friction: res.cd_friction,
              cd_pressure: res.cd_pressure,
              cm: res.cm,
              cp: res.cp,
              cp_x: res.cp_x,
              x_tr_upper: res.x_tr_upper,
              x_tr_lower: res.x_tr_lower,
              converged: res.converged,
              iterations: res.iterations,
              reynolds: res.reynolds,
              alpha: res.alpha,
            });
          } else {
            setError(res.error || 'Viscous analysis failed');
          }
        } else {
          const res = analyzeAirfoil(panels, targetAlpha);
          setResult(res);
          if (!res.success) {
            setError(res.error || 'Analysis failed');
          }
        }
      } else {
        // Iterate to find alpha for target CL
        let alpha = 0;
        let cl = 0;
        const maxIter = 20;
        const tol = 0.001;
        
        for (let i = 0; i < maxIter; i++) {
          const res = solverMode === 'viscous' 
            ? runViscousAnalysis(alpha)
            : analyzeAirfoil(panels, alpha);
            
          if (!res.success) {
            setError(res.error || 'Analysis failed during CL iteration');
            setIsRunning(false);
            return;
          }
          
          cl = res.cl;
          const clError = targetCl - cl;
          
          if (Math.abs(clError) < tol) {
            setResult(res);
            if (solverMode === 'viscous' && isViscousResult(res)) {
              setViscousSolution({
                cl: res.cl,
                cd: res.cd,
                cd_friction: res.cd_friction,
                cd_pressure: res.cd_pressure,
                cm: res.cm,
                cp: res.cp,
                cp_x: res.cp_x,
                x_tr_upper: res.x_tr_upper,
                x_tr_lower: res.x_tr_lower,
                converged: res.converged,
                iterations: res.iterations,
                reynolds: res.reynolds,
                alpha: res.alpha,
              });
            }
            break;
          }
          
          // Estimate cl_alpha (thin airfoil: ~2π per radian ≈ 0.11 per degree)
          const clAlpha = 0.11;
          alpha += clError / clAlpha;
          
          // Clamp to reasonable range
          alpha = Math.max(-20, Math.min(25, alpha));
          
          if (i === maxIter - 1) {
            // Last iteration - use whatever we got
            const finalRes = solverMode === 'viscous'
              ? runViscousAnalysis(alpha)
              : analyzeAirfoil(panels, alpha);
            setResult(finalRes);
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
  }, [panels, runMode, targetAlpha, targetCl, solverMode, reynolds, turbulentModel, nCrit, setViscousSolution]);

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
      
      if (solverMode === 'viscous') {
        // Use extended viscous analysis with selected model and nCrit
        // Loop through alphas to support model/nCrit selection
        for (let alpha = alphaStart; alpha <= alphaEnd + 0.001; alpha += alphaStep) {
          const res = analyzeViscousExtended(panels, alpha, reynolds, turbulentModel, nCrit);
          if (res.success) {
            points.push({
              alpha: res.alpha,
              cl: res.cl,
              cd: res.cd,
              cm: res.cm,
              x_tr_upper: res.x_tr_upper,
              x_tr_lower: res.x_tr_lower,
              converged: res.converged,
              reynolds,
            });
          }
        }
      } else {
        // Inviscid polar
        for (let alpha = alphaStart; alpha <= alphaEnd; alpha += alphaStep) {
          const res = analyzeAirfoil(panels, alpha);
          if (res.success) {
            points.push({
              alpha,
              cl: res.cl,
              cm: res.cm,
            });
          }
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
  }, [panels, alphaStart, alphaEnd, alphaStep, solverMode, reynolds, turbulentModel, nCrit, clearPolar, setPolarData]);

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

  // Format Reynolds number for display
  const formatRe = (re: number): string => {
    if (re >= 1e6) return `${(re / 1e6).toFixed(2)}M`;
    if (re >= 1e3) return `${(re / 1e3).toFixed(0)}k`;
    return re.toFixed(0);
  };

  return (
    <div className="panel">
      <div className="panel-header">Solve</div>
      <div className="panel-content">
        {/* Solver Mode */}
        <div className="form-group">
          <div className="form-label">Solver Mode</div>
          <div style={{ display: 'flex', gap: '4px' }}>
            <button
              onClick={() => setSolverMode('inviscid')}
              className={solverMode === 'inviscid' ? 'active' : ''}
              style={{ flex: 1 }}
            >
              Inviscid
            </button>
            <button
              onClick={() => setSolverMode('viscous')}
              className={solverMode === 'viscous' ? 'active' : ''}
              style={{ flex: 1 }}
            >
              Viscous
            </button>
          </div>
        </div>

        {/* Reynolds Number (only for viscous) */}
        {solverMode === 'viscous' && (
          <div className="form-group">
            <div className="form-label">Reynolds Number</div>
            <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
              <input
                type="number"
                value={reynolds}
                onChange={(e) => setReynolds(parseFloat(e.target.value) || 1e6)}
                step={1e5}
                min={1e3}
                max={1e9}
                style={{ flex: 1 }}
              />
              <span style={{ fontSize: '11px', color: 'var(--text-muted)', minWidth: '50px' }}>
                {formatRe(reynolds)}
              </span>
            </div>
            {/* Quick Re presets */}
            <div style={{ display: 'flex', gap: '4px', marginTop: '4px' }}>
              {[1e5, 5e5, 1e6, 3e6, 6e6].map((re) => (
                <button
                  key={re}
                  onClick={() => setReynolds(re)}
                  style={{ 
                    flex: 1, 
                    padding: '2px 4px', 
                    fontSize: '10px',
                    opacity: reynolds === re ? 1 : 0.6 
                  }}
                >
                  {formatRe(re)}
                </button>
              ))}
            </div>
          </div>
        )}

        {/* Turbulent Model (only for viscous) */}
        {solverMode === 'viscous' && (
          <div className="form-group">
            <div className="form-label">Turbulent Model</div>
            <select
              value={turbulentModel}
              onChange={(e) => setTurbulentModel(parseInt(e.target.value) as TurbulentModel)}
              style={{ width: '100%', padding: '6px 8px' }}
            >
              {([0, 1, 2] as TurbulentModel[]).map((m) => (
                <option key={m} value={m}>
                  {TurbulentModelNames[m]} - {TurbulentModelDescriptions[m]}
                </option>
              ))}
            </select>
          </div>
        )}

        {/* N_crit (only for viscous) */}
        {solverMode === 'viscous' && (
          <div className="form-group">
            <div className="form-label" style={{ display: 'flex', justifyContent: 'space-between' }}>
              <span>N_crit (Transition)</span>
              <span style={{ fontWeight: 'normal', color: 'var(--text-secondary)' }}>{nCrit.toFixed(1)}</span>
            </div>
            <input
              type="range"
              min={1}
              max={14}
              step={0.5}
              value={nCrit}
              onChange={(e) => setNCrit(parseFloat(e.target.value))}
              style={{ width: '100%' }}
            />
            <div style={{ 
              display: 'flex', 
              justifyContent: 'space-between', 
              fontSize: '9px', 
              color: 'var(--text-muted)',
              marginTop: '2px'
            }}>
              <span>Turbulent (1)</span>
              <span>Clean (9)</span>
              <span>Sailplane (14)</span>
            </div>
          </div>
        )}

        {/* Single Point Analysis */}
        <div className="form-group">
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
              gridTemplateColumns: isViscousResult(result) ? '1fr 1fr 1fr' : '1fr 1fr',
              gap: '8px',
            }}>
              <ResultCard label="CL" value={result.cl.toFixed(4)} />
              {isViscousResult(result) && (
                <ResultCard label="CD" value={result.cd.toFixed(5)} highlight />
              )}
              <ResultCard label="CM" value={result.cm.toFixed(4)} />
            </div>
            {isViscousResult(result) && (
              <div style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr',
                gap: '8px',
                marginTop: '8px',
              }}>
                <ResultCard 
                  label="L/D" 
                  value={result.cd > 0 ? (result.cl / result.cd).toFixed(1) : '-'} 
                />
                <ResultCard 
                  label="Xtr (U/L)" 
                  value={`${(result.x_tr_upper * 100).toFixed(0)}/${(result.x_tr_lower * 100).toFixed(0)}%`} 
                />
              </div>
            )}
          </div>
        )}

        {/* Polar Sweep */}
        <div className="form-group" style={{ marginTop: '16px', paddingTop: '16px', borderTop: '1px solid var(--border-color)' }}>
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
                gridTemplateColumns: solverMode === 'viscous' ? '1fr 1fr 1fr 1fr' : '1fr 1fr 1fr',
                gap: '4px',
                padding: '4px 8px',
                background: 'var(--bg-tertiary)',
                borderBottom: '1px solid var(--border-color)',
                position: 'sticky',
                top: 0,
              }}>
                <span>α (°)</span>
                <span>CL</span>
                {solverMode === 'viscous' && <span>CD</span>}
                <span>CM</span>
              </div>
              
              {/* Data rows */}
              {polar.map((p, i) => (
                <div
                  key={i}
                  style={{
                    display: 'grid',
                    gridTemplateColumns: solverMode === 'viscous' ? '1fr 1fr 1fr 1fr' : '1fr 1fr 1fr',
                    gap: '4px',
                    padding: '2px 8px',
                    borderBottom: '1px solid var(--border-color)',
                  }}
                >
                  <span>{p.alpha.toFixed(1)}</span>
                  <span style={{ color: p.cl >= 0 ? 'var(--accent-primary)' : 'var(--accent-secondary)' }}>
                    {p.cl.toFixed(4)}
                  </span>
                  {solverMode === 'viscous' && (
                    <span style={{ color: 'var(--text-secondary)' }}>
                      {(p.cd ?? 0).toFixed(5)}
                    </span>
                  )}
                  <span>{p.cm.toFixed(4)}</span>
                </div>
              ))}
            </div>
            
            {/* Export polar */}
            <button
              onClick={() => {
                const header = solverMode === 'viscous'
                  ? `# ${name} - Viscous Polar (Re=${formatRe(reynolds)})\n# Alpha(deg)  CL        CD         CM\n`
                  : `# ${name} - Inviscid Polar\n# Alpha(deg)  CL        CM\n`;
                const data = polar.map(p => 
                  solverMode === 'viscous'
                    ? `${p.alpha.toFixed(2).padStart(8)} ${p.cl.toFixed(6).padStart(10)} ${(p.cd ?? 0).toFixed(7).padStart(11)} ${p.cm.toFixed(6).padStart(10)}`
                    : `${p.alpha.toFixed(2).padStart(8)} ${p.cl.toFixed(6).padStart(10)} ${p.cm.toFixed(6).padStart(10)}`
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
            >
              Export Polar
            </button>
          </div>
        )}

        {/* Error display */}
        {error && (
          <div style={{
            padding: '8px',
            background: 'rgba(255, 100, 100, 0.1)',
            border: '1px solid rgba(255, 100, 100, 0.3)',
            borderRadius: '4px',
            fontSize: '11px',
            color: '#ff6b6b',
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
          <div>Solver: {solverMode === 'inviscid' 
            ? 'Linear Vorticity Panel' 
            : `Viscous (Re=${formatRe(reynolds)}, ${TurbulentModelNames[turbulentModel]})`}</div>
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
      background: highlight ? 'rgba(0, 200, 150, 0.1)' : 'var(--bg-tertiary)',
      borderRadius: '4px',
      textAlign: 'center',
      border: highlight ? '1px solid rgba(0, 200, 150, 0.3)' : 'none',
    }}>
      <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginBottom: '2px' }}>
        {label}
      </div>
      <div style={{ fontFamily: 'var(--font-mono)', fontWeight: 600, fontSize: '14px' }}>
        {value}
      </div>
    </div>
  );
}
