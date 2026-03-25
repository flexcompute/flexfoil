/**
 * CfdPanel - Controls for the WebGPU CFD solver.
 *
 * Provides configuration for grid size, physics mode, Mach/alpha/Re,
 * start/stop/reset controls, and convergence history plot.
 */

import { useState, useCallback, useRef, useEffect } from 'react';
import { useCfdStore } from '../../stores/cfdStore';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { CfdSolver } from '../../lib/webgpu/cfd/CfdSolver';
import type { CfdSolverResult } from '../../lib/webgpu/cfd/CfdSolver';
import type { PhysicsMode } from '../../lib/webgpu/cfd/CfdConfig';
import { cfd_generate_mesh, cfd_initial_conditions, cfd_boundary_types } from 'rustfoil-wasm';

export function CfdPanel() {
  const {
    config,
    isRunning,
    iteration,
    convergenceHistory,
    cl,
    cd,
    cm,
    residualL2,
    meshGenerated,
    setNi,
    setNj,
    setMach,
    setAlphaDeg,
    setReynolds,
    setCfl,
    setPhysics,
    setRunning,
    setMeshGenerated,
    appendResult,
    reset,
  } = useCfdStore();

  const coordinates = useAirfoilStore((s) => s.coordinates);

  const solverRef = useRef<CfdSolver | null>(null);
  const runningRef = useRef(false);
  const [error, setError] = useState<string | null>(null);

  // Sync running ref with state
  useEffect(() => {
    runningRef.current = isRunning;
  }, [isRunning]);

  const handleGenerateMesh = useCallback(async () => {
    if (!coordinates || coordinates.length < 6) {
      setError('Load an airfoil first');
      return;
    }

    try {
      setError(null);
      // Flatten [x,y] pairs into a flat Float64Array for WASM
      const flat = new Float64Array(coordinates.length * 2);
      for (let i = 0; i < coordinates.length; i++) {
        flat[i * 2] = coordinates[i].x;
        flat[i * 2 + 1] = coordinates[i].y;
      }
      const meshResult = cfd_generate_mesh(
        flat,
        config.ni,
        config.nj,
        config.farField,
        config.ds0,
      );

      if (!meshResult || !meshResult.x) {
        setError('Mesh generation failed');
        return;
      }

      // Get initial conditions
      const physicsMode = config.physics === 'euler' ? 0 : config.physics === 'laminar_ns' ? 1 : 2;
      const initialQ = cfd_initial_conditions(
        config.ni,
        config.nj,
        config.machInf,
        config.alphaDeg,
        config.gamma,
        physicsMode,
        config.reynolds,
      );

      const bcTypes = cfd_boundary_types(config.ni, config.nj);

      // Check for WebGPU
      if (!navigator.gpu) {
        setError('WebGPU not available in this browser');
        return;
      }

      const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
      if (!adapter) {
        setError('No WebGPU adapter found');
        return;
      }

      const device = await adapter.requestDevice({
        requiredLimits: {
          maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize,
          maxComputeWorkgroupsPerDimension: adapter.limits.maxComputeWorkgroupsPerDimension,
        },
      });

      // Destroy previous solver if any
      solverRef.current?.destroy();

      solverRef.current = CfdSolver.create(
        device,
        config,
        new Float32Array(meshResult.x),
        new Float32Array(meshResult.y),
        new Float32Array(initialQ),
        new Uint32Array(bcTypes),
      );

      setMeshGenerated(true);
      reset();
      setMeshGenerated(true); // re-set after reset clears it
    } catch (e) {
      setError(`Mesh generation error: ${e instanceof Error ? e.message : String(e)}`);
    }
  }, [coordinates, config, setMeshGenerated, reset]);

  const handleStart = useCallback(async () => {
    if (!solverRef.current) {
      setError('Generate mesh first');
      return;
    }

    setError(null);
    setRunning(true);
    runningRef.current = true;

    const BATCH_SIZE = 20;

    try {
      while (runningRef.current) {
        const result: CfdSolverResult = await solverRef.current.step(BATCH_SIZE);
        appendResult(result);

        // Check convergence (residual dropped 6 orders of magnitude)
        if (result.residualL2 < 1e-10) {
          setRunning(false);
          break;
        }

        // Yield to UI
        await new Promise((r) => setTimeout(r, 0));
      }
    } catch (e) {
      setError(`Solver error: ${e instanceof Error ? e.message : String(e)}`);
      setRunning(false);
    }
  }, [appendResult, setRunning]);

  const handleStop = useCallback(() => {
    setRunning(false);
    runningRef.current = false;
  }, [setRunning]);

  const handleReset = useCallback(() => {
    handleStop();
    reset();
    solverRef.current?.destroy();
    solverRef.current = null;
    setMeshGenerated(false);
  }, [handleStop, reset, setMeshGenerated]);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      solverRef.current?.destroy();
    };
  }, []);

  return (
    <div className="panel-content" style={{ padding: '8px', fontSize: '13px' }}>
      <h3 style={{ margin: '0 0 8px', fontSize: '14px' }}>CFD Solver (WebGPU)</h3>

      {error && (
        <div style={{ color: '#ff4444', marginBottom: 8, fontSize: '12px' }}>
          {error}
        </div>
      )}

      {/* Grid Settings */}
      <fieldset style={{ border: '1px solid #444', padding: '6px', marginBottom: '8px' }}>
        <legend style={{ fontSize: '12px', color: '#aaa' }}>Grid</legend>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '4px' }}>
          <label>
            ni: <input type="number" value={config.ni} onChange={(e) => setNi(+e.target.value)}
              min={32} max={1024} step={32} style={inputStyle} disabled={isRunning} />
          </label>
          <label>
            nj: <input type="number" value={config.nj} onChange={(e) => setNj(+e.target.value)}
              min={16} max={512} step={16} style={inputStyle} disabled={isRunning} />
          </label>
        </div>
      </fieldset>

      {/* Flow Conditions */}
      <fieldset style={{ border: '1px solid #444', padding: '6px', marginBottom: '8px' }}>
        <legend style={{ fontSize: '12px', color: '#aaa' }}>Flow Conditions</legend>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '4px' }}>
          <label>
            Mach: <input type="number" value={config.machInf} onChange={(e) => setMach(+e.target.value)}
              min={0.01} max={2.0} step={0.01} style={inputStyle} disabled={isRunning} />
          </label>
          <label>
            Alpha: <input type="number" value={config.alphaDeg} onChange={(e) => setAlphaDeg(+e.target.value)}
              min={-20} max={20} step={0.5} style={inputStyle} disabled={isRunning} />
          </label>
          <label>
            Re: <input type="number" value={config.reynolds} onChange={(e) => setReynolds(+e.target.value)}
              min={1e3} max={1e8} step={1e5} style={inputStyle} disabled={isRunning} />
          </label>
          <label>
            CFL: <input type="number" value={config.cfl} onChange={(e) => setCfl(+e.target.value)}
              min={0.01} max={10.0} step={0.1} style={inputStyle} disabled={isRunning} />
          </label>
        </div>
      </fieldset>

      {/* Physics Mode */}
      <fieldset style={{ border: '1px solid #444', padding: '6px', marginBottom: '8px' }}>
        <legend style={{ fontSize: '12px', color: '#aaa' }}>Physics</legend>
        <select
          value={config.physics}
          onChange={(e) => setPhysics(e.target.value as PhysicsMode)}
          style={inputStyle}
          disabled={isRunning}
        >
          <option value="euler">Euler (Inviscid)</option>
          <option value="laminar_ns">Laminar Navier-Stokes</option>
          <option value="rans_sa">RANS (Spalart-Allmaras)</option>
        </select>
      </fieldset>

      {/* Controls */}
      <div style={{ display: 'flex', gap: '4px', marginBottom: '8px' }}>
        <button
          onClick={handleGenerateMesh}
          disabled={isRunning || !coordinates?.length}
          style={btnStyle}
        >
          Generate Mesh
        </button>
        {!isRunning ? (
          <button onClick={handleStart} disabled={!meshGenerated} style={btnStyle}>
            Start
          </button>
        ) : (
          <button onClick={handleStop} style={{ ...btnStyle, background: '#aa4444' }}>
            Stop
          </button>
        )}
        <button onClick={handleReset} style={btnStyle}>
          Reset
        </button>
      </div>

      {/* Results */}
      <fieldset style={{ border: '1px solid #444', padding: '6px', marginBottom: '8px' }}>
        <legend style={{ fontSize: '12px', color: '#aaa' }}>Results</legend>
        <div style={{ fontFamily: 'monospace', fontSize: '12px', lineHeight: 1.5 }}>
          <div>Iteration: {iteration}</div>
          <div>Cl: {cl.toFixed(6)}</div>
          <div>Cd: {cd.toFixed(6)}</div>
          <div>Cm: {cm.toFixed(6)}</div>
          <div>Residual L2: {residualL2.toExponential(3)}</div>
        </div>
      </fieldset>

      {/* Convergence mini-plot (text-based for now) */}
      {convergenceHistory.length > 0 && (
        <fieldset style={{ border: '1px solid #444', padding: '6px' }}>
          <legend style={{ fontSize: '12px', color: '#aaa' }}>Convergence</legend>
          <div style={{ fontFamily: 'monospace', fontSize: '11px', maxHeight: '80px', overflow: 'auto' }}>
            {convergenceHistory.slice(-10).map((p) => (
              <div key={p.iteration}>
                [{p.iteration.toString().padStart(6)}] L2={p.residualL2.toExponential(3)}
              </div>
            ))}
          </div>
        </fieldset>
      )}
    </div>
  );
}

const inputStyle: React.CSSProperties = {
  width: '70px',
  background: '#333',
  color: '#eee',
  border: '1px solid #555',
  borderRadius: '3px',
  padding: '2px 4px',
  fontSize: '12px',
};

const btnStyle: React.CSSProperties = {
  flex: 1,
  padding: '4px 8px',
  background: '#555',
  color: '#eee',
  border: '1px solid #666',
  borderRadius: '3px',
  cursor: 'pointer',
  fontSize: '12px',
};
