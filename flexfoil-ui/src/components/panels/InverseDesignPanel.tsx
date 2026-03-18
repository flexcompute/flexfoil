/**
 * InverseDesignPanel - Interactive target Cp/velocity editor for inverse design (QDES).
 *
 * Users draw target distributions on upper/lower surfaces, then execute
 * the solver to get a modified airfoil that matches the target.
 */

import React, { useRef, useCallback, useMemo, useState, useEffect } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useVisualizationStore } from '../../stores/visualizationStore';
import { analyzeAirfoil, analyzeAirfoilInviscid, runInverseDesign, runFullInverseDesign, repanelXfoil, isWasmReady } from '../../lib/wasm';
import { useSolverJobStore } from '../../stores/solverJobStore';
import type { InverseDesignSurfaceTarget, AirfoilPoint } from '../../types';

type InverseMethod = 'qdes' | 'mdes';

const WIDTH = 500;
const HEIGHT = 340;
const MARGIN = { top: 20, right: 20, bottom: 40, left: 50 };
const PLOT_W = WIDTH - MARGIN.left - MARGIN.right;
const PLOT_H = HEIGHT - MARGIN.top - MARGIN.bottom;

const CP_MIN = -3.0;
const CP_MAX = 1.5;
const X_MIN = 0;
const X_MAX = 1;

function toSvgX(x: number): number {
  return MARGIN.left + ((x - X_MIN) / (X_MAX - X_MIN)) * PLOT_W;
}
function toSvgY(cp: number): number {
  return MARGIN.top + ((cp - CP_MAX) / (CP_MIN - CP_MAX)) * PLOT_H;
}
function toDataX(svgX: number): number {
  return X_MIN + ((svgX - MARGIN.left) / PLOT_W) * (X_MAX - X_MIN);
}
function toDataCp(svgY: number): number {
  return CP_MAX + ((svgY - MARGIN.top) / PLOT_H) * (CP_MIN - CP_MAX);
}

function getSvgCoords(svg: SVGSVGElement, clientX: number, clientY: number) {
  const ctm = svg.getScreenCTM();
  if (!ctm) {
    const rect = svg.getBoundingClientRect();
    const vb = svg.viewBox.baseVal;
    return {
      svgX: (clientX - rect.left) * vb.width / rect.width,
      svgY: (clientY - rect.top) * vb.height / rect.height,
    };
  }
  const pt = svg.createSVGPoint();
  pt.x = clientX;
  pt.y = clientY;
  const svgPt = pt.matrixTransform(ctm.inverse());
  return { svgX: svgPt.x, svgY: svgPt.y };
}

function ticks(min: number, max: number, step: number): number[] {
  const arr: number[] = [];
  for (let v = min; v <= max + step * 0.01; v += step) arr.push(Math.round(v * 1000) / 1000);
  return arr;
}

interface TargetCurveEditorProps {
  surface: 'upper' | 'lower';
  target: InverseDesignSurfaceTarget | null;
  achieved: InverseDesignSurfaceTarget | null;
  currentX: number[];
  currentValues: number[];
  onChange: (target: InverseDesignSurfaceTarget | null) => void;
  color: string;
  achievedColor: string;
}

function polyline(xs: number[], ys: number[]): string {
  if (xs.length === 0) return '';
  return xs.map((x, i) => `${toSvgX(x)},${toSvgY(ys[i])}`).join(' ');
}

const TargetCurveEditor: React.FC<TargetCurveEditorProps> = ({
  surface, target, achieved, currentX, currentValues, onChange, color, achievedColor,
}) => {
  if (!target) return null;

  return (
    <g>
      {/* Current distribution (thin, dashed) */}
      {currentX.length > 1 && (
        <polyline
          points={polyline(currentX, currentValues)}
          fill="none"
          stroke="var(--text-muted)"
          strokeWidth={1}
          strokeDasharray="4 3"
          opacity={0.6}
        />
      )}
      {/* Achieved distribution (if available) */}
      {achieved && achieved.x.length > 1 && (
        <polyline
          points={polyline(achieved.x, achieved.values)}
          fill="none"
          stroke={achievedColor}
          strokeWidth={1.5}
          strokeDasharray="6 3"
        />
      )}
      {/* Target distribution (solid, editable) */}
      {target.x.length > 1 && (
        <polyline
          points={polyline(target.x, target.values)}
          fill="none"
          stroke={color}
          strokeWidth={2}
        />
      )}
      {/* Draggable target knots */}
      {target.x.map((x, i) => (
        <circle
          key={`${surface}-${i}`}
          cx={toSvgX(x)}
          cy={toSvgY(target.values[i])}
          r={5}
          fill={color}
          stroke="var(--bg-primary)"
          strokeWidth={1.5}
          style={{ cursor: 'grab', pointerEvents: 'all' }}
        />
      ))}
    </g>
  );
};

export function InverseDesignPanel() {
  const svgRef = useRef<SVGSVGElement>(null);
  const [dragging, setDragging] = useState<{ surface: 'upper' | 'lower'; index: number } | null>(null);
  const [method, setMethod] = useState<InverseMethod>('qdes');
  const [mdesSymmetric, setMdesSymmetric] = useState(false);
  const [mdesFilter, setMdesFilter] = useState(0.0);
  const [mdesSolving, setMdesSolving] = useState(false);
  const [mdesResult, setMdesResultLocal] = useState<{ x: number[]; y: number[] } | null>(null);

  const {
    panels, displayAlpha, reynolds, mach, ncrit, maxIterations, solverMode,
    inverseDesign, setInverseDesignTarget, setInverseDesignTargetKind,
    setInverseDesignSolving, setInverseDesignResult, setInverseDesignMaxIterations,
    setInverseDesignDamping, applyInverseDesignResult, clearInverseDesign,
  } = useAirfoilStore();

  const inverseDesignRef = useRef(inverseDesign);
  inverseDesignRef.current = inverseDesign;

  // Compute current Cp distribution to show as reference
  const currentCp = useMemo(() => {
    if (panels.length < 5) return { upperX: [] as number[], upperCp: [] as number[], lowerX: [] as number[], lowerCp: [] as number[] };
    try {
      const result = solverMode === 'viscous'
        ? analyzeAirfoil(panels, displayAlpha, reynolds, mach, ncrit, maxIterations)
        : analyzeAirfoilInviscid(panels, displayAlpha);
      
      if (!result.success || !result.cp_x || result.cp_x.length === 0) {
        return { upperX: [], upperCp: [], lowerX: [], lowerCp: [] };
      }
      
      const leIdx = result.cp_x.reduce((minI, x, i, arr) => x < arr[minI] ? i : minI, 0);
      
      const upperX: number[] = [];
      const upperCp: number[] = [];
      const lowerX: number[] = [];
      const lowerCp: number[] = [];
      
      for (let i = 0; i < result.cp_x.length; i++) {
        if (i <= leIdx) {
          upperX.push(result.cp_x[i]);
          upperCp.push(result.cp[i]);
        } else {
          lowerX.push(result.cp_x[i]);
          lowerCp.push(result.cp[i]);
        }
      }
      
      return { upperX, upperCp, lowerX, lowerCp };
    } catch {
      return { upperX: [], upperCp: [], lowerX: [], lowerCp: [] };
    }
  }, [panels, displayAlpha, reynolds, mach, ncrit, maxIterations, solverMode]);

  const initFromCurrent = useCallback((surface: 'upper' | 'lower') => {
    const xs = surface === 'upper' ? currentCp.upperX : currentCp.lowerX;
    const vals = surface === 'upper' ? currentCp.upperCp : currentCp.lowerCp;
    if (xs.length < 2) return;
    
    // Sort by increasing x (upper surface arrives in TE→LE / decreasing-x order
    // from the solver, but targets must be monotonically increasing for QDES
    // validation and for the drag-constraint logic to work correctly)
    const pairs = xs.map((x, i) => ({ x, v: vals[i] }));
    pairs.sort((a, b) => a.x - b.x);
    
    // Subsample to ~15 editable points
    const step = Math.max(1, Math.floor(pairs.length / 15));
    const sampledX: number[] = [];
    const sampledV: number[] = [];
    for (let i = 0; i < pairs.length; i += step) {
      sampledX.push(pairs[i].x);
      sampledV.push(pairs[i].v);
    }
    if (sampledX[sampledX.length - 1] !== pairs[pairs.length - 1].x) {
      sampledX.push(pairs[pairs.length - 1].x);
      sampledV.push(pairs[pairs.length - 1].v);
    }
    
    setInverseDesignTarget(surface, { x: sampledX, values: sampledV });
  }, [currentCp, setInverseDesignTarget]);

  // Auto-initialize upper target from current Cp on first render if none exist
  const hasAutoInited = useRef(false);
  useEffect(() => {
    if (hasAutoInited.current) return;
    if (!inverseDesign.upperTarget && !inverseDesign.lowerTarget && currentCp.upperX.length > 2) {
      hasAutoInited.current = true;
      initFromCurrent('upper');
    }
  }, [inverseDesign.upperTarget, inverseDesign.lowerTarget, currentCp.upperX.length, initFromCurrent]);

  const handleMouseDown = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current) return;
    const { svgX, svgY } = getSvgCoords(svgRef.current, e.clientX, e.clientY);
    const design = inverseDesignRef.current;
    
    // Check if clicking on a target knot
    for (const surface of ['upper', 'lower'] as const) {
      const target = surface === 'upper' ? design.upperTarget : design.lowerTarget;
      if (!target) continue;
      for (let i = 0; i < target.x.length; i++) {
        const kx = toSvgX(target.x[i]);
        const ky = toSvgY(target.values[i]);
        if (Math.hypot(svgX - kx, svgY - ky) < 12) {
          if (e.shiftKey && target.x.length > 2) {
            const newX = [...target.x];
            const newV = [...target.values];
            newX.splice(i, 1);
            newV.splice(i, 1);
            setInverseDesignTarget(surface, { x: newX, values: newV });
          } else {
            setDragging({ surface, index: i });
          }
          e.preventDefault();
          return;
        }
      }
    }
    
    // Double-click to add a point
    if (e.detail === 2) {
      const x = toDataX(svgX);
      const cp = toDataCp(svgY);
      if (x < X_MIN || x > X_MAX) return;
      
      const surface = design.upperTarget ? 'upper' : design.lowerTarget ? 'lower' : null;
      if (!surface) return;
      
      const target = surface === 'upper' ? design.upperTarget : design.lowerTarget;
      if (!target) return;
      
      const newX = [...target.x, x];
      const newV = [...target.values, cp];
      const sorted = newX.map((xi, i) => ({ x: xi, v: newV[i] })).sort((a, b) => a.x - b.x);
      setInverseDesignTarget(surface, { x: sorted.map(s => s.x), values: sorted.map(s => s.v) });
    }
  }, [setInverseDesignTarget]);

  // Window-level drag listeners (fires even when mouse leaves SVG)
  useEffect(() => {
    if (!dragging) return;

    const onMove = (e: MouseEvent) => {
      if (!svgRef.current) return;
      const { svgX, svgY } = getSvgCoords(svgRef.current, e.clientX, e.clientY);
      const design = inverseDesignRef.current;
      const target = dragging.surface === 'upper' ? design.upperTarget : design.lowerTarget;
      if (!target) return;

      const newValues = [...target.values];
      newValues[dragging.index] = Math.max(CP_MIN, Math.min(CP_MAX, toDataCp(svgY)));

      const newX = [...target.x];
      if (dragging.index > 0 && dragging.index < target.x.length - 1) {
        const xVal = Math.max(
          target.x[dragging.index - 1] + 0.001,
          Math.min(target.x[dragging.index + 1] - 0.001, toDataX(svgX))
        );
        newX[dragging.index] = xVal;
      }

      setInverseDesignTarget(dragging.surface, { x: newX, values: newValues });
    };

    const onUp = () => setDragging(null);

    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
    return () => {
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseup', onUp);
    };
  }, [dragging, setInverseDesignTarget]);

  const jobDispatch = useSolverJobStore.getState().dispatch;
  const jobComplete = useSolverJobStore.getState().complete;

  const handleExecute = useCallback(() => {
    const design = inverseDesignRef.current;
    if (!design.upperTarget && !design.lowerTarget) return;
    
    const { id: jobId, signal } = jobDispatch('QDES Inverse Design');
    setInverseDesignSolving(true);
    
    setTimeout(() => {
      if (signal.aborted) { setInverseDesignSolving(false); return; }
      try {
        const result = runInverseDesign(panels, {
          alphaDeg: displayAlpha,
          reynolds,
          mach,
          ncrit,
          targetKind: design.targetKind,
          upper: design.upperTarget ? {
            x: design.upperTarget.x,
            values: design.upperTarget.values,
          } : undefined,
          lower: design.lowerTarget ? {
            x: design.lowerTarget.x,
            values: design.lowerTarget.values,
          } : undefined,
          maxDesignIterations: design.maxIterations,
          damping: design.damping,
        });
        
        if (result.success) {
          const points = [];
          for (let i = 0; i < result.output_coords.length; i += 2) {
            points.push({ x: result.output_coords[i], y: result.output_coords[i + 1] });
          }
          setInverseDesignResult({
            resultCoords: points,
            converged: result.converged,
            achievedUpper: result.achieved_upper_x.length > 0
              ? { x: result.achieved_upper_x, values: result.achieved_upper_values }
              : null,
            achievedLower: result.achieved_lower_x.length > 0
              ? { x: result.achieved_lower_x, values: result.achieved_lower_values }
              : null,
            history: result.history,
          });
          jobComplete(jobId);
        } else {
          setInverseDesignSolving(false);
          jobComplete(jobId, result.error || 'QDES failed');
        }
      } catch (e) {
        setInverseDesignSolving(false);
        jobComplete(jobId, e instanceof Error ? e.message : 'QDES error');
      }
    }, 10);
  }, [panels, displayAlpha, reynolds, mach, ncrit, setInverseDesignSolving, setInverseDesignResult, jobDispatch, jobComplete]);

  const handleExecuteMdes = useCallback(() => {
    const { id: jobId, signal } = jobDispatch('MDES Full-Inverse Design');
    setMdesSolving(true);
    setTimeout(() => {
      if (signal.aborted) { setMdesSolving(false); return; }
      try {
        const result = runFullInverseDesign(panels, {
          alphaDeg: displayAlpha,
          symmetric: mdesSymmetric,
          filterStrength: mdesFilter,
        });
        if (result.success && result.x.length > 0) {
          setMdesResultLocal({ x: result.x, y: result.y });
          const points = result.x.map((xi, i) => ({ x: xi, y: result.y[i] }));
          setInverseDesignResult({
            resultCoords: points,
            converged: true,
            achievedUpper: null,
            achievedLower: null,
            history: [],
          });
          jobComplete(jobId);
        } else {
          jobComplete(jobId, result.error || 'MDES failed');
        }
      } catch (e) {
        jobComplete(jobId, e instanceof Error ? e.message : 'MDES error');
      } finally {
        setMdesSolving(false);
      }
    }, 10);
  }, [panels, displayAlpha, mdesSymmetric, mdesFilter, setInverseDesignResult, jobDispatch, jobComplete]);

  const xTicks = ticks(0, 1, 0.2);
  const cpTicks = ticks(CP_MIN, CP_MAX, 0.5);

  return (
    <div className="inverse-design-panel" style={{ display: 'flex', flexDirection: 'column', flex: 1, minHeight: 0, overflow: 'auto' }}>
      {/* WIP banner */}
      <div style={{
        padding: '6px 10px',
        background: 'color-mix(in srgb, var(--warning-color, #f59e0b) 12%, transparent)',
        borderBottom: '1px solid color-mix(in srgb, var(--warning-color, #f59e0b) 30%, transparent)',
        fontSize: '11px',
        color: 'var(--text-secondary)',
        lineHeight: 1.4,
      }}>
        <strong style={{ color: 'var(--warning-color, #f59e0b)' }}>Experimental</strong> — QDES and MDES are works in progress and may produce unexpected results.
      </div>

      {/* Method selector */}
      <div style={{ padding: '6px 10px', display: 'flex', gap: '4px', borderBottom: '1px solid var(--border-color)' }}>
        <button
          className={`btn btn-xs ${method === 'qdes' ? 'btn-primary' : ''}`}
          onClick={() => setMethod('qdes')}
          title="Mixed-inverse: specify target Cp on a segment"
        >
          QDES (Mixed)
        </button>
        <button
          className={`btn btn-xs ${method === 'mdes' ? 'btn-primary' : ''}`}
          onClick={() => setMethod('mdes')}
          title="Full-inverse: circle-plane conformal mapping"
        >
          MDES (Full)
        </button>
      </div>
      
      {/* Controls */}
      <div style={{ padding: '8px 10px', display: 'flex', flexWrap: 'wrap', gap: '6px', alignItems: 'center', borderBottom: '1px solid var(--border-color)' }}>
        {method === 'qdes' && (
        <select
          value={inverseDesign.targetKind}
          onChange={(e) => setInverseDesignTargetKind(e.target.value as 'cp' | 'velocity')}
          style={{ fontSize: '11px', padding: '2px 4px' }}
        >
          <option value="cp">Target Cp</option>
          <option value="velocity">Target Ue</option>
        </select>
        )}
        
        <button
          className="btn btn-xs"
          onClick={() => initFromCurrent('upper')}
          title="Initialize upper target from current Cp"
        >
          + Upper
        </button>
        <button
          className="btn btn-xs"
          onClick={() => initFromCurrent('lower')}
          title="Initialize lower target from current Cp"
        >
          + Lower
        </button>
        
        {inverseDesign.upperTarget && (
          <button className="btn btn-xs" onClick={() => setInverseDesignTarget('upper', null)} title="Clear upper target">
            × Upper
          </button>
        )}
        {inverseDesign.lowerTarget && (
          <button className="btn btn-xs" onClick={() => setInverseDesignTarget('lower', null)} title="Clear lower target">
            × Lower
          </button>
        )}
      </div>
      
      {/* SVG Plot */}
      <div style={{ flex: 1, minHeight: 200, padding: '4px' }}>
        <svg
          ref={svgRef}
          viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
          style={{ width: '100%', height: '100%' }}
          onMouseDown={handleMouseDown}
        >
          {/* Background */}
          <rect x={MARGIN.left} y={MARGIN.top} width={PLOT_W} height={PLOT_H} fill="var(--bg-secondary)" rx={2} />
          
          {/* Grid */}
          {xTicks.map(x => (
            <line key={`gx-${x}`} x1={toSvgX(x)} y1={MARGIN.top} x2={toSvgX(x)} y2={MARGIN.top + PLOT_H}
              stroke="var(--border-color)" strokeWidth={0.5} opacity={0.4} />
          ))}
          {cpTicks.map(cp => (
            <line key={`gy-${cp}`} x1={MARGIN.left} y1={toSvgY(cp)} x2={MARGIN.left + PLOT_W} y2={toSvgY(cp)}
              stroke="var(--border-color)" strokeWidth={cp === 0 ? 1 : 0.5} opacity={cp === 0 ? 0.8 : 0.4} />
          ))}
          
          {/* Axes labels */}
          {xTicks.map(x => (
            <text key={`lx-${x}`} x={toSvgX(x)} y={HEIGHT - 5} textAnchor="middle" fontSize={10} fill="var(--text-muted)">
              {x.toFixed(1)}
            </text>
          ))}
          {cpTicks.map(cp => (
            <text key={`ly-${cp}`} x={MARGIN.left - 6} y={toSvgY(cp) + 3.5} textAnchor="end" fontSize={10} fill="var(--text-muted)">
              {cp.toFixed(1)}
            </text>
          ))}
          
          <text x={WIDTH / 2} y={HEIGHT - 22} textAnchor="middle" fontSize={11} fill="var(--text-secondary)">x/c</text>
          <text x={12} y={HEIGHT / 2} textAnchor="middle" fontSize={11} fill="var(--text-secondary)"
            transform={`rotate(-90, 12, ${HEIGHT / 2})`}>
            {inverseDesign.targetKind === 'cp' ? '-Cp' : 'Ue'}
          </text>
          
          {/* Empty state guidance */}
          {!inverseDesign.upperTarget && !inverseDesign.lowerTarget && (
            <>
              <text x={MARGIN.left + PLOT_W / 2} y={MARGIN.top + PLOT_H / 2 - 10}
                textAnchor="middle" fontSize={13} fill="var(--text-muted)">
                No target distribution set
              </text>
              <text x={MARGIN.left + PLOT_W / 2} y={MARGIN.top + PLOT_H / 2 + 12}
                textAnchor="middle" fontSize={11} fill="var(--text-muted)" opacity={0.7}>
                Click &quot;+ Upper&quot; or &quot;+ Lower&quot; above to begin
              </text>
            </>
          )}
          
          {/* Curves */}
          <TargetCurveEditor
            surface="upper"
            target={inverseDesign.upperTarget}
            achieved={inverseDesign.achievedUpper}
            currentX={currentCp.upperX}
            currentValues={currentCp.upperCp}
            onChange={(t) => setInverseDesignTarget('upper', t)}
            color="#ef4444"
            achievedColor="#f97316"
          />
          <TargetCurveEditor
            surface="lower"
            target={inverseDesign.lowerTarget}
            achieved={inverseDesign.achievedLower}
            currentX={currentCp.lowerX}
            currentValues={currentCp.lowerCp}
            onChange={(t) => setInverseDesignTarget('lower', t)}
            color="#3b82f6"
            achievedColor="#06b6d4"
          />
          
          {/* Legend */}
          <g transform={`translate(${MARGIN.left + 8}, ${MARGIN.top + 8})`}>
            <line x1={0} y1={0} x2={16} y2={0} stroke="var(--text-muted)" strokeDasharray="4 3" strokeWidth={1} />
            <text x={20} y={4} fontSize={9} fill="var(--text-muted)">Current</text>
            
            <line x1={0} y1={14} x2={16} y2={14} stroke="#ef4444" strokeWidth={2} />
            <text x={20} y={18} fontSize={9} fill="var(--text-secondary)">Upper target</text>
            
            <line x1={0} y1={28} x2={16} y2={28} stroke="#3b82f6" strokeWidth={2} />
            <text x={20} y={32} fontSize={9} fill="var(--text-secondary)">Lower target</text>
          </g>
        </svg>
      </div>
      
      {/* Solver settings & execute */}
      <div style={{ padding: '8px 10px', borderTop: '1px solid var(--border-color)' }}>
        {method === 'qdes' ? (
          <>
            <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px', fontSize: '11px' }}>
              <label style={{ color: 'var(--text-muted)' }}>Iterations:</label>
              <input
                type="number"
                min={1} max={20}
                value={inverseDesign.maxIterations}
                onChange={(e) => setInverseDesignMaxIterations(parseInt(e.target.value) || 6)}
                style={{ width: '40px', fontSize: '11px', padding: '1px 4px' }}
              />
              <label style={{ color: 'var(--text-muted)' }}>Damping:</label>
              <input
                type="number"
                min={0.1} max={1.0} step={0.1}
                value={inverseDesign.damping}
                onChange={(e) => setInverseDesignDamping(parseFloat(e.target.value) || 0.6)}
                style={{ width: '50px', fontSize: '11px', padding: '1px 4px' }}
              />
            </div>
            
            <div style={{ display: 'flex', gap: '6px' }}>
              <button
                className="btn btn-primary btn-sm"
                onClick={handleExecute}
                disabled={inverseDesign.solving || (!inverseDesign.upperTarget && !inverseDesign.lowerTarget)}
                style={{ flex: 1 }}
              >
                {inverseDesign.solving ? 'Solving...' : 'Execute QDES'}
              </button>
              
              {inverseDesign.resultCoords && (
                <button
                  className="btn btn-sm"
                  onClick={applyInverseDesignResult}
                  style={{ flex: 1 }}
                  title="Replace current airfoil with the inverse design result"
                >
                  Apply Result
                </button>
              )}
            </div>
          </>
        ) : (
          <>
            <div style={{ display: 'flex', gap: '8px', alignItems: 'center', marginBottom: '6px', fontSize: '11px' }}>
              <label style={{ display: 'flex', gap: '4px', alignItems: 'center', cursor: 'pointer' }}>
                <input
                  type="checkbox"
                  checked={mdesSymmetric}
                  onChange={(e) => setMdesSymmetric(e.target.checked)}
                />
                <span style={{ color: 'var(--text-muted)' }}>Symmetric</span>
              </label>
              <label style={{ color: 'var(--text-muted)' }}>Filter:</label>
              <input
                type="number"
                min={0} max={5} step={0.5}
                value={mdesFilter}
                onChange={(e) => setMdesFilter(parseFloat(e.target.value) || 0)}
                style={{ width: '45px', fontSize: '11px', padding: '1px 4px' }}
              />
            </div>
            
            <div style={{ display: 'flex', gap: '6px' }}>
              <button
                className="btn btn-primary btn-sm"
                onClick={handleExecuteMdes}
                disabled={mdesSolving}
                style={{ flex: 1 }}
              >
                {mdesSolving ? 'Solving...' : 'Execute MDES'}
              </button>
              
              {inverseDesign.resultCoords && (
                <button
                  className="btn btn-sm"
                  onClick={applyInverseDesignResult}
                  style={{ flex: 1 }}
                  title="Replace current airfoil with the MDES result"
                >
                  Apply Result
                </button>
              )}
            </div>
            
            <p style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '4px' }}>
              Full-inverse: designs entire airfoil from circle-plane mapping.
              Modify Cp targets above, then execute.
            </p>
          </>
        )}
        
        {/* Convergence info */}
        {inverseDesign.converged !== null && (
          <div style={{ fontSize: '10px', marginTop: '4px', color: 'var(--text-muted)' }}>
            {inverseDesign.converged ? 'Converged' : 'Not converged'}
            {inverseDesign.history.length > 0 && (
              <span> -- {inverseDesign.history.length} iter, RMS={inverseDesign.history[inverseDesign.history.length - 1].rms_error.toFixed(4)}</span>
            )}
          </div>
        )}
      </div>
    </div>
  );
}
