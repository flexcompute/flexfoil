/**
 * SolvePanel - Aerodynamic analysis controls with solver caching.
 *
 * Every solver evaluation is persisted to the run database (SQLite).
 * Before calling the WASM solver, the cache is checked — cache hits
 * skip the solver entirely and return the stored result.
 */

import { useState, useCallback, useMemo, useEffect, useRef } from 'react';
import { useAirfoilStore } from '../../stores/airfoilStore';
import { useRouteUiStore } from '../../stores/routeUiStore';
import { useRunStore } from '../../stores/runStore';
import { useCaseLogStore } from '../../stores/caseLogStore';
import { analyzeAirfoil, analyzeAirfoilInviscid, isWasmReady, type AnalysisResult } from '../../lib/wasm';
import type { PolarPoint } from '../../types';
import type { RunInsert } from '../../lib/runDatabase';

type RunMode = 'alpha' | 'cl';

type SolveOrCacheResult = {
  result: AnalysisResult | null;
  fromCache: boolean;
};

function roundLogNumber(value: number | null | undefined, digits = 4): number | null {
  if (value == null || !Number.isFinite(value)) return null;
  return Number(value.toFixed(digits));
}

function buildResultDetails(result: AnalysisResult): Record<string, unknown> {
  return {
    success: result.success,
    converged: result.converged,
    cl: roundLogNumber(result.cl),
    cd: roundLogNumber(result.cd, 6),
    cm: roundLogNumber(result.cm),
    iterations: result.iterations,
    residual: roundLogNumber(result.residual, 6),
    x_tr_upper: roundLogNumber(result.x_tr_upper),
    x_tr_lower: roundLogNumber(result.x_tr_lower),
    error: result.error ?? null,
  };
}
export function SolvePanel() {
  const {
    panels,
    name,
    polarData,
    displayAlpha,
    reynolds,
    mach,
    ncrit,
    maxIterations,
    solverMode,
    setDisplayAlpha,
    setReynolds,
    setMach,
    setNcrit,
    setMaxIterations,
    setSolverMode,
    upsertPolar,
    clearAllPolars,
  } = useAirfoilStore();

  const { lookup, addRun, hashPanels } = useRunStore();
  const startCase = useCaseLogStore((state) => state.startCase);
  const appendEvent = useCaseLogStore((state) => state.appendEvent);
  const finishCase = useCaseLogStore((state) => state.finishCase);

  const runMode = useRouteUiStore((state) => state.solveRunMode);
  const setRunMode = useRouteUiStore((state) => state.setSolveRunMode);
  const showAdvanced = useRouteUiStore((state) => state.solveShowAdvanced);
  const setShowAdvanced = useRouteUiStore((state) => state.setSolveShowAdvanced);
  const targetCl = useRouteUiStore((state) => state.solveTargetCl);
  const setTargetCl = useRouteUiStore((state) => state.setSolveTargetCl);
  const alphaStart = useRouteUiStore((state) => state.solvePolarStart);
  const setAlphaStart = useRouteUiStore((state) => state.setSolvePolarStart);
  const alphaEnd = useRouteUiStore((state) => state.solvePolarEnd);
  const setAlphaEnd = useRouteUiStore((state) => state.setSolvePolarEnd);
  const alphaStep = useRouteUiStore((state) => state.solvePolarStep);
  const setAlphaStep = useRouteUiStore((state) => state.setSolvePolarStep);

  // Single-point inputs
  const [targetAlpha, setTargetAlphaLocal] = useState(displayAlpha);

  const reynoldsDebounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const [reynoldsText, setReynoldsText] = useState(() => String(reynolds));
  const reynoldsFocusedRef = useRef(false);

  useEffect(() => {
    if (!reynoldsFocusedRef.current) {
      setReynoldsText(String(reynolds));
    }
  }, [reynolds]);

  // Sync local targetAlpha when store changes externally (e.g. URL load)
  useEffect(() => {
    setTargetAlphaLocal(displayAlpha);
  }, [displayAlpha]);

  // Commit alpha: updates both local state and the store (called on Enter/blur/arrows)
  const setTargetAlpha = useCallback((alpha: number) => {
    setTargetAlphaLocal(alpha);
    setDisplayAlpha(alpha);
  }, [setDisplayAlpha]);

  useEffect(() => () => {
    if (reynoldsDebounceRef.current) clearTimeout(reynoldsDebounceRef.current);
  }, []);

  const handleReynoldsChange = useCallback((raw: string) => {
    setReynoldsText(raw);
    if (reynoldsDebounceRef.current) clearTimeout(reynoldsDebounceRef.current);
    reynoldsDebounceRef.current = setTimeout(() => {
      const parsed = Number(raw);
      if (!isNaN(parsed) && parsed >= 1 && isFinite(parsed)) {
        setReynolds(parsed);
      }
    }, 800);
  }, [setReynolds]);

  const commitReynolds = useCallback(() => {
    reynoldsFocusedRef.current = false;
    if (reynoldsDebounceRef.current) clearTimeout(reynoldsDebounceRef.current);
    const parsed = Number(reynoldsText);
    if (!isNaN(parsed) && parsed >= 1 && isFinite(parsed)) {
      setReynolds(parsed);
    } else {
      setReynoldsText(String(reynolds));
    }
  }, [reynoldsText, reynolds, setReynolds]);
  // Results
  const [result, setResult] = useState<AnalysisResult | null>(null);
  const [resultAlpha, setResultAlpha] = useState<number | null>(null);
  const [isRunning, setIsRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [cacheStats, setCacheStats] = useState<{ hits: number; misses: number } | null>(null);

  const lastSeries = useMemo(
    () => (polarData.length > 0 ? polarData[polarData.length - 1] : null),
    [polarData],
  );
  const polar = useMemo(() => lastSeries?.points ?? [], [lastSeries]);

  const isViscous = solverMode === 'viscous';

  // --------------- helpers ---------------

  const runSolver = useCallback((coords: { x: number; y: number }[], alpha: number) => {
    if (isViscous) {
      return analyzeAirfoil(coords, alpha, reynolds, mach, ncrit, maxIterations);
    }
    return analyzeAirfoilInviscid(coords, alpha);
  }, [isViscous, reynolds, mach, ncrit, maxIterations]);

  /** Cache key uses Re=0/Mach=0/Ncrit=0/maxIter=0 for inviscid to separate caches. */
  const cacheRe = isViscous ? reynolds : 0;
  const cacheMach = isViscous ? mach : 0;
  const cacheNcrit = isViscous ? ncrit : 0;
  const cacheMaxIter = isViscous ? maxIterations : 0;

  const buildCaseMetadata = useCallback((overrides: Record<string, unknown> = {}) => ({
    airfoil: name,
    solver: isViscous ? 'viscous' : 'inviscid',
    reynolds: isViscous ? reynolds : 'n/a',
    mach: isViscous ? mach : 'n/a',
    ncrit: isViscous ? ncrit : 'n/a',
    max_iter: isViscous ? maxIterations : 'n/a',
    panels: panels.length - 1,
    ...overrides,
  }), [name, isViscous, reynolds, mach, ncrit, maxIterations, panels.length]);

  /** Run a single point — check cache first, fall back to solver. */
  const solveOrCache = useCallback(async (
    alpha: number,
    airfoilHash: string,
    nPanels: number,
    logContext?: { caseId: string; label: string },
  ): Promise<SolveOrCacheResult> => {
    const cached = lookup(airfoilHash, alpha, cacheRe, cacheMach, cacheNcrit, nPanels, cacheMaxIter);
    if (cached && cached.success) {
      const result = {
        cl: cached.cl ?? 0,
        cd: cached.cd ?? 0,
        cm: cached.cm ?? 0,
        cp: [],
        cp_x: [],
        gamma: [],
        psi_0: 0,
        converged: cached.converged,
        iterations: cached.iterations ?? 0,
        residual: cached.residual ?? 0,
        x_tr_upper: cached.x_tr_upper ?? 0,
        x_tr_lower: cached.x_tr_lower ?? 0,
        success: true,
      } as AnalysisResult;

      if (logContext) {
        appendEvent(logContext.caseId, {
          level: 'success',
          message: `${logContext.label}: cache hit`,
          details: {
            alpha: roundLogNumber(alpha),
            run_id: cached.id,
            ...buildResultDetails(result),
          },
        });
      }

      return { result, fromCache: true };
    }

    if (logContext) {
      appendEvent(logContext.caseId, {
        message: `${logContext.label}: cache miss`,
        details: {
          alpha: roundLogNumber(alpha),
          airfoil_hash: airfoilHash.slice(0, 12),
        },
      });
    }

    const res = runSolver(panels, alpha);

    if (logContext) {
      appendEvent(logContext.caseId, {
        level: res.success ? (res.converged ? 'success' : 'warning') : 'error',
        message: `${logContext.label}: solver completed`,
        details: {
          alpha: roundLogNumber(alpha),
          source: 'solver',
          ...buildResultDetails(res),
        },
      });
    }

    const run: RunInsert = {
      airfoil_name: name,
      airfoil_hash: airfoilHash,
      alpha,
      reynolds: cacheRe,
      mach: cacheMach,
      ncrit: cacheNcrit,
      n_panels: nPanels,
      max_iter: cacheMaxIter,
      cl: res.cl,
      cd: res.cd,
      cm: res.cm,
      converged: res.converged,
      iterations: res.iterations,
      residual: res.residual,
      x_tr_upper: res.x_tr_upper,
      x_tr_lower: res.x_tr_lower,
      success: res.success,
      error: res.error ?? null,
    };
    await addRun(run);

    if (logContext) {
      appendEvent(logContext.caseId, {
        message: `${logContext.label}: run stored`,
        details: {
          alpha: roundLogNumber(alpha),
          success: res.success,
          converged: res.converged,
        },
      });
    }

    return { result: res, fromCache: false };
  }, [panels, name, cacheRe, cacheMach, cacheNcrit, cacheMaxIter, runSolver, lookup, addRun, appendEvent]);

  // --------------- single-point ---------------

  const runAnalysis = useCallback(async () => {
    if (!isWasmReady() || panels.length < 3) {
      setError('WASM not ready or insufficient geometry');
      return;
    }

    setIsRunning(true);
    setError(null);
    setResultAlpha(null);
    setCacheStats(null);

    let caseId: string | null = null;

    try {
      caseId = startCase({
        kind: runMode === 'alpha' ? 'single-alpha' : 'single-cl',
        title: runMode === 'alpha'
          ? `Run to α ${targetAlpha.toFixed(2)}°`
          : `Run to CL ${targetCl.toFixed(3)}`,
        metadata: buildCaseMetadata(
          runMode === 'alpha'
            ? { target_alpha: roundLogNumber(targetAlpha) }
            : { target_cl: roundLogNumber(targetCl) },
        ),
      });

      appendEvent(caseId, {
        message: 'Starting analysis',
        details: {
          mode: runMode,
          solver: isViscous ? 'viscous' : 'inviscid',
        },
      });

      const airfoilHash = await hashPanels(panels);
      const nPanels = panels.length - 1;

      appendEvent(caseId, {
        message: 'Geometry hashed for cache lookup',
        details: {
          airfoil_hash: airfoilHash.slice(0, 12),
          n_panels: nPanels,
        },
      });

      if (runMode === 'alpha') {
        const { result: res, fromCache } = await solveOrCache(targetAlpha, airfoilHash, nPanels, {
          caseId,
          label: `Single-point α=${targetAlpha.toFixed(2)}°`,
        });
        if (!res) {
          setError('Analysis returned no result');
          finishCase(caseId, {
            status: 'error',
            summary: 'Analysis returned no result.',
          });
          return;
        }
        setResult(res);
        setResultAlpha(targetAlpha);
        if (!res.success) {
          setError(res.error || 'Analysis failed');
          finishCase(caseId, {
            status: 'error',
            summary: res.error || 'Single-point analysis failed.',
            metadata: { source: fromCache ? 'cache' : 'solver' },
          });
        } else if (isViscous && !res.converged) {
          setError(`Viscous solver did not converge (residual ${res.residual.toExponential(2)})`);
          finishCase(caseId, {
            status: 'warning',
            summary: `Single-point result returned from ${fromCache ? 'cache' : 'solver'}, but the viscous solve did not converge.`,
            metadata: { source: fromCache ? 'cache' : 'solver' },
          });
        } else {
          finishCase(caseId, {
            status: 'success',
            summary: `Single-point analysis completed from ${fromCache ? 'cache' : 'solver'}.`,
            metadata: { source: fromCache ? 'cache' : 'solver' },
          });
        }
      } else {
        let alpha = 0;
        const clMaxIter = 20;
        const tol = 0.001;

        appendEvent(caseId, {
          message: 'Starting CL targeting loop',
          details: {
            initial_alpha: 0,
            tolerance: tol,
            max_iterations: clMaxIter,
          },
        });

        for (let i = 0; i < clMaxIter; i++) {
          const iterationAlpha = alpha;
          const { result: res, fromCache } = await solveOrCache(iterationAlpha, airfoilHash, nPanels, {
            caseId,
            label: `CL iteration ${i + 1} at α=${iterationAlpha.toFixed(2)}°`,
          });
          if (!res || !res.success) {
            setError(res?.error || 'Analysis failed during CL iteration');
            finishCase(caseId, {
              status: 'error',
              summary: res?.error || 'CL targeting failed during iteration.',
              metadata: {
                failed_iteration: i + 1,
                source: fromCache ? 'cache' : 'solver',
              },
            });
            return;
          }

          const clError = targetCl - res.cl;
          appendEvent(caseId, {
            level: Math.abs(clError) < tol ? 'success' : 'info',
            message: `CL iteration ${i + 1} residual evaluated`,
            details: {
              alpha: roundLogNumber(iterationAlpha),
              achieved_cl: roundLogNumber(res.cl),
              target_cl: roundLogNumber(targetCl),
              cl_error: roundLogNumber(clError, 6),
              source: fromCache ? 'cache' : 'solver',
            },
          });

          if (Math.abs(clError) < tol) {
            setResult(res);
            setResultAlpha(iterationAlpha);
            setDisplayAlpha(iterationAlpha);
            finishCase(caseId, {
              status: 'success',
              summary: `CL targeting converged in ${i + 1} iteration${i === 0 ? '' : 's'}.`,
              metadata: {
                converged_alpha: roundLogNumber(iterationAlpha),
                achieved_cl: roundLogNumber(res.cl),
              },
            });
            break;
          }

          alpha += clError / 0.11;
          alpha = Math.max(-20, Math.min(25, alpha));

          appendEvent(caseId, {
            message: `Adjusting α for iteration ${i + 2}`,
            details: {
              previous_alpha: roundLogNumber(iterationAlpha),
              next_alpha: roundLogNumber(alpha),
            },
          });

          if (i === clMaxIter - 1) {
            const finalRes = await solveOrCache(alpha, airfoilHash, nPanels, {
              caseId,
              label: 'Final CL targeting evaluation',
            });
            setResult(finalRes.result);
            setResultAlpha(alpha);
            setDisplayAlpha(alpha);
            if (finalRes.result) {
              appendEvent(caseId, {
                level: Math.abs(targetCl - finalRes.result.cl) > 0.01 ? 'warning' : 'success',
                message: 'Final CL targeting evaluation completed',
                details: {
                  alpha: roundLogNumber(alpha),
                  achieved_cl: roundLogNumber(finalRes.result.cl),
                  source: finalRes.fromCache ? 'cache' : 'solver',
                },
              });
            }

            if (finalRes.result && Math.abs(targetCl - finalRes.result.cl) > 0.01) {
              setError(`Could not converge to CL=${targetCl.toFixed(3)}. Got CL=${finalRes.result.cl.toFixed(3)} at α=${alpha.toFixed(2)}°`);
              finishCase(caseId, {
                status: 'warning',
                summary: `CL targeting stopped at α=${alpha.toFixed(2)}° without reaching the requested tolerance.`,
                metadata: {
                  final_cl: roundLogNumber(finalRes.result.cl),
                  final_alpha: roundLogNumber(alpha),
                },
              });
            } else {
              finishCase(caseId, {
                status: 'success',
                summary: `CL targeting completed after ${clMaxIter} iterations.`,
                metadata: {
                  final_alpha: roundLogNumber(alpha),
                  final_cl: roundLogNumber(finalRes.result?.cl ?? null),
                },
              });
            }
          }
        }
      }
    } catch (e) {
      const message = e instanceof Error ? e.message : 'Unknown error';
      setError(message);
      if (caseId) {
        appendEvent(caseId, {
          level: 'error',
          message: 'Analysis threw an exception',
          details: { error: message },
        });
        finishCase(caseId, {
          status: 'error',
          summary: message,
        });
      }
    } finally {
      setIsRunning(false);
    }
  }, [
    panels,
    runMode,
    targetAlpha,
    targetCl,
    isViscous,
    setDisplayAlpha,
    hashPanels,
    solveOrCache,
    startCase,
    appendEvent,
    finishCase,
    buildCaseMetadata,
  ]);

  // --------------- polar sweep ---------------

  const runPolar = useCallback(async () => {
    if (!isWasmReady() || panels.length < 3) {
      setError('WASM not ready or insufficient geometry');
      return;
    }

    setIsRunning(true);
    setError(null);
    setCacheStats(null);

    let caseId: string | null = null;

    try {
      caseId = startCase({
        kind: 'polar',
        title: `Alpha polar ${alphaStart.toFixed(1)}° to ${alphaEnd.toFixed(1)}°`,
        metadata: buildCaseMetadata({
          alpha_start: roundLogNumber(alphaStart),
          alpha_end: roundLogNumber(alphaEnd),
          alpha_step: roundLogNumber(alphaStep),
        }),
      });

      appendEvent(caseId, {
        message: 'Starting polar sweep',
        details: {
          alpha_start: roundLogNumber(alphaStart),
          alpha_end: roundLogNumber(alphaEnd),
          alpha_step: roundLogNumber(alphaStep),
        },
      });

      const airfoilHash = await hashPanels(panels);
      const nPanels = panels.length - 1;

      appendEvent(caseId, {
        message: 'Geometry hashed for polar sweep',
        details: {
          airfoil_hash: airfoilHash.slice(0, 12),
          n_panels: nPanels,
        },
      });

      const polarKey = `${airfoilHash.slice(0, 12)}_${cacheRe}_${cacheMach}_${cacheNcrit}_${nPanels}_${cacheMaxIter}`;
      const polarLabel = isViscous
        ? `${name} Re=${cacheRe.toExponential(1)} Nc=${cacheNcrit}`
        : `${name} (inviscid)`;

      const points: PolarPoint[] = [];
      let hits = 0;
      let misses = 0;
      let failures = 0;
      const totalPoints = Math.max(0, Math.floor(((alphaEnd - alphaStart) / alphaStep) + 1 + 1e-9));

      for (let alpha = alphaStart, pointIndex = 0; alpha <= alphaEnd + 1e-9; alpha += alphaStep, pointIndex++) {
        const roundedAlpha = Math.round(alpha * 1e6) / 1e6;

        const { result: res, fromCache } = await solveOrCache(roundedAlpha, airfoilHash, nPanels, {
          caseId,
          label: `Polar point ${pointIndex + 1}/${totalPoints || '?'}`,
        });

        if (res && fromCache) {
          points.push({ alpha: roundedAlpha, cl: res.cl ?? 0, cd: res.cd ?? 0, cm: res.cm ?? 0 });
          hits++;
          continue;
        }

        if (res) {
          misses++;
          if (res.success) {
            points.push({ alpha: roundedAlpha, cl: res.cl, cd: res.cd, cm: res.cm });
          } else {
            failures++;
            appendEvent(caseId, {
              level: 'warning',
              message: `Polar point ${pointIndex + 1}/${totalPoints || '?'} skipped`,
              details: {
                alpha: roundLogNumber(roundedAlpha),
                error: res.error ?? 'Solver returned an unsuccessful result.',
              },
            });
          }
        }
      }

      upsertPolar({ key: polarKey, label: polarLabel, points });
      setCacheStats({ hits, misses });

      if (points.length === 0) {
        setError('No valid points in polar');
        finishCase(caseId, {
          status: 'error',
          summary: 'Polar sweep finished without any valid points.',
          metadata: { cache_hits: hits, computed_points: misses, failed_points: failures },
        });
      } else if (failures > 0) {
        finishCase(caseId, {
          status: 'warning',
          summary: `Polar sweep completed with ${points.length} valid points, ${failures} failed point${failures === 1 ? '' : 's'}, ${hits} cache hit${hits === 1 ? '' : 's'}, and ${misses} computed point${misses === 1 ? '' : 's'}.`,
          metadata: { cache_hits: hits, computed_points: misses, failed_points: failures },
        });
      } else {
        finishCase(caseId, {
          status: 'success',
          summary: `Polar sweep completed with ${points.length} valid points, ${hits} cache hit${hits === 1 ? '' : 's'}, and ${misses} computed point${misses === 1 ? '' : 's'}.`,
          metadata: { cache_hits: hits, computed_points: misses, failed_points: failures },
        });
      }
    } catch (e) {
      const message = e instanceof Error ? e.message : 'Unknown error';
      setError(message);
      if (caseId) {
        appendEvent(caseId, {
          level: 'error',
          message: 'Polar sweep threw an exception',
          details: { error: message },
        });
        finishCase(caseId, {
          status: 'error',
          summary: message,
        });
      }
    } finally {
      setIsRunning(false);
    }
  }, [panels, name, alphaStart, alphaEnd, alphaStep, cacheRe, cacheMach, cacheNcrit, cacheMaxIter,
      isViscous, upsertPolar, hashPanels, solveOrCache, startCase, appendEvent, finishCase, buildCaseMetadata]);

  // --------------- derived ---------------

  const clAlpha = useMemo(() => {
    if (polar.length < 2) return null;
    const linearPoints = polar.filter(p => p.alpha >= -5 && p.alpha <= 10);
    if (linearPoints.length < 2) return null;
    const n = linearPoints.length;
    const sumX = linearPoints.reduce((s, p) => s + p.alpha, 0);
    const sumY = linearPoints.reduce((s, p) => s + p.cl, 0);
    const sumXY = linearPoints.reduce((s, p) => s + p.alpha * p.cl, 0);
    const sumX2 = linearPoints.reduce((s, p) => s + p.alpha * p.alpha, 0);
    return (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
  }, [polar]);

  // --------------- render ---------------

  return (
    <div className="panel">
      <div className="panel-header">Solve</div>
      <div className="panel-content">
        <div className="form-group" data-tour="solve-mode">
          <div className="form-label">Solver</div>
          <div style={{ display: 'flex', gap: '4px', marginBottom: '4px' }}>
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
          <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '2px' }}>
            {isViscous
              ? 'XFOIL viscous solver with boundary layer coupling'
              : 'Linear-vorticity panel method (CD = 0)'}
          </div>
        </div>

        {isViscous && (
          <div className="form-group">
            <div className="form-label">Reynolds Number</div>
            <input
              type="text"
              inputMode="numeric"
              value={reynoldsText}
              onChange={(e) => handleReynoldsChange(e.target.value)}
              onFocus={() => { reynoldsFocusedRef.current = true; }}
              onBlur={commitReynolds}
              onKeyDown={(e) => { if (e.key === 'Enter') commitReynolds(); }}
              placeholder="e.g. 6e6, 1000000"
            />
          </div>
        )}

        {isViscous && (
          <div className="form-group">
            <button
              onClick={() => setShowAdvanced(!showAdvanced)}
              style={{
                width: '100%',
                padding: '4px 8px',
                fontSize: '11px',
                background: 'transparent',
                border: '1px solid var(--border-color)',
                borderRadius: '4px',
                color: 'var(--text-secondary)',
                cursor: 'pointer',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
              }}
            >
              <span>Advanced Settings</span>
              <span style={{ fontSize: '10px' }}>{showAdvanced ? '▲' : '▼'}</span>
            </button>
            {showAdvanced && (
              <div style={{ display: 'flex', flexDirection: 'column', gap: '8px', marginTop: '8px' }}>
                <div>
                  <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Mach Number</label>
                  <input
                    type="number"
                    value={mach}
                    onChange={(e) => setMach(Math.max(0, Math.min(0.8, parseFloat(e.target.value) || 0)))}
                    step={0.01}
                    min={0}
                    max={0.8}
                  />
                </div>
                <div>
                  <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>
                    Ncrit (e<sup>N</sup> transition)
                  </label>
                  <input
                    type="number"
                    value={ncrit}
                    onChange={(e) => setNcrit(Math.max(1, Math.min(14, parseFloat(e.target.value) || 9)))}
                    step={0.5}
                    min={1}
                    max={14}
                  />
                </div>
                <div>
                  <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Max Iterations</label>
                  <input
                    type="number"
                    value={maxIterations}
                    onChange={(e) => setMaxIterations(Math.max(10, Math.min(500, parseInt(e.target.value) || 100)))}
                    step={10}
                    min={10}
                    max={500}
                  />
                </div>
              </div>
            )}
          </div>
        )}

        {/* Single Point Analysis */}
        <div className="form-group" data-tour="solve-alpha">
          <div className="form-label">Single Point</div>

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

          <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
            <label style={{ fontSize: '12px', minWidth: '60px' }}>
              {runMode === 'alpha' ? 'Alpha (°):' : 'Target CL:'}
            </label>
            <input
              type="number"
              value={runMode === 'alpha' ? targetAlpha : targetCl}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (runMode === 'alpha') setTargetAlpha(val);
                else setTargetCl(val);
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
              gridTemplateColumns: `repeat(${(runMode === 'cl' ? 1 : 0) + (isViscous ? 3 : 2)}, 1fr)`,
              gap: '8px',
            }}>
              {runMode === 'cl' && resultAlpha !== null && (
                <ResultCard label="α (°)" value={resultAlpha.toFixed(2)} highlight />
              )}
              <ResultCard label="CL" value={result.cl.toFixed(4)} />
              {isViscous && <ResultCard label="CD" value={result.cd.toFixed(5)} />}
              <ResultCard label="CM" value={result.cm.toFixed(4)} />
            </div>
            <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '8px' }}>
              {isViscous
                ? (result.converged
                    ? `Converged in ${result.iterations} iterations`
                    : `Not converged after ${result.iterations} iterations (residual ${result.residual.toExponential(2)})`)
                : 'Direct panel solve'}
            </div>
          </div>
        )}

        {/* Polar Sweep */}
        <div className="form-group" style={{ marginTop: '16px', paddingTop: '16px', borderTop: '1px solid var(--border-color)' }} data-tour="solve-polar">
          <div className="form-label">Alpha Polar</div>

          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '8px', marginBottom: '8px' }}>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Start (°)</label>
              <input type="number" value={alphaStart} onChange={(e) => setAlphaStart(parseFloat(e.target.value))} step={1} />
            </div>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>End (°)</label>
              <input type="number" value={alphaEnd} onChange={(e) => setAlphaEnd(parseFloat(e.target.value))} step={1} />
            </div>
            <div>
              <label style={{ fontSize: '10px', color: 'var(--text-muted)' }}>Step (°)</label>
              <input type="number" value={alphaStep} onChange={(e) => setAlphaStep(parseFloat(e.target.value))} step={0.5} min={0.1} />
            </div>
          </div>

          <div style={{ display: 'flex', gap: '4px' }}>
            <button
              onClick={runPolar}
              disabled={isRunning || !isWasmReady()}
              style={{ flex: 1 }}
            >
              {isRunning ? 'Running...' : 'Generate Polar'}
            </button>
            {polarData.length > 0 && (
              <button
                onClick={clearAllPolars}
                disabled={isRunning}
                title="Clear all polar series"
                style={{ padding: '4px 8px', fontSize: '11px' }}
              >
                Clear
              </button>
            )}
          </div>
          {polarData.length > 1 && (
            <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '4px' }}>
              {polarData.length} polar series overlaid — see Polar panel
            </div>
          )}
        </div>

        {/* Cache stats */}
        {cacheStats && (
          <div style={{ fontSize: '10px', color: 'var(--text-muted)', padding: '4px 0' }}>
            Cache: {cacheStats.hits} hits, {cacheStats.misses} computed
          </div>
        )}

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

            <div style={{
              maxHeight: '200px',
              overflow: 'auto',
              border: '1px solid var(--border-color)',
              borderRadius: '4px',
              fontFamily: 'var(--font-mono)',
              fontSize: '10px',
            }}>
              <div style={{
                display: 'grid',
                gridTemplateColumns: isViscous ? '1fr 1fr 1fr 1fr' : '1fr 1fr 1fr',
                gap: '4px',
                padding: '4px 8px',
                background: 'var(--bg-tertiary)',
                borderBottom: '1px solid var(--border-color)',
                position: 'sticky',
                top: 0,
              }}>
                <span>α (°)</span>
                <span>CL</span>
                {isViscous && <span>CD</span>}
                <span>CM</span>
              </div>

              {polar.map((p, i) => (
                <div key={i} style={{
                  display: 'grid',
                  gridTemplateColumns: isViscous ? '1fr 1fr 1fr 1fr' : '1fr 1fr 1fr',
                  gap: '4px',
                  padding: '2px 8px',
                  borderBottom: '1px solid var(--border-color)',
                }}>
                  <span>{p.alpha.toFixed(1)}</span>
                  <span style={{ color: p.cl >= 0 ? 'var(--accent-primary)' : 'var(--accent-secondary)' }}>
                    {p.cl.toFixed(4)}
                  </span>
                  {isViscous && <span>{p.cd !== undefined ? p.cd.toFixed(5) : '—'}</span>}
                  <span>{p.cm.toFixed(4)}</span>
                </div>
              ))}
            </div>

            <button
              onClick={() => {
                const header = isViscous
                  ? `# ${lastSeries?.label ?? name} - Viscous Polar (Re=${reynolds})\n# Alpha(deg)  CL        CD         CM\n`
                  : `# ${lastSeries?.label ?? name} - Inviscid Polar\n# Alpha(deg)  CL        CM\n`;
                const data = polar.map(p =>
                  isViscous
                    ? `${p.alpha.toFixed(2).padStart(8)} ${p.cl.toFixed(6).padStart(10)} ${(p.cd ?? 0).toFixed(6).padStart(10)} ${p.cm.toFixed(6).padStart(10)}`
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
          <div>Solver: {isViscous ? 'XFOIL viscous' : 'Inviscid panel method'}</div>
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
