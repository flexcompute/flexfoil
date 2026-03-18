/**
 * Multi-parameter sweep engine.
 *
 * Runs single or matrix sweeps over any combination of solver parameters
 * (alpha, Re, Mach, Ncrit) and geometry parameters (flap deflection,
 * flap hinge x/c).
 */

import type {
  AirfoilPoint,
  FlapDefinition,
  PolarPoint,
  SolverMode,
  SweepAxis,
  SweepParam,
} from '../types';
import { analyzeAirfoil, analyzeAirfoilInviscid, isWasmReady } from './wasm';
import { applyFlapsToBase } from './flapGeometry';
import { computeAirfoilHash } from './airfoilHash';

export interface SweepConfig {
  primary: SweepAxis;
  secondary: SweepAxis | null;
  alpha: number;
  reynolds: number;
  mach: number;
  ncrit: number;
  maxIterations: number;
  solverMode: SolverMode;
  baseCoordinates: AirfoilPoint[];
  panels: AirfoilPoint[];
  flaps: FlapDefinition[];
  nPanels: number;
  name: string;
}

export interface SweepSeriesResult {
  key: string;
  label: string;
  points: PolarPoint[];
}

/** Per-point result passed to onRun for database persistence */
export interface SweepRunData {
  airfoilHash: string;
  alpha: number;
  reynolds: number;
  mach: number;
  ncrit: number;
  nPanels: number;
  maxIterations: number;
  solverMode: SolverMode;
  cl: number;
  cd: number;
  cm: number;
  converged: boolean;
  iterations: number;
  residual: number;
  x_tr_upper: number;
  x_tr_lower: number;
  success: boolean;
  error: string | null;
  panels: AirfoilPoint[];
  flaps: FlapDefinition[];
  name: string;
}

export interface SweepCallbacks {
  onPoint: (series: SweepSeriesResult) => void;
  onRun: (data: SweepRunData) => Promise<void> | void;
  onProgress: (completed: number, total: number) => Promise<void> | void;
  /** Fires after all primary-axis points for one secondary value are done. */
  onSeriesComplete?: () => Promise<void> | void;
  signal: AbortSignal;
}

const GEOMETRY_PARAMS: SweepParam[] = ['flapDeflection', 'flapHingeX'];

function isGeometryParam(p: SweepParam): boolean {
  return GEOMETRY_PARAMS.includes(p);
}

function generateValues(axis: SweepAxis): number[] {
  const values: number[] = [];
  const step = Math.abs(axis.step) || 1;
  const direction = axis.end >= axis.start ? 1 : -1;
  for (let v = axis.start; direction > 0 ? v <= axis.end + 1e-9 : v >= axis.end - 1e-9; v += step * direction) {
    values.push(Math.round(v * 1e8) / 1e8);
  }
  return values;
}

function paramLabel(param: SweepParam): string {
  switch (param) {
    case 'alpha': return 'α';
    case 'reynolds': return 'Re';
    case 'mach': return 'M';
    case 'ncrit': return 'Ncrit';
    case 'flapDeflection': return 'δ_flap';
    case 'flapHingeX': return 'x_hinge';
  }
}

function formatValue(param: SweepParam, value: number): string {
  switch (param) {
    case 'alpha': return `${value.toFixed(1)}°`;
    case 'reynolds': return value.toExponential(1);
    case 'mach': return value.toFixed(2);
    case 'ncrit': return value.toFixed(1);
    case 'flapDeflection': return `${value.toFixed(1)}°`;
    case 'flapHingeX': return value.toFixed(3);
  }
}

function applyParamToFlaps(
  flaps: FlapDefinition[],
  param: SweepParam,
  flapId: string | undefined,
  value: number,
): FlapDefinition[] {
  if (param === 'flapDeflection') {
    return flaps.map(f => f.id === flapId ? { ...f, deflection: value } : f);
  }
  if (param === 'flapHingeX') {
    return flaps.map(f => f.id === flapId ? { ...f, hingeX: value } : f);
  }
  return flaps;
}

/**
 * Run a single- or matrix-parameter sweep.
 *
 * For a single sweep, the primary axis is iterated and one polar series is produced.
 * For a matrix sweep, the secondary axis is the outer loop (one series per value)
 * and the primary axis is the inner loop.
 */
export async function runSweep(
  config: SweepConfig,
  callbacks: SweepCallbacks,
): Promise<SweepSeriesResult[]> {
  if (!isWasmReady()) throw new Error('WASM not initialized');

  const { primary, secondary } = config;
  const isViscous = config.solverMode === 'viscous';

  const primaryValues = generateValues(primary);
  const secondaryValues = secondary ? generateValues(secondary) : [null];
  const totalPoints = primaryValues.length * secondaryValues.length;

  const allSeries: SweepSeriesResult[] = [];
  let completed = 0;

  for (let secIdx = 0; secIdx < secondaryValues.length; secIdx++) {
    const secVal = secondaryValues[secIdx];
    if (callbacks.signal.aborted) break;

    const seriesPoints: PolarPoint[] = [];

    // Build a stable series key per secondary value (not per point)
    const seriesKey = secondary && secVal != null
      ? `sweep_${config.name}_${secondary.param}_${secIdx}`
      : `sweep_${config.name}_single`;

    let seriesLabel = config.name;
    if (isViscous) {
      const reForLabel = secondary?.param === 'reynolds' && secVal != null ? secVal : config.reynolds;
      const ncritForLabel = secondary?.param === 'ncrit' && secVal != null ? secVal : config.ncrit;
      seriesLabel += ` Re=${reForLabel.toExponential(1)} Nc=${ncritForLabel}`;
    } else {
      seriesLabel += ' (inviscid)';
    }

    if (secVal != null && secondary) {
      seriesLabel += ` ${paramLabel(secondary.param)}=${formatValue(secondary.param, secVal)}`;
    }

    // Pre-apply secondary geometry param if it's a geometry param
    let secFlaps = config.flaps;
    let secPanels = config.panels;
    let secHash: string | null = null;

    if (secVal != null && secondary && isGeometryParam(secondary.param)) {
      secFlaps = applyParamToFlaps(config.flaps, secondary.param, secondary.flapId, secVal);
      const geom = applyFlapsToBase(config.baseCoordinates, secFlaps, config.nPanels);
      if (!geom) continue;
      secPanels = geom.panels;
      secHash = await computeAirfoilHash(secPanels);
    }

    for (const priVal of primaryValues) {
      if (callbacks.signal.aborted) break;

      let alpha = config.alpha;
      let reynolds = config.reynolds;
      let mach = config.mach;
      let ncrit = config.ncrit;
      let panels = secPanels;
      let airfoilHash = secHash;
      let currentFlaps = secFlaps;

      // Apply primary value
      switch (primary.param) {
        case 'alpha': alpha = priVal; break;
        case 'reynolds': reynolds = priVal; break;
        case 'mach': mach = priVal; break;
        case 'ncrit': ncrit = priVal; break;
        case 'flapDeflection':
        case 'flapHingeX': {
          currentFlaps = applyParamToFlaps(currentFlaps, primary.param, primary.flapId, priVal);
          const geom = applyFlapsToBase(config.baseCoordinates, currentFlaps, config.nPanels);
          if (!geom) { completed++; continue; }
          panels = geom.panels;
          airfoilHash = await computeAirfoilHash(panels);
          break;
        }
      }

      // Apply secondary solver param
      if (secVal != null && secondary && !isGeometryParam(secondary.param)) {
        switch (secondary.param) {
          case 'alpha': alpha = secVal; break;
          case 'reynolds': reynolds = secVal; break;
          case 'mach': mach = secVal; break;
          case 'ncrit': ncrit = secVal; break;
        }
      }

      if (!airfoilHash) {
        airfoilHash = await computeAirfoilHash(panels);
      }

      // Solve
      try {
        const res = isViscous
          ? analyzeAirfoil(panels, alpha, reynolds, mach, ncrit, config.maxIterations)
          : analyzeAirfoilInviscid(panels, alpha);

        if (res.success) {
          const point: PolarPoint = {
            alpha,
            cl: res.cl,
            cd: res.cd,
            cm: res.cm,
          };
          if (primary.param === 'reynolds' || secondary?.param === 'reynolds') point.reynolds = reynolds;
          if (primary.param === 'mach' || secondary?.param === 'mach') point.mach = mach;
          if (primary.param === 'ncrit' || secondary?.param === 'ncrit') point.ncrit = ncrit;
          if (primary.param === 'flapDeflection' || secondary?.param === 'flapDeflection') {
            const fid = primary.param === 'flapDeflection' ? primary.flapId : secondary?.flapId;
            point.flapDeflection = currentFlaps.find(f => f.id === fid)?.deflection ?? priVal;
          }
          if (primary.param === 'flapHingeX' || secondary?.param === 'flapHingeX') {
            const fid = primary.param === 'flapHingeX' ? primary.flapId : secondary?.flapId;
            point.flapHingeX = currentFlaps.find(f => f.id === fid)?.hingeX ?? priVal;
          }

          seriesPoints.push(point);
        }

        // Queue run data for database persistence (caller handles actual insert)
        callbacks.onRun({
          airfoilHash,
          alpha,
          reynolds: isViscous ? reynolds : 0,
          mach: isViscous ? mach : 0,
          ncrit: isViscous ? ncrit : 0,
          nPanels: config.nPanels,
          maxIterations: isViscous ? config.maxIterations : 0,
          solverMode: config.solverMode,
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
          panels,
          flaps: currentFlaps,
          name: config.name,
        });
      } catch (err) {
        console.warn('[sweep] solver error at point', completed + 1, '/', totalPoints, err);
      }

      completed++;
      await callbacks.onProgress(completed, totalPoints);

      callbacks.onPoint({ key: seriesKey, label: seriesLabel, points: [...seriesPoints] });

      await new Promise(r => setTimeout(r, 0));
    }

    if (seriesPoints.length > 0) {
      allSeries.push({ key: seriesKey, label: seriesLabel, points: seriesPoints });
    }

    await callbacks.onSeriesComplete?.();
  }

  return allSeries;
}
