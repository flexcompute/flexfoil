import { compressToEncodedURIComponent, decompressFromEncodedURIComponent } from 'lz-string';
import type { IJsonModel } from 'flexlayout-react';
import type {
  AirfoilPoint,
  AirfoilState,
  CamberControlPoint,
  ControlMode,
  FlapDefinition,
  RunRow,
  SpacingKnot,
  ThicknessControlPoint,
  VisualizationState,
} from '../types';
import { isPanelId, type PanelId } from '../layoutConfig';
import {
  AUTO_GROUP_KEY,
  DEFAULT_ROUTE_UI_STATE,
  type RouteUiSnapshot,
  type RouteViewportState,
} from '../stores/routeUiStore';
import { loadFromUrl, parseNacaFromName } from './urlState';

export type RouteTheme = 'dark' | 'light';

interface HeavyRouteData {
  foil?: SerializedAirfoilData;
  spacing?: SpacingKnot[];
  layout?: IJsonModel;
  filters?: unknown;
  flaps?: FlapDefinition[];
}

interface SerializedAirfoilData {
  name: string;
  coordinates: AirfoilPoint[];
  panels: AirfoilPoint[];
  baseCoordinates: AirfoilPoint[];
  camberControlPoints: CamberControlPoint[];
  thicknessControlPoints: ThicknessControlPoint[];
}

export interface RouteStateSnapshot {
  basePath: string;
  panel: PanelId;
  theme: RouteTheme | null;
  airfoil: Partial<AirfoilState> & {
    nacaCode?: string;
    nacaPoints?: number;
    exactGeometry?: SerializedAirfoilData;
    spacingKnots?: SpacingKnot[];
  };
  visualization: Partial<VisualizationState>;
  ui: Partial<RouteUiSnapshot>;
  layoutJson: IJsonModel | null;
}

export interface CanonicalRouteStateSnapshot extends RouteStateSnapshot {
  theme: RouteTheme;
  ui: RouteUiSnapshot;
}

interface BuildRouteStateArgs {
  theme: RouteTheme;
  airfoil: Pick<
    AirfoilState,
    | 'name'
    | 'coordinates'
    | 'panels'
    | 'baseCoordinates'
    | 'camberControlPoints'
    | 'thicknessControlPoints'
    | 'controlMode'
    | 'displayAlpha'
    | 'reynolds'
    | 'mach'
    | 'ncrit'
    | 'maxIterations'
    | 'solverMode'
    | 'thicknessScale'
    | 'camberScale'
    | 'nPanels'
    | 'spacingKnots'
    | 'curvatureWeight'
    | 'spacingPanelMode'
    | 'sspInterpolation'
    | 'sspVisualization'
    | 'geometryDesign'
  >;
  visualization: Pick<
    VisualizationState,
    | 'showGrid'
    | 'showCurve'
    | 'showPanels'
    | 'showPoints'
    | 'showControls'
    | 'showStreamlines'
    | 'showSmoke'
    | 'showPsiContours'
    | 'showCp'
    | 'showForces'
    | 'showBoundaryLayer'
    | 'showWake'
    | 'showDisplacementThickness'
    | 'enableMorphing'
    | 'morphDuration'
    | 'streamlineDensity'
    | 'adaptiveStreamlines'
    | 'smokeDensity'
    | 'smokeParticlesPerBlob'
    | 'smokeWaveSpacing'
    | 'flowSpeed'
    | 'cpDisplayMode'
    | 'cpBarScale'
    | 'forceScale'
    | 'blThicknessScale'
    | 'useGPU'
  >;
  ui: RouteUiSnapshot;
}

function encodeCompressed<T>(value: T): string {
  return compressToEncodedURIComponent(JSON.stringify(value));
}

function decodeCompressed<T>(value: string | null): T | null {
  if (!value) return null;
  try {
    const decompressed = decompressFromEncodedURIComponent(value);
    return decompressed ? (JSON.parse(decompressed) as T) : null;
  } catch (error) {
    console.warn('Failed to decode compressed route state:', error);
    return null;
  }
}

function parseNumber(value: string | null): number | undefined {
  if (value == null || value === '') return undefined;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : undefined;
}

function parseBoolean(value: string | null): boolean | undefined {
  if (value == null || value === '') return undefined;
  if (value === '1' || value === 'true') return true;
  if (value === '0' || value === 'false') return false;
  return undefined;
}

function setNumberParam(params: URLSearchParams, key: string, value: number, defaultValue?: number): void {
  if (defaultValue === undefined || value !== defaultValue) {
    params.set(key, String(value));
  }
}

export function deriveBasePath(pathname: string): string {
  const trimmed = pathname.replace(/\/+$/, '') || '/';
  const segments = trimmed.split('/').filter(Boolean);
  const maybePanel = segments.at(-1);
  if (isPanelId(maybePanel)) {
    const baseSegments = segments.slice(0, -1);
    return baseSegments.length > 0 ? `/${baseSegments.join('/')}` : '/';
  }
  return trimmed || '/';
}

export function extractPanelFromPath(pathname: string): PanelId | null {
  const trimmed = pathname.replace(/\/+$/, '') || '/';
  const segments = trimmed.split('/').filter(Boolean);
  const maybePanel = segments.at(-1);
  return isPanelId(maybePanel) ? maybePanel : null;
}

function buildPathname(basePath: string, panel: PanelId): string {
  const trimmedBase = basePath === '/' ? '' : basePath.replace(/\/+$/, '');
  return `${trimmedBase}/${panel}` || `/${panel}`;
}

function toExactGeometry(airfoil: BuildRouteStateArgs['airfoil']): SerializedAirfoilData {
  return {
    name: airfoil.name,
    coordinates: airfoil.coordinates,
    panels: airfoil.panels,
    baseCoordinates: airfoil.baseCoordinates,
    camberControlPoints: airfoil.camberControlPoints,
    thicknessControlPoints: airfoil.thicknessControlPoints,
  };
}

export function buildRouteStateSnapshot({
  theme,
  airfoil,
  visualization,
  ui,
}: BuildRouteStateArgs): CanonicalRouteStateSnapshot {
  const nacaCode = parseNacaFromName(airfoil.name) ?? ui.libraryNacaCode;
  const spacingKnots = airfoil.spacingKnots;

  return {
    basePath: '/',
    panel: ui.activePanel,
    theme,
    airfoil: {
      nacaCode,
      nacaPoints: ui.libraryNPoints,
      exactGeometry: toExactGeometry(airfoil),
      controlMode: airfoil.controlMode,
      displayAlpha: airfoil.displayAlpha,
      reynolds: airfoil.reynolds,
      mach: airfoil.mach,
      ncrit: airfoil.ncrit,
      maxIterations: airfoil.maxIterations,
      solverMode: airfoil.solverMode,
      thicknessScale: airfoil.thicknessScale,
      camberScale: airfoil.camberScale,
      nPanels: airfoil.nPanels,
      spacingKnots,
      curvatureWeight: airfoil.curvatureWeight,
      spacingPanelMode: airfoil.spacingPanelMode,
      sspInterpolation: airfoil.sspInterpolation,
      sspVisualization: airfoil.sspVisualization,
      geometryDesign: airfoil.geometryDesign,
    },
    visualization: { ...visualization },
    ui,
    layoutJson: ui.layoutJson,
  };
}

export function serializeRouteState(snapshot: CanonicalRouteStateSnapshot, basePath: string): string {
  const params = new URLSearchParams();
  const { airfoil } = snapshot;

  if (airfoil.nacaCode) params.set('naca', airfoil.nacaCode);
  if (airfoil.nacaPoints !== undefined) setNumberParam(params, 'npts', airfoil.nacaPoints, DEFAULT_ROUTE_UI_STATE.libraryNPoints);
  if (airfoil.controlMode) params.set('mode', airfoil.controlMode);
  if (airfoil.displayAlpha !== undefined) setNumberParam(params, 'alpha', airfoil.displayAlpha, 0);
  if (airfoil.reynolds !== undefined) setNumberParam(params, 're', airfoil.reynolds, 1e6);
  if (airfoil.solverMode) params.set('solver', airfoil.solverMode);
  if (airfoil.mach !== undefined) setNumberParam(params, 'mach', airfoil.mach, 0);
  if (airfoil.ncrit !== undefined) setNumberParam(params, 'ncrit', airfoil.ncrit, 9);
  if (airfoil.maxIterations !== undefined) setNumberParam(params, 'iter', airfoil.maxIterations, 100);
  if (airfoil.thicknessScale !== undefined) setNumberParam(params, 'thickness', airfoil.thicknessScale, 1);
  if (airfoil.camberScale !== undefined) setNumberParam(params, 'camber', airfoil.camberScale, 1);
  if (airfoil.nPanels !== undefined) setNumberParam(params, 'npan', airfoil.nPanels, 160);
  if (airfoil.curvatureWeight !== undefined) setNumberParam(params, 'curvature', airfoil.curvatureWeight, 0);

  const heavy: HeavyRouteData = {};
  if (airfoil.spacingKnots?.length) heavy.spacing = airfoil.spacingKnots;
  if (airfoil.exactGeometry) heavy.foil = airfoil.exactGeometry;
  const flaps = airfoil.geometryDesign?.flaps;
  if (flaps && flaps.length > 0) {
    heavy.flaps = flaps;
    console.log('[routeState] Serializing', flaps.length, 'flap(s) into URL');
  }

  const pathname = buildPathname(basePath, snapshot.panel);
  const search = params.toString();
  const hasHeavy = Object.keys(heavy).length > 0;
  const fragment = hasHeavy ? `#h=${encodeCompressed(heavy)}` : '';

  if (search) return `${pathname}?${search}${fragment}`;
  return `${pathname}${fragment}`;
}

function buildLegacyFallback(basePath: string, hash?: string): RouteStateSnapshot | null {
  const legacy =
    typeof window !== 'undefined'
      ? loadFromUrl()
      : hash && hash.startsWith('#s=')
        ? null
        : null;
  if (!legacy) return null;

  return {
    basePath,
    panel: DEFAULT_ROUTE_UI_STATE.activePanel,
    theme: legacy.theme ?? 'dark',
    airfoil: {
      nacaCode: legacy.naca,
      controlMode: legacy.mode,
      displayAlpha: legacy.alpha,
      thicknessScale: legacy.thicknessScale,
      camberScale: legacy.camberScale,
      nPanels: legacy.nPanels,
      spacingKnots: legacy.spacing,
      curvatureWeight: legacy.curvatureWeight,
      spacingPanelMode: legacy.spacingPanelMode,
      sspInterpolation: legacy.sspInterpolation,
      sspVisualization: legacy.sspVisualization,
    },
    visualization: {
      showGrid: legacy.showGrid,
      showCurve: legacy.showCurve,
      showPanels: legacy.showPanels,
      showPoints: legacy.showPoints,
      showControls: legacy.showControls,
      showStreamlines: legacy.showStreamlines,
      showSmoke: legacy.showSmoke,
      showPsiContours: legacy.showPsiContours,
      showCp: legacy.showCp,
      showForces: legacy.showForces,
      enableMorphing: legacy.enableMorphing,
      morphDuration: legacy.morphDuration,
      streamlineDensity: legacy.streamlineDensity,
      adaptiveStreamlines: legacy.adaptiveStreamlines,
      smokeDensity: legacy.smokeDensity,
      smokeParticlesPerBlob: legacy.smokeParticlesPerBlob,
      flowSpeed: legacy.flowSpeed,
      cpDisplayMode: legacy.cpDisplayMode,
      cpBarScale: legacy.cpBarScale,
      forceScale: legacy.forceScale,
    },
    ui: {
      viewport: {
        centerX: legacy.viewCenterX ?? DEFAULT_ROUTE_UI_STATE.viewport.centerX,
        centerY: legacy.viewCenterY ?? DEFAULT_ROUTE_UI_STATE.viewport.centerY,
        zoom: legacy.viewZoom ?? DEFAULT_ROUTE_UI_STATE.viewport.zoom,
      },
      solvePolarStart: legacy.polarStart ?? DEFAULT_ROUTE_UI_STATE.solvePolarStart,
      solvePolarEnd: legacy.polarEnd ?? DEFAULT_ROUTE_UI_STATE.solvePolarEnd,
      solvePolarStep: legacy.polarStep ?? DEFAULT_ROUTE_UI_STATE.solvePolarStep,
    },
    layoutJson: (legacy.layout as IJsonModel | undefined) ?? null,
  };
}

export function parseRouteStateFromLocation(
  locationLike: Pick<Location, 'pathname' | 'search' | 'hash'>,
): RouteStateSnapshot | null {
  const basePath = deriveBasePath(locationLike.pathname);
  const panel = extractPanelFromPath(locationLike.pathname) ?? DEFAULT_ROUTE_UI_STATE.activePanel;
  const params = new URLSearchParams(locationLike.search);

  // Decode heavy data from fragment (new format) or fall back to query params (old URLs)
  let heavy: HeavyRouteData | null = null;
  const hash = locationLike.hash;
  if (hash.startsWith('#h=')) {
    heavy = decodeCompressed<HeavyRouteData>(hash.slice(3));
  }

  if ([...params.keys()].length === 0 && !heavy) {
    return (
      buildLegacyFallback(basePath, hash) ?? {
        basePath,
        panel,
        theme: null,
        airfoil: {},
        visualization: {},
        ui: {
          activePanel: panel,
        },
        layoutJson: null,
      }
    );
  }

  const spacingKnots = heavy?.spacing ?? decodeCompressed<SpacingKnot[]>(params.get('spacing')) ?? undefined;
  const exactGeometry = heavy?.foil ?? decodeCompressed<SerializedAirfoilData>(params.get('foil')) ?? undefined;
  const restoredFlaps = heavy?.flaps;
  if (restoredFlaps?.length) console.log('[routeState] Decoded', restoredFlaps.length, 'flap(s) from URL');
  const layoutJson = heavy?.layout ?? decodeCompressed<IJsonModel>(params.get('layout'));
  const filterModel = heavy?.filters ?? decodeCompressed(params.get('filters'));

  const plotGroup = params.get('plotGroup');
  const dataView = params.get('dataView');
  const splom = params.get('splom');
  const themeParam = params.get('theme');
  const theme: RouteTheme | null = themeParam === 'light' ? 'light' : themeParam === 'dark' ? 'dark' : null;

  const viewport: Partial<RouteViewportState> = {};
  const centerX = parseNumber(params.get('cx'));
  const centerY = parseNumber(params.get('cy'));
  const zoom = parseNumber(params.get('zoom'));
  if (centerX !== undefined) viewport.centerX = centerX;
  if (centerY !== undefined) viewport.centerY = centerY;
  if (zoom !== undefined) viewport.zoom = zoom;

  return {
    basePath,
    panel,
    theme,
    airfoil: {
      nacaCode: params.get('naca') ?? undefined,
      nacaPoints: parseNumber(params.get('npts')),
      exactGeometry,
      controlMode: params.get('mode') as ControlMode | null ?? undefined,
      displayAlpha: parseNumber(params.get('alpha')),
      reynolds: parseNumber(params.get('re')),
      solverMode: (params.get('solver') as AirfoilState['solverMode'] | null) ?? undefined,
      mach: parseNumber(params.get('mach')),
      ncrit: parseNumber(params.get('ncrit')),
      maxIterations: parseNumber(params.get('iter')),
      thicknessScale: parseNumber(params.get('thickness')),
      camberScale: parseNumber(params.get('camber')),
      nPanels: parseNumber(params.get('npan')),
      spacingKnots,
      curvatureWeight: parseNumber(params.get('curvature')),
      spacingPanelMode: (params.get('spacingMode') as AirfoilState['spacingPanelMode'] | null) ?? undefined,
      sspInterpolation: (params.get('sspInterp') as AirfoilState['sspInterpolation'] | null) ?? undefined,
      sspVisualization: (params.get('sspViz') as AirfoilState['sspVisualization'] | null) ?? undefined,
      ...(restoredFlaps && restoredFlaps.length > 0
        ? { geometryDesign: { flaps: restoredFlaps, teGap: 0, teGapBlend: 0.8, leRadiusFactor: 1.0 } }
        : {}),
    },
    visualization: {
      showGrid: parseBoolean(params.get('grid')),
      showCurve: parseBoolean(params.get('curve')),
      showPanels: parseBoolean(params.get('panels')),
      showPoints: parseBoolean(params.get('points')),
      showControls: parseBoolean(params.get('controls')),
      showStreamlines: parseBoolean(params.get('streamlines')),
      showSmoke: parseBoolean(params.get('smoke')),
      showPsiContours: parseBoolean(params.get('psi')),
      showCp: parseBoolean(params.get('cp')),
      showForces: parseBoolean(params.get('forces')),
      showBoundaryLayer: parseBoolean(params.get('bl')),
      showWake: parseBoolean(params.get('wake')),
      showDisplacementThickness: parseBoolean(params.get('deltaStar')),
      enableMorphing: parseBoolean(params.get('morph')),
      adaptiveStreamlines: parseBoolean(params.get('adaptive')),
      useGPU: parseBoolean(params.get('gpu')),
      morphDuration: parseNumber(params.get('morphMs')),
      streamlineDensity: parseNumber(params.get('streamDensity')),
      smokeDensity: parseNumber(params.get('smokeDensity')),
      smokeParticlesPerBlob: parseNumber(params.get('smokeBlob')),
      smokeWaveSpacing: parseNumber(params.get('smokeWave')),
      flowSpeed: parseNumber(params.get('flow')),
      cpDisplayMode: (params.get('cpMode') as VisualizationState['cpDisplayMode'] | null) ?? undefined,
      cpBarScale: parseNumber(params.get('cpScale')),
      forceScale: parseNumber(params.get('forceScale')),
      blThicknessScale: parseNumber(params.get('blScale')),
    },
    ui: {
      libraryNacaCode: params.get('naca') ?? undefined,
      libraryNPoints: parseNumber(params.get('npts')),
      spacingLiveUpdate: parseBoolean(params.get('live')),
      solveRunMode: (params.get('runMode') as RouteUiSnapshot['solveRunMode'] | null) ?? undefined,
      solveShowAdvanced: parseBoolean(params.get('adv')),
      solveTargetCl: parseNumber(params.get('targetCl')),
      solvePolarStart: parseNumber(params.get('polarStart')),
      solvePolarEnd: parseNumber(params.get('polarEnd')),
      solvePolarStep: parseNumber(params.get('polarStep')),
      polarXAxis: (params.get('polarX') as RouteUiSnapshot['polarXAxis'] | null) ?? undefined,
      polarYAxis: (params.get('polarY') as RouteUiSnapshot['polarYAxis'] | null) ?? undefined,
      plotDataSource: (params.get('plotData') as RouteUiSnapshot['plotDataSource'] | null) ?? undefined,
      plotChartType: (params.get('plotType') as RouteUiSnapshot['plotChartType'] | null) ?? undefined,
      plotXField: (params.get('plotX') as keyof RunRow | null) ?? undefined,
      plotYField: (params.get('plotY') as keyof RunRow | null) ?? undefined,
      plotGroupBy: plotGroup
        ? ((plotGroup === AUTO_GROUP_KEY ? AUTO_GROUP_KEY : plotGroup) as RouteUiSnapshot['plotGroupBy'])
        : undefined,
      plotXScale: (params.get('plotXScale') as RouteUiSnapshot['plotXScale'] | null) ?? undefined,
      plotYScale: (params.get('plotYScale') as RouteUiSnapshot['plotYScale'] | null) ?? undefined,
      dataExplorerView: (dataView as RouteUiSnapshot['dataExplorerView'] | null) ?? undefined,
      dataExplorerSplomKeys: splom ? (splom.split(',').filter(Boolean) as (keyof RunRow)[]) : undefined,
      dataExplorerColorBy: (params.get('colorBy') as keyof RunRow | null) ?? undefined,
      dataExplorerFilterModel: filterModel,
      activePanel: panel,
      viewport: Object.keys(viewport).length > 0 ? (viewport as RouteViewportState) : undefined,
      layoutJson: layoutJson ?? undefined,
    },
    layoutJson,
  };
}

export function writeRouteState(
  snapshot: CanonicalRouteStateSnapshot,
  options: { basePath: string; mode?: 'replace' | 'push' },
): void {
  const href = serializeRouteState(snapshot, options.basePath);
  const method = options.mode === 'push' ? 'pushState' : 'replaceState';
  window.history[method](null, '', href);
}
