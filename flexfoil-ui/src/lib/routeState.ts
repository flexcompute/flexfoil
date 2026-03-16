import { compressToEncodedURIComponent, decompressFromEncodedURIComponent } from 'lz-string';
import type { IJsonModel } from 'flexlayout-react';
import type {
  AirfoilPoint,
  AirfoilState,
  CamberControlPoint,
  ControlMode,
  RunRow,
  SpacingKnot,
  ThicknessControlPoint,
  VisualizationState,
} from '../types';
import { defaultLayoutJson, isPanelId, type PanelId } from '../layoutConfig';
import {
  AUTO_GROUP_KEY,
  DEFAULT_ROUTE_UI_STATE,
  type RouteUiSnapshot,
  type RouteViewportState,
} from '../stores/routeUiStore';
import { DEFAULT_VISUALIZATION_STATE } from '../stores/visualizationStore';
import { loadFromUrl, parseNacaFromName } from './urlState';

export type RouteTheme = 'dark' | 'light';

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
  theme: RouteTheme;
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

function setBooleanParam(
  params: URLSearchParams,
  key: string,
  value: boolean,
  defaultValue: boolean,
): void {
  if (value !== defaultValue) {
    params.set(key, value ? '1' : '0');
  }
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

function isDefaultLayout(layoutJson: IJsonModel | null): boolean {
  if (!layoutJson) return true;
  return JSON.stringify(layoutJson) === JSON.stringify(defaultLayoutJson);
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
    },
    visualization: { ...visualization },
    ui,
    layoutJson: ui.layoutJson,
  };
}

export function serializeRouteState(snapshot: CanonicalRouteStateSnapshot, basePath: string): string {
  const params = new URLSearchParams();
  const { airfoil, visualization, ui } = snapshot;

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
  if (airfoil.spacingPanelMode) params.set('spacingMode', airfoil.spacingPanelMode);
  if (airfoil.sspInterpolation) params.set('sspInterp', airfoil.sspInterpolation);
  if (airfoil.sspVisualization) params.set('sspViz', airfoil.sspVisualization);
  if (airfoil.spacingKnots?.length) params.set('spacing', encodeCompressed(airfoil.spacingKnots));
  if (airfoil.exactGeometry) params.set('foil', encodeCompressed(airfoil.exactGeometry));

  setBooleanParam(params, 'grid', visualization.showGrid ?? DEFAULT_VISUALIZATION_STATE.showGrid, DEFAULT_VISUALIZATION_STATE.showGrid);
  setBooleanParam(params, 'curve', visualization.showCurve ?? DEFAULT_VISUALIZATION_STATE.showCurve, DEFAULT_VISUALIZATION_STATE.showCurve);
  setBooleanParam(params, 'panels', visualization.showPanels ?? DEFAULT_VISUALIZATION_STATE.showPanels, DEFAULT_VISUALIZATION_STATE.showPanels);
  setBooleanParam(params, 'points', visualization.showPoints ?? DEFAULT_VISUALIZATION_STATE.showPoints, DEFAULT_VISUALIZATION_STATE.showPoints);
  setBooleanParam(params, 'controls', visualization.showControls ?? DEFAULT_VISUALIZATION_STATE.showControls, DEFAULT_VISUALIZATION_STATE.showControls);
  setBooleanParam(params, 'streamlines', visualization.showStreamlines ?? DEFAULT_VISUALIZATION_STATE.showStreamlines, DEFAULT_VISUALIZATION_STATE.showStreamlines);
  setBooleanParam(params, 'smoke', visualization.showSmoke ?? DEFAULT_VISUALIZATION_STATE.showSmoke, DEFAULT_VISUALIZATION_STATE.showSmoke);
  setBooleanParam(params, 'psi', visualization.showPsiContours ?? DEFAULT_VISUALIZATION_STATE.showPsiContours, DEFAULT_VISUALIZATION_STATE.showPsiContours);
  setBooleanParam(params, 'cp', visualization.showCp ?? DEFAULT_VISUALIZATION_STATE.showCp, DEFAULT_VISUALIZATION_STATE.showCp);
  setBooleanParam(params, 'forces', visualization.showForces ?? DEFAULT_VISUALIZATION_STATE.showForces, DEFAULT_VISUALIZATION_STATE.showForces);
  setBooleanParam(params, 'bl', visualization.showBoundaryLayer ?? DEFAULT_VISUALIZATION_STATE.showBoundaryLayer, DEFAULT_VISUALIZATION_STATE.showBoundaryLayer);
  setBooleanParam(params, 'wake', visualization.showWake ?? DEFAULT_VISUALIZATION_STATE.showWake, DEFAULT_VISUALIZATION_STATE.showWake);
  setBooleanParam(
    params,
    'deltaStar',
    visualization.showDisplacementThickness ?? DEFAULT_VISUALIZATION_STATE.showDisplacementThickness,
    DEFAULT_VISUALIZATION_STATE.showDisplacementThickness,
  );
  setBooleanParam(params, 'morph', visualization.enableMorphing ?? DEFAULT_VISUALIZATION_STATE.enableMorphing, DEFAULT_VISUALIZATION_STATE.enableMorphing);
  setBooleanParam(params, 'adaptive', visualization.adaptiveStreamlines ?? DEFAULT_VISUALIZATION_STATE.adaptiveStreamlines, DEFAULT_VISUALIZATION_STATE.adaptiveStreamlines);
  setBooleanParam(params, 'gpu', visualization.useGPU ?? DEFAULT_VISUALIZATION_STATE.useGPU, DEFAULT_VISUALIZATION_STATE.useGPU);
  setNumberParam(params, 'morphMs', visualization.morphDuration ?? DEFAULT_VISUALIZATION_STATE.morphDuration, DEFAULT_VISUALIZATION_STATE.morphDuration);
  setNumberParam(params, 'streamDensity', visualization.streamlineDensity ?? DEFAULT_VISUALIZATION_STATE.streamlineDensity, DEFAULT_VISUALIZATION_STATE.streamlineDensity);
  setNumberParam(params, 'smokeDensity', visualization.smokeDensity ?? DEFAULT_VISUALIZATION_STATE.smokeDensity, DEFAULT_VISUALIZATION_STATE.smokeDensity);
  setNumberParam(
    params,
    'smokeBlob',
    visualization.smokeParticlesPerBlob ?? DEFAULT_VISUALIZATION_STATE.smokeParticlesPerBlob,
    DEFAULT_VISUALIZATION_STATE.smokeParticlesPerBlob,
  );
  setNumberParam(
    params,
    'smokeWave',
    visualization.smokeWaveSpacing ?? DEFAULT_VISUALIZATION_STATE.smokeWaveSpacing,
    DEFAULT_VISUALIZATION_STATE.smokeWaveSpacing,
  );
  setNumberParam(params, 'flow', visualization.flowSpeed ?? DEFAULT_VISUALIZATION_STATE.flowSpeed, DEFAULT_VISUALIZATION_STATE.flowSpeed);
  if (visualization.cpDisplayMode && visualization.cpDisplayMode !== DEFAULT_VISUALIZATION_STATE.cpDisplayMode) {
    params.set('cpMode', visualization.cpDisplayMode);
  }
  setNumberParam(params, 'cpScale', visualization.cpBarScale ?? DEFAULT_VISUALIZATION_STATE.cpBarScale, DEFAULT_VISUALIZATION_STATE.cpBarScale);
  setNumberParam(params, 'forceScale', visualization.forceScale ?? DEFAULT_VISUALIZATION_STATE.forceScale, DEFAULT_VISUALIZATION_STATE.forceScale);
  setNumberParam(params, 'blScale', visualization.blThicknessScale ?? DEFAULT_VISUALIZATION_STATE.blThicknessScale, DEFAULT_VISUALIZATION_STATE.blThicknessScale);

  if (snapshot.theme !== 'dark') params.set('theme', snapshot.theme);
  setBooleanParam(params, 'live', ui.spacingLiveUpdate ?? DEFAULT_ROUTE_UI_STATE.spacingLiveUpdate, DEFAULT_ROUTE_UI_STATE.spacingLiveUpdate);
  if ((ui.solveRunMode ?? DEFAULT_ROUTE_UI_STATE.solveRunMode) !== DEFAULT_ROUTE_UI_STATE.solveRunMode) {
    params.set('runMode', ui.solveRunMode ?? DEFAULT_ROUTE_UI_STATE.solveRunMode);
  }
  setBooleanParam(params, 'adv', ui.solveShowAdvanced ?? DEFAULT_ROUTE_UI_STATE.solveShowAdvanced, DEFAULT_ROUTE_UI_STATE.solveShowAdvanced);
  setNumberParam(params, 'targetCl', ui.solveTargetCl ?? DEFAULT_ROUTE_UI_STATE.solveTargetCl, DEFAULT_ROUTE_UI_STATE.solveTargetCl);
  setNumberParam(params, 'polarStart', ui.solvePolarStart ?? DEFAULT_ROUTE_UI_STATE.solvePolarStart, DEFAULT_ROUTE_UI_STATE.solvePolarStart);
  setNumberParam(params, 'polarEnd', ui.solvePolarEnd ?? DEFAULT_ROUTE_UI_STATE.solvePolarEnd, DEFAULT_ROUTE_UI_STATE.solvePolarEnd);
  setNumberParam(params, 'polarStep', ui.solvePolarStep ?? DEFAULT_ROUTE_UI_STATE.solvePolarStep, DEFAULT_ROUTE_UI_STATE.solvePolarStep);
  if ((ui.polarXAxis ?? DEFAULT_ROUTE_UI_STATE.polarXAxis) !== DEFAULT_ROUTE_UI_STATE.polarXAxis) params.set('polarX', ui.polarXAxis ?? DEFAULT_ROUTE_UI_STATE.polarXAxis);
  if ((ui.polarYAxis ?? DEFAULT_ROUTE_UI_STATE.polarYAxis) !== DEFAULT_ROUTE_UI_STATE.polarYAxis) params.set('polarY', ui.polarYAxis ?? DEFAULT_ROUTE_UI_STATE.polarYAxis);
  if ((ui.plotDataSource ?? DEFAULT_ROUTE_UI_STATE.plotDataSource) !== DEFAULT_ROUTE_UI_STATE.plotDataSource) params.set('plotData', ui.plotDataSource ?? DEFAULT_ROUTE_UI_STATE.plotDataSource);
  if ((ui.plotChartType ?? DEFAULT_ROUTE_UI_STATE.plotChartType) !== DEFAULT_ROUTE_UI_STATE.plotChartType) params.set('plotType', ui.plotChartType ?? DEFAULT_ROUTE_UI_STATE.plotChartType);
  if ((ui.plotXField ?? DEFAULT_ROUTE_UI_STATE.plotXField) !== DEFAULT_ROUTE_UI_STATE.plotXField) params.set('plotX', String(ui.plotXField ?? DEFAULT_ROUTE_UI_STATE.plotXField));
  if ((ui.plotYField ?? DEFAULT_ROUTE_UI_STATE.plotYField) !== DEFAULT_ROUTE_UI_STATE.plotYField) params.set('plotY', String(ui.plotYField ?? DEFAULT_ROUTE_UI_STATE.plotYField));
  if ((ui.plotGroupBy ?? DEFAULT_ROUTE_UI_STATE.plotGroupBy) !== DEFAULT_ROUTE_UI_STATE.plotGroupBy) params.set('plotGroup', String(ui.plotGroupBy ?? DEFAULT_ROUTE_UI_STATE.plotGroupBy));
  if ((ui.plotXScale ?? DEFAULT_ROUTE_UI_STATE.plotXScale) !== DEFAULT_ROUTE_UI_STATE.plotXScale) params.set('plotXScale', ui.plotXScale ?? DEFAULT_ROUTE_UI_STATE.plotXScale);
  if ((ui.plotYScale ?? DEFAULT_ROUTE_UI_STATE.plotYScale) !== DEFAULT_ROUTE_UI_STATE.plotYScale) params.set('plotYScale', ui.plotYScale ?? DEFAULT_ROUTE_UI_STATE.plotYScale);
  if ((ui.dataExplorerView ?? DEFAULT_ROUTE_UI_STATE.dataExplorerView) !== DEFAULT_ROUTE_UI_STATE.dataExplorerView) params.set('dataView', ui.dataExplorerView ?? DEFAULT_ROUTE_UI_STATE.dataExplorerView);
  if ((ui.dataExplorerSplomKeys ?? DEFAULT_ROUTE_UI_STATE.dataExplorerSplomKeys).join(',') !== DEFAULT_ROUTE_UI_STATE.dataExplorerSplomKeys.join(',')) {
    params.set('splom', (ui.dataExplorerSplomKeys ?? DEFAULT_ROUTE_UI_STATE.dataExplorerSplomKeys).join(','));
  }
  if (ui.dataExplorerColorBy) params.set('colorBy', String(ui.dataExplorerColorBy));
  if (ui.dataExplorerFilterModel) params.set('filters', encodeCompressed(ui.dataExplorerFilterModel));
  setNumberParam(params, 'cx', ui.viewport?.centerX ?? DEFAULT_ROUTE_UI_STATE.viewport.centerX, DEFAULT_ROUTE_UI_STATE.viewport.centerX);
  setNumberParam(params, 'cy', ui.viewport?.centerY ?? DEFAULT_ROUTE_UI_STATE.viewport.centerY, DEFAULT_ROUTE_UI_STATE.viewport.centerY);
  setNumberParam(params, 'zoom', ui.viewport?.zoom ?? DEFAULT_ROUTE_UI_STATE.viewport.zoom, DEFAULT_ROUTE_UI_STATE.viewport.zoom);
  if (!isDefaultLayout(ui.layoutJson ?? null)) {
    params.set('layout', encodeCompressed(ui.layoutJson ?? defaultLayoutJson));
  }

  const pathname = buildPathname(basePath, snapshot.panel);
  const search = params.toString();
  return search ? `${pathname}?${search}` : pathname;
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

  if ([...params.keys()].length === 0) {
    return (
      buildLegacyFallback(basePath, locationLike.hash) ?? {
        basePath,
        panel,
        theme: 'dark',
        airfoil: {},
        visualization: {},
        ui: {
          activePanel: panel,
        },
        layoutJson: null,
      }
    );
  }

  const spacingKnots = decodeCompressed<SpacingKnot[]>(params.get('spacing')) ?? undefined;
  const exactGeometry = decodeCompressed<SerializedAirfoilData>(params.get('foil')) ?? undefined;
  const layoutJson = decodeCompressed<IJsonModel>(params.get('layout'));
  const filterModel = decodeCompressed(params.get('filters'));

  const plotGroup = params.get('plotGroup');
  const dataView = params.get('dataView');
  const splom = params.get('splom');
  const theme = params.get('theme') === 'light' ? 'light' : 'dark';

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
