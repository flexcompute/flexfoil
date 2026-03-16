import { describe, expect, it } from 'vitest';
import { defaultLayoutJson } from '../layoutConfig';
import { DEFAULT_ROUTE_UI_STATE } from '../stores/routeUiStore';
import { DEFAULT_VISUALIZATION_STATE } from '../stores/visualizationStore';
import {
  buildRouteStateSnapshot,
  parseRouteStateFromLocation,
  serializeRouteState,
} from './routeState';

describe('routeState', () => {
  it('round-trips a fully populated route snapshot', () => {
    const snapshot = buildRouteStateSnapshot({
      theme: 'light',
      airfoil: {
        name: 'Custom Foil',
        coordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        panels: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        baseCoordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        camberControlPoints: [{ id: 'c1', x: 0.5, y: 0.02 }],
        thicknessControlPoints: [{ id: 't1', x: 0.5, t: 0.08 }],
        controlMode: 'camber-spline',
        displayAlpha: 4,
        reynolds: 2_000_000,
        mach: 0.1,
        ncrit: 7,
        maxIterations: 120,
        solverMode: 'viscous',
        thicknessScale: 1.1,
        camberScale: 0.9,
        nPanels: 180,
        spacingKnots: [{ S: 0, F: 1 }, { S: 1, F: 1 }],
        curvatureWeight: 0.4,
        spacingPanelMode: 'advanced',
        sspInterpolation: 'linear',
        sspVisualization: 'foil',
      },
      visualization: {
        ...DEFAULT_VISUALIZATION_STATE,
        showCp: true,
        showBoundaryLayer: true,
        showWake: true,
        blThicknessScale: 12,
        smokeWaveSpacing: 0.8,
      },
      ui: {
        ...DEFAULT_ROUTE_UI_STATE,
        activePanel: 'solve',
        solveShowAdvanced: true,
        solveTargetCl: 0.8,
        solvePolarStart: -3,
        solvePolarEnd: 18,
        solvePolarStep: 0.5,
        polarXAxis: 'cd',
        polarYAxis: 'cl',
        plotChartType: 'line',
        plotXField: 'cd',
        plotYField: 'cl',
        plotGroupBy: 'airfoil_name',
        dataExplorerView: 'correlogram',
        dataExplorerSplomKeys: ['alpha', 'cl', 'cd'],
        dataExplorerColorBy: 'reynolds',
        dataExplorerFilterModel: { alpha: { filterType: 'number', type: 'greaterThan', filter: 2 } },
        viewport: { centerX: 0.4, centerY: 0.1, zoom: 550 },
        layoutJson: defaultLayoutJson,
      },
    });

    const href = serializeRouteState(snapshot, '/app');
    const [pathname, search = ''] = href.split('?');
    const parsed = parseRouteStateFromLocation({
      pathname,
      search: search ? `?${search}` : '',
      hash: '',
    });

    expect(parsed?.panel).toBe('solve');
    expect(parsed?.theme).toBe('light');
    expect(parsed?.airfoil.controlMode).toBe('camber-spline');
    expect(parsed?.airfoil.maxIterations).toBe(120);
    expect(parsed?.visualization.showCp).toBe(true);
    expect(parsed?.visualization.showBoundaryLayer).toBe(true);
    expect(parsed?.ui.solveShowAdvanced).toBe(true);
    expect(parsed?.ui.polarXAxis).toBe('cd');
    expect(parsed?.ui.dataExplorerView).toBe('correlogram');
    expect(parsed?.ui.dataExplorerSplomKeys).toEqual(['alpha', 'cl', 'cd']);
    expect(parsed?.ui.dataExplorerFilterModel).toEqual({
      alpha: { filterType: 'number', type: 'greaterThan', filter: 2 },
    });
    expect(parsed?.ui.viewport).toEqual({ centerX: 0.4, centerY: 0.1, zoom: 550 });
    expect(parsed?.airfoil.exactGeometry?.camberControlPoints).toEqual([{ id: 'c1', x: 0.5, y: 0.02 }]);
  });

  it('supports path-only panel routes', () => {
    const parsed = parseRouteStateFromLocation({
      pathname: '/app/visualization',
      search: '',
      hash: '',
    });

    expect(parsed?.basePath).toBe('/app');
    expect(parsed?.panel).toBe('visualization');
    expect(parsed?.ui.activePanel).toBe('visualization');
  });

  it('omits default values from serialized search params', () => {
    const snapshot = buildRouteStateSnapshot({
      theme: 'dark',
      airfoil: {
        name: 'NACA 0012',
        coordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        panels: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        baseCoordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        camberControlPoints: [],
        thicknessControlPoints: [],
        controlMode: 'parameters',
        displayAlpha: 0,
        reynolds: 1e6,
        mach: 0,
        ncrit: 9,
        maxIterations: 100,
        solverMode: 'viscous',
        thicknessScale: 1,
        camberScale: 1,
        nPanels: 160,
        spacingKnots: [{ S: 0, F: 1.5 }, { S: 0.25, F: 0.4 }, { S: 0.5, F: 1.5 }, { S: 0.75, F: 0.4 }, { S: 1, F: 1.5 }],
        curvatureWeight: 0,
        spacingPanelMode: 'simple',
        sspInterpolation: 'linear',
        sspVisualization: 'plot',
      },
      visualization: DEFAULT_VISUALIZATION_STATE,
      ui: DEFAULT_ROUTE_UI_STATE,
    });

    const href = serializeRouteState(snapshot, '/');
    expect(href.startsWith('/canvas')).toBe(true);
    expect(href).not.toContain('theme=');
    expect(href).not.toContain('runMode=');
  });
});
