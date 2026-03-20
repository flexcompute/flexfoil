import { describe, expect, it } from 'vitest';
import { DEFAULT_ROUTE_UI_STATE } from '../stores/routeUiStore';
import { DEFAULT_VISUALIZATION_STATE } from '../stores/visualizationStore';
import {
  buildRouteStateSnapshot,
  parseRouteStateFromLocation,
  serializeRouteState,
} from './routeState';

function parseHref(href: string, _basePath?: string) {
  const url = new URL(href, 'http://localhost');
  return parseRouteStateFromLocation({
    pathname: url.pathname,
    search: url.search,
    hash: url.hash,
  });
}

describe('routeState', () => {
  it('round-trips essential state through URL', () => {
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
        geometryDesign: { flaps: [], teGap: 0, teGapBlend: 0.8, leRadiusFactor: 1.0 },
      },
      visualization: {
        ...DEFAULT_VISUALIZATION_STATE,
        showCp: true,
        showBoundaryLayer: true,
      },
      ui: {
        ...DEFAULT_ROUTE_UI_STATE,
        activePanel: 'solve',
        solveShowAdvanced: true,
        viewport: { centerX: 0.4, centerY: 0.1, zoom: 550 },
      },
    });

    const href = serializeRouteState(snapshot, '/app');
    const parsed = parseHref(href, '/app');

    expect(parsed?.panel).toBe('solve');
    expect(parsed?.airfoil.controlMode).toBe('camber-spline');
    expect(parsed?.airfoil.displayAlpha).toBe(4);
    expect(parsed?.airfoil.reynolds).toBe(2_000_000);
    expect(parsed?.airfoil.mach).toBe(0.1);
    expect(parsed?.airfoil.ncrit).toBe(7);
    expect(parsed?.airfoil.maxIterations).toBe(120);
    expect(parsed?.airfoil.thicknessScale).toBe(1.1);
    expect(parsed?.airfoil.camberScale).toBe(0.9);
    expect(parsed?.airfoil.nPanels).toBe(180);
    expect(parsed?.airfoil.curvatureWeight).toBe(0.4);
    expect(parsed?.airfoil.exactGeometry?.camberControlPoints).toEqual([{ id: 'c1', x: 0.5, y: 0.02 }]);
    expect(parsed?.airfoil.spacingKnots).toEqual([{ S: 0, F: 1 }, { S: 1, F: 1 }]);
  });

  it('does not serialize ephemeral state to URL', () => {
    const snapshot = buildRouteStateSnapshot({
      theme: 'light',
      airfoil: {
        name: 'NACA 2412',
        coordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        panels: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        baseCoordinates: [{ x: 0, y: 0 }, { x: 1, y: 0 }],
        camberControlPoints: [],
        thicknessControlPoints: [],
        controlMode: 'parameters',
        displayAlpha: 4,
        reynolds: 1e6,
        mach: 0,
        ncrit: 9,
        maxIterations: 100,
        solverMode: 'viscous',
        thicknessScale: 1,
        camberScale: 1,
        nPanels: 160,
        spacingKnots: [],
        curvatureWeight: 0,
        spacingPanelMode: 'simple',
        sspInterpolation: 'linear',
        sspVisualization: 'plot',
        geometryDesign: { flaps: [], teGap: 0, teGapBlend: 0.8, leRadiusFactor: 1.0 },
      },
      visualization: {
        ...DEFAULT_VISUALIZATION_STATE,
        showCp: true,
        showStreamlines: true,
        streamlineDensity: 80,
      },
      ui: {
        ...DEFAULT_ROUTE_UI_STATE,
        activePanel: 'canvas',
        solveShowAdvanced: true,
        polarXAxis: 'cd',
        plotChartType: 'line',
        dataExplorerView: 'correlogram',
        viewport: { centerX: 0.4, centerY: 0.1, zoom: 550 },
      },
    });

    const href = serializeRouteState(snapshot, '/');
    expect(href).not.toContain('theme=');
    expect(href).not.toContain('grid=');
    expect(href).not.toContain('cp=');
    expect(href).not.toContain('streamlines=');
    expect(href).not.toContain('streamDensity=');
    expect(href).not.toContain('adv=');
    expect(href).not.toContain('polarX=');
    expect(href).not.toContain('plotType=');
    expect(href).not.toContain('dataView=');
    expect(href).not.toContain('cx=');
    expect(href).not.toContain('zoom=');
    expect(href).not.toContain('runMode=');
    expect(href).toContain('naca=2412');
    expect(href).toContain('alpha=4');
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

  it('still parses legacy URLs with ephemeral params (backward compat)', () => {
    const parsed = parseRouteStateFromLocation({
      pathname: '/solve',
      search: '?naca=2412&alpha=4&re=2000000&cp=1&bl=1&adv=1&polarX=cd&cx=0.4&zoom=550&theme=light',
      hash: '',
    });

    expect(parsed?.airfoil.nacaCode).toBe('2412');
    expect(parsed?.airfoil.displayAlpha).toBe(4);
    expect(parsed?.airfoil.reynolds).toBe(2_000_000);
    expect(parsed?.visualization.showCp).toBe(true);
    expect(parsed?.visualization.showBoundaryLayer).toBe(true);
    expect(parsed?.ui.solveShowAdvanced).toBe(true);
    expect(parsed?.ui.polarXAxis).toBe('cd');
    expect(parsed?.ui.viewport).toEqual({ centerX: 0.4, zoom: 550 });
    expect(parsed?.theme).toBe('light');
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
        geometryDesign: { flaps: [], teGap: 0, teGapBlend: 0.8, leRadiusFactor: 1.0 },
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
