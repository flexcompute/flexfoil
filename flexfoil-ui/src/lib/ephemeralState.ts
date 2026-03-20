import type { IJsonModel } from 'flexlayout-react';
import type { AirfoilState, VisualizationState } from '../types';
import type { RouteUiSnapshot, RouteViewportState } from '../stores/routeUiStore';
import type { RouteStateSnapshot, RouteTheme } from './routeState';

const STORAGE_KEY = 'flexfoil-ephemeral-v1';

interface EphemeralPayload {
  theme?: RouteTheme;
  visualization?: Partial<VisualizationState>;
  ui?: Partial<Omit<RouteUiSnapshot, 'activePanel'>>;
  airfoilUi?: {
    spacingPanelMode?: AirfoilState['spacingPanelMode'];
    sspInterpolation?: AirfoilState['sspInterpolation'];
    sspVisualization?: AirfoilState['sspVisualization'];
  };
}

export function saveEphemeralState(snapshot: {
  theme: RouteTheme;
  visualization: Partial<VisualizationState>;
  ui: Partial<RouteUiSnapshot> & { viewport?: RouteViewportState; layoutJson?: IJsonModel | null };
  airfoil?: Partial<AirfoilState>;
  layoutJson?: IJsonModel | null;
}): void {
  try {
    const payload: EphemeralPayload = {
      theme: snapshot.theme,
      visualization: snapshot.visualization,
      ui: {
        spacingLiveUpdate: snapshot.ui.spacingLiveUpdate,
        solveRunMode: snapshot.ui.solveRunMode,
        solveShowAdvanced: snapshot.ui.solveShowAdvanced,
        solveTargetCl: snapshot.ui.solveTargetCl,
        solvePolarStart: snapshot.ui.solvePolarStart,
        solvePolarEnd: snapshot.ui.solvePolarEnd,
        solvePolarStep: snapshot.ui.solvePolarStep,
        polarXAxis: snapshot.ui.polarXAxis,
        polarYAxis: snapshot.ui.polarYAxis,
        plotDataSource: snapshot.ui.plotDataSource,
        plotChartType: snapshot.ui.plotChartType,
        plotXField: snapshot.ui.plotXField,
        plotYField: snapshot.ui.plotYField,
        plotGroupBy: snapshot.ui.plotGroupBy,
        plotXScale: snapshot.ui.plotXScale,
        plotYScale: snapshot.ui.plotYScale,
        dataExplorerView: snapshot.ui.dataExplorerView,
        dataExplorerSplomKeys: snapshot.ui.dataExplorerSplomKeys,
        dataExplorerColorBy: snapshot.ui.dataExplorerColorBy,
        dataExplorerFilterModel: snapshot.ui.dataExplorerFilterModel,
        dataExplorerSmartGroup: snapshot.ui.dataExplorerSmartGroup,
        dataExplorerDataSource: snapshot.ui.dataExplorerDataSource,
        layoutJson: snapshot.layoutJson ?? snapshot.ui.layoutJson,
        viewport: snapshot.ui.viewport,
      },
      airfoilUi: snapshot.airfoil
        ? {
            spacingPanelMode: snapshot.airfoil.spacingPanelMode,
            sspInterpolation: snapshot.airfoil.sspInterpolation,
            sspVisualization: snapshot.airfoil.sspVisualization,
          }
        : undefined,
    };
    localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
  } catch {
    // localStorage unavailable or quota exceeded
  }
}

export function loadEphemeralState(): EphemeralPayload | null {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return null;
    return JSON.parse(raw) as EphemeralPayload;
  } catch {
    return null;
  }
}

/**
 * Merge ephemeral (localStorage) state into a URL-parsed snapshot.
 * URL-defined fields take priority; ephemeral fills undefined gaps.
 */
export function mergeEphemeralIntoSnapshot(
  snapshot: RouteStateSnapshot,
  ephemeral: EphemeralPayload | null,
): RouteStateSnapshot {
  if (!ephemeral) return { ...snapshot, theme: snapshot.theme ?? 'dark' };

  const mergedViz: Partial<VisualizationState> = { ...ephemeral.visualization };
  for (const [k, v] of Object.entries(snapshot.visualization)) {
    if (v !== undefined) (mergedViz as Record<string, unknown>)[k] = v;
  }

  const mergedUi: Partial<RouteUiSnapshot> = {};
  if (ephemeral.ui) {
    for (const [k, v] of Object.entries(ephemeral.ui)) {
      if (v !== undefined) (mergedUi as Record<string, unknown>)[k] = v;
    }
  }
  for (const [k, v] of Object.entries(snapshot.ui)) {
    if (v !== undefined) (mergedUi as Record<string, unknown>)[k] = v;
  }

  const mergedAirfoil = { ...snapshot.airfoil };
  if (ephemeral.airfoilUi) {
    if (mergedAirfoil.spacingPanelMode === undefined && ephemeral.airfoilUi.spacingPanelMode) {
      mergedAirfoil.spacingPanelMode = ephemeral.airfoilUi.spacingPanelMode;
    }
    if (mergedAirfoil.sspInterpolation === undefined && ephemeral.airfoilUi.sspInterpolation) {
      mergedAirfoil.sspInterpolation = ephemeral.airfoilUi.sspInterpolation;
    }
    if (mergedAirfoil.sspVisualization === undefined && ephemeral.airfoilUi.sspVisualization) {
      mergedAirfoil.sspVisualization = ephemeral.airfoilUi.sspVisualization;
    }
  }

  const mergedLayoutJson =
    snapshot.layoutJson ?? (ephemeral.ui?.layoutJson !== undefined ? ephemeral.ui.layoutJson : null);

  return {
    ...snapshot,
    theme: snapshot.theme ?? ephemeral.theme ?? 'dark',
    airfoil: mergedAirfoil,
    visualization: mergedViz,
    ui: mergedUi,
    layoutJson: mergedLayoutJson,
  };
}
