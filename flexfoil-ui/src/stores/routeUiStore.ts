import { create } from 'zustand';
import type { IJsonModel } from 'flexlayout-react';
import type {
  AxisScale,
  AxisVariable,
  ChartType,
  DataExplorerViewMode,
  DataSource,
  RunMode,
  RunRow,
  SweepAxis,
} from '../types';
import type { PanelId } from '../layoutConfig';

export const AUTO_GROUP_KEY = '__auto__' as const;

export interface RouteViewportState {
  centerX: number;
  centerY: number;
  zoom: number;
}

export interface RouteUiSnapshot {
  libraryNacaCode: string;
  libraryNPoints: number;
  spacingLiveUpdate: boolean;
  solveRunMode: RunMode;
  solveShowAdvanced: boolean;
  solveTargetCl: number;
  solvePolarStart: number;
  solvePolarEnd: number;
  solvePolarStep: number;
  sweepPrimary: SweepAxis;
  sweepSecondary: SweepAxis | null;
  polarXAxis: AxisVariable;
  polarYAxis: AxisVariable;
  plotDataSource: DataSource;
  plotChartType: ChartType;
  plotXField: keyof RunRow;
  plotYField: keyof RunRow;
  plotGroupBy: keyof RunRow | '' | typeof AUTO_GROUP_KEY;
  plotXScale: AxisScale;
  plotYScale: AxisScale;
  dataExplorerView: DataExplorerViewMode;
  dataExplorerSplomKeys: (keyof RunRow)[];
  dataExplorerColorBy: keyof RunRow | '';
  dataExplorerFilterModel: unknown | null;
  outlierFilterEnabled: boolean;
  activePanel: PanelId;
  layoutJson: IJsonModel | null;
  viewport: RouteViewportState;
}

interface RouteUiStore extends RouteUiSnapshot {
  viewportRevision: number;
  layoutRevision: number;
  activePanelRevision: number;
  setLibraryNacaCode: (value: string) => void;
  setLibraryNPoints: (value: number) => void;
  setSpacingLiveUpdate: (value: boolean) => void;
  setSolveRunMode: (value: RunMode) => void;
  setSolveShowAdvanced: (value: boolean) => void;
  setSolveTargetCl: (value: number) => void;
  setSolvePolarStart: (value: number) => void;
  setSolvePolarEnd: (value: number) => void;
  setSolvePolarStep: (value: number) => void;
  setSweepPrimary: (value: SweepAxis) => void;
  updateSweepPrimary: (partial: Partial<SweepAxis>) => void;
  setSweepSecondary: (value: SweepAxis | null) => void;
  updateSweepSecondary: (partial: Partial<SweepAxis>) => void;
  setPolarXAxis: (value: AxisVariable) => void;
  setPolarYAxis: (value: AxisVariable) => void;
  setPlotDataSource: (value: DataSource) => void;
  setPlotChartType: (value: ChartType) => void;
  setPlotXField: (value: keyof RunRow) => void;
  setPlotYField: (value: keyof RunRow) => void;
  setPlotGroupBy: (value: keyof RunRow | '' | typeof AUTO_GROUP_KEY) => void;
  setPlotXScale: (value: AxisScale) => void;
  setPlotYScale: (value: AxisScale) => void;
  setDataExplorerView: (value: DataExplorerViewMode) => void;
  setDataExplorerSplomKeys: (value: (keyof RunRow)[]) => void;
  setDataExplorerColorBy: (value: keyof RunRow | '') => void;
  setDataExplorerFilterModel: (value: unknown | null) => void;
  setOutlierFilterEnabled: (value: boolean) => void;
  setViewport: (value: RouteViewportState) => void;
  applyRouteViewport: (value: Partial<RouteViewportState>) => void;
  setLayoutJson: (value: IJsonModel | null) => void;
  applyRouteLayout: (value: IJsonModel | null) => void;
  setActivePanel: (value: PanelId) => void;
  applyRouteActivePanel: (value: PanelId) => void;
  hydrateFromRoute: (value: Partial<RouteUiSnapshot> & { viewport?: Partial<RouteViewportState> }) => void;
}

export const DEFAULT_ROUTE_UI_STATE: RouteUiSnapshot = {
  libraryNacaCode: '0012',
  libraryNPoints: 50,
  spacingLiveUpdate: true,
  solveRunMode: 'alpha',
  solveShowAdvanced: false,
  solveTargetCl: 0.5,
  solvePolarStart: -5,
  solvePolarEnd: 15,
  solvePolarStep: 1,
  sweepPrimary: { param: 'alpha', start: -5, end: 15, step: 1 },
  sweepSecondary: null,
  polarXAxis: 'alpha',
  polarYAxis: 'cl',
  plotDataSource: 'full',
  plotChartType: 'scatter',
  plotXField: 'alpha',
  plotYField: 'cl',
  plotGroupBy: AUTO_GROUP_KEY,
  plotXScale: 'linear',
  plotYScale: 'linear',
  dataExplorerView: 'table',
  dataExplorerSplomKeys: ['alpha', 'cl', 'cd', 'cm'],
  dataExplorerColorBy: '',
  dataExplorerFilterModel: null,
  outlierFilterEnabled: false,
  activePanel: 'canvas',
  layoutJson: null,
  viewport: {
    centerX: 0.5,
    centerY: 0,
    zoom: 400,
  },
};

export const useRouteUiStore = create<RouteUiStore>((set) => ({
  ...DEFAULT_ROUTE_UI_STATE,
  viewportRevision: 0,
  layoutRevision: 0,
  activePanelRevision: 0,
  setLibraryNacaCode: (libraryNacaCode) => set({ libraryNacaCode }),
  setLibraryNPoints: (libraryNPoints) => set({ libraryNPoints: Math.max(10, Math.min(500, libraryNPoints)) }),
  setSpacingLiveUpdate: (spacingLiveUpdate) => set({ spacingLiveUpdate }),
  setSolveRunMode: (solveRunMode) => set({ solveRunMode }),
  setSolveShowAdvanced: (solveShowAdvanced) => set({ solveShowAdvanced }),
  setSolveTargetCl: (solveTargetCl) => set({ solveTargetCl }),
  setSolvePolarStart: (solvePolarStart) => set({ solvePolarStart }),
  setSolvePolarEnd: (solvePolarEnd) => set({ solvePolarEnd }),
  setSolvePolarStep: (solvePolarStep) => set({ solvePolarStep: Math.max(0.01, solvePolarStep) }),
  setSweepPrimary: (sweepPrimary) => set({ sweepPrimary }),
  updateSweepPrimary: (partial) => set((s) => ({ sweepPrimary: { ...s.sweepPrimary, ...partial } })),
  setSweepSecondary: (sweepSecondary) => set({ sweepSecondary }),
  updateSweepSecondary: (partial) => set((s) => ({
    sweepSecondary: s.sweepSecondary ? { ...s.sweepSecondary, ...partial } : null,
  })),
  setPolarXAxis: (polarXAxis) => set({ polarXAxis }),
  setPolarYAxis: (polarYAxis) => set({ polarYAxis }),
  setPlotDataSource: (plotDataSource) => set({ plotDataSource }),
  setPlotChartType: (plotChartType) => set({ plotChartType }),
  setPlotXField: (plotXField) => set({ plotXField }),
  setPlotYField: (plotYField) => set({ plotYField }),
  setPlotGroupBy: (plotGroupBy) => set({ plotGroupBy }),
  setPlotXScale: (plotXScale) => set({ plotXScale }),
  setPlotYScale: (plotYScale) => set({ plotYScale }),
  setDataExplorerView: (dataExplorerView) => set({ dataExplorerView }),
  setDataExplorerSplomKeys: (dataExplorerSplomKeys) => set({ dataExplorerSplomKeys }),
  setDataExplorerColorBy: (dataExplorerColorBy) => set({ dataExplorerColorBy }),
  setDataExplorerFilterModel: (dataExplorerFilterModel) => set({ dataExplorerFilterModel }),
  setOutlierFilterEnabled: (outlierFilterEnabled) => set({ outlierFilterEnabled }),
  setViewport: (viewport) => set({ viewport }),
  applyRouteViewport: (value) =>
    set((state) => ({
      viewport: {
        centerX: value.centerX ?? state.viewport.centerX,
        centerY: value.centerY ?? state.viewport.centerY,
        zoom: value.zoom ?? state.viewport.zoom,
      },
      viewportRevision: state.viewportRevision + 1,
    })),
  setLayoutJson: (layoutJson) => set({ layoutJson }),
  applyRouteLayout: (layoutJson) =>
    set((state) => ({
      layoutJson,
      layoutRevision: state.layoutRevision + 1,
    })),
  setActivePanel: (activePanel) => set({ activePanel }),
  applyRouteActivePanel: (activePanel) =>
    set((state) => ({
      activePanel,
      activePanelRevision: state.activePanelRevision + 1,
    })),
  hydrateFromRoute: (value) =>
    set((state) => {
      const definedRouteState = Object.fromEntries(
        Object.entries(value).filter(([key, entryValue]) => key !== 'viewport' && entryValue !== undefined),
      ) as Partial<RouteUiSnapshot>;

      return {
        ...state,
        ...definedRouteState,
        viewport: value.viewport
          ? {
              centerX: value.viewport.centerX ?? state.viewport.centerX,
              centerY: value.viewport.centerY ?? state.viewport.centerY,
              zoom: value.viewport.zoom ?? state.viewport.zoom,
            }
          : state.viewport,
        layoutJson: value.layoutJson !== undefined ? value.layoutJson : state.layoutJson,
        viewportRevision: state.viewportRevision + (value.viewport ? 1 : 0),
        layoutRevision: state.layoutRevision + (value.layoutJson !== undefined ? 1 : 0),
        activePanelRevision: state.activePanelRevision + (value.activePanel ? 1 : 0),
      };
    }),
}));
