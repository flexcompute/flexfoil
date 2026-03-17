import { create } from 'zustand';
import type { AxisScale, ChartType, DataSource, RunRow } from '../types';
import { AUTO_GROUP_KEY } from './routeUiStore';

const STORAGE_KEY = 'flexfoil-saved-plots-v1';

export interface SavedPlot {
  id: string;
  title: string;
  chartType: ChartType;
  xField: keyof RunRow;
  yField: keyof RunRow;
  groupBy: keyof RunRow | '' | typeof AUTO_GROUP_KEY;
  colorBy: keyof RunRow | '';
  sizeBy: keyof RunRow | '';
  symbolBy: keyof RunRow | '';
  xScale: AxisScale;
  yScale: AxisScale;
  dataSource: DataSource;
  createdAt: string;
}

interface SavedPlotStore {
  savedPlots: SavedPlot[];
  addPlot: (plot: Omit<SavedPlot, 'id' | 'createdAt'>) => void;
  deletePlot: (id: string) => void;
}

function loadFromStorage(): SavedPlot[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) return JSON.parse(raw);
  } catch {
    // ignore
  }
  return [];
}

function saveToStorage(plots: SavedPlot[]) {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(plots));
  } catch {
    // ignore
  }
}

export const useSavedPlotStore = create<SavedPlotStore>((set) => ({
  savedPlots: loadFromStorage(),

  addPlot: (plot) =>
    set((state) => {
      const newPlot: SavedPlot = {
        ...plot,
        id: `sp_${Date.now()}_${Math.random().toString(36).slice(2, 8)}`,
        createdAt: new Date().toISOString(),
      };
      const next = [newPlot, ...state.savedPlots];
      saveToStorage(next);
      return { savedPlots: next };
    }),

  deletePlot: (id) =>
    set((state) => {
      const next = state.savedPlots.filter((p) => p.id !== id);
      saveToStorage(next);
      return { savedPlots: next };
    }),
}));
