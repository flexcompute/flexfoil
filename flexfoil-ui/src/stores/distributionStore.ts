/**
 * Store for surface distribution panel state.
 *
 * Tracks which runs are pinned for overlay and the axis/surface configuration.
 */

import { create } from 'zustand';
import type { DistributionQuantity, SurfaceCoordinate } from '../types';

interface DistributionStore {
  pinnedRunIds: Set<number>;
  distributionXAxis: SurfaceCoordinate;
  distributionYAxis: DistributionQuantity;
  showUpper: boolean;
  showLower: boolean;

  togglePinned: (id: number) => void;
  setPinned: (ids: number[]) => void;
  clearPinned: () => void;
  setDistributionXAxis: (axis: SurfaceCoordinate) => void;
  setDistributionYAxis: (axis: DistributionQuantity) => void;
  setShowUpper: (v: boolean) => void;
  setShowLower: (v: boolean) => void;
}

export const useDistributionStore = create<DistributionStore>((set) => ({
  pinnedRunIds: new Set(),
  distributionXAxis: 'x',
  distributionYAxis: 'cp',
  showUpper: true,
  showLower: true,

  togglePinned: (id) =>
    set((state) => {
      const next = new Set(state.pinnedRunIds);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return { pinnedRunIds: next };
    }),

  setPinned: (ids) => set({ pinnedRunIds: new Set(ids) }),
  clearPinned: () => set({ pinnedRunIds: new Set() }),

  setDistributionXAxis: (distributionXAxis) => set({ distributionXAxis }),
  setDistributionYAxis: (distributionYAxis) => set({ distributionYAxis }),
  setShowUpper: (showUpper) => set({ showUpper }),
  setShowLower: (showLower) => set({ showLower }),
}));
