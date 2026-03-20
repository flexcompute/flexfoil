/**
 * In-memory store for manually flagged polar-plot outliers.
 *
 * PolarPoints are ephemeral (not persisted in the DB), so their outlier
 * flags live here and reset when the page reloads — matching the lifetime
 * of the polar data itself.
 */

import { create } from 'zustand';

function makeKey(seriesKey: string, pointIndex: number): string {
  return `${seriesKey}:${pointIndex}`;
}

interface PolarOutlierState {
  /** Composite keys of flagged points: "seriesKey:pointIndex" */
  flaggedKeys: Set<string>;
  toggleFlag: (seriesKey: string, pointIndex: number) => void;
  isFlagged: (seriesKey: string, pointIndex: number) => boolean;
  clearAll: () => void;
  /** Revision counter so React selectors re-render on mutation */
  revision: number;
}

export const usePolarOutlierStore = create<PolarOutlierState>()((set, get) => ({
  flaggedKeys: new Set(),
  revision: 0,

  toggleFlag: (seriesKey, pointIndex) => {
    const key = makeKey(seriesKey, pointIndex);
    const next = new Set(get().flaggedKeys);
    if (next.has(key)) next.delete(key);
    else next.add(key);
    set({ flaggedKeys: next, revision: get().revision + 1 });
  },

  isFlagged: (seriesKey, pointIndex) => {
    return get().flaggedKeys.has(makeKey(seriesKey, pointIndex));
  },

  clearAll: () => set({ flaggedKeys: new Set(), revision: get().revision + 1 }),
}));
