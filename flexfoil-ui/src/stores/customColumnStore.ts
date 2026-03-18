import { create } from 'zustand';
import { compileExpression, validateExpression } from '../lib/expressionEngine';
import type { RunRow } from '../types';

const STORAGE_KEY = 'flexfoil-custom-columns-v1';

export interface CustomColumn {
  id: string;
  name: string;
  expression: string;
}

interface CustomColumnStore {
  columns: CustomColumn[];
  addColumn: (name: string, expression: string) => string | null;
  removeColumn: (id: string) => void;
  updateColumn: (id: string, name: string, expression: string) => string | null;
}

function loadFromStorage(): CustomColumn[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) return JSON.parse(raw);
  } catch {
    // ignore
  }
  return [];
}

function saveToStorage(columns: CustomColumn[]) {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(columns));
  } catch {
    // ignore
  }
}

export const useCustomColumnStore = create<CustomColumnStore>((set) => ({
  columns: loadFromStorage(),

  addColumn: (name, expression) => {
    const validation = validateExpression(expression);
    if (!validation.valid) return validation.error ?? 'Invalid expression';
    const col: CustomColumn = {
      id: `cc_${Date.now()}_${Math.random().toString(36).slice(2, 8)}`,
      name: name.trim(),
      expression: expression.trim(),
    };
    set((state) => {
      const next = [...state.columns, col];
      saveToStorage(next);
      return { columns: next };
    });
    return null;
  },

  removeColumn: (id) =>
    set((state) => {
      const next = state.columns.filter((c) => c.id !== id);
      saveToStorage(next);
      return { columns: next };
    }),

  updateColumn: (id, name, expression) => {
    const validation = validateExpression(expression);
    if (!validation.valid) return validation.error ?? 'Invalid expression';
    set((state) => {
      const next = state.columns.map((c) =>
        c.id === id ? { ...c, name: name.trim(), expression: expression.trim() } : c,
      );
      saveToStorage(next);
      return { columns: next };
    });
    return null;
  },
}));

/**
 * Build a compiled evaluator map for all custom columns.
 * Returns a map of columnId -> evaluator function.
 * Columns that fail to compile are silently skipped.
 */
export function compileCustomColumns(
  columns: CustomColumn[],
): Map<string, (row: RunRow) => number | null> {
  const map = new Map<string, (row: RunRow) => number | null>();
  for (const col of columns) {
    try {
      map.set(col.id, compileExpression(col.expression));
    } catch {
      map.set(col.id, () => null);
    }
  }
  return map;
}
