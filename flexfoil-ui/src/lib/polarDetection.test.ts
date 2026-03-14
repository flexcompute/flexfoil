import { describe, it, expect } from 'vitest';
import {
  detectPolarGroups,
  splitByAlphaGap,
  formatReynolds,
} from './polarDetection';
import type { RunRow } from '../types';

/** Minimal RunRow factory — only fields used by polar detection. */
function makeRow(overrides: Partial<RunRow> = {}): RunRow {
  return {
    id: 1,
    airfoil_name: 'naca0012',
    airfoil_hash: 'abc123',
    alpha: 0,
    reynolds: 1e6,
    mach: 0,
    ncrit: 9,
    n_panels: 160,
    max_iter: 100,
    cl: 0.5,
    cd: 0.01,
    cm: -0.02,
    converged: true,
    iterations: 30,
    residual: 1e-6,
    x_tr_upper: 0.3,
    x_tr_lower: 0.6,
    success: true,
    error: null,
    created_at: '2025-01-01T00:00:00',
    session_id: 's_abc',
    ...overrides,
  };
}

// ─── formatReynolds ─────────────────────────────────────────────────

describe('formatReynolds', () => {
  it('formats millions exactly', () => {
    expect(formatReynolds(1e6)).toBe('1M');
    expect(formatReynolds(3e6)).toBe('3M');
  });

  it('formats fractional millions', () => {
    expect(formatReynolds(1.5e6)).toBe('1.5M');
  });

  it('formats thousands exactly', () => {
    expect(formatReynolds(500_000)).toBe('500k');
    expect(formatReynolds(100_000)).toBe('100k');
  });

  it('formats fractional thousands', () => {
    expect(formatReynolds(250_500)).toBe('250.5k');
  });

  it('formats small values as plain numbers', () => {
    expect(formatReynolds(500)).toBe('500');
  });
});

// ─── splitByAlphaGap ────────────────────────────────────────────────

describe('splitByAlphaGap', () => {
  it('returns a single segment for fewer than 3 rows', () => {
    const rows = [makeRow({ alpha: 0 }), makeRow({ alpha: 1 })];
    const result = splitByAlphaGap(rows);
    expect(result).toHaveLength(1);
    expect(result[0]).toHaveLength(2);
  });

  it('returns a single segment for evenly spaced alphas', () => {
    const rows = Array.from({ length: 10 }, (_, i) => makeRow({ alpha: i }));
    const result = splitByAlphaGap(rows);
    expect(result).toHaveLength(1);
  });

  it('splits on a large gap', () => {
    const rows = [
      ...Array.from({ length: 5 }, (_, i) => makeRow({ alpha: i })),
      ...Array.from({ length: 5 }, (_, i) => makeRow({ alpha: 20 + i })),
    ];
    const result = splitByAlphaGap(rows);
    expect(result).toHaveLength(2);
    expect(result[0]).toHaveLength(5);
    expect(result[1]).toHaveLength(5);
  });

  it('sorts by alpha regardless of input order', () => {
    const rows = [makeRow({ alpha: 5 }), makeRow({ alpha: 1 }), makeRow({ alpha: 3 })];
    const result = splitByAlphaGap(rows);
    expect(result[0].map(r => r.alpha)).toEqual([1, 3, 5]);
  });
});

// ─── detectPolarGroups ──────────────────────────────────────────────

describe('detectPolarGroups', () => {
  it('returns empty array for empty input', () => {
    expect(detectPolarGroups([])).toEqual([]);
  });

  it('groups a single polar into one group', () => {
    const rows = Array.from({ length: 5 }, (_, i) =>
      makeRow({ id: i, alpha: -2 + i }),
    );
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(1);
    expect(groups[0].rows).toHaveLength(5);
    expect(groups[0].label).toBe('naca0012');
    expect(groups[0].isSinglePoint).toBe(false);
  });

  it('separates polars by Reynolds number', () => {
    const rows = [
      ...Array.from({ length: 3 }, (_, i) =>
        makeRow({ id: i, alpha: i, reynolds: 1e6 }),
      ),
      ...Array.from({ length: 3 }, (_, i) =>
        makeRow({ id: 10 + i, alpha: i, reynolds: 2e6 }),
      ),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    expect(groups.map(g => g.label).sort()).toEqual(['Re=1M', 'Re=2M']);
  });

  it('separates polars by airfoil', () => {
    const rows = [
      ...Array.from({ length: 3 }, (_, i) =>
        makeRow({ id: i, alpha: i, airfoil_name: 'naca0012', airfoil_hash: 'aaa' }),
      ),
      ...Array.from({ length: 3 }, (_, i) =>
        makeRow({ id: 10 + i, alpha: i, airfoil_name: 'naca2412', airfoil_hash: 'bbb' }),
      ),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    expect(groups.map(g => g.label).sort()).toEqual(['naca0012', 'naca2412']);
  });

  it('shows both airfoil and Re when both vary', () => {
    const rows = [
      makeRow({ id: 1, alpha: 0, airfoil_name: 'naca0012', airfoil_hash: 'aaa', reynolds: 1e6 }),
      makeRow({ id: 2, alpha: 1, airfoil_name: 'naca0012', airfoil_hash: 'aaa', reynolds: 1e6 }),
      makeRow({ id: 3, alpha: 0, airfoil_name: 'naca2412', airfoil_hash: 'bbb', reynolds: 2e6 }),
      makeRow({ id: 4, alpha: 1, airfoil_name: 'naca2412', airfoil_hash: 'bbb', reynolds: 2e6 }),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    const labels = groups.map(g => g.label).sort();
    expect(labels[0]).toContain('naca0012');
    expect(labels[0]).toContain('Re=1M');
    expect(labels[1]).toContain('naca2412');
    expect(labels[1]).toContain('Re=2M');
  });

  it('marks single-point groups', () => {
    const rows = [makeRow({ id: 1, alpha: 5 })];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(1);
    expect(groups[0].isSinglePoint).toBe(true);
  });

  it('splits on alpha gaps within one polar', () => {
    const rows = [
      ...Array.from({ length: 5 }, (_, i) => makeRow({ id: i, alpha: i })),
      ...Array.from({ length: 5 }, (_, i) => makeRow({ id: 10 + i, alpha: 20 + i })),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    expect(groups[0].rows.map(r => r.alpha)).toEqual([0, 1, 2, 3, 4]);
    expect(groups[1].rows.map(r => r.alpha)).toEqual([20, 21, 22, 23, 24]);
  });

  it('appends sub-index to labels for gap-split groups', () => {
    const rows = [
      ...Array.from({ length: 5 }, (_, i) => makeRow({ id: i, alpha: i })),
      ...Array.from({ length: 5 }, (_, i) => makeRow({ id: 10 + i, alpha: 20 + i })),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups[0].label).toContain('(1)');
    expect(groups[1].label).toContain('(2)');
  });

  it('sorts rows within each group by the sort field', () => {
    const rows = [
      makeRow({ id: 1, alpha: 5 }),
      makeRow({ id: 2, alpha: -2 }),
      makeRow({ id: 3, alpha: 10 }),
      makeRow({ id: 4, alpha: 0 }),
    ];
    const groups = detectPolarGroups(rows, 'alpha');
    expect(groups[0].rows.map(r => r.alpha)).toEqual([-2, 0, 5, 10]);
  });

  it('separates by Mach number', () => {
    const rows = [
      makeRow({ id: 1, alpha: 0, mach: 0 }),
      makeRow({ id: 2, alpha: 1, mach: 0 }),
      makeRow({ id: 3, alpha: 0, mach: 0.3 }),
      makeRow({ id: 4, alpha: 1, mach: 0.3 }),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    expect(groups.map(g => g.label).sort()).toEqual(['M=0', 'M=0.3']);
  });

  it('separates by Ncrit', () => {
    const rows = [
      makeRow({ id: 1, alpha: 0, ncrit: 9 }),
      makeRow({ id: 2, alpha: 1, ncrit: 9 }),
      makeRow({ id: 3, alpha: 0, ncrit: 4 }),
      makeRow({ id: 4, alpha: 1, ncrit: 4 }),
    ];
    const groups = detectPolarGroups(rows);
    expect(groups).toHaveLength(2);
    expect(groups.map(g => g.label).sort()).toEqual(['Nc=4', 'Nc=9']);
  });
});
