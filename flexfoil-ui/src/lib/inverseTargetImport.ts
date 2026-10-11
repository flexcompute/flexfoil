import type { InverseDesignSurfaceTarget } from '../types';

/** Import one surface, without sorting or changing the prescribed distribution. */
export function parseInverseTarget(text: string, kind: 'cp' | 'velocity'): InverseDesignSurfaceTarget {
  const x: number[] = [];
  const values: number[] = [];
  for (const [index, line] of text.split(/\r?\n/).entries()) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;
    const parts = trimmed.includes(',') ? trimmed.split(',').map((part) => part.trim()) : trimmed.split(/\s+/);
    if (x.length === 0 && /^(x|x\/c)$/i.test(parts[0])) {
      const expected = kind === 'cp' ? 'cp' : 'ue';
      if (parts.length !== 2 || parts[1].toLowerCase() !== expected) {
        throw new Error(`Line ${index + 1}: expected x,${expected} for the selected target type.`);
      }
      continue;
    }
    const [position, value] = parts.map(Number);
    if (parts.length !== 2 || parts.some((part) => part === '') || !Number.isFinite(position) || !Number.isFinite(value)) {
      throw new Error(`Line ${index + 1}: expected two finite numbers (x/c and target value).`);
    }
    if (position < 0 || position > 1 || (x.length > 0 && position <= x[x.length - 1])) {
      throw new Error(`Line ${index + 1}: x/c must increase strictly within 0–1.`);
    }
    x.push(position);
    values.push(value);
  }
  if (x.length < 2) throw new Error('Expected at least two target points.');
  return { x, values };
}
