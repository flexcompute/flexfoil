import { describe, it, expect } from 'vitest';
import { parseSweepValues, formatSweepValues } from './parseSweepValues';

describe('parseSweepValues', () => {
  it('parses start:step:end range', () => {
    expect(parseSweepValues('-5:1:5')).toEqual([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]);
  });

  it('parses start:end (step defaults to 1)', () => {
    expect(parseSweepValues('0:3')).toEqual([0, 1, 2, 3]);
  });

  it('parses bracketed explicit values', () => {
    expect(parseSweepValues('[5e5, 1e6, 3e6]')).toEqual([5e5, 1e6, 3e6]);
  });

  it('parses bare CSV (no brackets)', () => {
    expect(parseSweepValues('5e5, 1e6, 3e6')).toEqual([5e5, 1e6, 3e6]);
  });

  it('parses mixed range + explicit', () => {
    expect(parseSweepValues('-5:1:5, 7, 10, 15')).toEqual([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 7, 10, 15]);
  });

  it('deduplicates values', () => {
    expect(parseSweepValues('1, 2, 3, 2, 1')).toEqual([1, 2, 3]);
  });

  it('handles scientific notation', () => {
    expect(parseSweepValues('1e5, 5e5, 1e6')).toEqual([1e5, 5e5, 1e6]);
  });

  it('handles negative values', () => {
    expect(parseSweepValues('-10, -5, 0, 5')).toEqual([-10, -5, 0, 5]);
  });

  it('handles single value', () => {
    expect(parseSweepValues('42')).toEqual([42]);
  });

  it('returns empty for empty string', () => {
    expect(parseSweepValues('')).toEqual([]);
  });

  it('handles descending range', () => {
    expect(parseSweepValues('5:1:-5')).toEqual([5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]);
  });

  it('handles fractional step', () => {
    const result = parseSweepValues('0:0.5:2');
    expect(result).toEqual([0, 0.5, 1, 1.5, 2]);
  });

  it('handles zero step gracefully', () => {
    expect(parseSweepValues('5:0:10')).toEqual([5]);
  });

  it('ignores invalid tokens', () => {
    expect(parseSweepValues('1, abc, 3')).toEqual([1, 3]);
  });

  it('handles whitespace gracefully', () => {
    expect(parseSweepValues('  1 , 2 ,  3  ')).toEqual([1, 2, 3]);
  });
});

describe('formatSweepValues', () => {
  it('formats uniform range as start:step:end', () => {
    expect(formatSweepValues([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])).toBe('-5:1:5');
  });

  it('formats non-uniform as CSV', () => {
    expect(formatSweepValues([5e5, 1e6, 3e6])).toBe('5e+5, 1e+6, 3e+6');
  });

  it('formats single value', () => {
    expect(formatSweepValues([42])).toBe('42');
  });

  it('formats empty array', () => {
    expect(formatSweepValues([])).toBe('');
  });

  it('round-trips a range', () => {
    const text = '-5:1:15';
    const values = parseSweepValues(text);
    expect(formatSweepValues(values)).toBe('-5:1:15');
  });

  it('round-trips explicit values', () => {
    const text = '5e5, 1e6, 3e6';
    const values = parseSweepValues(text);
    const formatted = formatSweepValues(values);
    expect(parseSweepValues(formatted)).toEqual(values);
  });
});
