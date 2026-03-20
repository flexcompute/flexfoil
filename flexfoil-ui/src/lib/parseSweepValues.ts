/**
 * Parse a sweep text input into an array of numeric values.
 *
 * Accepted formats (freely combinable, comma-separated):
 *   "start:step:end"    → generates range  e.g. "-5:1:15"
 *   "[v1, v2, v3]"      → explicit values  e.g. "[5e5, 1e6, 3e6]"
 *   "v1, v2, v3"        → bare CSV         e.g. "5e5, 1e6, 3e6"
 *   "-5:1:5, 7, 10, 15" → mixed range + explicit values
 */

function generateRange(start: number, step: number, end: number): number[] {
  if (!Number.isFinite(start) || !Number.isFinite(step) || !Number.isFinite(end)) return [];
  if (step === 0) return [start];

  const direction = end >= start ? 1 : -1;
  const absStep = Math.abs(step) * direction;
  const values: number[] = [];
  for (let v = start; direction > 0 ? v <= end + 1e-9 : v >= end - 1e-9; v += absStep) {
    values.push(Math.round(v * 1e8) / 1e8);
  }
  return values;
}

function isRangeToken(token: string): boolean {
  return token.includes(':') && token.split(':').length >= 2;
}

function parseRangeToken(token: string): number[] {
  const parts = token.split(':').map((s) => s.trim());
  if (parts.length === 3) {
    const [start, step, end] = parts.map(Number);
    return generateRange(start, step, end);
  }
  if (parts.length === 2) {
    const [start, end] = parts.map(Number);
    return generateRange(start, 1, end);
  }
  return [];
}

export function parseSweepValues(text: string): number[] {
  const stripped = text.replace(/[\[\]]/g, '').trim();
  if (!stripped) return [];

  const values: number[] = [];

  // Split on commas, but re-join colon-containing tokens that were broken
  // (ranges like -5:1:15 won't contain commas, so splitting on comma is safe)
  const tokens = stripped.split(',').map((s) => s.trim()).filter(Boolean);

  for (const token of tokens) {
    if (isRangeToken(token)) {
      values.push(...parseRangeToken(token));
    } else {
      const num = Number(token);
      if (Number.isFinite(num)) {
        values.push(num);
      }
    }
  }

  // Deduplicate while preserving order
  const seen = new Set<number>();
  return values.filter((v) => {
    if (seen.has(v)) return false;
    seen.add(v);
    return true;
  });
}

/**
 * Format a values array back into display text.
 * If the array looks like a uniform range, formats as start:step:end.
 * Otherwise, formats as comma-separated values.
 */
export function formatSweepValues(values: number[]): string {
  if (values.length === 0) return '';
  if (values.length === 1) return formatNum(values[0]);

  // Check for uniform spacing
  const step = values[1] - values[0];
  if (values.length >= 3 && Math.abs(step) > 1e-12) {
    let isUniform = true;
    for (let i = 2; i < values.length; i++) {
      if (Math.abs((values[i] - values[i - 1]) - step) > Math.abs(step) * 1e-6) {
        isUniform = false;
        break;
      }
    }
    if (isUniform) {
      return `${formatNum(values[0])}:${formatNum(step)}:${formatNum(values[values.length - 1])}`;
    }
  }

  return values.map(formatNum).join(', ');
}

function formatNum(v: number): string {
  if (v === 0) return '0';
  const abs = Math.abs(v);
  if (abs >= 1e5 || (abs > 0 && abs < 0.01)) return v.toExponential();
  return String(v);
}
