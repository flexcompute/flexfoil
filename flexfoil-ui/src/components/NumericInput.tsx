/**
 * NumericInput - text-based numeric input that commits on Enter/blur/arrows.
 *
 * Solves three problems with native <input type="number">:
 *  1. Scientific notation (e.g. "6e6") is rejected by some browsers
 *  2. Typing a leading minus sign produces NaN mid-keystroke
 *  3. Immediate onChange triggers re-computation before the user finishes typing
 *
 * The component keeps a local *string* while the user types and only
 * pushes a parsed number outward on commit (Enter, blur, arrow-key step).
 */

import { useState, useEffect, useRef, useCallback, type CSSProperties } from 'react';

interface NumericInputProps {
  value: number;
  onChange: (value: number) => void;
  step?: number;
  min?: number;
  max?: number;
  style?: CSSProperties;
  /** Format function for display when the input is not focused */
  format?: (v: number) => string;
}

function clamp(v: number, min?: number, max?: number): number {
  if (min !== undefined && v < min) return min;
  if (max !== undefined && v > max) return max;
  return v;
}

export function NumericInput({
  value,
  onChange,
  step = 1,
  min,
  max,
  style,
  format,
}: NumericInputProps) {
  const [localText, setLocalText] = useState(() =>
    format ? format(value) : String(value),
  );
  const [isFocused, setIsFocused] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);

  // Sync external value → local text when the input is NOT focused
  useEffect(() => {
    if (!isFocused) {
      setLocalText(format ? format(value) : String(value));
    }
  }, [value, isFocused, format]);

  const commit = useCallback(
    (text: string) => {
      const parsed = Number(text);
      if (!Number.isFinite(parsed)) return; // ignore garbage
      const clamped = clamp(parsed, min, max);
      onChange(clamped);
      setLocalText(format ? format(clamped) : String(clamped));
    },
    [onChange, min, max, format],
  );

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent<HTMLInputElement>) => {
      if (e.key === 'Enter') {
        commit(localText);
        inputRef.current?.blur();
        return;
      }

      if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
        e.preventDefault();
        const base = Number(localText);
        const current = Number.isFinite(base) ? base : value;
        const delta = e.key === 'ArrowUp' ? step : -step;
        const next = clamp(current + delta, min, max);
        onChange(next);
        setLocalText(format ? format(next) : String(next));
      }
    },
    [localText, value, step, min, max, onChange, commit, format],
  );

  return (
    <input
      ref={inputRef}
      type="text"
      inputMode="decimal"
      value={localText}
      onChange={(e) => setLocalText(e.target.value)}
      onFocus={(e) => {
        setIsFocused(true);
        // Show raw number on focus so the user can edit freely
        setLocalText(String(value));
        // Select all for easy replacement
        requestAnimationFrame(() => e.target.select());
      }}
      onBlur={() => {
        setIsFocused(false);
        commit(localText);
      }}
      onKeyDown={handleKeyDown}
      style={style}
    />
  );
}
