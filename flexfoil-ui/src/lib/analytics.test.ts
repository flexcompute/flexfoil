import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { trackEvent } from './analytics';

// Tests run in the node environment, so stand up the `window` the module targets.
const globals = globalThis as { window?: { gtag?: (...args: unknown[]) => void } };

beforeEach(() => {
  globals.window = {};
});

afterEach(() => {
  delete globals.window;
});

describe('trackEvent', () => {
  it('forwards the event name and params to gtag', () => {
    const gtag = vi.fn();
    globals.window = { gtag };

    trackEvent('solve_run', { solve_mode: 'polar', solver_mode: 'viscous', n_panels: 160 });

    expect(gtag).toHaveBeenCalledWith('event', 'solve_run', {
      solve_mode: 'polar',
      solver_mode: 'viscous',
      n_panels: 160,
    });
  });

  it('does not throw when the tag is blocked', () => {
    expect(() => trackEvent('solve_run')).not.toThrow();
  });
});
