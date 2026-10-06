import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { trackEvent } from './analytics';

// Tests run in the node environment, so stand up the `window` the module targets.
const globals = globalThis as { window?: { gtag?: (...args: unknown[]) => void } };

beforeEach(() => {
  globals.window = {};
});

afterEach(() => {
  delete globals.window;
  vi.unstubAllEnvs();
  vi.resetModules();
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

describe('solve analytics opt-in', () => {
  it.each([undefined, 'false', '1'])('is off without the exact explicit true value (%s)', async value => {
    vi.stubEnv('VITE_SOLVE_RUN_ANALYTICS', value);
    vi.resetModules();
    expect((await import('./analytics')).SOLVE_RUN_ANALYTICS_ENABLED).toBe(false);
  });
  it('is enabled by the explicit build setting', async () => {
    vi.stubEnv('VITE_SOLVE_RUN_ANALYTICS', 'true');
    vi.resetModules();
    expect((await import('./analytics')).SOLVE_RUN_ANALYTICS_ENABLED).toBe(true);
  });
});
