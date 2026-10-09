import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

const values = new Map<string, string>();
const gtag = vi.fn();

beforeEach(() => {
  values.clear();
  values.set('ff_cookie_consent', 'granted');
  gtag.mockReset();
  vi.stubEnv('VITE_SOLVE_RUN_ANALYTICS', 'true');
  vi.stubEnv('BASE_URL', '/flexfoil/');
  vi.stubGlobal('window', { gtag, location: { origin: 'https://foil.flexcompute.com' } });
  vi.stubGlobal('document', { referrer: 'https://example.org/private?email=secret#airfoil' });
  vi.stubGlobal('localStorage', {
    getItem: (key: string) => values.get(key) ?? null,
    setItem: (key: string, value: string) => values.set(key, value),
  });
  vi.resetModules();
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.unstubAllEnvs();
  vi.resetModules();
});

describe('consent-gated feature events', () => {
  it('sends coarse parameters with safe page and referrer URLs', async () => {
    const { trackEvent } = await import('./analytics');
    trackEvent('solve_run', { solve_mode: 'polar', solver_mode: 'viscous', n_panels: 160,
      page_location: 'https://private.example/?geometry=secret' });
    expect(gtag).toHaveBeenCalledExactlyOnceWith('event', 'solve_run', {
      solve_mode: 'polar', solver_mode: 'viscous', n_panels: 160,
      page_location: 'https://foil.flexcompute.com/flexfoil/', page_referrer: 'https://example.org',
    });
  });

  it.each([null, 'denied', 'invalid'])('does not send without consent (%s)', async consent => {
    values.delete('ff_cookie_consent');
    if (consent) values.set('ff_cookie_consent', consent);
    (await import('./analytics')).trackEvent('feature_use', { feature: 'export_dat' });
    expect(gtag).not.toHaveBeenCalled();
  });

  it.each([undefined, 'false', '1'])('does not send without build opt-in (%s)', async value => {
    vi.stubEnv('VITE_SOLVE_RUN_ANALYTICS', value);
    const { trackEvent, SOLVE_RUN_ANALYTICS_ENABLED } = await import('./analytics');
    expect(SOLVE_RUN_ANALYTICS_ENABLED).toBe(false);
    trackEvent('feature_use', { feature: 'export_dat' });
    expect(gtag).not.toHaveBeenCalled();
  });

  it('tolerates blocked tags and unavailable storage', async () => {
    const { trackEvent } = await import('./analytics');
    vi.stubGlobal('window', {});
    expect(() => trackEvent('solve_run')).not.toThrow();
    vi.stubGlobal('localStorage', { getItem: () => { throw new Error('blocked'); } });
    expect(() => trackEvent('solve_run')).not.toThrow();
  });

  it('configures safe page metadata only after acceptance and stops after revocation', async () => {
    const { updateGtagConsent, trackEvent } = await import('./analytics');
    updateGtagConsent('granted');
    expect(gtag).toHaveBeenCalledWith('config', 'G-065GK6XBSR', {
      page_location: 'https://foil.flexcompute.com/flexfoil/', page_referrer: 'https://example.org',
      allow_google_signals: false, allow_ad_personalization_signals: false,
    });
    updateGtagConsent('denied');
    gtag.mockClear();
    trackEvent('feature_use', { feature: 'panel_select', panel_id: 'solve' });
    expect(gtag).not.toHaveBeenCalled();
  });
});
