// Explicit build-time opt-in; ordinary solve paths do not construct event parameters.
export const SOLVE_RUN_ANALYTICS_ENABLED = import.meta.env.VITE_SOLVE_RUN_ANALYTICS === 'true';

const CONSENT_KEY = 'ff_cookie_consent';
const GA_ID = 'G-065GK6XBSR';

export type ConsentStatus = 'granted' | 'denied' | null;

export function getStoredConsent(): ConsentStatus {
  try {
    const value = localStorage.getItem(CONSENT_KEY);
    if (value === 'granted' || value === 'denied') return value;
  } catch {
    // localStorage may be unavailable in some contexts
  }
  return null;
}

export function setStoredConsent(status: 'granted' | 'denied') {
  try {
    localStorage.setItem(CONSENT_KEY, status);
  } catch {
    // silently fail
  }
}

/**
 * Update Google Consent Mode and, if granted, begin measurement.
 * Safe to call before gtag is loaded -- commands queue in dataLayer.
 */
export function updateGtagConsent(status: 'granted' | 'denied') {
  setStoredConsent(status);

  window.gtag?.('consent', 'update', {
    analytics_storage: status,
  });

  if (status === 'granted') {
    window.gtag?.('config', GA_ID, {
      page_location: new URL(import.meta.env.BASE_URL, window.location.origin).href,
      page_referrer: document.referrer ? new URL(document.referrer).origin : '',
      allow_google_signals: false,
      allow_ad_personalization_signals: false,
    });
  }
}

/**
 * Boot-time initialization: set consent defaults, then apply any stored choice.
 * Called once from index.html inline script via the dataLayer default command,
 * and again from React to apply stored consent.
 */
export function initAnalyticsConsent() {
  const stored = getStoredConsent();
  if (stored) {
    updateGtagConsent(stored);
  }
}

/**
 * Send coarse feature metadata only for opted-in builds and consenting visitors.
 * Override page URLs so shared airfoil state in queries/fragments is not sent.
 */
export function trackEvent(name: string, params?: Record<string, string | number | boolean>) {
  if (!SOLVE_RUN_ANALYTICS_ENABLED || getStoredConsent() !== 'granted') return;
  window.gtag?.('event', name, {
    ...params,
    page_location: new URL(import.meta.env.BASE_URL, window.location.origin).href,
    page_referrer: document.referrer ? new URL(document.referrer).origin : '',
  });
}

declare global {
  interface Window {
    gtag?: (...args: unknown[]) => void;
    dataLayer?: unknown[];
  }
}
