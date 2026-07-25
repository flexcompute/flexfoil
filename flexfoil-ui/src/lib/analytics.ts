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
    window.gtag?.('config', GA_ID);
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
 * Send a custom GA4 event. Consent Mode gates delivery, so callers do not need
 * to check consent first. gtag is absent when the tag is blocked, hence optional.
 */
export function trackEvent(name: string, params?: Record<string, string | number | boolean>) {
  window.gtag?.('event', name, params);
}

declare global {
  interface Window {
    gtag?: (...args: unknown[]) => void;
    dataLayer?: unknown[];
  }
}
