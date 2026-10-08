import { useCallback, useEffect, useState } from 'react';
import { getStoredConsent, updateGtagConsent } from '../lib/analytics';

export function CookieConsent() {
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const show = () => setVisible(true);
    window.addEventListener('flexfoil:analytics-preferences', show);
    const stored = getStoredConsent();
    const timer = stored === null ? setTimeout(() => setVisible(true), 1200) : undefined;
    return () => {
      window.removeEventListener('flexfoil:analytics-preferences', show);
      clearTimeout(timer);
    };
  }, []);

  const handleAccept = useCallback(() => {
    updateGtagConsent('granted');
    setVisible(false);
  }, []);

  const handleReject = useCallback(() => {
    updateGtagConsent('denied');
    setVisible(false);
  }, []);

  if (!visible) return null;

  return (
    <div className="cookie-banner" role="dialog" aria-label="Cookie consent">
      <div className="cookie-banner__body">
        <p className="cookie-banner__text">
          With your permission, Google Analytics measures feature use and approximate
          location to help us improve Flexfoil. We do not send airfoil geometry,
          feedback text, or contact details as analytics events.{' '}
          <a
            href="https://policies.google.com/technologies/partner-sites"
            target="_blank"
            rel="noopener noreferrer"
            className="cookie-banner__link"
          >
            Learn more
          </a>
        </p>
        <div className="cookie-banner__actions">
          <button className="cookie-banner__btn cookie-banner__btn--reject" onClick={handleReject}>
            Reject
          </button>
          <button className="cookie-banner__btn cookie-banner__btn--accept" onClick={handleAccept}>
            Accept
          </button>
        </div>
      </div>
    </div>
  );
}
