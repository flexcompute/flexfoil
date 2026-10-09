import { useCallback, useRef, useState } from 'react';
import { trackEvent } from '../lib/analytics';

const REQUEST_TRACKER_URL = 'https://github.com/flexcompute/flexfoil/issues?q=is%3Aissue+label%3Aenhancement';

const GOOGLE_SHEET_URL: string | undefined = import.meta.env.VITE_FEEDBACK_SHEET_URL;

type FeedbackType = 'bug' | 'feature' | 'general';

const TYPE_LABELS: Record<FeedbackType, { label: string; icon: string; template: string }> = {
  bug: { label: 'Bug', icon: '🐛', template: 'bug-report.md' },
  feature: { label: 'Feature', icon: '💡', template: 'feature-request.md' },
  general: { label: 'General', icon: '💬', template: 'question.md' },
};

type SubmitState = 'idle' | 'sending' | 'success' | 'github' | 'error';

export function FeedbackWidget() {
  const [open, setOpen] = useState(false);
  const [type, setType] = useState<FeedbackType>('general');
  const [message, setMessage] = useState('');
  const [contact, setContact] = useState('');
  const [submitState, setSubmitState] = useState<SubmitState>('idle');
  const formRef = useRef<HTMLFormElement>(null);

  const reset = useCallback(() => {
    setType('general');
    setMessage('');
    setContact('');
    setSubmitState('idle');
  }, []);

  const handleClose = useCallback(() => {
    setOpen(false);
    if (submitState === 'success') reset();
    else if (submitState === 'github') setSubmitState('idle');
  }, [submitState, reset]);

  const handleSubmit = useCallback(
    async (e: React.FormEvent) => {
      e.preventDefault();
      if (!message.trim()) return;

      if (!GOOGLE_SHEET_URL) {
        const params = new URLSearchParams({
          title: message.trim().split('\n')[0].slice(0, 80),
          body: message.trim(),
          template: TYPE_LABELS[type].template,
        });
        window.open(`https://github.com/flexcompute/flexfoil/issues/new?${params}`, '_blank', 'noopener,noreferrer');
        trackEvent('feedback_handoff', { feedback_type: type });
        setSubmitState('github');
        return;
      }

      setSubmitState('sending');

      const payload = {
        type,
        message: message.trim(),
        contact: contact.trim() || undefined,
        timestamp: new Date().toISOString(),
        url: window.location.origin + window.location.pathname,
        userAgent: navigator.userAgent,
      };

      try {
        await fetch(GOOGLE_SHEET_URL, {
          method: 'POST',
          // ponytail: opaque responses cannot confirm persistence; use a CORS-aware service for receipts.
          mode: 'no-cors',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(payload),
        });
        setSubmitState('success');
        trackEvent('feedback_sent', { feedback_type: type });
      } catch (err) {
        console.error('Feedback submission failed:', err);
        setSubmitState('error');
      }
    },
    [type, message, contact],
  );

  return (
    <>
      <button
        className="feedback-trigger"
        onClick={() => {
          setOpen((prev) => !prev);
          if (submitState === 'success') reset();
          else if (submitState === 'github') setSubmitState('idle');
          if (!open) trackEvent('feature_use', { feature: 'feedback_open' });
        }}
        aria-label="Send feedback"
        title="Send feedback"
      >
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" aria-hidden>
          <path
            d="M22 2L11 13M22 2l-7 20-4-9-9-4 20-7z"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
        Feedback & requests
      </button>

      {open && (
        <div className="feedback-panel" role="dialog" aria-label="Feedback form">
          {submitState === 'success' || submitState === 'github' ? (
            <div className="feedback-panel__success">
              <span className="feedback-panel__success-icon">✓</span>
              <p className="feedback-panel__success-title">{submitState === 'github' ? 'Finish your submission on GitHub' : 'Feedback sent'}</p>
              <p className="feedback-panel__success-sub">{submitState === 'github' ? 'Submit your draft on GitHub to add it to the tracker. Opening a draft does not submit it.' : 'The request was sent to our feedback service; delivery cannot be confirmed here.'}</p>
              <button className="feedback-panel__btn feedback-panel__btn--primary" onClick={handleClose}>
                Close
              </button>
            </div>
          ) : (
            <form ref={formRef} className="feedback-panel__form" onSubmit={handleSubmit}>
              <div className="feedback-panel__header">
                <span className="feedback-panel__title">Feedback & feature requests</span>
              </div>

              <div className="feedback-panel__type-row">
                {(Object.keys(TYPE_LABELS) as FeedbackType[]).map((t) => (
                  <button
                    key={t}
                    type="button"
                    className={`feedback-panel__type-btn${type === t ? ' feedback-panel__type-btn--active' : ''}`}
                    aria-label={TYPE_LABELS[t].label}
                    aria-pressed={type === t}
                    onClick={() => setType(t)}
                  >
                    <span>{TYPE_LABELS[t].icon}</span>
                    <span>{TYPE_LABELS[t].label}</span>
                  </button>
                ))}
              </div>

              <textarea
                className="feedback-panel__textarea"
                aria-label="Feedback message"
                value={message}
                onChange={(e) => setMessage(e.target.value)}
                placeholder={
                  type === 'bug'
                    ? 'What happened? What did you expect?'
                    : type === 'feature'
                      ? "What would you like to see? What problem does it solve?"
                      : 'Tell us what you think...'
                }
                rows={4}
                required
                autoFocus
              />

              {!GOOGLE_SHEET_URL && (
                <p>Continue on GitHub to submit a public issue. A GitHub account is required; do not include private information.</p>
              )}
              <a href={REQUEST_TRACKER_URL} target="_blank" rel="noopener noreferrer"
                onClick={() => trackEvent('feature_use', { feature: 'request_tracker' })}>View feature request tracker</a>

              {GOOGLE_SHEET_URL && <input
                className="feedback-panel__input"
                type="text"
                value={contact}
                onChange={(e) => setContact(e.target.value)}
                placeholder="Name or email (optional, for follow-up)"
              />}

              {submitState === 'error' && (
                <p className="feedback-panel__error">
                  Something went wrong. Please try again.
                </p>
              )}

              <div className="feedback-panel__actions">
                <button
                  type="button"
                  className="feedback-panel__btn feedback-panel__btn--ghost"
                  onClick={handleClose}
                >
                  Cancel
                </button>
                <button
                  type="submit"
                  className="feedback-panel__btn feedback-panel__btn--primary"
                  disabled={!message.trim() || submitState === 'sending'}
                >
                  {submitState === 'sending' ? 'Sending...' : GOOGLE_SHEET_URL ? 'Send' : 'Continue on GitHub'}
                </button>
              </div>
            </form>
          )}
        </div>
      )}
    </>
  );
}
