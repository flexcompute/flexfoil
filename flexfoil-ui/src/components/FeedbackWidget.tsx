import { useCallback, useRef, useState } from 'react';

const GOOGLE_SHEET_URL: string | undefined = import.meta.env.VITE_FEEDBACK_SHEET_URL;

type FeedbackType = 'bug' | 'feature' | 'general';

const TYPE_LABELS: Record<FeedbackType, { label: string; icon: string }> = {
  bug: { label: 'Bug', icon: '🐛' },
  feature: { label: 'Feature', icon: '💡' },
  general: { label: 'General', icon: '💬' },
};

type SubmitState = 'idle' | 'sending' | 'success' | 'error';

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
  }, [submitState, reset]);

  const handleSubmit = useCallback(
    async (e: React.FormEvent) => {
      e.preventDefault();
      if (!message.trim()) return;

      setSubmitState('sending');

      const payload = {
        type,
        message: message.trim(),
        contact: contact.trim() || undefined,
        timestamp: new Date().toISOString(),
        url: window.location.href,
        userAgent: navigator.userAgent,
      };

      if (!GOOGLE_SHEET_URL) {
        console.info('[Feedback] No VITE_FEEDBACK_SHEET_URL configured. Payload:', payload);
        setSubmitState('success');
        return;
      }

      try {
        await fetch(GOOGLE_SHEET_URL, {
          method: 'POST',
          mode: 'no-cors',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(payload),
        });
        setSubmitState('success');
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
        Feedback
      </button>

      {open && (
        <div className="feedback-panel" role="dialog" aria-label="Feedback form">
          {submitState === 'success' ? (
            <div className="feedback-panel__success">
              <span className="feedback-panel__success-icon">✓</span>
              <p className="feedback-panel__success-title">Thanks for your feedback!</p>
              <p className="feedback-panel__success-sub">We read every submission.</p>
              <button className="feedback-panel__btn feedback-panel__btn--primary" onClick={handleClose}>
                Close
              </button>
            </div>
          ) : (
            <form ref={formRef} className="feedback-panel__form" onSubmit={handleSubmit}>
              <div className="feedback-panel__header">
                <span className="feedback-panel__title">Send Feedback</span>
              </div>

              <div className="feedback-panel__type-row">
                {(Object.keys(TYPE_LABELS) as FeedbackType[]).map((t) => (
                  <button
                    key={t}
                    type="button"
                    className={`feedback-panel__type-btn${type === t ? ' feedback-panel__type-btn--active' : ''}`}
                    onClick={() => setType(t)}
                  >
                    <span>{TYPE_LABELS[t].icon}</span>
                    <span>{TYPE_LABELS[t].label}</span>
                  </button>
                ))}
              </div>

              <textarea
                className="feedback-panel__textarea"
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

              <input
                className="feedback-panel__input"
                type="text"
                value={contact}
                onChange={(e) => setContact(e.target.value)}
                placeholder="Name or email (optional, for follow-up)"
              />

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
                  {submitState === 'sending' ? 'Sending...' : 'Send'}
                </button>
              </div>
            </form>
          )}
        </div>
      )}
    </>
  );
}
