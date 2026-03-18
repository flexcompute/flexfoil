import { useState, useEffect, useCallback } from 'react';
import { CHANGELOG, type ChangeCategory, type TourSlide } from '../lib/version';

const CHANGELOG_SEEN_KEY = 'flexfoil-changelog-seen';

export function getLastSeenChangelogVersion(): string | null {
  return localStorage.getItem(CHANGELOG_SEEN_KEY);
}

export function setLastSeenChangelogVersion(version: string): void {
  localStorage.setItem(CHANGELOG_SEEN_KEY, version);
}

interface ChangelogDialogProps {
  open: boolean;
  onClose: () => void;
  onNavigateToPanel?: (panelId: string) => void;
}

const CATEGORY_STYLES: Record<ChangeCategory, { label: string; color: string; bg: string }> = {
  added: { label: 'Added', color: '#22c55e', bg: 'rgba(34,197,94,0.12)' },
  changed: { label: 'Changed', color: '#3b82f6', bg: 'rgba(59,130,246,0.12)' },
  fixed: { label: 'Fixed', color: '#f59e0b', bg: 'rgba(245,158,11,0.12)' },
};

export function ChangelogDialog({ open, onClose, onNavigateToPanel }: ChangelogDialogProps) {
  const [mode, setMode] = useState<'tour' | 'list'>('tour');
  const [slideIndex, setSlideIndex] = useState(0);

  const latestEntry = CHANGELOG[0];
  const slides = latestEntry?.tourSlides;
  const hasTour = slides && slides.length > 0;

  useEffect(() => {
    if (open) {
      setMode(hasTour ? 'tour' : 'list');
      setSlideIndex(0);
    }
  }, [open, hasTour]);

  const handleClose = useCallback(() => {
    if (latestEntry) setLastSeenChangelogVersion(latestEntry.version);
    onClose();
  }, [latestEntry, onClose]);

  const handleGoTo = useCallback((panelId: string) => {
    if (latestEntry) setLastSeenChangelogVersion(latestEntry.version);
    onClose();
    onNavigateToPanel?.(panelId);
  }, [latestEntry, onClose, onNavigateToPanel]);

  useEffect(() => {
    if (!open) return;
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        handleClose();
      } else if (mode === 'tour' && slides) {
        if (e.key === 'ArrowRight' || e.key === 'ArrowDown') {
          e.preventDefault();
          setSlideIndex((i) => Math.min(i + 1, slides.length - 1));
        } else if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') {
          e.preventDefault();
          setSlideIndex((i) => Math.max(i - 1, 0));
        }
      }
    };
    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [open, mode, slides, handleClose]);

  if (!open) return null;

  return (
    <div
      onClick={handleClose}
      style={{
        position: 'fixed',
        inset: 0,
        background: 'rgba(0,0,0,0.45)',
        zIndex: 9999,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          background: 'var(--bg-secondary)',
          border: '1px solid var(--border-color)',
          borderRadius: '12px',
          boxShadow: '0 24px 64px rgba(0,0,0,0.3)',
          width: '520px',
          maxWidth: '90vw',
          maxHeight: '80vh',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
        }}
      >
        {mode === 'tour' && hasTour ? (
          <TourView
            slides={slides!}
            version={latestEntry.version}
            slideIndex={slideIndex}
            onSlideChange={setSlideIndex}
            onClose={handleClose}
            onViewAll={() => setMode('list')}
            onGoTo={onNavigateToPanel ? handleGoTo : undefined}
          />
        ) : (
          <ListView onClose={handleClose} onViewTour={hasTour ? () => setMode('tour') : undefined} />
        )}
      </div>
    </div>
  );
}

/* ---------- Tour View ---------- */

function TourView({
  slides,
  version,
  slideIndex,
  onSlideChange,
  onClose,
  onViewAll,
  onGoTo,
}: {
  slides: TourSlide[];
  version: string;
  slideIndex: number;
  onSlideChange: (i: number) => void;
  onClose: () => void;
  onViewAll: () => void;
  onGoTo?: (panelId: string) => void;
}) {
  const slide = slides[slideIndex];
  const isLast = slideIndex === slides.length - 1;
  const isFirst = slideIndex === 0;

  return (
    <>
      {/* Header */}
      <div
        style={{
          padding: '20px 24px 0',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          flexShrink: 0,
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
          <span style={{ fontSize: '14px', fontWeight: 700, color: 'var(--text-primary)' }}>
            What's New
          </span>
          <span
            style={{
              padding: '2px 8px',
              borderRadius: '10px',
              background: 'var(--accent-primary)',
              color: 'var(--bg-primary)',
              fontSize: '10px',
              fontWeight: 700,
              fontFamily: 'var(--font-mono)',
            }}
          >
            v{version}
          </span>
        </div>
        <button onClick={onClose} style={closeButtonStyle}>
          &times;
        </button>
      </div>

      {/* Slide content */}
      <div style={{ flex: 1, padding: '24px', overflow: 'hidden' }}>
        <div
          key={slideIndex}
          style={{
            animation: 'changelogFadeIn 0.2s ease-out',
          }}
        >
          {/* Icon + title */}
          <div style={{ display: 'flex', alignItems: 'center', gap: '14px', marginBottom: '8px' }}>
            <div
              style={{
                width: '44px',
                height: '44px',
                borderRadius: '12px',
                background: 'var(--accent-primary)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontSize: '22px',
                flexShrink: 0,
              }}
            >
              {slide.icon}
            </div>
            <div>
              <div
                style={{
                  fontSize: '16px',
                  fontWeight: 700,
                  color: 'var(--text-primary)',
                  lineHeight: 1.3,
                }}
              >
                {slide.title}
              </div>
              <div
                style={{
                  fontSize: '12px',
                  color: 'var(--text-muted)',
                  marginTop: '2px',
                }}
              >
                {slide.description}
              </div>
            </div>
          </div>

          {/* Items */}
          <div
            style={{
              marginTop: '16px',
              display: 'flex',
              flexDirection: 'column',
              gap: '10px',
            }}
          >
            {slide.items.map((item, i) => (
              <div
                key={i}
                style={{
                  display: 'flex',
                  alignItems: 'flex-start',
                  gap: '10px',
                  fontSize: '13px',
                  lineHeight: 1.5,
                  color: 'var(--text-secondary)',
                }}
              >
                <span
                  style={{
                    flexShrink: 0,
                    width: '6px',
                    height: '6px',
                    borderRadius: '50%',
                    background: 'var(--accent-primary)',
                    marginTop: '6px',
                  }}
                />
                {item}
              </div>
            ))}
          </div>

          {/* "Go to" action */}
          {onGoTo && slide.goTo && (
            <button
              onClick={() => onGoTo(slide.goTo!.panel)}
              style={{
                marginTop: '16px',
                display: 'inline-flex',
                alignItems: 'center',
                gap: '6px',
                background: 'var(--bg-tertiary)',
                border: '1px solid var(--border-color)',
                borderRadius: '6px',
                padding: '6px 12px',
                fontSize: '12px',
                fontWeight: 600,
                color: 'var(--accent-primary)',
                cursor: 'pointer',
              }}
            >
              {slide.goTo.label} →
            </button>
          )}
        </div>
      </div>

      {/* Footer: dots + nav */}
      <div
        style={{
          padding: '0 24px 20px',
          display: 'flex',
          flexDirection: 'column',
          gap: '14px',
          flexShrink: 0,
        }}
      >
        {/* Navigation row */}
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <button
            onClick={() => onSlideChange(slideIndex - 1)}
            disabled={isFirst}
            style={{
              ...navButtonStyle,
              opacity: isFirst ? 0.3 : 1,
              cursor: isFirst ? 'default' : 'pointer',
            }}
          >
            ← Back
          </button>

          {/* Dot indicators */}
          <div style={{ display: 'flex', gap: '6px', alignItems: 'center' }}>
            {slides.map((_, i) => (
              <button
                key={i}
                onClick={() => onSlideChange(i)}
                style={{
                  width: i === slideIndex ? '18px' : '7px',
                  height: '7px',
                  borderRadius: '4px',
                  border: 'none',
                  background:
                    i === slideIndex ? 'var(--accent-primary)' : 'var(--text-muted)',
                  opacity: i === slideIndex ? 1 : 0.3,
                  cursor: 'pointer',
                  padding: 0,
                  transition: 'all 0.2s ease',
                }}
                aria-label={`Go to slide ${i + 1}`}
              />
            ))}
          </div>

          {isLast ? (
            <button onClick={onClose} style={{ ...navButtonStyle, ...primaryButtonStyle }}>
              Got it
            </button>
          ) : (
            <button onClick={() => onSlideChange(slideIndex + 1)} style={navButtonStyle}>
              Next →
            </button>
          )}
        </div>

        {/* View all link */}
        <div style={{ textAlign: 'center' }}>
          <button onClick={onViewAll} style={linkButtonStyle}>
            View full changelog
          </button>
        </div>
      </div>

      {/* Inline keyframes */}
      <style>{`
        @keyframes changelogFadeIn {
          from { opacity: 0; transform: translateY(6px); }
          to   { opacity: 1; transform: translateY(0); }
        }
      `}</style>
    </>
  );
}

/* ---------- List View ---------- */

function ListView({
  onClose,
  onViewTour,
}: {
  onClose: () => void;
  onViewTour?: () => void;
}) {
  return (
    <>
      {/* Header */}
      <div
        style={{
          padding: '20px 24px 12px',
          borderBottom: '1px solid var(--border-color)',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          flexShrink: 0,
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <span style={{ fontSize: '16px', fontWeight: 700, color: 'var(--text-primary)' }}>
            Changelog
          </span>
          {onViewTour && (
            <button onClick={onViewTour} style={linkButtonStyle}>
              ← Feature tour
            </button>
          )}
        </div>
        <button onClick={onClose} style={closeButtonStyle}>
          &times;
        </button>
      </div>

      {/* Scrollable content */}
      <div style={{ flex: 1, overflow: 'auto', padding: '16px 24px 24px' }}>
        {CHANGELOG.map((entry, idx) => (
          <div key={entry.version} style={{ marginBottom: idx < CHANGELOG.length - 1 ? '24px' : 0 }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '10px' }}>
              <span
                style={{
                  display: 'inline-block',
                  padding: '2px 10px',
                  borderRadius: '12px',
                  background: idx === 0 ? 'var(--accent-primary)' : 'var(--bg-tertiary)',
                  color: idx === 0 ? 'var(--bg-primary)' : 'var(--text-secondary)',
                  fontSize: '12px',
                  fontWeight: 700,
                  fontFamily: 'var(--font-mono)',
                }}
              >
                v{entry.version}
              </span>
              <span style={{ fontSize: '11px', color: 'var(--text-muted)' }}>{entry.date}</span>
            </div>

            <div style={{ display: 'flex', flexDirection: 'column', gap: '6px' }}>
              {entry.items.map((item, i) => {
                const style = CATEGORY_STYLES[item.category];
                return (
                  <div
                    key={i}
                    style={{
                      display: 'flex',
                      alignItems: 'flex-start',
                      gap: '8px',
                      fontSize: '12px',
                      lineHeight: 1.5,
                    }}
                  >
                    <span
                      style={{
                        flexShrink: 0,
                        padding: '1px 6px',
                        borderRadius: '4px',
                        fontSize: '10px',
                        fontWeight: 600,
                        background: style.bg,
                        color: style.color,
                        marginTop: '1px',
                      }}
                    >
                      {style.label}
                    </span>
                    <span style={{ color: 'var(--text-secondary)' }}>{item.text}</span>
                  </div>
                );
              })}
            </div>
          </div>
        ))}
      </div>
    </>
  );
}

/* ---------- Shared styles ---------- */

const closeButtonStyle: React.CSSProperties = {
  background: 'none',
  border: 'none',
  fontSize: '18px',
  color: 'var(--text-muted)',
  cursor: 'pointer',
  padding: '4px 8px',
  borderRadius: '4px',
  lineHeight: 1,
};

const navButtonStyle: React.CSSProperties = {
  background: 'none',
  border: '1px solid var(--border-color)',
  borderRadius: '6px',
  padding: '6px 14px',
  fontSize: '12px',
  fontWeight: 600,
  color: 'var(--text-secondary)',
  cursor: 'pointer',
  minWidth: '72px',
  textAlign: 'center',
};

const primaryButtonStyle: React.CSSProperties = {
  background: 'var(--accent-primary)',
  color: 'var(--bg-primary)',
  border: '1px solid var(--accent-primary)',
};

const linkButtonStyle: React.CSSProperties = {
  background: 'none',
  border: 'none',
  color: 'var(--text-muted)',
  fontSize: '11px',
  cursor: 'pointer',
  padding: 0,
  textDecoration: 'underline',
  textUnderlineOffset: '2px',
};
