import { CHANGELOG, type ChangeCategory } from '../lib/version';

interface ChangelogDialogProps {
  open: boolean;
  onClose: () => void;
}

const CATEGORY_STYLES: Record<ChangeCategory, { label: string; color: string; bg: string }> = {
  added: { label: 'Added', color: '#22c55e', bg: 'rgba(34,197,94,0.12)' },
  changed: { label: 'Changed', color: '#3b82f6', bg: 'rgba(59,130,246,0.12)' },
  fixed: { label: 'Fixed', color: '#f59e0b', bg: 'rgba(245,158,11,0.12)' },
};

export function ChangelogDialog({ open, onClose }: ChangelogDialogProps) {
  if (!open) return null;

  return (
    <div
      onClick={onClose}
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
        }}
      >
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
          <div style={{ fontSize: '16px', fontWeight: 700, color: 'var(--text-primary)' }}>
            What's New
          </div>
          <button
            onClick={onClose}
            style={{
              background: 'none',
              border: 'none',
              fontSize: '16px',
              color: 'var(--text-muted)',
              cursor: 'pointer',
              padding: '4px 8px',
              borderRadius: '4px',
            }}
          >
            &times;
          </button>
        </div>

        {/* Scrollable content */}
        <div style={{ flex: 1, overflow: 'auto', padding: '16px 24px 24px' }}>
          {CHANGELOG.map((entry, idx) => (
            <div key={entry.version} style={{ marginBottom: idx < CHANGELOG.length - 1 ? '24px' : 0 }}>
              {/* Version header */}
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
                <span style={{ fontSize: '11px', color: 'var(--text-muted)' }}>
                  {entry.date}
                </span>
              </div>

              {/* Items */}
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
      </div>
    </div>
  );
}
