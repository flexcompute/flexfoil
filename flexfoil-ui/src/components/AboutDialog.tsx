import { APP_VERSION } from '../lib/version';

interface AboutDialogProps {
  open: boolean;
  onClose: () => void;
  onOpenChangelog: () => void;
}

export function AboutDialog({ open, onClose, onOpenChangelog }: AboutDialogProps) {
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
          padding: '32px',
          width: '360px',
          maxWidth: '90vw',
          textAlign: 'center',
        }}
      >
        <div
          style={{
            fontSize: '22px',
            fontWeight: 700,
            letterSpacing: '0.06em',
            textTransform: 'uppercase',
            color: 'var(--brand-primary)',
            marginBottom: '4px',
          }}
        >
          FlexFoil
        </div>
        <div style={{ fontSize: '11px', color: 'var(--text-muted)', marginBottom: '16px' }}>
          by Flexcompute
        </div>

        <div
          style={{
            display: 'inline-block',
            padding: '4px 12px',
            borderRadius: '20px',
            background: 'var(--bg-tertiary)',
            border: '1px solid var(--border-color)',
            fontSize: '13px',
            fontFamily: 'var(--font-mono)',
            fontWeight: 600,
            color: 'var(--text-primary)',
            marginBottom: '16px',
          }}
        >
          v{APP_VERSION}
        </div>

        <div style={{ fontSize: '12px', color: 'var(--text-secondary)', lineHeight: 1.6, marginBottom: '20px' }}>
          Real-time 2D airfoil analysis powered by a Rust/WASM panel-method and boundary-layer solver.
        </div>

        <div style={{ display: 'flex', gap: '8px', justifyContent: 'center' }}>
          <button
            onClick={() => {
              onClose();
              onOpenChangelog();
            }}
            style={{
              padding: '6px 16px',
              fontSize: '12px',
              fontWeight: 600,
              background: 'transparent',
              border: '1px solid var(--border-color)',
              borderRadius: '6px',
              color: 'var(--text-secondary)',
              cursor: 'pointer',
            }}
          >
            What's New
          </button>
          <button
            onClick={onClose}
            style={{
              padding: '6px 16px',
              fontSize: '12px',
              fontWeight: 600,
              background: 'var(--accent-primary)',
              border: '1px solid var(--accent-primary)',
              borderRadius: '6px',
              color: 'var(--bg-primary)',
              cursor: 'pointer',
            }}
          >
            Close
          </button>
        </div>
      </div>
    </div>
  );
}
