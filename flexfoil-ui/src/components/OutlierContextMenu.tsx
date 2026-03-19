/**
 * Shared context menu for flagging data points as outliers.
 *
 * Used by PolarPanel, PlotBuilderPanel, and DataExplorerPanel.
 * Positioned absolutely relative to the nearest positioned ancestor.
 */

import { useEffect, useRef } from 'react';
import { useTheme } from '../contexts/ThemeContext';

export interface OutlierMenuExtra {
  label: string;
  detail?: string;
  onClick: () => void;
}

interface Props {
  x: number;
  y: number;
  isFlagged: boolean;
  onToggle: () => void;
  onDismiss: () => void;
  extras?: OutlierMenuExtra[];
}

export function OutlierContextMenu({ x, y, isFlagged, onToggle, onDismiss, extras }: Props) {
  const { isDark } = useTheme();
  const menuRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const dismiss = (e: MouseEvent | KeyboardEvent) => {
      if (e instanceof KeyboardEvent && e.key !== 'Escape') return;
      onDismiss();
    };
    window.addEventListener('click', dismiss);
    window.addEventListener('keydown', dismiss);
    window.addEventListener('scroll', onDismiss, true);
    return () => {
      window.removeEventListener('click', dismiss);
      window.removeEventListener('keydown', dismiss);
      window.removeEventListener('scroll', onDismiss, true);
    };
  }, [onDismiss]);

  const itemStyle: React.CSSProperties = {
    display: 'block',
    width: '100%',
    padding: '7px 14px',
    background: 'transparent',
    border: 'none',
    borderRadius: 0,
    textAlign: 'left',
    fontSize: '12px',
    color: isDark ? '#f0f2f0' : '#111813',
    cursor: 'pointer',
    whiteSpace: 'nowrap',
  };

  const handleHover = (e: React.MouseEvent<HTMLButtonElement>) => {
    e.currentTarget.style.background = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.06)';
  };
  const handleLeave = (e: React.MouseEvent<HTMLButtonElement>) => {
    e.currentTarget.style.background = 'transparent';
  };

  return (
    <div
      ref={menuRef}
      style={{
        position: 'absolute',
        left: x,
        top: y,
        zIndex: 1000,
        background: isDark ? '#1c241f' : '#fff',
        border: `1px solid ${isDark ? '#3a463d' : '#d2dbd0'}`,
        borderRadius: '6px',
        boxShadow: isDark ? '0 4px 16px rgba(0,0,0,0.4)' : '0 4px 16px rgba(0,0,0,0.12)',
        padding: '4px 0',
        minWidth: '180px',
      }}
      onClick={(e) => e.stopPropagation()}
    >
      <button
        onClick={(e) => { e.stopPropagation(); onToggle(); onDismiss(); }}
        style={itemStyle}
        onMouseEnter={handleHover}
        onMouseLeave={handleLeave}
      >
        {isFlagged ? 'Remove Outlier Flag' : 'Flag as Outlier'}
      </button>

      {extras?.map((extra, i) => (
        <button
          key={i}
          onClick={(e) => { e.stopPropagation(); extra.onClick(); onDismiss(); }}
          style={itemStyle}
          onMouseEnter={handleHover}
          onMouseLeave={handleLeave}
        >
          {extra.label}
          {extra.detail && (
            <span style={{ fontSize: '10px', color: 'var(--text-muted)', marginLeft: '8px' }}>
              {extra.detail}
            </span>
          )}
        </button>
      ))}
    </div>
  );
}
