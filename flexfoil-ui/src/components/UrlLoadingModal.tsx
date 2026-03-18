import { useEffect, useMemo, useState } from 'react';
import type { RouteStateSnapshot } from '../lib/routeState';

const MODE_LABELS: Record<string, string> = {
  parameters: 'Parameters',
  'camber-spline': 'Camber Spline',
  'thickness-spline': 'Thickness Spline',
  'inverse-design': 'Inverse Design',
  'geometry-design': 'Geometry Design',
};

interface DetailItem {
  label: string;
  value: string;
}

function extractDetails(snapshot: RouteStateSnapshot): DetailItem[] {
  const items: DetailItem[] = [];
  const { airfoil, visualization, ui } = snapshot;

  if (airfoil.exactGeometry?.name) {
    items.push({ label: 'Airfoil', value: airfoil.exactGeometry.name });
  } else if (airfoil.nacaCode) {
    items.push({ label: 'Airfoil', value: `NACA ${airfoil.nacaCode}` });
  }

  if (airfoil.controlMode) {
    items.push({ label: 'Mode', value: MODE_LABELS[airfoil.controlMode] ?? airfoil.controlMode });
  }

  if (airfoil.displayAlpha !== undefined) {
    items.push({ label: 'Angle of attack', value: `${airfoil.displayAlpha}°` });
  }

  if (airfoil.reynolds !== undefined) {
    const re = airfoil.reynolds;
    const formatted = re >= 1e6 ? `${(re / 1e6).toFixed(re % 1e6 === 0 ? 0 : 1)}M` : re.toLocaleString();
    items.push({ label: 'Reynolds', value: formatted });
  }

  if (airfoil.mach !== undefined && airfoil.mach > 0) {
    items.push({ label: 'Mach', value: String(airfoil.mach) });
  }

  if (airfoil.nPanels !== undefined) {
    items.push({ label: 'Panels', value: String(airfoil.nPanels) });
  }

  if (airfoil.spacingKnots && airfoil.spacingKnots.length > 0) {
    items.push({ label: 'Spacing', value: `${airfoil.spacingKnots.length} knots` });
  }

  if (airfoil.geometryDesign?.flaps && airfoil.geometryDesign.flaps.length > 0) {
    const n = airfoil.geometryDesign.flaps.length;
    items.push({ label: 'Flaps', value: `${n} flap${n > 1 ? 's' : ''} defined` });
  }

  const visFlags: string[] = [];
  if (visualization.showStreamlines) visFlags.push('Streamlines');
  if (visualization.showSmoke) visFlags.push('Smoke');
  if (visualization.showCp) visFlags.push('Cp');
  if (visualization.showForces) visFlags.push('Forces');
  if (visualization.showPsiContours) visFlags.push('Psi contours');
  if (visualization.showBoundaryLayer) visFlags.push('Boundary layer');
  if (visualization.showPanels) visFlags.push('Panels');
  if (visualization.showGrid) visFlags.push('Grid');
  if (visFlags.length > 0) {
    items.push({ label: 'Visualizations', value: visFlags.join(', ') });
  }

  if (ui.viewport && (ui.viewport.centerX !== undefined || ui.viewport.zoom !== undefined)) {
    items.push({ label: 'Viewport', value: 'Custom view restored' });
  }

  if (snapshot.layoutJson) {
    items.push({ label: 'Layout', value: 'Custom panel layout' });
  }

  return items;
}

interface UrlLoadingModalProps {
  snapshot: RouteStateSnapshot;
  wasmReady: boolean;
  hydrated: boolean;
}

export function UrlLoadingModal({ snapshot, wasmReady, hydrated }: UrlLoadingModalProps) {
  const details = useMemo(() => extractDetails(snapshot), [snapshot]);
  const [visibleCount, setVisibleCount] = useState(0);
  const [dismissed, setDismissed] = useState(false);

  useEffect(() => {
    if (details.length === 0) return;
    let i = 0;
    const interval = setInterval(() => {
      i++;
      setVisibleCount(i);
      if (i >= details.length) clearInterval(interval);
    }, 90);
    return () => clearInterval(interval);
  }, [details.length]);

  useEffect(() => {
    if (hydrated) {
      const timer = setTimeout(() => setDismissed(true), 600);
      return () => clearTimeout(timer);
    }
  }, [hydrated]);

  if (dismissed || details.length === 0) return null;

  const allRevealed = visibleCount >= details.length;
  const phase = !wasmReady ? 'Initializing solver...' : !hydrated ? 'Applying state...' : 'Ready';

  return (
    <div style={{
      position: 'fixed',
      inset: 0,
      background: 'rgba(0,0,0,0.55)',
      zIndex: 9998,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      animation: hydrated ? 'urlLoadFadeOut 0.5s ease forwards' : undefined,
    }}>
      <div style={{
        background: 'var(--bg-secondary)',
        border: '1px solid var(--border-color)',
        borderRadius: '12px',
        boxShadow: '0 24px 64px rgba(0,0,0,0.4)',
        width: '380px',
        maxWidth: '90vw',
        overflow: 'hidden',
        animation: 'urlLoadSlideIn 0.3s ease-out',
      }}>
        {/* Header */}
        <div style={{
          padding: '20px 24px 12px',
          display: 'flex',
          alignItems: 'center',
          gap: '12px',
        }}>
          <div style={{
            width: '36px',
            height: '36px',
            borderRadius: '10px',
            background: 'var(--accent-primary)',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            flexShrink: 0,
          }}>
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#fff" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
              <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71" />
              <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71" />
            </svg>
          </div>
          <div>
            <div style={{
              fontSize: '14px',
              fontWeight: 700,
              color: 'var(--text-primary)',
              lineHeight: 1.3,
            }}>
              Restoring shared state
            </div>
            <div style={{
              fontSize: '11px',
              color: 'var(--text-muted)',
              fontFamily: 'var(--font-mono)',
              marginTop: '2px',
              display: 'flex',
              alignItems: 'center',
              gap: '6px',
            }}>
              {!hydrated && (
                <span style={{
                  display: 'inline-block',
                  width: '6px',
                  height: '6px',
                  borderRadius: '50%',
                  background: 'var(--accent-warning)',
                  animation: 'urlLoadPulse 1s ease-in-out infinite',
                }} />
              )}
              {phase}
            </div>
          </div>
        </div>

        {/* Detail rows */}
        <div style={{
          padding: '8px 24px 20px',
          display: 'flex',
          flexDirection: 'column',
          gap: '0',
        }}>
          {details.map((item, i) => {
            const visible = i < visibleCount;
            return (
              <div
                key={i}
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  padding: '7px 0',
                  borderBottom: i < details.length - 1 ? '1px solid var(--border-color)' : 'none',
                  opacity: visible ? 1 : 0,
                  transform: visible ? 'translateX(0)' : 'translateX(-8px)',
                  transition: 'opacity 0.25s ease, transform 0.25s ease',
                }}
              >
                <span style={{
                  fontSize: '12px',
                  color: 'var(--text-muted)',
                  fontWeight: 500,
                  flexShrink: 0,
                }}>
                  {item.label}
                </span>
                <span style={{
                  fontSize: '12px',
                  color: 'var(--text-primary)',
                  fontFamily: 'var(--font-mono)',
                  fontWeight: 600,
                  textAlign: 'right',
                  marginLeft: '12px',
                  overflow: 'hidden',
                  textOverflow: 'ellipsis',
                  whiteSpace: 'nowrap',
                }}>
                  {item.value}
                </span>
              </div>
            );
          })}
        </div>

        {/* Progress bar */}
        <div style={{
          height: '3px',
          background: 'var(--bg-tertiary)',
          overflow: 'hidden',
        }}>
          <div style={{
            height: '100%',
            background: hydrated
              ? 'var(--accent-primary)'
              : 'linear-gradient(90deg, var(--accent-primary), var(--accent-secondary))',
            width: hydrated ? '100%' : allRevealed ? '85%' : `${Math.min(90, (visibleCount / details.length) * 70)}%`,
            transition: 'width 0.4s ease',
            ...(hydrated ? {} : {
              animation: allRevealed ? 'urlLoadShimmer 1.5s ease-in-out infinite' : undefined,
            }),
          }} />
        </div>
      </div>

      <style>{`
        @keyframes urlLoadSlideIn {
          from { opacity: 0; transform: translateY(12px) scale(0.97); }
          to   { opacity: 1; transform: translateY(0) scale(1); }
        }
        @keyframes urlLoadFadeOut {
          from { opacity: 1; }
          to   { opacity: 0; pointer-events: none; }
        }
        @keyframes urlLoadPulse {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.3; }
        }
        @keyframes urlLoadShimmer {
          0%   { opacity: 1; }
          50%  { opacity: 0.7; }
          100% { opacity: 1; }
        }
      `}</style>
    </div>
  );
}
