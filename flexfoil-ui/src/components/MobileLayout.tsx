import { useState } from 'react';
import { AirfoilCanvas } from './AirfoilCanvas';
import { SolvePanel } from './panels/SolvePanel';
import { AirfoilLibraryPanel } from './panels/AirfoilLibraryPanel';
import { PolarPanel } from './panels/PolarPanel';
import { DarkModeToggle } from './DarkModeToggle';
import flexcomputeLogo from '../assets/flexcompute-logo.svg';

type MobileTab = 'canvas' | 'solve' | 'library' | 'polar';

const TABS: { id: MobileTab; label: string; shortLabel: string }[] = [
  { id: 'canvas', label: 'Canvas', shortLabel: 'Canvas' },
  { id: 'solve', label: 'Solve', shortLabel: 'Solve' },
  { id: 'library', label: 'Library', shortLabel: 'Library' },
  { id: 'polar', label: 'Polar', shortLabel: 'Polar' },
];

interface MobileLayoutProps {
  wasmStatus: 'loading' | 'ready' | 'error';
}

export function MobileLayout({ wasmStatus }: MobileLayoutProps) {
  const [activeTab, setActiveTab] = useState<MobileTab>('canvas');
  const [bannerDismissed, setBannerDismissed] = useState(false);

  return (
    <div className="mobile-layout">
      {/* Compact header */}
      <header className="mobile-header">
        <a
          href="https://www.flexcompute.com/"
          target="_blank"
          rel="noopener noreferrer"
          className="mobile-header__brand"
        >
          <img src={flexcomputeLogo} alt="Flexcompute" className="mobile-header__logo" />
          <span className="mobile-header__title">FlexFoil</span>
        </a>
        <div className="mobile-header__right">
          <span className={`mobile-header__status mobile-header__status--${wasmStatus}`}>
            {wasmStatus === 'loading' && '…'}
            {wasmStatus === 'ready' && '●'}
            {wasmStatus === 'error' && '✕'}
          </span>
          <DarkModeToggle />
        </div>
      </header>

      {/* Mobile info banner */}
      {!bannerDismissed && (
        <div className="mobile-banner">
          <p className="mobile-banner__text">
            Solver works really well on mobile, but UI on desktop.
          </p>
          <button
            className="mobile-banner__dismiss"
            onClick={() => setBannerDismissed(true)}
            aria-label="Dismiss"
          >
            ✕
          </button>
        </div>
      )}

      {/* Tab bar */}
      <nav className="mobile-tabs">
        {TABS.map((tab) => (
          <button
            key={tab.id}
            className={`mobile-tabs__btn${activeTab === tab.id ? ' mobile-tabs__btn--active' : ''}`}
            onClick={() => setActiveTab(tab.id)}
          >
            {tab.shortLabel}
          </button>
        ))}
      </nav>

      {/* Content area - only render active panel */}
      <div className="mobile-content">
        {activeTab === 'canvas' && (
          <div className="mobile-card mobile-card--canvas">
            <AirfoilCanvas />
          </div>
        )}
        {activeTab === 'solve' && (
          <div className="mobile-card">
            <SolvePanel />
          </div>
        )}
        {activeTab === 'library' && (
          <div className="mobile-card">
            <AirfoilLibraryPanel />
          </div>
        )}
        {activeTab === 'polar' && (
          <div className="mobile-card">
            <PolarPanel />
          </div>
        )}
      </div>

      {/* Compact footer */}
      <footer className="mobile-footer">
        <span>Powered by Flexcompute Thread</span>
      </footer>
    </div>
  );
}
