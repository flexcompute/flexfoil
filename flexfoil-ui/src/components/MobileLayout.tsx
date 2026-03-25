import { useState, useCallback, useMemo } from 'react';
import { AirfoilCanvas } from './AirfoilCanvas';
import { SolvePanel } from './panels/SolvePanel';
import { AirfoilLibraryPanel } from './panels/AirfoilLibraryPanel';
import { PolarPanel } from './panels/PolarPanel';
import { ControlModePanel } from './panels/ControlModePanel';
import { SpacingPanel } from './panels/SpacingPanel';
import { PropertiesPanel } from './panels/PropertiesPanel';
import { VisualizationPanel } from './panels/VisualizationPanel';
import { DataExplorerPanel } from './panels/DataExplorerPanel';
import { PlotBuilderPanel } from './panels/PlotBuilderPanel';
import { CaseLogsPanel } from './panels/CaseLogsPanel';
import { DistributionPanel } from './panels/DistributionPanel';
import { CfdPanel } from './panels/CfdPanel';
import { DarkModeToggle } from './DarkModeToggle';
import { LayoutContext } from '../contexts/LayoutContext';
import { useOnboarding } from '../onboarding';
import type { PanelId } from '../layoutConfig';
import { FlexcomputeLogo } from './FlexcomputeLogo';

const TABS: { id: PanelId; label: string }[] = [
  { id: 'canvas', label: 'Canvas' },
  { id: 'solve', label: 'Solve' },
  { id: 'library', label: 'Library' },
  { id: 'polar', label: 'Polar' },
  { id: 'control', label: 'Controls' },
  { id: 'spacing', label: 'Spacing' },
  { id: 'properties', label: 'Properties' },
  { id: 'visualization', label: 'Vis' },
  { id: 'data-explorer', label: 'Explorer' },
  { id: 'plot-builder', label: 'Plots' },
  { id: 'case-logs', label: 'Logs' },
  { id: 'distributions', label: 'Distrib.' },
  { id: 'cfd', label: 'CFD' },
];

const PANEL_COMPONENTS: Record<PanelId, React.FC> = {
  'canvas': AirfoilCanvas,
  'solve': SolvePanel,
  'library': AirfoilLibraryPanel,
  'polar': PolarPanel,
  'control': ControlModePanel,
  'spacing': SpacingPanel,
  'properties': PropertiesPanel,
  'visualization': VisualizationPanel,
  'data-explorer': DataExplorerPanel,
  'plot-builder': PlotBuilderPanel,
  'case-logs': CaseLogsPanel,
  'distributions': DistributionPanel,
  'cfd': CfdPanel,
};

interface MobileLayoutProps {
  wasmStatus: 'loading' | 'ready' | 'error';
}

export function MobileLayout({ wasmStatus }: MobileLayoutProps) {
  const [activeTab, setActiveTab] = useState<PanelId>('canvas');
  const [bannerDismissed, setBannerDismissed] = useState(false);
  const { startTour, hasStartedTour, isActive } = useOnboarding();

  const openPanel = useCallback((panelId: string) => {
    const match = TABS.find((t) => t.id === panelId);
    if (match) setActiveTab(match.id);
  }, []);

  const layoutCtx = useMemo(() => ({ openPanel }), [openPanel]);

  const ActivePanel = PANEL_COMPONENTS[activeTab];
  const isCanvas = activeTab === 'canvas';

  return (
    <LayoutContext.Provider value={layoutCtx}>
      <div className="mobile-layout">
        {/* Compact header */}
        <header className="mobile-header">
          <a
            href="https://www.flexcompute.com/"
            target="_blank"
            rel="noopener noreferrer"
            className="mobile-header__brand"
          >
            <FlexcomputeLogo height={18} />
            <span className="mobile-header__title">FlexFoil</span>
          </a>
          <div className="mobile-header__right">
            {hasStartedTour('welcome') && !isActive && (
              <button
                className="mobile-header__tutorial"
                onClick={() => startTour('welcome')}
              >
                ? Tour
              </button>
            )}
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
              The solver runs great on mobile — for the full editing experience, try desktop.
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

        {/* Swipeable tab ribbon */}
        <nav className="mobile-tabs">
          {TABS.map((tab) => (
            <button
              key={tab.id}
              className={`mobile-tabs__btn${activeTab === tab.id ? ' mobile-tabs__btn--active' : ''}`}
              onClick={() => setActiveTab(tab.id)}
            >
              {tab.label}
            </button>
          ))}
        </nav>

        {/* Content area */}
        <div className="mobile-content">
          <div
            className={`mobile-card${isCanvas ? ' mobile-card--canvas' : ''}`}
            data-tour={`panel-${activeTab}`}
          >
            <ActivePanel />
          </div>
        </div>

        {/* Compact footer */}
        <footer className="mobile-footer">
          <span>Powered by Flexcompute Thread</span>
        </footer>
      </div>
    </LayoutContext.Provider>
  );
}
