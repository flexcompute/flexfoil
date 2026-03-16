/**
 * DockingLayout - FlexLayout-based dockable panel system
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { Layout, Model, Actions, DockLocation } from 'flexlayout-react';

// Panel components
import { AirfoilCanvas } from './AirfoilCanvas';
import { AirfoilLibraryPanel } from './panels/AirfoilLibraryPanel';
import { ControlModePanel } from './panels/ControlModePanel';
import { SpacingPanel } from './panels/SpacingPanel';
import { PropertiesPanel } from './panels/PropertiesPanel';
import { SolvePanel } from './panels/SolvePanel';
import { PolarPanel } from './panels/PolarPanel';
import { VisualizationPanel } from './panels/VisualizationPanel';
import { DataExplorerPanel } from './panels/DataExplorerPanel';
import { PlotBuilderPanel } from './panels/PlotBuilderPanel';
import { CaseLogsPanel } from './panels/CaseLogsPanel';
import { MenuBar } from './MenuBar';
import { MobileLayout } from './MobileLayout';
import { LayoutProvider } from '../contexts/LayoutContext';
import { defaultLayoutJson, PANELS } from '../layoutConfig';
import { useRouteUiStore } from '../stores/routeUiStore';
import { useIsMobile } from '../hooks/useIsMobile';

// Storage keys
const LAYOUT_STORAGE_KEY = 'flexfoil-layout-v4';

interface DockingLayoutProps {
  wasmStatus: 'loading' | 'ready' | 'error';
}

function BrandFooter() {
  return (
    <footer className="brand-footer">
      <span className="brand-footer__label">Powered by Flexcompute Thread</span>
      <a
        className="brand-footer__link"
        href="https://www.flexcompute.com/"
        target="_blank"
        rel="noopener noreferrer"
      >
        Flexcompute.com
      </a>
    </footer>
  );
}

function findSelectedPanelFromJson(node: any): string | null {
  if (!node || typeof node !== 'object') return null;

  if (node.type === 'tabset' && Array.isArray(node.children) && typeof node.selected === 'number') {
    const selected = node.children[node.selected];
    if (selected?.component) {
      return selected.component;
    }
  }

  if (Array.isArray(node.children)) {
    for (const child of node.children) {
      const match = findSelectedPanelFromJson(child);
      if (match) return match;
    }
  }

  return null;
}

export function DockingLayout({ wasmStatus }: DockingLayoutProps) {
  const isMobile = useIsMobile();

  if (isMobile) {
    return <MobileLayout wasmStatus={wasmStatus} />;
  }

  return <DesktopLayout wasmStatus={wasmStatus} />;
}

function DesktopLayout({ wasmStatus }: DockingLayoutProps) {
  const layoutJson = useRouteUiStore((state) => state.layoutJson);
  const layoutRevision = useRouteUiStore((state) => state.layoutRevision);
  const activePanel = useRouteUiStore((state) => state.activePanel);
  const activePanelRevision = useRouteUiStore((state) => state.activePanelRevision);
  const setLayoutJson = useRouteUiStore((state) => state.setLayoutJson);
  const setActivePanel = useRouteUiStore((state) => state.setActivePanel);

  // Initialize model from localStorage or default
  const [model, setModel] = useState(() => {
    if (layoutJson) {
      return Model.fromJson(layoutJson);
    }
    try {
      const savedLayout = localStorage.getItem(LAYOUT_STORAGE_KEY);
      if (savedLayout) {
        return Model.fromJson(JSON.parse(savedLayout));
      }
    } catch (error) {
      console.warn('Failed to load saved layout, using default:', error);
    }
    return Model.fromJson(defaultLayoutJson);
  });

  const layoutRef = useRef<Layout>(null);
  const [closedPanels, setClosedPanels] = useState<Set<string>>(new Set());
  const lastAppliedLayoutRevision = useRef(layoutRevision);
  const lastAppliedActivePanelRevision = useRef(activePanelRevision);

  const applyClosedPanelsFromModel = useCallback((nextModel: Model) => {
    const visiblePanels = new Set<string>();
    const traverse = (node: any) => {
      if (node.getType?.() === 'tab') {
        const component = node.getComponent?.();
        if (component && PANELS.some((panel) => panel.id === component)) {
          visiblePanels.add(component);
        }
      }
      if (node.getChildren) {
        node.getChildren().forEach(traverse);
      }
    };
    traverse(nextModel.getRoot());

    const nextClosed = new Set<string>();
    PANELS.forEach((panel) => {
      if (!visiblePanels.has(panel.id)) {
        nextClosed.add(panel.id);
      }
    });
    setClosedPanels(nextClosed);
  }, []);

  // Save layout to localStorage when it changes
  const handleModelChange = useCallback((newModel: Model) => {
    const json = newModel.toJson();
    try {
      localStorage.setItem(LAYOUT_STORAGE_KEY, JSON.stringify(json));
    } catch (error) {
      console.warn('Failed to save layout:', error);
    }

    setLayoutJson(json);
    applyClosedPanelsFromModel(newModel);

    const selected = findSelectedPanelFromJson(json.layout);
    if (selected && PANELS.some((panel) => panel.id === selected)) {
      setActivePanel(selected as typeof activePanel);
    }
  }, [activePanel, applyClosedPanelsFromModel, setActivePanel, setLayoutJson]);

  // Reset layout to default
  const handleResetLayout = useCallback(() => {
    const newModel = Model.fromJson(defaultLayoutJson);
    setModel(newModel);
    setClosedPanels(new Set());
    setLayoutJson(null);
    setActivePanel('canvas');
    localStorage.removeItem(LAYOUT_STORAGE_KEY);
  }, [setActivePanel, setLayoutJson]);

  // Focus/select an existing panel (bring to front)
  const handleOpenPanel = useCallback(
    (panelId: string) => {
      const root = model.getRoot();
      let existingTabId: string | null = null;
      
      const findTab = (node: any) => {
        if (node.getType?.() === 'tab' && node.getComponent?.() === panelId) {
          existingTabId = node.getId();
          return true;
        }
        if (node.getChildren) {
          for (const child of node.getChildren()) {
            if (findTab(child)) return true;
          }
        }
        return false;
      };
      findTab(root);

      if (existingTabId) {
        try {
          model.doAction(Actions.selectTab(existingTabId));
          setActivePanel(panelId as typeof activePanel);
        } catch (e) {
          console.warn('Failed to select tab:', e);
        }
      }
    },
    [activePanel, model, setActivePanel]
  );

  // Restore a closed panel
  const handleRestorePanel = useCallback(
    (panelId: string) => {
      const panelInfo = PANELS.find((p) => p.id === panelId);
      if (!panelInfo) return;

      // Find a tabset to add it to
      const root = model.getRoot();
      const tabsets: any[] = [];
      const findTabsets = (node: any) => {
        if (node.getType?.() === 'tabset') {
          tabsets.push(node);
        }
        if (node.getChildren) {
          node.getChildren().forEach(findTabsets);
        }
      };
      findTabsets(root);

      if (tabsets.length > 0) {
        try {
          model.doAction(
            Actions.addNode(
              {
                type: 'tab',
                id: `${panelId}-${Date.now()}`,
                name: panelInfo.name,
                component: panelId,
              },
              tabsets[0].getId(),
              DockLocation.CENTER,
              -1
            )
          );
          setClosedPanels((prev) => {
            const next = new Set(prev);
            next.delete(panelId);
            return next;
          });
          setActivePanel(panelId as typeof activePanel);
        } catch (e) {
          console.warn('Failed to restore panel:', e);
        }
      }
    },
    [activePanel, model, setActivePanel]
  );

  useEffect(() => {
    applyClosedPanelsFromModel(model);
  }, [applyClosedPanelsFromModel, model]);

  useEffect(() => {
    if (layoutRevision === lastAppliedLayoutRevision.current) {
      return;
    }
    lastAppliedLayoutRevision.current = layoutRevision;

    if (layoutJson) {
      try {
        const nextModel = Model.fromJson(layoutJson);
        setModel(nextModel);
        applyClosedPanelsFromModel(nextModel);
      } catch (error) {
        console.warn('Failed to apply route layout:', error);
      }
      return;
    }

    const nextModel = Model.fromJson(defaultLayoutJson);
    setModel(nextModel);
    applyClosedPanelsFromModel(nextModel);
  }, [applyClosedPanelsFromModel, layoutJson, layoutRevision]);

  useEffect(() => {
    if (activePanelRevision === lastAppliedActivePanelRevision.current) {
      return;
    }
    lastAppliedActivePanelRevision.current = activePanelRevision;

    if (closedPanels.has(activePanel)) {
      handleRestorePanel(activePanel);
    } else {
      handleOpenPanel(activePanel);
    }
  }, [activePanel, activePanelRevision, closedPanels, handleOpenPanel, handleRestorePanel]);

  // Factory function to create panel components
  // Each panel is wrapped with data-tour attribute for onboarding targeting
  const factory = useCallback((node: any) => {
    const component = node.getComponent();

    switch (component) {
      case 'canvas':
        return (
          <div data-tour="panel-canvas" style={{ width: '100%', height: '100%' }}>
            <AirfoilCanvas />
          </div>
        );
      case 'library':
        return (
          <div data-tour="panel-library" style={{ width: '100%', height: '100%' }}>
            <AirfoilLibraryPanel />
          </div>
        );
      case 'control':
        return (
          <div data-tour="panel-control" style={{ width: '100%', height: '100%' }}>
            <ControlModePanel />
          </div>
        );
      case 'spacing':
        return (
          <div data-tour="panel-spacing" style={{ width: '100%', height: '100%' }}>
            <SpacingPanel />
          </div>
        );
      case 'properties':
        return (
          <div data-tour="panel-properties" style={{ width: '100%', height: '100%' }}>
            <PropertiesPanel />
          </div>
        );
      case 'solve':
        return (
          <div data-tour="panel-solve" style={{ width: '100%', height: '100%' }}>
            <SolvePanel />
          </div>
        );
      case 'polar':
        return (
          <div data-tour="panel-polar" style={{ width: '100%', height: '100%' }}>
            <PolarPanel />
          </div>
        );
      case 'visualization':
        return (
          <div data-tour="panel-visualization" style={{ width: '100%', height: '100%' }}>
            <VisualizationPanel />
          </div>
        );
      case 'data-explorer':
        return (
          <div data-tour="panel-data-explorer" style={{ width: '100%', height: '100%' }}>
            <DataExplorerPanel />
          </div>
        );
      case 'plot-builder':
        return (
          <div data-tour="panel-plot-builder" style={{ width: '100%', height: '100%' }}>
            <PlotBuilderPanel />
          </div>
        );
      case 'case-logs':
        return (
          <div style={{ width: '100%', height: '100%' }}>
            <CaseLogsPanel />
          </div>
        );
      default:
        return <div className="panel">Unknown panel: {component}</div>;
    }
  }, []);

  const layoutPanels = useMemo(() => [...PANELS], []);

  // Custom tabset header buttons
  const onRenderTabSet = useCallback(
    (tabsetNode: any, renderValues: any) => {
      // Add minimize button
      renderValues.buttons.push(
        <button
          key="close"
          onClick={(e) => {
            e.stopPropagation();
            const selectedTab = tabsetNode.getSelectedNode();
            if (selectedTab) {
              const panelId = selectedTab.getComponent();
              setClosedPanels((prev) => new Set([...prev, panelId]));
              model.doAction(Actions.deleteTab(selectedTab.getId()));
            }
          }}
          className="flexlayout__tab_toolbar_button"
          title="Close panel"
          style={{ padding: '2px 6px' }}
        >
          ✕
        </button>
      );
    },
    [model]
  );

  return (
    <LayoutProvider model={model} panels={layoutPanels}>
      <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
        {/* Menu bar */}
        <MenuBar
          panels={layoutPanels}
          closedPanels={closedPanels}
          onRestorePanel={handleRestorePanel}
          onOpenPanel={handleOpenPanel}
          onResetLayout={handleResetLayout}
          wasmStatus={wasmStatus}
        />

        {/* FlexLayout container */}
        <div style={{ flex: 1, position: 'relative' }}>
          <Layout
            ref={layoutRef}
            model={model}
            factory={factory}
            onModelChange={handleModelChange}
            onRenderTabSet={onRenderTabSet}
          />
        </div>
        <BrandFooter />
      </div>
    </LayoutProvider>
  );
}
