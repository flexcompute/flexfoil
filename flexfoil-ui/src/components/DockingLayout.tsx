/**
 * DockingLayout - FlexLayout-based dockable panel system
 */

import { useCallback, useRef, useState } from 'react';
import { Layout, Model, Actions, DockLocation } from 'flexlayout-react';
import type { IJsonModel } from 'flexlayout-react';

// Panel components
import { AirfoilCanvas } from './AirfoilCanvas';
import { AirfoilLibraryPanel } from './panels/AirfoilLibraryPanel';
import { ControlModePanel } from './panels/ControlModePanel';
import { SpacingPanel } from './panels/SpacingPanel';
import { PropertiesPanel } from './panels/PropertiesPanel';
import { SolvePanel } from './panels/SolvePanel';
import { MenuBar } from './MenuBar';

// Storage keys
const LAYOUT_STORAGE_KEY = 'flexfoil-layout-v1';

// Panel definitions
export const PANELS = [
  { id: 'canvas', name: 'Airfoil Canvas', component: 'canvas' },
  { id: 'library', name: 'Airfoil Library', component: 'library' },
  { id: 'control', name: 'Control Mode', component: 'control' },
  { id: 'spacing', name: 'Spacing', component: 'spacing' },
  { id: 'properties', name: 'Properties', component: 'properties' },
  { id: 'solve', name: 'Solve', component: 'solve' },
];

// Default layout configuration
const defaultLayoutJson: IJsonModel = {
  global: {
    tabSetMinWidth: 200,
    tabSetMinHeight: 150,
    borderMinSize: 100,
    splitterSize: 6,
  },
  borders: [],
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'row',
        weight: 25,
        children: [
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'library', name: 'Airfoil Library', component: 'library' },
            ],
          },
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'spacing', name: 'Spacing', component: 'spacing' },
            ],
          },
        ],
      },
      {
        type: 'tabset',
        weight: 50,
        children: [
          { type: 'tab', id: 'canvas', name: 'Airfoil Canvas', component: 'canvas' },
        ],
      },
      {
        type: 'row',
        weight: 25,
        children: [
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'properties', name: 'Properties', component: 'properties' },
              { type: 'tab', id: 'solve', name: 'Solve', component: 'solve' },
            ],
          },
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'control', name: 'Control Mode', component: 'control' },
            ],
          },
        ],
      },
    ],
  },
};

interface DockingLayoutProps {
  wasmStatus: 'loading' | 'ready' | 'error';
}

export function DockingLayout({ wasmStatus }: DockingLayoutProps) {
  // Initialize model from localStorage or default
  const [model, setModel] = useState(() => {
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

  // Save layout to localStorage when it changes
  const handleModelChange = useCallback((newModel: Model) => {
    try {
      localStorage.setItem(LAYOUT_STORAGE_KEY, JSON.stringify(newModel.toJson()));
    } catch (error) {
      console.warn('Failed to save layout:', error);
    }

    // Track closed panels
    const visiblePanels = new Set<string>();
    const traverse = (node: any) => {
      if (node.getType?.() === 'tab') {
        const component = node.getComponent?.();
        if (component && PANELS.some((p) => p.id === component)) {
          visiblePanels.add(component);
        }
      }
      if (node.getChildren) {
        node.getChildren().forEach(traverse);
      }
    };
    traverse(newModel.getRoot());

    const newClosed = new Set<string>();
    PANELS.forEach((p) => {
      if (!visiblePanels.has(p.id)) {
        newClosed.add(p.id);
      }
    });
    setClosedPanels(newClosed);
  }, []);

  // Reset layout to default
  const handleResetLayout = useCallback(() => {
    const newModel = Model.fromJson(defaultLayoutJson);
    setModel(newModel);
    setClosedPanels(new Set());
    localStorage.removeItem(LAYOUT_STORAGE_KEY);
  }, []);

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
        } catch (e) {
          console.warn('Failed to restore panel:', e);
        }
      }
    },
    [model]
  );

  // Factory function to create panel components
  const factory = useCallback((node: any) => {
    const component = node.getComponent();

    switch (component) {
      case 'canvas':
        return <AirfoilCanvas />;
      case 'library':
        return <AirfoilLibraryPanel />;
      case 'control':
        return <ControlModePanel />;
      case 'spacing':
        return <SpacingPanel />;
      case 'properties':
        return <PropertiesPanel />;
      case 'solve':
        return <SolvePanel />;
      default:
        return <div className="panel">Unknown panel: {component}</div>;
    }
  }, []);

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
    <div style={{ width: '100%', height: '100%', display: 'flex', flexDirection: 'column' }}>
      {/* Menu bar */}
      <MenuBar
        panels={PANELS}
        closedPanels={closedPanels}
        onRestorePanel={handleRestorePanel}
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
    </div>
  );
}
