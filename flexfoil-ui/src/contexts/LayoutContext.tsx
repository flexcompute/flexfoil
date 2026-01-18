/**
 * LayoutContext - Provides panel control actions to child components
 */

import { createContext, useContext, useCallback, type ReactNode } from 'react';
import type { Model } from 'flexlayout-react';
import { Actions, DockLocation } from 'flexlayout-react';

interface LayoutContextValue {
  openPanel: (panelId: string) => void;
}

const LayoutContext = createContext<LayoutContextValue | null>(null);

interface LayoutProviderProps {
  children: ReactNode;
  model: Model;
  panels: Array<{ id: string; name: string; component: string }>;
}

export function LayoutProvider({ children, model, panels }: LayoutProviderProps) {
  // Open/focus a panel by ID
  const openPanel = useCallback((panelId: string) => {
    const panelInfo = panels.find((p) => p.id === panelId);
    if (!panelInfo) return;

    // First, try to find an existing tab with this component
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
      // Select the existing tab
      try {
        model.doAction(Actions.selectTab(existingTabId));
      } catch (e) {
        console.warn('Failed to select tab:', e);
      }
    } else {
      // Add a new tab
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
        } catch (e) {
          console.warn('Failed to add panel:', e);
        }
      }
    }
  }, [model, panels]);

  return (
    <LayoutContext.Provider value={{ openPanel }}>
      {children}
    </LayoutContext.Provider>
  );
}

export function useLayout() {
  const context = useContext(LayoutContext);
  if (!context) {
    // Return a no-op if not within provider (for standalone testing)
    return { openPanel: () => {} };
  }
  return context;
}
