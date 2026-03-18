/**
 * Welcome Tour - Auto-triggered on first visit
 * 
 * Walks users through all main panels and features.
 * Includes interactive challenges and visualization demos.
 * 
 * Uses focusPanel to bring tabbed panels to the front when needed.
 */

import type { TourStep } from './index';
import { useVisualizationStore } from '../../stores/visualizationStore';

export const welcomeTour: TourStep[] = [
  {
    popover: {
      title: 'Welcome to FlexFoil',
      description: `Interactive airfoil analysis in your browser. Let's take a quick tour of the interface!
        <div style="margin-top: 12px; padding: 10px; background: var(--bg-tertiary); border-radius: 6px; border-left: 3px solid var(--accent-secondary);">
          <strong>Tip:</strong> If you accidentally dismiss this tutorial, look for the <strong>Tutorial</strong> button in the top-right corner of the screen.
        </div>`,
      side: 'over',
      align: 'center',
    },
  },
  // Docking layout introduction
  {
    popover: {
      title: 'Flexible Layout',
      description: 'FlexFoil uses a docking layout system. All panels can be <strong>resized</strong> by dragging the dividers between them, <strong>moved</strong> by dragging their tabs, and <strong>rearranged</strong> to suit your workflow.',
      side: 'over',
      align: 'center',
    },
  },
  {
    // NOTE: We intentionally don't highlight .flexlayout__splitter directly because
    // dragging it causes layout changes that can break the tour. Instead, we show
    // a centered popover explaining the concept.
    popover: {
      title: 'Resize Panels',
      description: `The thin lines between panels are <strong>splitters</strong>. You can drag them to resize any panel.
        <div style="margin-top: 12px; padding: 10px; background: var(--bg-tertiary); border-radius: 6px; font-size: 13px; color: var(--text-secondary);">
          <strong>Try it:</strong> After this tour, drag any divider to resize panels. Make the canvas larger for detailed editing!
        </div>`,
      side: 'over',
      align: 'center',
    },
  },
  {
    element: '.flexlayout__tab_button',
    popover: {
      title: 'Move & Dock Tabs',
      description: 'Drag any tab to move it. Drop it on another panel to dock alongside existing tabs, or drop it on an edge to create a new split. Double-click a tab to maximize that panel.',
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="menu-window"]',
    popover: {
      title: 'Manage Panels',
      description: 'Use the <strong>Window</strong> menu to restore closed panels or reset the layout to default. You can also close panels using the ✕ button on each tab group.',
      side: 'right',
      align: 'start',
    },
  },
  // Airfoil library
  {
    element: '[data-tour="panel-library"]',
    focusPanel: 'library',
    popover: {
      title: 'Airfoil Library',
      description: 'Generate NACA airfoils by entering a 4-digit code, or use the quick presets to load common airfoils instantly. Try loading a different airfoil now!',
      side: 'right',
      align: 'start',
    },
    challengeId: 'change-airfoil',
  },
  {
    element: '[data-tour="panel-canvas"]',
    focusPanel: 'canvas',
    popover: {
      title: 'Main Canvas',
      description: `Your airfoil visualization workspace. <strong>Pan</strong> with middle-mouse, <strong>zoom</strong> with scroll wheel.
        <div style="margin-top: 10px; font-size: 13px; color: var(--text-secondary);">
          Notice the <strong>Vis Options</strong> button in the bottom-right corner — it's a quick shortcut to the Visualization panel we'll cover later.
        </div>`,
      side: 'left',
      align: 'center',
    },
  },
  {
    element: '[data-tour="panel-control"]',
    focusPanel: 'control',
    popover: {
      title: 'Geometry Control',
      description: 'The <strong>Parameters</strong> tab lets you quickly adjust thickness and camber using sliders. The Camber and Thickness tabs provide finer control with draggable spline points.',
      side: 'left',
      align: 'start',
    },
  },
  // Parameter editing exercise
  {
    element: '[data-tour="thickness-slider"]',
    focusPanel: 'control',
    popover: {
      title: 'Try: Adjust Thickness',
      description: 'Drag the <strong>Thickness Scale</strong> slider to make the airfoil thicker or thinner. Watch the shape update in real-time on the canvas!',
      side: 'left',
      align: 'start',
    },
    challengeId: 'adjust-thickness',
  },
  {
    element: '[data-tour="camber-slider"]',
    focusPanel: 'control',
    popover: {
      title: 'Try: Adjust Camber',
      description: `Now try the <strong>Camber Scale</strong> slider. Increasing camber curves the airfoil for more lift. Setting it to 0% makes a symmetric airfoil.
        <div style="margin-top: 10px; font-size: 13px; color: var(--text-secondary);">
          <strong>Tip:</strong> Click "Apply & Reset Sliders" to save your changes as the new base shape.
        </div>`,
      side: 'left',
      align: 'start',
    },
    challengeId: 'adjust-camber',
  },
  {
    element: '[data-tour="panel-spacing"]',
    focusPanel: 'spacing',
    popover: {
      title: 'Panel Spacing',
      description: 'Control how points are distributed along the airfoil surface. This affects analysis accuracy and mesh quality.',
      side: 'right',
      align: 'end',
    },
  },
  // Panels can be hidden - need to ensure Visualization is visible
  {
    element: '[data-tour="menu-window"]',
    popover: {
      title: 'Not All Panels Are Visible',
      description: 'Some panels may be hidden or closed. You can show, hide, and restore panels from the <strong>Window</strong> menu. Before we continue, make sure the <strong>Visualization</strong> panel is visible.',
      side: 'right',
      align: 'start',
    },
    challengeId: 'open-visualization-panel',
  },
  // Visualization panel and flow features
  {
    element: '[data-tour="panel-visualization"]',
    focusPanel: 'visualization',
    popover: {
      title: 'Visualization Options',
      description: 'Control what\'s displayed on the canvas. Let\'s explore the flow visualization features!',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="viz-streamlines"]',
    focusPanel: 'visualization',
    popover: {
      title: 'Streamlines',
      description: 'Streamlines show the path fluid particles follow. Enable them to visualize the flow pattern around the airfoil.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'enable-streamlines',
  },
  {
    element: '[data-tour="viz-psi"]',
    focusPanel: 'visualization',
    popover: {
      title: 'Stream Function (ψ)',
      description: 'The stream function creates a colored contour map of the flow. Blue shows flow going under the airfoil, red shows flow going over. The dividing streamline marks the boundary.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'enable-psi',
  },
  {
    element: '[data-tour="viz-smoke"]',
    focusPanel: 'visualization',
    popover: {
      title: 'Smoke Visualization',
      description: 'Animated smoke particles show the flow in real-time. Watch how the flow accelerates over the upper surface and separates at the trailing edge!',
      side: 'left',
      align: 'start',
    },
    challengeId: 'enable-smoke',
    onBeforeHighlight: () => {
      // Turn off psi contours so smoke is more visible
      const viz = useVisualizationStore.getState();
      if (viz.showPsiContours) {
        viz.setShowPsiContours(false);
      }
    },
  },
  // Solver panel
  {
    element: '[data-tour="panel-solve"]',
    focusPanel: 'solve',
    popover: {
      title: 'Solver Panel',
      description: 'Run aerodynamic analysis to compute lift (Cl), drag (Cd), and moment (Cm) coefficients. Set the angle of attack and click Run.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-polar"]',
    focusPanel: 'solve',
    popover: {
      title: 'Generate a Polar',
      description: 'A polar sweep analyzes multiple angles of attack at once, creating a complete lift curve. Set the alpha range and click Generate Polar.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'generate-polar',
  },
  // Ensure Solve panel is visible before showing Polar (for context)
  {
    element: '[data-tour="menu-window"]',
    popover: {
      title: 'Keep Solve Panel Visible',
      description: 'Before we look at the polar plot, make sure the <strong>Solve</strong> panel is visible so you can see both panels together.',
      side: 'right',
      align: 'start',
    },
    challengeId: 'open-solve-panel',
  },
  {
    element: '[data-tour="panel-polar"]',
    focusPanel: 'polar',
    popover: {
      title: 'Polar Plot',
      description: 'View your polar data here! You can change the X and Y axes to see Cl vs Alpha, Cd vs Alpha, Cl vs Cm, or other combinations. The data can be exported from the Solve panel.',
      side: 'left',
      align: 'center',
    },
  },
  // Dark mode toggle
  {
    element: '[data-tour="dark-mode-toggle"]',
    popover: {
      title: 'Toggle Dark Mode',
      description: 'For if you\'re one of those dark mode weirdos, like the developer of this code. Click to switch between light and dark themes.',
      side: 'bottom',
      align: 'end',
    },
  },
  // URL sharing
  {
    popover: {
      title: 'Shareable URLs',
      description: 'Every setting you change updates the URL in your browser. Copy and share the URL to give others the exact same view - airfoil, angle of attack, visualization options, viewport position, and more. Great for collaboration or bookmarking interesting configurations!',
      side: 'over',
      align: 'center',
    },
  },
  {
    element: '[data-tour="menu-help"]',
    popover: {
      title: 'Need Help?',
      description: 'The Help menu has additional tutorials for specific features like editing and solving. You can also reset your tutorial progress here if you want to start fresh.',
      side: 'bottom',
      align: 'end',
    },
  },
];
