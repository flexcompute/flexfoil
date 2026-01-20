/**
 * Welcome Tour - Auto-triggered on first visit
 * 
 * Walks users through all main panels and features.
 * Includes a simple challenge to try loading a different airfoil.
 */

import type { TourStep } from './index';

export const welcomeTour: TourStep[] = [
  {
    popover: {
      title: 'Welcome to FlexFoil',
      description: 'Interactive airfoil analysis in your browser. Let\'s take a quick tour of the interface.',
      side: 'over',
      align: 'center',
    },
  },
  {
    element: '[data-tour="panel-library"]',
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
    popover: {
      title: 'Main Canvas',
      description: 'Your airfoil visualization workspace. See the pressure distribution (Cp) and interact with control points to modify the shape in real-time.',
      side: 'left',
      align: 'center',
    },
  },
  {
    element: '[data-tour="panel-control"]',
    popover: {
      title: 'Control Modes',
      description: 'Switch between different editing modes: Surface Points for direct manipulation, Bezier for smooth curves, or B-Spline for maximum flexibility.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="panel-spacing"]',
    popover: {
      title: 'Panel Spacing',
      description: 'Control how points are distributed along the airfoil surface. This affects analysis accuracy and mesh quality.',
      side: 'right',
      align: 'end',
    },
  },
  {
    element: '[data-tour="panel-properties"]',
    popover: {
      title: 'Properties',
      description: 'View computed geometric properties like area, thickness, and camber. These update automatically as you modify the airfoil.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="panel-solve"]',
    popover: {
      title: 'Solver',
      description: 'Run aerodynamic analysis to compute lift (CL) and moment (CM) coefficients. Generate polar sweeps across multiple angles of attack.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="menu-help"]',
    popover: {
      title: 'Need Help?',
      description: 'Access this tutorial anytime from Help → Tutorial. You can also find guides for specific features like editing and solving.',
      side: 'bottom',
      align: 'end',
    },
  },
];
