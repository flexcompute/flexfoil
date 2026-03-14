/**
 * Solving Tour - On-demand guide for aerodynamic analysis
 * 
 * Includes interactive challenges where users perform actual analysis.
 */

import type { TourStep } from './index';

export const solvingTour: TourStep[] = [
  {
    element: '[data-tour="solve-mode"]',
    focusPanel: 'solve',
    popover: {
      title: 'Solver Mode',
      description: 'FlexFoil now uses the faithful Rust reproduction of XFOIL\'s viscous solver path for aerodynamic analysis.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-alpha"]',
    focusPanel: 'solve',
    popover: {
      title: 'Set Angle of Attack',
      description: 'The angle of attack (α) determines how the airfoil meets the oncoming flow. Let\'s try setting it to 10°.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'set-alpha-10',
  },
  {
    element: '[data-tour="solve-run"]',
    focusPanel: 'solve',
    popover: {
      title: 'Run Analysis',
      description: 'Now click the Run button to compute the aerodynamic coefficients. You\'ll see the lift (CL) and moment (CM) results appear below.',
      side: 'left',
      align: 'start',
    },
    // Note: No challenge here since analysis results aren't stored in global state
  },
  {
    element: '[data-tour="solve-polar"]',
    focusPanel: 'solve',
    popover: {
      title: 'Generate Polar',
      description: 'A polar sweep analyzes multiple angles at once, creating a complete lift curve. Try generating one now!',
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
      title: 'View Polar Plot',
      description: 'Your polar data is displayed here. Change X/Y axes to view different relationships (Cl vs α, Cl vs Cm, etc.).',
      side: 'left',
      align: 'center',
    },
  },
  {
    element: '[data-tour="solve-export"]',
    focusPanel: 'solve',
    popover: {
      title: 'Export Data',
      description: 'Export your polar data as a text file for use in other tools, spreadsheets, or further analysis.',
      side: 'left',
      align: 'end',
    },
  },
];
