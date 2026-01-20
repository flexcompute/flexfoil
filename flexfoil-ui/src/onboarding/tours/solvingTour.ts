/**
 * Solving Tour - On-demand guide for aerodynamic analysis
 * 
 * Includes interactive challenges where users perform actual analysis.
 */

import type { TourStep } from './index';

export const solvingTour: TourStep[] = [
  {
    element: '[data-tour="solve-mode"]',
    popover: {
      title: 'Solver Mode',
      description: 'Currently using the inviscid panel method for fast, accurate potential flow analysis. Viscous solver coming in a future update!',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-alpha"]',
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
    popover: {
      title: 'Generate Polar',
      description: 'A polar sweep analyzes multiple angles at once, creating a complete lift curve. Try generating one now!',
      side: 'left',
      align: 'start',
    },
    challengeId: 'generate-polar',
  },
  {
    element: '[data-tour="solve-export"]',
    popover: {
      title: 'Export Data',
      description: 'Export your polar data as a text file for use in other tools, spreadsheets, or further analysis.',
      side: 'left',
      align: 'end',
    },
  },
];
