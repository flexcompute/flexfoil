/**
 * Show Me: Reynolds Type Modes (1/2/3)
 *
 * Short guided tour showing the Mode 1/2/3 selector in the Solve panel.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const reTypeModesShowMe: TourStep[] = [
  {
    element: '[data-tour="solve-mode"]',
    focusPanel: 'solve',
    popover: {
      title: 'Viscous Solver',
      description:
        'Make sure the <strong>Viscous</strong> solver is active — Mode selection only applies to viscous analysis.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-re-type"]',
    focusPanel: 'solve',
    popover: {
      title: 'Reynolds Type Modes',
      description:
        '<strong>Mode 1</strong> (default): constant Re. ' +
        '<strong>Mode 2</strong>: fixed Re·√CL — Re varies with speed as CL changes (aircraft cruise). ' +
        '<strong>Mode 3</strong>: fixed Re·CL — for propeller or turbomachinery blades.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-alpha"]',
    focusPanel: 'solve',
    popover: {
      title: 'Try It',
      description:
        'Set a mode, enter an angle of attack, and click <strong>Run</strong>. ' +
        'In Mode 2 or 3, the effective Reynolds number will vary with CL — you\'ll see different drag predictions compared to Mode 1.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'run-solve',
  },
];
