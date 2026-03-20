/**
 * Show Me: Multi-Parameter Sweeps
 *
 * Short guided tour showing the sweep configuration and matrix sweep feature.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const multiSweepShowMe: TourStep[] = [
  {
    element: '[data-tour="solve-polar"]',
    focusPanel: 'solve',
    popover: {
      title: 'Polar Sweep Section',
      description:
        'The sweep engine lives here. You can sweep alpha, Reynolds, Mach, Ncrit, flap deflection, or flap hinge x/c — all from this panel.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="solve-polar"]',
    focusPanel: 'solve',
    popover: {
      title: 'Generate a Polar',
      description:
        'Click <strong>Generate Polar</strong> to run a sweep. The solver caches previously computed points — repeated sweeps skip what\'s already done.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'generate-polar',
  },
  {
    element: '[data-tour="solver-status"]',
    popover: {
      title: 'Solver Queue',
      description:
        'The status indicator turns yellow while a sweep is running. Click it to see the queue popover with per-job progress bars, elapsed time, and cancel buttons.',
      side: 'top',
      align: 'end',
    },
  },
  {
    element: '[data-tour="panel-polar"]',
    focusPanel: 'polar',
    popover: {
      title: 'Polar Plot',
      description:
        'Results appear here. The axis dropdowns now include all swept parameters (Re, Mach, Ncrit, flap deflection, hinge x/c) so you can explore multi-dimensional data.',
      side: 'left',
      align: 'center',
    },
  },
];
