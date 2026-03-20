/**
 * Show Me: Solver Queue & Status Bar
 *
 * Short guided tour showing the solver status indicator and job management.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const solverQueueShowMe: TourStep[] = [
  {
    element: '[data-tour="solver-status"]',
    popover: {
      title: 'Solver Status Indicator',
      description:
        'This dot shows the solver state at a glance: <strong style="color:#22c55e">green</strong> when ready, <strong style="color:#eab308">pulsing yellow</strong> when a job is running.',
      side: 'top',
      align: 'end',
    },
  },
  {
    element: '[data-tour="solve-polar"]',
    focusPanel: 'solve',
    popover: {
      title: 'Start a Sweep',
      description:
        'Run a polar sweep to see the queue in action. Click <strong>Generate Polar</strong> to queue a job.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'run-sweep',
  },
  {
    element: '[data-tour="solver-status"]',
    popover: {
      title: 'Queue Popover',
      description:
        'Click the status indicator to open the queue popover. It shows recent jobs with progress bars, elapsed time, and cancel buttons. Works for polar sweeps, multi-sweeps, and inverse design jobs.',
      side: 'top',
      align: 'end',
    },
  },
];
