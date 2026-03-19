/**
 * Show Me: Enhanced Data Analysis
 *
 * Short guided tour highlighting L/D ratio and custom computed columns.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const dataAnalysisShowMe: TourStep[] = [
  {
    element: '[data-tour="de-tab-table"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Data Explorer Table',
      description:
        'The table now includes an automatic <strong>L/D</strong> column computed from Cl and Cd for every run. You can sort, filter, and export it like any other column.',
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="de-column-chips"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Custom Computed Columns',
      description:
        'Click <strong>Add Column</strong> to define algebraic expressions over your data fields — for example <code>cl^2 / cd</code> or <code>cm / cl</code>. These columns appear in the table, correlogram, and Plot Builder.',
      side: 'bottom',
      align: 'start',
    },
    challengeId: 'add-computed-column',
  },
  {
    element: '[data-tour="de-tab-correlogram"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Correlogram',
      description:
        'Switch to the Correlogram view to see pairwise scatter plots of all your data columns — including custom columns. Great for spotting trends across sweep parameters.',
      side: 'bottom',
      align: 'start',
    },
  },
];
