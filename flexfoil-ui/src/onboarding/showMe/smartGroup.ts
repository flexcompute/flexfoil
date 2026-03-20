/**
 * Show Me: Smart Group & Aggregation
 *
 * Short guided tour demonstrating row grouping, Smart Group, and aggregated statistics.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const smartGroupShowMe: TourStep[] = [
  {
    element: '[data-tour="de-tab-table"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Data Explorer Table',
      description:
        'The table now supports <strong>row grouping</strong> with built-in aggregation — group your runs to see summary statistics like CL_max, L/D_max, and α_stall per group.',
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="de-smart-group"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Smart Group',
      description:
        'One click groups by polar configuration <strong>(airfoil + Re + Mach + Ncrit + flap)</strong> and switches to aggregated statistics. Try it now!',
      side: 'bottom',
      align: 'start',
    },
    challengeId: 'smart-group',
  },
  {
    element: '[data-tour="de-table-view"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Aggregated Columns',
      description:
        'When grouping is active, summary columns appear: <strong>CL_max</strong>, <strong>CD_min</strong>, <strong>L/D_max</strong>, <strong>α_stall</strong>, and <strong>α @ L/D_max</strong>. These use argmax/argmin aggregation to extract the right operating point.',
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="de-tab-correlogram"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Aggregated Data in Plots',
      description:
        'The correlogram and Polar Plot overlay can also use aggregated data — plot group-level statistics like L/D_max vs Re across all your configurations.',
      side: 'bottom',
      align: 'start',
    },
  },
];
