/**
 * Show Me: Selig Airfoil Database
 *
 * Short guided tour demonstrating the ~1,600 airfoil database browser and Random Foil.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const seligDatabaseShowMe: TourStep[] = [
  {
    element: '[data-tour="panel-library"]',
    focusPanel: 'library',
    popover: {
      title: 'Airfoil Library',
      description:
        'The Library panel now includes the full <strong>UIUC Selig database</strong> — over 1,600 real-world airfoils you can search and load instantly.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="library-selig"]',
    focusPanel: 'library',
    popover: {
      title: 'Search the Database',
      description:
        'Type at least two characters to search by name or designator — try <strong>e387</strong>, <strong>clark</strong>, or <strong>naca</strong>. Click any result to load it.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="library-random-foil"]',
    focusPanel: 'library',
    popover: {
      title: 'Random Foil',
      description:
        'Feeling adventurous? Click <strong>Random Foil</strong> to load a surprise airfoil from the database. Great for exploration and benchmarking.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'random-foil',
  },
  {
    element: '[data-tour="panel-canvas"]',
    focusPanel: 'canvas',
    popover: {
      title: 'Instant Preview',
      description:
        'The loaded airfoil appears immediately on the canvas, ready to analyze. Run a solve or polar sweep right away.',
      side: 'left',
      align: 'center',
    },
  },
];
