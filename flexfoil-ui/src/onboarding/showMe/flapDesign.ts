/**
 * Show Me: Flap Design System
 *
 * Short guided tour demonstrating flap creation and deflection.
 * Launched from the "What's New" changelog dialog.
 */

import type { TourStep } from '../tours';

export const flapDesignShowMe: TourStep[] = [
  {
    element: '[data-tour="control-mode-gdes"]',
    focusPanel: 'control',
    popover: {
      title: 'Switch to Geometry Mode',
      description:
        'The flap editor lives inside the <strong>Geometry</strong> control mode. Click it to open the flap and trailing-edge controls.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'switch-to-gdes',
  },
  {
    element: '[data-tour="gdes-add-flap"]',
    focusPanel: 'control',
    popover: {
      title: 'Add a Flap',
      description:
        'Click <strong>Add Flap</strong> to create a new flap. It starts at 75% chord with zero deflection.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'add-flap',
  },
  {
    element: '[data-tour="gdes-flaps"]',
    focusPanel: 'control',
    popover: {
      title: 'Deflect the Flap',
      description:
        'Change the deflection angle to see the airfoil shape update in real time. XFOIL-style fold trimming keeps the hinge surface clean.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'deflect-flap',
  },
  {
    element: '[data-tour="panel-canvas"]',
    focusPanel: 'canvas',
    popover: {
      title: 'Instant Geometry Update',
      description:
        'The canvas reflects every deflection change instantly — no Apply button needed. Flap definitions are saved with your run data and included in shareable URLs.',
      side: 'left',
      align: 'center',
    },
  },
];
