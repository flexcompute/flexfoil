/**
 * Airfoil Editing Tour - On-demand guide for canvas interaction
 * 
 * Includes challenges to try different editing modes.
 */

import type { TourStep } from './index';

export const airfoilEditingTour: TourStep[] = [
  {
    element: '[data-tour="control-mode-parameters"]',
    popover: {
      title: 'Parameters Mode',
      description: 'Adjust thickness and camber using sliders. Great for quick, proportional scaling of the airfoil shape.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="control-mode-camber"]',
    popover: {
      title: 'Camber Mode',
      description: 'Drag control points to reshape the camber line. Try switching to Camber mode now!',
      side: 'left',
      align: 'start',
    },
    challengeId: 'switch-to-camber',
  },
  {
    element: '[data-tour="control-mode-thickness"]',
    popover: {
      title: 'Thickness Mode',
      description: 'Modify the thickness distribution along the chord. Now try switching to Thickness mode.',
      side: 'left',
      align: 'start',
    },
    challengeId: 'switch-to-thickness',
  },
  {
    element: '[data-tour="panel-canvas"]',
    popover: {
      title: 'Interactive Canvas',
      description: 'Click and drag control points to modify the airfoil. Changes update the pressure distribution in real-time.',
      side: 'left',
      align: 'center',
    },
  },
  {
    element: '[data-tour="menu-edit"]',
    popover: {
      title: 'Undo & Redo',
      description: 'Made a mistake? Use Edit → Undo (Cmd+Z) to go back, or Edit → Redo (Cmd+Shift+Z) to restore changes.',
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="library-reset"]',
    popover: {
      title: 'Reset to Default',
      description: 'Click this button to restore the airfoil to its original shape. Useful when experimenting goes too far!',
      side: 'right',
      align: 'end',
    },
  },
];
