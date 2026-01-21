/**
 * Airfoil Editing Tour - On-demand guide for canvas interaction
 * 
 * Includes challenges to try different editing modes.
 * Uses focusPanel to bring tabbed panels to the front.
 */

import type { TourStep } from './index';

export const airfoilEditingTour: TourStep[] = [
  {
    element: '[data-tour="control-mode-parameters"]',
    focusPanel: 'control',
    popover: {
      title: 'Parameters Mode',
      description: 'The Parameters tab lets you quickly scale the airfoil shape using simple sliders.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="thickness-slider"]',
    focusPanel: 'control',
    popover: {
      title: 'Thickness Slider',
      description: 'Drag this slider to make the airfoil thicker or thinner. 100% is the original thickness.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="camber-slider"]',
    focusPanel: 'control',
    popover: {
      title: 'Camber Slider',
      description: 'Adjust the camber (curvature) of the airfoil. Set to 0% for a symmetric airfoil, or increase for more lift.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="panel-canvas"]',
    focusPanel: 'canvas',
    popover: {
      title: 'Real-time Preview',
      description: 'Watch the airfoil shape and pressure distribution update instantly as you adjust the sliders.',
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
    focusPanel: 'library',
    popover: {
      title: 'Reset to Default',
      description: 'Click this button to restore the airfoil to its original shape. Useful when experimenting goes too far!',
      side: 'right',
      align: 'end',
    },
  },
];
