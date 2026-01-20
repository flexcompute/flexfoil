/**
 * Tour registry - exports all available tours
 */

import type { DriveStep } from 'driver.js';
import { welcomeTour } from './welcomeTour';
import { airfoilEditingTour } from './airfoilEditingTour';
import { solvingTour } from './solvingTour';

export type TourId = 'welcome' | 'airfoilEditing' | 'solving';

/**
 * Extended tour step that can include a challenge
 */
export interface TourStep extends DriveStep {
  /** Optional challenge ID - if set, user must complete the challenge to proceed */
  challengeId?: string;
}

export const tours: Record<TourId, TourStep[]> = {
  welcome: welcomeTour,
  airfoilEditing: airfoilEditingTour,
  solving: solvingTour,
};

export { welcomeTour, airfoilEditingTour, solvingTour };
