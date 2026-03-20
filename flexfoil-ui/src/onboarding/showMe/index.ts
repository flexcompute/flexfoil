/**
 * Show Me micro-tours — short, feature-focused tutorials launched from the changelog.
 *
 * Each tour is a small TourStep[] that demonstrates a single new feature.
 * The same step definitions are consumed by:
 *   1. driver.js (in-app "Show Me" tutorial via OnboardingProvider)
 *   2. Playwright E2E tests (via e2e/helpers/tourRunner.ts)
 */

import type { TourStep } from '../tours';
import { flapDesignShowMe } from './flapDesign';
import { multiSweepShowMe } from './multiSweep';
import { dataAnalysisShowMe } from './dataAnalysis';
import { solverQueueShowMe } from './solverQueue';
import { seligDatabaseShowMe } from './seligDatabase';
import { smartGroupShowMe } from './smartGroup';

export const showMeTours: Record<string, TourStep[]> = {
  'showMe:flapDesign': flapDesignShowMe,
  'showMe:multiSweep': multiSweepShowMe,
  'showMe:dataAnalysis': dataAnalysisShowMe,
  'showMe:solverQueue': solverQueueShowMe,
  'showMe:seligDatabase': seligDatabaseShowMe,
  'showMe:smartGroup': smartGroupShowMe,
};

export {
  flapDesignShowMe,
  multiSweepShowMe,
  dataAnalysisShowMe,
  solverQueueShowMe,
  seligDatabaseShowMe,
  smartGroupShowMe,
};
