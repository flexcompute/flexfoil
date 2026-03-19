import { test, expect } from '@playwright/test';
import { runTourSteps } from '../helpers/tourRunner';
import { solverQueueShowMe } from '../../flexfoil-ui/src/onboarding/showMe/solverQueue';

test.describe('Show Me: Solver Queue & Status Bar', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('/');
    await page.waitForSelector('[data-tour="panel-canvas"]', { timeout: 15000 });
  });

  test('all tour elements are visible when their panels are focused', async ({ page }) => {
    await runTourSteps(page, solverQueueShowMe, { skipChallenges: true });
  });

  test('full interactive walkthrough completes', async ({ page }) => {
    await runTourSteps(page, solverQueueShowMe);
    const statusIndicator = page.locator('[data-tour="solver-status"]');
    await expect(statusIndicator).toBeVisible();
  });
});
