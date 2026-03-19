import { test, expect } from '@playwright/test';
import { runTourSteps } from '../helpers/tourRunner';
import { multiSweepShowMe } from '../../flexfoil-ui/src/onboarding/showMe/multiSweep';

test.describe('Show Me: Multi-Parameter Sweeps', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('/');
    await page.waitForSelector('[data-tour="panel-canvas"]', { timeout: 15000 });
  });

  test('all tour elements are visible when their panels are focused', async ({ page }) => {
    await runTourSteps(page, multiSweepShowMe, { skipChallenges: true });
  });

  test('full interactive walkthrough completes', async ({ page }) => {
    await runTourSteps(page, multiSweepShowMe);
    const polarPanel = page.locator('[data-tour="panel-polar"]');
    await expect(polarPanel).toBeVisible();
  });
});
