import { test, expect } from '@playwright/test';
import { runTourSteps } from '../helpers/tourRunner';
import { dataAnalysisShowMe } from '../../flexfoil-ui/src/onboarding/showMe/dataAnalysis';

test.describe('Show Me: Enhanced Data Analysis', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('/');
    await page.waitForSelector('[data-tour="panel-canvas"]', { timeout: 15000 });
  });

  test('all tour elements are visible when their panels are focused', async ({ page }) => {
    await runTourSteps(page, dataAnalysisShowMe, { skipChallenges: true });
  });

  test('full interactive walkthrough completes', async ({ page }) => {
    await runTourSteps(page, dataAnalysisShowMe);
    const correlogramTab = page.locator('[data-tour="de-tab-correlogram"]');
    await expect(correlogramTab).toBeVisible();
  });
});
