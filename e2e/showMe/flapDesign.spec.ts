import { test, expect } from '@playwright/test';
import { runTourSteps } from '../helpers/tourRunner';
import { flapDesignShowMe } from '../../flexfoil-ui/src/onboarding/showMe/flapDesign';

test.describe('Show Me: Flap Design', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('/');
    await page.waitForSelector('[data-tour="panel-canvas"]', { timeout: 15000 });
  });

  test('all tour elements are visible when their panels are focused', async ({ page }) => {
    await runTourSteps(page, flapDesignShowMe, { skipChallenges: true });
  });

  test('full interactive walkthrough completes', async ({ page }) => {
    await runTourSteps(page, flapDesignShowMe);
    const flaps = page.locator('[data-tour="gdes-flaps"]');
    await expect(flaps).toContainText('Flap');
  });
});
