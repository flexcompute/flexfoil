/**
 * tourRunner — Converts TourStep[] definitions into Playwright test actions.
 *
 * This is the bridge between in-app "Show Me" tutorials and E2E tests.
 * Both consume the same TourStep[] from onboarding/showMe/*.ts, ensuring
 * the tutorial script and the test script stay in sync.
 *
 * For each step the runner:
 *   1. Focuses the required panel (clicks its tab)
 *   2. Waits for the target element to be visible
 *   3. Executes the challenge action (if the step has a challengeId)
 *   4. Asserts the expected state after the action
 */

import { type Page, expect } from '@playwright/test';

export interface TourStepLike {
  element?: string;
  focusPanel?: string;
  challengeId?: string;
  popover?: { title?: string; description?: string };
}

/**
 * Maps challengeId -> a function that performs the user action in Playwright
 * and asserts the expected post-condition.
 */
type ChallengeAction = (page: Page) => Promise<void>;

const challengeActions: Record<string, ChallengeAction> = {
  'switch-to-gdes': async (page) => {
    await page.locator('[data-tour="control-mode-gdes"]').click();
    await expect(page.locator('[data-tour="control-mode-gdes"].active')).toBeVisible();
  },

  'add-flap': async (page) => {
    await page.locator('[data-tour="gdes-add-flap"]').click();
    await expect(page.locator('[data-tour="gdes-flaps"]')).toContainText('Flap');
  },

  'deflect-flap': async (page) => {
    const deflectionInput = page.locator('[data-tour="gdes-flaps"] input[type="number"]').first();
    await deflectionInput.fill('10');
    await deflectionInput.press('Enter');
  },

  'generate-polar': async (page) => {
    const generateBtn = page.locator('[data-tour="solve-polar"] button', { hasText: /generate/i });
    await generateBtn.click();
    await page.waitForTimeout(2000);
  },

  'run-sweep': async (page) => {
    const generateBtn = page.locator('[data-tour="solve-polar"] button', { hasText: /generate/i });
    await generateBtn.click();
    await page.waitForTimeout(1000);
  },

  'set-alpha-10': async (page) => {
    const alphaInput = page.locator('[data-tour="solve-alpha"] input[type="number"]');
    await alphaInput.fill('10');
    await alphaInput.press('Enter');
  },

  'change-airfoil': async (page) => {
    const nacaInput = page.locator('[data-tour="panel-library"] input').first();
    await nacaInput.fill('4412');
    const generateBtn = page.locator('[data-tour="panel-library"] button', { hasText: /generate/i });
    await generateBtn.click();
    await page.waitForTimeout(500);
  },

  'add-computed-column': async (page) => {
    const addBtn = page.locator('[data-tour="de-column-chips"] button', { hasText: /add/i });
    if (await addBtn.isVisible()) {
      await addBtn.click();
      await page.waitForTimeout(500);
    }
  },

  'open-visualization-panel': async (page) => {
    await focusPanelByMenu(page, 'Visualization');
  },

  'open-solve-panel': async (page) => {
    await focusPanelByMenu(page, 'Solve');
  },

  'enable-streamlines': async (page) => {
    const toggle = page.locator('[data-tour="viz-streamlines"] input[type="checkbox"]');
    if (!(await toggle.isChecked())) await toggle.click();
  },

  'enable-psi': async (page) => {
    const toggle = page.locator('[data-tour="viz-psi"] input[type="checkbox"]');
    if (!(await toggle.isChecked())) await toggle.click();
  },

  'enable-smoke': async (page) => {
    const toggle = page.locator('[data-tour="viz-smoke"] input[type="checkbox"]');
    if (!(await toggle.isChecked())) await toggle.click();
  },

  'adjust-thickness': async (page) => {
    const slider = page.locator('[data-tour="thickness-slider"] input[type="range"]');
    await slider.fill('1.2');
  },

  'adjust-camber': async (page) => {
    const slider = page.locator('[data-tour="camber-slider"] input[type="range"]');
    await slider.fill('1.3');
  },
};

async function focusPanelByMenu(page: Page, panelName: string) {
  await page.locator('[data-tour="menu-window"]').click();
  await page.locator(`text=${panelName}`).click();
  await page.waitForTimeout(300);
}

async function focusPanel(page: Page, panelId: string) {
  const tabButton = page.locator(`.flexlayout__tab_button`, { hasText: new RegExp(panelId, 'i') });
  if (await tabButton.isVisible()) {
    await tabButton.click();
    await page.waitForTimeout(200);
    return;
  }
  const panelNames: Record<string, string> = {
    'control': 'Control',
    'solve': 'Solve',
    'library': 'Library',
    'visualization': 'Visualization',
    'polar': 'Polar',
    'canvas': 'Canvas',
    'spacing': 'Spacing',
    'data-explorer': 'Data Explorer',
    'plot-builder': 'Plot Builder',
    'properties': 'Properties',
  };
  const name = panelNames[panelId] ?? panelId;
  await focusPanelByMenu(page, name);
}

/**
 * Run a tour's steps as Playwright actions.
 *
 * @param page - Playwright Page object
 * @param steps - TourStep[] from a showMe micro-tour
 * @param options.skipChallenges - if true, only verify element visibility (no actions)
 */
export async function runTourSteps(
  page: Page,
  steps: TourStepLike[],
  options: { skipChallenges?: boolean } = {},
) {
  for (const step of steps) {
    if (step.focusPanel) {
      await focusPanel(page, step.focusPanel);
    }

    if (step.element) {
      const locator = page.locator(step.element).first();
      await expect(locator).toBeVisible({ timeout: 5000 });
    }

    if (step.challengeId && !options.skipChallenges) {
      const action = challengeActions[step.challengeId];
      if (action) {
        await action(page);
      }
    }
  }
}
