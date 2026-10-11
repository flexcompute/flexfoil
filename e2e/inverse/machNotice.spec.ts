import { test, expect } from '@playwright/test';
import { version } from '../../flexfoil-ui/package.json';

test('viscous Mach notice follows the selected operating condition', async ({ page }) => {
  await page.addInitScript((appVersion) => {
    localStorage.setItem('flexfoil-onboarding', JSON.stringify({ completedTours: ['welcome'], tourProgress: {}, lastSeenVersion: '1.0.0' }));
    localStorage.setItem('flexfoil-changelog-seen', appVersion);
  }, version);
  await page.goto('/');
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  await page.getByRole('button', { name: 'Viscous', exact: true }).click();
  await page.getByRole('button', { name: /Advanced Settings/ }).click();
  const mach = page.getByText('Mach Number', { exact: true }).locator('..').locator('input');
  const notice = page.getByRole('note').filter({ hasText: 'Nonzero-Mach viscous analysis' });
  await expect(notice).toHaveCount(0);
  await mach.fill('0.25');
  await expect(notice).toBeVisible();
  await expect(notice.getByRole('link')).toHaveAttribute('href', 'https://github.com/flexcompute/flexfoil/issues/21');
  await page.getByRole('button', { name: /Advanced Settings/ }).click();
  await expect(notice).toBeVisible();
  await page.getByRole('button', { name: /Advanced Settings/ }).click();
  await mach.fill('0');
  await expect(notice).toHaveCount(0);
  await page.locator('[data-tour="solve-polar"] select').first().selectOption('mach');
  await expect(notice).toBeVisible();
  await page.getByRole('button', { name: 'Inviscid', exact: true }).click();
  await expect(notice).toHaveCount(0);
});
