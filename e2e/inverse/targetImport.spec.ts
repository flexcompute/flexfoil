import { test, expect } from '@playwright/test';
import { version } from '../../flexfoil-ui/package.json';

test('QDES imports exact targets and retains them after a malformed file', async ({ page }) => {
  await page.addInitScript((appVersion) => {
    localStorage.setItem('flexfoil-onboarding', JSON.stringify({ completedTours: ['welcome'], tourProgress: {}, lastSeenVersion: '1.0.0' }));
    localStorage.setItem('flexfoil-changelog-seen', appVersion);
  }, version);
  await page.goto('/');
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  await page.locator('.flexlayout__tab_button', { hasText: 'Geometry Control' }).click();
  await page.locator('[data-tour="control-mode-inverse"]').click();
  const input = page.getByLabel('Import upper target');
  await input.setInputFiles({ name: 'target.csv', mimeType: 'text/csv', buffer: Buffer.from('x,cp\n.2,-1.2\n.8,-.5') });
  const knots = page.locator('.inverse-design-panel svg circle');
  await expect(knots).toHaveCount(2);
  const before = await knots.first().getAttribute('cy');
  await input.setInputFiles({ name: 'broken.csv', mimeType: 'text/csv', buffer: Buffer.from('x,cp\n.2,NaN\n.8,-.5') });
  await expect(page.getByRole('alert')).toContainText('Line 2');
  await expect(knots).toHaveCount(2);
  await expect(knots.first()).toHaveAttribute('cy', before!);
});

test('a successful lower import retains an upper import error', async ({ page }) => {
  await page.addInitScript((appVersion) => {
    localStorage.setItem('flexfoil-onboarding', JSON.stringify({ completedTours: ['welcome'], tourProgress: {}, lastSeenVersion: '1.0.0' }));
    localStorage.setItem('flexfoil-changelog-seen', appVersion);
    const read = File.prototype.text;
    File.prototype.text = async function () {
      if (this.name === 'lower.csv') await new Promise(resolve => setTimeout(resolve, 500));
      return read.call(this);
    };
  }, version);
  await page.goto('/');
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  await page.locator('.flexlayout__tab_button', { hasText: 'Geometry Control' }).click();
  await page.locator('[data-tour="control-mode-inverse"]').click();
  await page.getByLabel('Import lower target').setInputFiles({ name: 'lower.csv', mimeType: 'text/csv', buffer: Buffer.from('x,cp\n.2,.1\n.8,.2') });
  await page.getByLabel('Import upper target').setInputFiles({ name: 'broken.csv', mimeType: 'text/csv', buffer: Buffer.from('x,cp\n.2,NaN\n.8,-.5') });
  await expect(page.getByRole('alert')).toContainText('Line 2');
  await expect(page.getByRole('button', { name: '× Lower', exact: true })).toBeVisible();
  await expect(page.getByRole('alert')).toContainText('Line 2');
  await page.getByLabel('Import upper target').setInputFiles({ name: 'upper.csv', mimeType: 'text/csv', buffer: Buffer.from('x,cp\n.2,-1\n.8,-.5') });
  await expect(page.getByRole('alert')).toHaveCount(0);
});
