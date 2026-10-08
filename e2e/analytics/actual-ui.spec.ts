import { test, expect } from '@playwright/test';

test.beforeEach(async ({ page }) => {
  await page.route('https://**', route => route.abort());
  await page.addInitScript(() => {
    localStorage.setItem('flexfoil-onboarding', JSON.stringify({
      lastSeenVersion: '1.0.0', completedTours: ['welcome'], tourProgress: {},
    }));
    localStorage.setItem('ff_cookie_consent', 'granted');
  });
});

test('desktop panel selection and export use the real layout and WASM app', async ({ page }, info) => {
  await page.goto('/');
  await page.locator('.flexlayout__tab_button').filter({ hasText: 'Data Explorer' }).click();
  await expect.poll(() => page.evaluate(() => window.dataLayer?.some(value => {
    const event = Array.from(value as ArrayLike<unknown>);
    return event[0] === 'event' && event[1] === 'feature_use'
      && (event[2] as any).panel_id === 'data-explorer';
  }))).toBe(true);
  await page.getByRole('button', { name: 'File', exact: true }).click();
  const download = page.waitForEvent('download');
  await page.getByText('Export .dat...', { exact: true }).click();
  await download;
  await expect.poll(() => page.evaluate(() => window.dataLayer?.some(value => {
    const event = Array.from(value as ArrayLike<unknown>);
    return event[0] === 'event' && (event[2] as any)?.feature === 'export_dat';
  }))).toBe(true);
  await page.getByRole('button', { name: 'Help', exact: true }).click();
  await page.getByText('Analytics preferences', { exact: true }).click();
  await expect(page.getByRole('button', { name: 'Reject', exact: true })).toBeVisible();
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  await page.getByRole('button', { name: 'Send feedback' }).click();
  await expect(page.getByRole('link', { name: 'View feature request tracker' })).toBeVisible();
  await page.screenshot({ path: info.outputPath('feedback.png'), animations: 'disabled' });
  await info.attach('desktop-feedback', { path: info.outputPath('feedback.png'), contentType: 'image/png' });
});

test('mobile panels, feedback and consent controls are reachable', async ({ page }, info) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/');
  await page.locator('.mobile-tabs').getByRole('button', { name: 'Solve', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.dataLayer?.some(value => {
    const event = Array.from(value as ArrayLike<unknown>);
    return event[0] === 'event' && (event[2] as any)?.panel_id === 'solve';
  }))).toBe(true);
  await page.getByRole('button', { name: 'Send feedback' }).click();
  await expect(page.getByRole('textbox', { name: 'Feedback message' })).toBeVisible();
  await expect(page.getByRole('link', { name: 'View feature request tracker' })).toBeVisible();
  await page.screenshot({ path: info.outputPath('feedback.png'), animations: 'disabled' });
  await info.attach('mobile-feedback', { path: info.outputPath('feedback.png'), contentType: 'image/png' });
  await page.setViewportSize({ width: 320, height: 568 });
  const dialog = page.getByRole('dialog', { name: 'Feedback form' });
  const box = await dialog.boundingBox();
  expect(box?.x).toBeGreaterThanOrEqual(0);
  expect((box?.x ?? 0) + (box?.width ?? 0)).toBeLessThanOrEqual(320);
  await page.getByRole('button', { name: 'Cancel', exact: true }).click();
  await page.getByRole('button', { name: 'Analytics preferences', exact: true }).click();
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  expect(await page.evaluate(() => localStorage.getItem('ff_cookie_consent'))).toBe('denied');
});
