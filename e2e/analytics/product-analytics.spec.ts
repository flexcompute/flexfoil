import { test, expect } from '@playwright/test';

for (const consent of ['missing', 'denied']) {
  test(`solving works without ${consent} analytics consent`, async ({ page }) => {
    await page.route('https://**', route => route.abort());
    await page.goto(`/e2e/analytics/fixture.html?consent=${consent}`);
    await page.getByRole('button', { name: 'Run', exact: true }).click();
    await expect.poll(() => page.evaluate(() => (window as any).__solves.length)).toBe(1);
    expect(await page.evaluate(() => (window as any).__events)).toEqual([]);
  });
}

test('feature request handoff is public, correctly categorized and never claims submission', async ({ page }, info) => {
  test.skip(info.project.name === 'feedback-service');
  await page.route('https://**', route => route.abort());
  await page.goto('/e2e/analytics/fixture.html?geometry=private#private-foil');
  await page.evaluate(() => { window.open = url => { (window as any).__draftUrl = String(url); return null; }; });
  await page.getByRole('button', { name: 'Send feedback' }).click();
  await page.getByRole('button', { name: 'Feature', exact: true }).click();
  await page.getByRole('textbox', { name: 'Feedback message' }).fill('Please add a comparison plot\nPrivate example for the request only');
  const tracker = page.getByRole('link', { name: 'View feature request tracker' });
  expect(new URL(await tracker.getAttribute('href') as string).searchParams.get('q')).toBe('is:issue label:enhancement');
  await expect(page.getByText('A GitHub account is required', { exact: false })).toBeVisible();
  await page.getByRole('button', { name: 'Continue on GitHub' }).click();
  const url = new URL(await page.evaluate(() => (window as any).__draftUrl));
  expect(url.origin).toBe('https://github.com');
  expect(url.pathname).toBe('/flexcompute/flexfoil/issues/new');
  expect(url.searchParams.get('template')).toBe('feature-request.md');
  expect(url.searchParams.get('labels')).toBe('enhancement');
  expect(url.searchParams.get('body')).toContain('Please add a comparison plot');
  await expect(page.getByText('Opening a draft does not submit it.', { exact: false })).toBeVisible();
  const events = await page.evaluate(() => (window as any).__events);
  if (info.project.name === 'explicit-on') {
    expect(events.map((event: any[]) => event[1])).toEqual(['feature_use', 'feedback_handoff']);
    expect(events[1][2].feedback_type).toBe('feature');
    expect(events[1][2].page_location).toBe(new URL('/', info.project.use.baseURL).href);
  } else expect(events).toEqual([]);
  expect(JSON.stringify(events)).not.toContain('private');
  expect(JSON.stringify(events)).not.toContain('comparison plot');
  await info.attach('request-handoff', { body: await page.screenshot(), contentType: 'image/png' });
  await page.getByRole('button', { name: 'Close', exact: true }).click();
  await page.getByRole('button', { name: 'Send feedback' }).click();
  await expect(page.getByRole('textbox', { name: 'Feedback message' })).toHaveValue('Please add a comparison plot\nPrivate example for the request only');
});

test('acceptance and revocation control subsequent custom events', async ({ page }, info) => {
  await page.route('https://**', route => route.abort());
  await page.goto('/e2e/analytics/fixture.html?consent=missing');
  await page.getByRole('button', { name: 'Accept', exact: true }).click();
  await page.getByRole('button', { name: 'Run', exact: true }).click();
  await expect.poll(() => page.evaluate(() => (window as any).__solves.length)).toBe(1);
  const events = await page.evaluate(() => (window as any).__events.filter((event: any[]) => event[0] === 'event'));
  expect(events.length).toBe(info.project.name === 'default-off' ? 0 : 1);
  await page.getByRole('button', { name: 'Analytics preferences', exact: true }).click();
  await page.getByRole('button', { name: 'Reject', exact: true }).click();
  await page.evaluate(() => { (window as any).__events = []; });
  await page.getByRole('button', { name: 'Run', exact: true }).click();
  await expect.poll(() => page.evaluate(() => (window as any).__solves.length)).toBe(2);
  expect(await page.evaluate(() => (window as any).__events)).toEqual([]);
});

test('configured feedback handles transport failure and distinguishes sending from persistence', async ({ page }, info) => {
  test.skip(info.project.name !== 'feedback-service');
  await page.route('https://**', route => route.abort());
  await page.route('**/feedback-service', route => route.abort());
  await page.goto('/e2e/analytics/fixture.html?geometry=private#private-foil');
  await page.getByRole('button', { name: 'Send feedback' }).click();
  await page.getByRole('textbox', { name: 'Feedback message' }).fill('A feature suggestion');
  await page.getByRole('button', { name: 'Send', exact: true }).click();
  await expect(page.getByText('Something went wrong. Please try again.')).toBeVisible();
  await page.unroute('**/feedback-service');
  let body: any;
  await page.route('**/feedback-service', async route => {
    body = route.request().postDataJSON();
    await route.fulfill({ status: 200, body: '{}' });
  });
  await page.getByRole('button', { name: 'Send', exact: true }).click();
  await expect(page.getByText('delivery cannot be confirmed here.', { exact: false })).toBeVisible();
  expect(body.url).not.toContain('private');
  expect(body.message).toBe('A feature suggestion');
  const events = await page.evaluate(() => (window as any).__events);
  expect(events.at(-1)[1]).toBe('feedback_sent');
  expect(JSON.stringify(events)).not.toContain('A feature suggestion');
});
