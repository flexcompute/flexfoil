import {test, expect} from '@playwright/test';
for (const mode of ['single_alpha', 'single_cl', 'polar', 'sweep_1d', 'sweep_2d']) {
  test(`${mode}: analytics opt-in preserves solve dispatch`, async ({page}, info) => {
    await page.route('https://**', route => route.abort());
    await page.goto(`/e2e/analytics/fixture.html?mode=${mode}`);
    const button = page.getByRole('button', {name: mode === 'polar' ? 'α Polar' : mode.startsWith('sweep') ? 'Generate Sweep' : 'Run', exact: true});
    await expect(button).toBeEnabled();
    await button.click();
    await expect.poll(() => page.evaluate(() => (window as any).__solves.length)).toBeGreaterThan(0);
    await expect(button).toBeEnabled();
    const events = await page.evaluate(() => (window as any).__events);
    const on = info.project.name === 'explicit-on';
    if (on) {
      expect(events).toEqual([['event', 'solve_run', {solve_mode: mode, solver_mode: 'inviscid', n_panels: 24}]]);
    } else expect(events).toEqual([]);
    const solves = await page.evaluate(() => (window as any).__solves);
    if (mode === 'single_cl') expect(solves.at(-1).cl).toBeCloseTo(0.5, 2);
    else expect(solves).toEqual((mode === 'single_alpha' ? [0] : mode === 'sweep_2d' ? [0, 1, 0, 1] : [0, 1])
      .map(alpha => ({alpha, cl: alpha * 0.1})));
    await info.attach('solve-dispatch', {body: JSON.stringify({mode, on, solves}), contentType: 'application/json'});
  });
}
