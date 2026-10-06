import {defineConfig} from '@playwright/test';
const root = new URL('../../', import.meta.url).pathname;
export default defineConfig({testDir: '.', workers: 1, reporter: 'list',
  outputDir: 'test-results/solve-analytics', use: {channel: 'chrome'},
  projects: [{name: 'default-off', use: {baseURL: 'http://127.0.0.1:18997'}},
    {name: 'explicit-on', use: {baseURL: 'http://127.0.0.1:18998'}}],
  webServer: [false, true].map(on => ({cwd: root,
    command: `flexfoil-ui/node_modules/.bin/vite --config flexfoil-ui/e2e/analytics.vite.config.ts --host 127.0.0.1 --port ${on ? 18998 : 18997} --strictPort`,
    env: {VITE_SOLVE_RUN_ANALYTICS: on ? 'true' : 'false'},
    url: `http://127.0.0.1:${on ? 18998 : 18997}/e2e/analytics/fixture.html`, reuseExistingServer: false})),
});
