import {defineConfig} from '@playwright/test';
const root = new URL('../../', import.meta.url).pathname;
export default defineConfig({testDir: '.', workers: 1, reporter: 'list',
  outputDir: 'test-results/solve-analytics', use: {channel: 'chrome'},
  projects: [{name: 'default-off', testIgnore: 'actual-ui.spec.ts', use: {baseURL: 'http://127.0.0.1:18997'}},
    {name: 'explicit-on', testIgnore: 'actual-ui.spec.ts', use: {baseURL: 'http://127.0.0.1:18998'}},
    {name: 'feedback-service', testIgnore: 'actual-ui.spec.ts', use: {baseURL: 'http://127.0.0.1:18999'}},
    {name: 'actual-app', testMatch: 'actual-ui.spec.ts', use: {baseURL: 'http://127.0.0.1:19000'}}],
  webServer: [...[18997, 18998, 18999].map(port => ({cwd: root,
    command: `flexfoil-ui/node_modules/.bin/vite --config flexfoil-ui/e2e/analytics.vite.config.ts --host 127.0.0.1 --port ${port} --strictPort`,
    env: {VITE_SOLVE_RUN_ANALYTICS: port === 18997 ? 'false' : 'true',
      VITE_FEEDBACK_SHEET_URL: port === 18999 ? '/feedback-service' : ''},
    url: `http://127.0.0.1:${port}/e2e/analytics/fixture.html`, reuseExistingServer: false})), {
      command: 'npm run dev -- --host 127.0.0.1 --port 19000 --strictPort',
      cwd: `${root}flexfoil-ui`, env: {VITE_SOLVE_RUN_ANALYTICS: 'true', VITE_FEEDBACK_SHEET_URL: ''},
      url: 'http://127.0.0.1:19000', reuseExistingServer: false,
    }],
});
