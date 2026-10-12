/** @file Browser and security checks against the actual packaged static site. */
import { defineConfig, devices } from '@playwright/test';

export default defineConfig({
  testDir: './tests/web',
  testMatch: '**/*.spec.mjs',
  timeout: 35000,
  expect: { timeout: 8000 },
  workers: 2,
  reporter: 'list',
  use: {
    baseURL: process.env.GEODE_TEST_URL || 'http://127.0.0.1:4173',
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
  },
  projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
  webServer: process.env.GEODE_TEST_URL ? undefined : {
    command: 'python3 -m http.server 4173 --bind 127.0.0.1 --directory out/web',
    url: 'http://127.0.0.1:4173',
    reuseExistingServer: false,
    timeout: 10000,
  },
});
