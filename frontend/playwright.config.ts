import { defineConfig } from '@playwright/test'
export default defineConfig({
  testDir: './tests',
  timeout: 120_000,
  expect: { timeout: 30_000 },
  retries: 0,
  workers: 1,
  use: { baseURL: 'http://127.0.0.1:5173', viewport: { width: 1400, height: 1000 }, screenshot: 'only-on-failure' },
  webServer: { command: 'npx vite --port 5173 --strictPort --host 127.0.0.1', url: 'http://127.0.0.1:5173', reuseExistingServer: true, timeout: 60_000 },
  reporter: [['list'], ['html', { open: 'never', outputFolder: 'playwright-report' }]],
})
