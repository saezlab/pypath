import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './test/browser',
  fullyParallel: true,
  use: {
    baseURL: 'http://127.0.0.1:4183',
    headless: true,
    launchOptions: { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE || undefined },
  },
  webServer: {
    command: 'pnpm dev --host 127.0.0.1 --port 4183',
    url: 'http://127.0.0.1:4183/explore',
    env: { API_SERVICE_URL: 'http://127.0.0.1:9' },
    reuseExistingServer: !process.env.CI,
  },
});
