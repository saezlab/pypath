import { test, expect } from '@playwright/test';
import { defaultLayout } from '../../src/lib/components/workspace/layout';

for (const savedLayout of [false, true]) {
  test(`taxonomy stays visible and filters searches with ${savedLayout ? 'an earlier saved' : 'a fresh'} layout`, async ({
    page,
  }) => {
    if (savedLayout) {
      await page.addInitScript(
        (layout) => {
          localStorage.setItem('omnipath-workspace-v3-entities', JSON.stringify(layout));
        },
        defaultLayout(['results', 'entity_types', 'sources', 'taxonomy_ids']),
      );
    }
    let requestedTaxa: string[] = [];
    await page.route('**/app-api/**', async (route) => {
      const url = new URL(route.request().url());
      if (url.pathname.endsWith('/entities/scoped-facets')) {
        await route.fulfill({
          json: [
            { facetName: 'taxonomy_id', facetValue: '9606', facetLabel: 'Human', scopedCount: 7 },
          ],
        });
      } else if (url.pathname.endsWith('/entities/search')) {
        const filters = JSON.parse(url.searchParams.get('filters') || '{}');
        requestedTaxa = filters.taxonomy_ids || filters.ncbi_tax_id || [];
        await route.fulfill({ json: { entities: [], nextCursor: null } });
      } else {
        await route.fulfill({ json: { entities: [], nextCursor: null } });
      }
    });
    await page.goto('/explore');
    const panel = page.getByRole('region', { name: 'Taxonomy panel', exact: true });
    await expect(panel).toBeVisible();
    await panel.getByRole('checkbox', { name: /Human/ }).check();
    await expect.poll(() => requestedTaxa).toEqual(['9606']);
    await expect(panel.getByRole('checkbox', { name: /Human/ })).toBeChecked();
    await page.reload();
    await expect(panel).toBeVisible();
  });
}
