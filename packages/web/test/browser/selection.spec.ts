import { test, expect, type Page } from '@playwright/test';

const key = 'a'.repeat(64);
async function mockApi(page: Page, scopes: string[][] = []) {
  await page.route('**/app-api/**', async (route) => {
    const path = new URL(route.request().url()).pathname;
    const body = route.request().postDataJSON() ?? {};
    let response: unknown = [];
    if (path.endsWith('/selection/scope')) {
      scopes.push(body.entityPks ?? []);
      response = {
        entityPks: body.entityPks ?? [],
        seedEntityPks: body.entityPks ?? [],
        termEntityPks: [],
        ontologyTermIds: body.annotationTermIds ?? [],
        criteriaCount: 1,
        expandedEntityCount: 0,
      };
    } else if (path.includes('/entities/')) {
      response = path.includes('facets') ? [] : { entities: [], nextCursor: null };
    } else if (path.endsWith('/relations/search'))
      response = { rows: [], relations: [], total: 0, nextCursor: null };
    await route.fulfill({ json: response });
  });
}

test('a fresh shared SHA-256 selection hydrates scope without cached metadata', async ({
  page,
}) => {
  const scopes: string[][] = [];
  await mockApi(page, scopes);
  await page.goto(`/selection?entities=${key}&annotations=`);
  await expect(page.getByRole('button', { name: /Open Selection/ })).toBeVisible();
  await expect.poll(() => scopes.some((ids) => ids.includes(key))).toBe(true);
});

test('clear removes both lists atomically, persists, and survives reload', async ({ page }) => {
  await mockApi(page);
  await page.goto(`/selection?entities=${key}&annotations=GO:0001`);
  await page.getByRole('button', { name: /Open Selection/ }).click();
  await page.getByRole('button', { name: 'Clear all', exact: true }).click();
  await expect.poll(() => new URL(page.url()).searchParams.get('entities')).toBe('');
  await expect.poll(() => new URL(page.url()).searchParams.get('annotations')).toBe('');
  await expect(page.getByRole('button', { name: /Open Selection/ })).toHaveCount(0);
  await page.reload();
  await expect(page.getByRole('button', { name: /Open Selection/ })).toHaveCount(0);
  expect(
    await page.evaluate(() => JSON.parse(localStorage.getItem('omnipath-selection-ids')!)),
  ).toEqual([]);
  expect(
    await page.evaluate(() =>
      JSON.parse(localStorage.getItem('omnipath-selection-annotation-ids')!),
    ),
  ).toEqual([]);
});

test('explicit empty URL lists override unrelated browser selections', async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.setItem('omnipath-selection-ids', JSON.stringify(['b'.repeat(64)]));
    localStorage.setItem('omnipath-selection-annotation-ids', JSON.stringify(['STALE:1']));
  });
  await mockApi(page);
  await page.goto('/selection?entities=&annotations=');
  await expect(page.getByRole('button', { name: /Open Selection/ })).toHaveCount(0);
});

test('relation text search includes entity pages beyond the first 200 matches', async ({
  page,
}) => {
  let requestedRelations: string[] = [];
  await mockApi(page);
  await page.route('**/app-api/entities/search**', async (route) => {
    const cursor = new URL(route.request().url()).searchParams.get('cursor');
    const start = cursor ? 200 : 0;
    const count = cursor ? 3 : 200;
    await route.fulfill({
      json: {
        entities: Array.from({ length: count }, (_, index) => ({
          entityPk: String(start + index),
        })),
        nextCursor: cursor ? null : { entityPk: '199', relationCount: 0 },
      },
    });
  });
  await page.route('**/app-api/relations/search', async (route) => {
    requestedRelations = route.request().postDataJSON().filters.entity_ids || [];
    await route.fulfill({ json: { rows: [], relations: [], total: 0, nextCursor: null } });
  });
  await page.goto('/explore?tab=relations&q=protein');
  await expect.poll(() => requestedRelations.length).toBe(203);
  expect(requestedRelations).toContain('202');
});

test('annotation navigation cancels an earlier request and keeps the newest results', async ({
  page,
}) => {
  await mockApi(page);
  let started = false;
  let releaseOld!: () => void;
  const old = new Promise<void>((resolve) => {
    releaseOld = resolve;
  });
  await page.route('**/app-api/ontology/scoped-search', async (route) => {
    const query = route.request().postDataJSON().query;
    if (query === 'old') {
      started = true;
      await old;
    }
    const label = query === 'old' ? 'Outdated term' : 'Newest term';
    await route
      .fulfill({
        json: [
          {
            entityPk: `${query}-key`,
            termId: `${query}:1`,
            label,
            ontologyPrefix: query,
            ontologyId: query,
            definition: null,
            synonyms: [],
            sources: [],
            annotatedEntityCount: 1,
            annotatedRelationCount: 1,
          },
        ],
      })
      .catch(() => {});
  });
  await page.goto(`/selection?tab=ontology&entities=${key}&annotations=&q=old`);
  await expect.poll(() => started).toBe(true);
  const search = page.getByPlaceholder('Search associated ontology terms…');
  await search.fill('new');
  await search.press('Enter');
  await expect(page.getByText('Newest term', { exact: true })).toBeVisible();
  releaseOld();
  await expect(page.getByText('Outdated term', { exact: true })).toHaveCount(0);
  await expect(page.getByText('Newest term', { exact: true })).toBeVisible();
});
