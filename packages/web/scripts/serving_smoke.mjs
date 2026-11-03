/**
 * Real browser smoke against the existing capped migration SIGNOR fixture.
 * Start the API/web stack first; no mock routes, source builds or downloads run.
 */
import assert from "node:assert/strict";
import { createRequire } from "node:module";

const require = createRequire(
  new URL("../packages/web/package.json", import.meta.url),
);
const { chromium, expect } = require("@playwright/test");
const base = process.env.OMNIPATH_WEB_URL || "http://127.0.0.1:5173";
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE || undefined,
  headless: true,
});
try {
  const page = await browser.newPage({
    viewport: { width: 1440, height: 1000 },
  });
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  const catalog = await page.request.get(
    base + "/app-api/resources?shape=svelte",
  );
  assert.equal(catalog.status(), 200);
  const resource = (await catalog.json()).resources.find(
    (row) => row.resource_id === "signor",
  );
  assert.ok(resource, "The capped SIGNOR fixture must be served");

  await page.goto(base + "/explore?q=TP53&source=signor");
  await page.getByRole("button", { name: /^TP53 UniProt P04637/ }).click();
  const details = page.getByRole("dialog");
  await expect(
    details.getByRole("heading", { name: "TP53", exact: true }),
  ).toBeVisible();
  await details.getByRole("tab", { name: /Identifiers/ }).click();
  await expect(
    details.getByText("P04637", { exact: true }).first(),
  ).toBeVisible();
  await details.getByRole("tab", { name: /Relationships/ }).click();
  // This binary-interaction fixture has no participant/composition relationships.
  await expect(
    details.getByText("No participant or composition relationships recorded."),
  ).toBeVisible();

  await page.goto(base + "/explore?q=TP53&source=signor");
  const taxonomy = page.getByRole("region", {
    name: "Taxonomy panel",
    exact: true,
  });
  const filteredSearch = page.waitForResponse((response) => {
    const url = new URL(response.url());
    if (!url.pathname.endsWith("/entities/search")) return false;
    const filters = JSON.parse(url.searchParams.get("filters") || "{}");
    return (filters.ncbi_tax_id || filters.taxonomy_ids || []).includes("9606");
  });
  await taxonomy.getByRole("checkbox", { name: /9606/ }).check();
  assert.equal((await filteredSearch).status(), 200);
  await expect(
    page.getByRole("button", { name: /^TP53 UniProt P04637/ }),
  ).toBeVisible();

  await page.goto(base + "/explore?tab=relations&source=signor");
  await expect(page.locator("tbody tr")).toHaveCount(2);
  const evidenceResponse = page.waitForResponse((response) =>
    new URL(response.url()).pathname.endsWith("/evidence"),
  );
  await page.locator("tbody tr").first().click();
  const evidence = await (await evidenceResponse).json();
  assert.ok(evidence.evidence.length > 0 && evidence.evidence.length <= 20);
  await expect(
    page
      .getByRole("dialog")
      .getByText("Observation 1 of " + evidence.evidence.length),
  ).toBeVisible();

  const exported = await page.request.post(base + "/app-api/export", {
    data: { resources: ["signor"], filters: {} },
  });
  assert.equal(exported.status(), 200);
  // Exports are Parquet: the file starts with its magic bytes.
  assert.equal((await exported.body()).subarray(0, 4).toString(), "PAR1");
  assert.deepEqual(errors, []);
  if (process.env.OMNIPATH_SMOKE_SCREENSHOT) {
    await page.screenshot({
      path: process.env.OMNIPATH_SMOKE_SCREENSHOT,
      fullPage: true,
    });
  }
  console.log(
    JSON.stringify(
      {
        resource: resource.resource_id,
        version: resource.version,
        exportedRelations: rows.length,
        sourceEvidence: 20,
        checks: [
          "search",
          "identifiers",
          "relationships",
          "taxonomy_filter",
          "relation_evidence",
          "proxy_export",
        ],
      },
      null,
      2,
    ),
  );
} finally {
  await browser.close();
}
