# Parquet migration: API and explorer

Integrated the prototype API and web explorer from revision
`86f978e39bd1f7a40180714acf03747f3426e87c` into the migration repository.

## Behavior preserved

Resources publish independently. API and explorer Latest selects the newest
published version of each resource, including newly added resources. Existing
optional named Parquet snapshots stay pinned. PostgreSQL continues to load its
own fixed monthly snapshot; resource publication neither changes that database
nor changes the API default to a PostgreSQL release.

The serving backend remains the prototype's DuckDB implementation. The imported
API covers resource inspection, entities, relations, evidence, facets, ontology,
selection and exports. The SvelteKit explorer proxies requests to that API.
Read-only serving needs neither pypath, the build package nor a resolver.
Administration and the existing separate worker remain explicitly configured.

## Integration work

- Added the API to the Python workspace and locked its serving dependencies.
- Imported the web package with its existing frozen pnpm lockfile.
- Carried over OpenAPI, generated web types and raw documentation assets.
- Added `make api`, `make web`, `make test-api` and a local `run_all.sh` launcher.
- Added schema-3 query/export checks and an independent-update regression that
  compares all PostgreSQL base rows and complete release metadata before/after.
- Renamed the imported API release test module to avoid colliding with the
  existing PostgreSQL release test module in combined collection.
- Added a repeatable real-browser check in `scripts/serving_smoke.mjs`.

## Verification

- Combined PostgreSQL, subsets, core and API suite: **324 passed**.
- Web unit tests: **26 passed**; mocked browser tests: **7 passed**.
- Svelte type checks: **0 errors, 0 warnings**.
- Python/web lint and formatting, production web build and OpenAPI contract:
  passed.
- Standalone API wheel: imports and read-only requests passed without builder,
  resolver or pypath installed; both new API integration tests passed against it.
- Real production explorer: search, identifiers, lazy composition relationships,
  taxonomy filtering, relation evidence and JSON export passed through the
  actual proxy and API. The existing offline SIGNOR fixture exports two
  qualified relations with 20 source evidence records.
- Local development launcher and proxy checked; temporary servers shut down.

No resource builds or upstream data downloads ran. New migration fixtures use
one or two original records per resource; inherited synthetic pagination tests
retain enough rows to exercise their page boundaries.

The production build retains the prototype's large JavaScript chunk warnings.
The tests also emit existing dependency deprecation and missing taxonomy-reference
warnings for tiny fixtures. These do not establish production performance or
full-release compatibility.

## Run locally

```sh
uv sync --frozen --all-packages
pnpm --dir packages/web install --frozen-lockfile
OMNIPATH_DATA_ROOT=data/migration-smoke ./run_all.sh
```

Open `http://127.0.0.1:5173/explore`. This reads the existing capped sample.
For another data root, change `OMNIPATH_DATA_ROOT`. The two optional script
arguments select API and web ports.

With a Playwright browser available, run the real fixture smoke in another
terminal:

```sh
# Set PLAYWRIGHT_CHROMIUM_EXECUTABLE if using an existing Chrome installation.
node scripts/serving_smoke.mjs
```

Use `OMNIPATH_WEB_URL` for another local web port. This check expects the
20-record SIGNOR migration fixture, not arbitrary production datasets.

PostgreSQL product endpoints, the Python client, full product-output comparisons,
deployment and meaningful backend performance benchmarks remain later work.
