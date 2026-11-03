# OmniPath Web Explorer (`packages/web/`)

SvelteKit explorer talking to the DuckDB Parquet API in `../api`.

Latest follows independent Parquet resource updates. The existing release
selector can still choose a named Parquet snapshot; PostgreSQL's monthly release
schedule does not control this explorer.

## Developing

```sh
pnpm install
API_SERVICE_URL=http://127.0.0.1:8085 pnpm run dev
```

The SvelteKit proxy forwards browser `/api/*` and `/app-api/*` requests to the configured API origin.

## Production preview

```sh
pnpm run build
pnpm run preview --port 5173
```

Or start API + web together from the repo root:

```sh
./run_all.sh 8085 5173
```

## Code organization

- `src/lib/api/types.generated.ts` is generated from the repository OpenAPI document. `contracts.ts` names transport schemas; `adapters.ts` normalizes optional values into UI records.
- `src/lib/domain/` owns entity identity, vocabulary, presentation, and UI enrichment types. `src/lib/types/` contains small compatibility/UI aliases.
- `src/lib/features/selection/` owns URL codecs, persistence, exact-key restoration, navigation batching, and scope state. The old `stores/` and `navigation/` paths only re-export its public interface.
- `src/lib/features/explorer/` owns cancellation and paging; replaced queries cannot publish initial results or append stale pages.
- `src/lib/features/admin/` owns admin transport, shared job types, and polling. Polls wait for completion before scheduling another run.
- `src/lib/server/proxy-core.ts` forwards only to the configured backend origin and rejects external redirects. Development and production use the same proxy.

A URL containing `selection`, `entities`, or `annotations` owns both selection lists. Missing or empty lists in such a URL mean empty selection. URLs without these parameters may restore local selection; generated share links always include explicit ownership.

## Validation

```sh
pnpm test
pnpm check
pnpm lint
pnpm format:check
pnpm build
pnpm exec playwright install --with-deps chromium
pnpm test:browser
```

Browser tests start an isolated server on port 4183 and mock API responses; no local dataset is required. To use an already installed Chromium executable, set `PLAYWRIGHT_CHROMIUM_EXECUTABLE`.

Regenerate transport types after updating the root OpenAPI document with `pnpm generate:types`. Generated types and Biolink files are excluded from formatting; use `pnpm format` for maintained source files.
