# Current architecture

Resource builds resolve entities before publication. Each immutable resource
version contains entities, relations, source payloads and a build manifest.
Named Parquet releases pin resource versions; API Latest follows independent
resource updates. PostgreSQL chooses its own explicit release on its monthly schedule.

| Owner | Responsibility |
| --- | --- |
| core | Published schemas, vocabulary, identities, manifest contracts and version validation |
| resolver | Runtime matching, normalization, identity library lookup, entity records and the native policy kernel |
| build | Inputs, hub export and the identity library, resource assembly and atomic publication |
| postgres | Pinned local/HTTPS input, DuckDB projection, COPY, relational schema/indexes and shared derivations |
| subsets | MetSigDB, network views, COSMOS and their shared transaction/checkpoint runner |
| api | Inventory snapshots, bounded DuckDB queries, exports and optional job adapters |
| web | Explorer, generated HTTP/vocabulary types and shipped documentation |
| client | Python queries over local or remote published Parquets |

The package dependency direction is postgres -> subsets -> core; subsets does
not import PostgreSQL build internals. Build depends on resolver and core.
Serving requires core and API only; a separately installed worker also uses build.
The Python client runs independently of the API, builder and resolver. The
serving tables are written by the build and scoped to one resource version; gene
and chemical grouping across resources happens at query time.
Public Python names remain `omnipath_*`; outer project folders use short names.
Only Python's extra `src` level was removed; Rust and Svelte use their usual layout.

PostgreSQL never repeats entity resolution. It stages the resolved normalized
tables in DuckDB, projects the established relational tables, streams COPY,
restores constraints/indexes, then commits the base. Shared derivations and each
complete subset have durable checkpoints. The standalone subsets command explicitly
rebuilds selected products. PostgreSQL finish resumes unfinished phases/products.
Both use one runner with release identity checks and a schema lock across product
commits. Products publish their tables and metadata atomically.

Frozen legacy source is an independent test oracle. The earlier record-layout
PostgreSQL loader has been removed; PostgreSQL and subsets contain only the
current implementation. Normal commands use the current layout.
Developer fixture tests and wheel checks are separate from production execution;
there is no mandatory full-release rebuild/rollback verification phase.

Artifacts can live on another server: PostgreSQL accepts an explicit HTTPS release
and file root with byte ranges. It hashes artifacts before loading and again before
the base checkpoint; projection reads add traffic. The API currently inventories a
local mounted artifact directory. Public download URLs need not be on that server.
The Parquet contract is the boundary between systems, not a shared database or cache.

Start from [the root workflow](../README.md), [entity resolution](../core_documentation/resolution.md),
[PostgreSQL](../packages/postgres/README.md) or [serving operations](../deploy/README.md).
The numbered [old pipeline notes](pipeline/README.md) are historical.
