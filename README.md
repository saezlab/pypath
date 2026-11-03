# OmniPath build

Build resolved, independently versioned Parquet resources from pypath `inputs_v2`.
Serve those files through the API and explorer, or load an explicitly pinned
release into PostgreSQL. Entity resolution happens before Parquet publication.
API Latest and PostgreSQL's monthly snapshots have independent schedules.

## Setup

For development, install Python 3.11+, uv, Rust, Node.js 22 and pnpm, then run:

```sh
make setup
make help
```

Setup uses the locked dependencies and the committed pypath submodule revision; the
resolver's Rust extension is built by `uv sync`. A serving container needs no resolver,
reference rebuild or sibling prototype checkout.

## Run existing data

```sh
make dev DATA_ROOT=/absolute/path/to/existing/data
# API: http://127.0.0.1:8085; explorer: http://127.0.0.1:5173
```

For containers, set an external artifact directory in an environment file:

```sh
cp deploy/serving.env.example deploy/serving.env
# Edit DATA_DIR and the public URLs in deploy/serving.env.
make serving-build COMPOSE_ENV=deploy/serving.env
make serving-up COMPOSE_ENV=deploy/serving.env
make serving-status COMPOSE_ENV=deploy/serving.env
```

These commands start only the read-only API, web app and file service. Existing
Parquets, references and caches can remain outside the checkout.
See [deployment](deploy/README.md) for HTTPS routing, parallel cutover and rollback.

## Build and publish resources

For a small offline example, use `make sample`; it processes at most 20 fixture
records under `data/migration-smoke/`. It does not download source datasets.

A real build uses an existing identity library:

```sh
OMNIPATH_LIBRARY_DIR=/path/to/reference/identity/<fingerprint> \
  make build SOURCE=signor VERSION=1.0.0 MAX_RECORDS=20 DATA_ROOT=data
```

The limit applies per dataset; a source parser may still download its full input.
Versions are immutable. Prepare hubs and the identity library only when needed:

```sh
make hubs HUB_ARGS="--max-records 0"   # every hub, uncapped
make reference                         # hub indexes, kv stores, identity library
```

`make reference` indexes every hub in `HUBS_DIR` (default `$(DATA_ROOT)/reference/hubs`)
into `HUB_INDEX_ROOT` (`.../hub-index`), reusing an index while its hub export is
unchanged, then writes the identity library to `IDENTITY_DIR/<fingerprint>`
(`.../identity`). `IDENTITY_HUBS` restricts the hubs; `HUB_INDEX_ARGS`, `HUB_KV_ARGS`
and `IDENTITY_ARGS` pass memory and thread limits through. See the
[identity layer guide](packages/build/REFERENCE.md), the
[resource pipeline](packages/build/README.md),
[entity resolution](core_documentation/resolution.md) and
[orchestration guide](packages/build/ORCHESTRATION.md).

A named snapshot lists exact resource versions:

```json
{"schema_version": 1, "version": "2026.09", "resources": {"signor": "1.0.0"}}
```

Publish it from existing local artifacts with
`make publish-release RELEASE_MANIFEST=release.json DATA_ROOT=data`.
Publication validates availability and adds a cached taxonomy reference when
provisioned. It runs no source build and needs no separate index step: the
serving tables are part of each resource version. API Latest continues selecting
each resource's newest version.

## PostgreSQL and products

Use PostgreSQL 18 with `roaringbitmap` and `pg_trgm`; the image definition is in
[postgres/](postgres/). Set `OMNIPATH_DATABASE_URL` for the intended database.
Load a new schema from a local directory or an HTTPS endpoint:

```sh
make load-postgres RELEASE_MANIFEST=release.json DATA_ROOT=/path/to/data \
  POSTGRES_SCHEMA=release_2026_09

make load-postgres \
  RELEASE_MANIFEST=https://data.example.org/releases/2026.09.json \
  DATA_ROOT=https://data.example.org POSTGRES_SCHEMA=release_2026_09
```

The release is explicit; the loader never substitutes Latest or repeats entity
resolution. HTTP input requires byte-range support. Checksums, schemas, sizes
and row counts are checked. DuckDB stages entity/relation inputs locally once,
then projects relational tables and streams CSV bytes through PostgreSQL COPY.
No persistent downloaded-Parquet cache is kept. Working storage can be selected
with `POSTGRES_ARGS="--temp-directory /path/to/spool --memory-limit 4GB"`.

The validated base commits before downstream work. Completed phases and products
are durable; retrying skips them. Raw source payloads and nested record JSON stay
in Parquet. Main's relational tables and indexes hold the scientific data.

```sh
# Optional: import only the base first.
make load-postgres RELEASE_MANIFEST=release.json DATA_ROOT=data \
  POSTGRES_SCHEMA=release_2026_09 POSTGRES_ARGS=--base-only
# Finish remaining derivations and products without reading Parquet or repeating COPY.
make finish-postgres POSTGRES_SCHEMA=release_2026_09
make build-subsets POSTGRES_SCHEMA=release_2026_09 SUBSET_PRODUCTS=cosmos
```

See the [PostgreSQL guide](packages/postgres/README.md) for schema,
checkpoint and product details. The [subsets package](packages/subsets/README.md)
owns MetSigDB, network views and COSMOS. Both commands use its shared runner:
`make build-subsets` resumes unfinished products; `make rebuild-subsets` explicitly
rebuilds the selected products, committing each complete product. The earlier
record-layout PostgreSQL loader has been removed; subsets contains only the
current product implementations.

## Packages and checks

| Package | Responsibility |
| --- | --- |
| [core](packages/core/README.md) | Schemas, vocabulary, IDs and versions |
| [resolver](packages/resolver/README.md) | Reference indexes and identifier matching |
| [build](packages/build/README.md) | Reference preparation and resource publication |
| [postgres](packages/postgres/README.md) | Parquet projection, relational layout, indexes and shared derivations |
| [subsets](packages/subsets/README.md) | MetSigDB, network views, COSMOS and durable product publication |
| [api](packages/api/README.md) | REST API, inspection and administration; its Parquet queries are `omnipath_api.queries` |
| [web](packages/web/README.md) | Explorer and API proxy |
| [query](packages/query/README.md) | The query layer: SQL macros (`sql/`) and the Python and R clients that install resources and query them locally |

The root Makefile delegates to package commands; implementation stays inside
packages. `make check`, `make test`, `make check-web` and `make test-api` cover
local checks. `make test-native` checks Rust policy; `make check-generated` checks
generated contracts; `make test-wheels` checks isolated installations. PostgreSQL
fixtures require explicit `make test-postgres` and create a disposable local cluster.
CI runs developer checks without imposing a production verification phase.
The [versioning policy](VERSIONING.md) describes independent publication schedules.
The Python client queries the same published Parquets independently of the server.
`legacy/` is excluded from the active workspace.

## Architecture and dependency policy

Python projects use `packages/<short-name>/omnipath_<name>/`, with explicit
package discovery and wheel checks. Rust and Svelte keep their own `src` folders.
See [current architecture](docs/architecture.md) and [refactoring log](docs/reviews/refactoring-progress-20261003.md).

Only pypath is a runtime source submodule. The lockfile selects published cachedir
and dlmachine releases; their optional checkouts and omnipath-utils are not installed
or copied into images. The cachedir workaround remains scoped to the pinned release
until its upstream fix is available in a selected version. Proxy settings are honored
as supplied by the launcher; library calls do not rewrite the process environment.

Role-specific setup: `make setup-serving`, `make setup-build`, or
`make setup-postgres`. The full `make setup` remains available for maintainers.
For a local database, copy `deploy/postgres-dev.env.example` to a private
`deploy/postgres.env`, set a password, then run
`make postgres-up POSTGRES_ENV=deploy/postgres.env`. The default bind is loopback.
Large server examples are separate and must be sized to that host.
