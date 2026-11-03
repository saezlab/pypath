# OmniPath

OmniPath builds resolved, independently versioned Parquet resources from its resource
input modules, serves them through an API and explorer, queries them locally through
Python and R clients, and loads pinned releases into PostgreSQL. Entity resolution
happens before Parquet publication. API Latest and PostgreSQL's monthly snapshots have
independent schedules ([versioning](docs/versioning.md)).

## Packages

| Package | What it does |
| --- | --- |
| [resources](packages/resources/README.rst) | Resource input modules (`omnipath_resources`, formerly pypath): download and parse each source into typed entities and relations |
| [core](packages/core/README.md) | Schemas, vocabulary, identifiers and versions shared by every package |
| [resolver](packages/resolver/README.md) | Reference indexes and identifier matching (Rust extension) |
| [build](packages/build/README.md) | Builds the identity library and each resource version; publishes releases |
| [api](packages/api/README.md) | REST API over the published Parquets; the serving deployment |
| [web](packages/web/README.md) | The explorer and its API proxy |
| [query](packages/query/README.md) | SQL macros and the Python and R clients that install resources and query them locally |
| [postgres](packages/postgres/README.md) | Loads a pinned release into PostgreSQL: relational layout, indexes, shared derivations |
| [subsets](packages/subsets/README.md) | Products built in PostgreSQL: MetSigDB, network views, COSMOS |
| [utils](packages/utils/README.md) | ID translation, taxonomy, orthology and reference lists (`omnipath_utils`) |

Data flows from resources through build (with resolver and core) to published
Parquet, which api, web and query read and postgres and subsets load. The
[documentation index](docs/README.md) lists each package's guides.

## Layout

Each package keeps what belongs to it: its README and code, `tests/`, and where it
has them `docs/` (guides, plans, dated reports, history), `scripts/` (tools run from
the repository root) and `deploy/` (containers and services). The root holds the
workspace: `pyproject.toml` and `uv.lock`, the `Makefile`, [`docs/`](docs/README.md)
for what spans packages (overview, glossary, architecture, versioning), and
`scripts/` for the two workspace-wide tools: the isolated wheel check and the build
constraints it shares with the container images.

## Setup

For development, install Python 3.11+, uv, Rust, Node.js 22 and pnpm, then run:

```sh
make setup
make help
```

Setup uses the locked dependencies; every Python package is a workspace member and
the resolver's Rust extension is built by `uv sync`. Role-specific setup:
`make setup-serving`, `make setup-build` or `make setup-postgres`.

## Common tasks

| Task | Command | Details |
| --- | --- | --- |
| Run the API and explorer on existing data | `make dev DATA_ROOT=…` | [api](packages/api/README.md) |
| Serve in containers | `make serving-up COMPOSE_ENV=…` | [deployment](packages/api/deploy/README.md) |
| Build resources, references and releases | `make sample`, `make build`, `make reference`, `make publish-release` | [build](packages/build/README.md) |
| Load PostgreSQL and build products | `make load-postgres`, `make finish-postgres`, `make build-subsets` | [postgres](packages/postgres/README.md), [subsets](packages/subsets/README.md) |
| Run a local build database | `make postgres-up` | [postgres](packages/postgres/README.md) |
| Read the core documentation | `make docs` | [docs](docs/README.md) |

## Checks

`make check`, `make test`, `make check-web` and `make test-api` cover local checks.
`make test-resources` and `make test-utils` run those packages' suites with their
own settings. `make test-native` checks Rust policy, `make check-generated` checks
generated contracts and `make test-wheels` checks isolated installations. PostgreSQL
fixtures require `make test-postgres` and create a disposable local cluster. CI runs
the developer checks.

## Dependency policy

Python packages use `packages/<short-name>/omnipath_<name>/`, with explicit package
discovery and wheel checks; Rust and Svelte keep their own `src` folders. See the
[architecture](docs/architecture.md). The lockfile selects published cachedir and
dlmachine releases; the cachedir workaround remains scoped to the pinned release
until its upstream fix is available in a selected version. omnipath-utils is not
copied into images. Proxy settings are honored as supplied by the launcher; library
calls do not rewrite the process environment.
