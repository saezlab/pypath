# Documentation

This folder holds what spans packages. Everything else lives with its package:
each package's README, its `docs/` folder (guides, plans, dated reports and
history) and, where it has one, its `deploy/` folder.

## Core documentation

A short, current description of the system: what the data means, how entities are
resolved, how PostgreSQL and the subsets are built, and the decisions behind them.

| Page | Read it to learn |
| --- | --- |
| [Overview](overview.md) | The stages, versions and package ownership |
| [Glossary](glossary.md) | The terms we use, key terms first |
| [Entities and relations](../packages/core/docs/entities.md) | What each Parquet row and field means |
| [Entity resolution](../packages/resolver/docs/resolution.md) | How the reference library is built and how matching decides |
| [PostgreSQL](../packages/postgres/docs/postgres.md) | How a release is loaded and mapped to tables |
| [Subsets](../packages/subsets/docs/subsets.md) | MetSigDB, network views and COSMOS |
| [Architecture](architecture.md) | Package boundaries and dependencies |
| [Versioning](versioning.md) | Resource versions and OmniPath releases |

Decisions are written as **Decision** blocks on the page they belong to, each with
its reason. The Markdown renders directly on GitHub. For the HTML viewer with
navigation and glossary pop-ups, run `make docs` and open
<http://127.0.0.1:8090/docs/>.

These pages describe the system as it is now. A change to a schema, a resolution
rule, the PostgreSQL mapping or a product rule updates the matching page; a changed
decision is edited in place in its **Decision** block with the new date, and a new
term goes into the [glossary](glossary.md). Keep each page short and link to the
package guides for operational detail.

## Package guides

| Package | Guides |
| --- | --- |
| [resources](../packages/resources/README.rst) | [Contributing](../packages/resources/CONTRIBUTING.md), [inputs_v2](../packages/resources/omnipath_resources/inputs_v2/README.md), [docs](../packages/resources/docs/index.md), [master changes not on this branch](../packages/resources/docs/pypath-master-since-branch.md) |
| [core](../packages/core/README.md) | [Entities and relations](../packages/core/docs/entities.md), [taxonomy](../packages/core/TAXONOMY.md) |
| [resolver](../packages/resolver/README.md) | [Entity resolution](../packages/resolver/docs/resolution.md), [resolution and molecular forms](../packages/resolver/docs/entity-resolution-and-molecular-forms.md); the earlier compact reference: [methods](../packages/resolver/docs/reference-resolver.md), [compiler semantics](../packages/resolver/docs/reference-compiler-semantics.md) |
| [build](../packages/build/README.md) | [Identity layer](../packages/build/REFERENCE.md), [identity layer spec](../packages/build/docs/identity-layer-spec.md), [orchestration](../packages/build/ORCHESTRATION.md), [replaying published resources](../packages/build/docs/replaying-published-resources.md), [build log](../packages/build/docs/build-log.md), [history](../packages/build/docs/history/README.md) |
| [api](../packages/api/README.md) | [Serving deployment](../packages/api/deploy/README.md), [nicesrv](../packages/api/deploy/nicesrv.md), [serving milestone](../packages/api/docs/parquet-serving-milestone.md) |
| [web](../packages/web/README.md) | [Explorer prototypes](../packages/web/docs/prototypes/) |
| [query](../packages/query/README.md) | [Python client](../packages/query/python/README.md) |
| [postgres](../packages/postgres/README.md) | [PostgreSQL](../packages/postgres/docs/postgres.md), [column map](../packages/postgres/docs/postgres-parquet-column-map.md), [parity with main](../packages/postgres/docs/postgres-main-parity-checklist.md), [build tuning](../packages/postgres/docs/build-tuning.md), [promotion](../packages/postgres/deploy/PROMOTE.md) |
| [subsets](../packages/subsets/README.md) | [Subsets](../packages/subsets/docs/subsets.md), [subsets milestone](../packages/subsets/docs/parquet-subsets-milestone.md) |
| [utils](../packages/utils/README.md) | [docs](../packages/utils/docs/index.md) |

Dated reports (`docs/reports/` in a package) record specific experiments and
deployments; their paths, producer revisions and numbers describe those runs.
