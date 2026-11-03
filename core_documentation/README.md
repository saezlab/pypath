# Core documentation

A short, current description of the core of OmniPath build: what the data
means, how entities are resolved, how PostgreSQL and the subsets are built,
and the decisions behind them.

| Page | Read it to learn |
| --- | --- |
| [Overview](overview.md) | The stages, versions and package ownership |
| [Glossary](glossary.md) | The terms we use, key terms first |
| [Entities and relations](entities.md) | What each Parquet row and field means |
| [Entity resolution](resolution.md) | How the reference library is built and how matching decides |
| [PostgreSQL](postgres.md) | How a release is loaded and mapped to tables |
| [Subsets](subsets.md) | MetSigDB, network views and COSMOS |

Decisions are written as **Decision** blocks on the page they belong to, each
with its reason.

## Reading it

The Markdown files are the source and render directly on GitHub. For the HTML
viewer with navigation and glossary pop-ups, run:

```bash
make docs
```

and open <http://127.0.0.1:8090>.

## Keeping it current

- These pages describe the system as it is now. History, experiments and dated
  measurements stay in [`docs/`](../docs/README.md).
- A pull request that changes a schema, a resolution rule, the PostgreSQL
  mapping or a product rule updates the matching page here.
- A changed decision is edited in place in its **Decision** block, with the
  new date; a new term goes into [glossary.md](glossary.md).
- Keep each page short. Link to package READMEs for operational detail.
