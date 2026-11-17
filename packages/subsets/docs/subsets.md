# Subsets

Each subset is a curated product of exactly one published release. It uses the
published identities. It runs no resource builds and no identifier resolution. There are
four products. We build two as **Parquet subsets**, once, from the published release.
We build two in a loaded PostgreSQL schema.

```sh
make build-subsets RELEASE=2026.11                                              # the Parquet subsets
make rebuild-subsets POSTGRES_SCHEMA=release_2026_09 SUBSET_PRODUCTS=cosmos     # PostgreSQL products
```

Order: `make build`, `make publish-release`, `make build-subsets`, `make release-db`
(copies the subsets into the serving database), `make load-postgres` (loads them from
Parquet). The layout is `<data_root>/subsets/<release>/<product>/` with a
`subset_manifest.json` that pins the release. A subset is immutable.

## Products

| Product | Built as | What it is | Key rules |
| --- | --- | --- | --- |
| **Exclusion metabolites** | Parquet | Ions, atoms, small molecules, xenobiotics and cofactors that mass spectrometry does not measure | Curated seeds plus ChEBI derivation. Variants that share a structure key get rows of their own. A repository file adds or vetoes members. Runs first |
| **MetSigDB** | Parquet | Metabolite memberships in Reactome, WikiPathways, KEGG, MACdb and ClassyFire sets | Uses published chemical entities, including native fallbacks |
| **Network views** | PostgreSQL | The MetaLinksDB, LIANA and Reactions preset definitions | Registers recipes, not materialized networks. Keeps qualified statement keys separate. LIANA roles need matching payload IDs; transport needs catalytic and compartment evidence |
| **COSMOS** | PostgreSQL | Metabolite–enzyme edges from the shared reaction tables | GeneID catalysts stay GeneID. No protein is picked from a gene's UniProt aliases. Reversible reactions get both halves; unknown direction does not |

> **Decision: exclusion is a read-time choice.** The exclusion table lists
> entities and removes no membership row. A reader leaves listed entities out
> with an anti-join. An unfiltered read does not change, and COSMOS output does
> not change. *Why:* what counts as unmeasurable depends on the question. Only a
> rebuild can restore a removed row.

> **Decision: MetSigDB uses published chemical entities.** It includes
> native fallback identities. *Why:* the Parquets do not carry the legacy
> matched flag, and we chose not to re-resolve to recover it.

> **Decision: COSMOS keeps gene catalysts and separate source events.** GeneID
> catalysts keep their GeneID label and `entrez` namespace; no protein is picked
> from UniProt aliases. Reactions from different sources stay separate events.
> Unknown direction does not create reverse edges, and gene–protein–reaction
> rules do not become catalysis. *Why:* no invented participants, as in
> [resolution](../../resolver/docs/resolution.md#gene-and-product-level).

## How a product build works

`omnipath-subsets build` writes the Parquet subsets. It follows the layout rules: write to
a temporary directory, validate, then move into place. The digest of the release
manifest must match. `omnipath-postgres` loads the subsets. The PostgreSQL products share
one runner (`omnipath_subsets.runner`). Both `omnipath-postgres --finish` and
`omnipath-subsets build-postgres` use it.

- The target schema must contain exactly one release. Each product checks the
  release ID and manifest digest before and after its work.
- A product's tables and its `subset_build_metadata` row commit together.
  A failed product rolls back only itself.
- `finish` skips products already committed for this release. An explicit
  `build-postgres` replaces the selected built products.
- The release is complete when every loaded and built product is committed.
- The same schema advisory lock as the loader prevents overlapping work.

## Where the code is

| Product | Package |
| --- | --- |
| Exclusion metabolites | `packages/subsets/omnipath_subsets/exclusion_metabolites/` (DuckDB SQL in `duckdb_sql/`) |
| MetSigDB | `packages/subsets/omnipath_subsets/metsigdb/` (DuckDB SQL in `duckdb_sql/`) |
| Subset layout | `packages/core/omnipath_core/subset_layout.py` |
| Release checks | `packages/subsets/omnipath_subsets/checks.py` |
| Loading into PostgreSQL | `packages/postgres/omnipath_postgres/subset_loading.py` |
| Network views | `packages/subsets/omnipath_subsets/network_views/` |
| COSMOS | `packages/subsets/omnipath_subsets/cosmos/` |
| Shared runner | `packages/subsets/omnipath_subsets/runner.py` |
