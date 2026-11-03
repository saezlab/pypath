# Subsets

Subsets are curated products built inside a loaded PostgreSQL schema. They read
the normalized tables of exactly one release and use its published identities.
They run no resource builds and no identifier resolution.

```sh
make build-subsets POSTGRES_SCHEMA=release_2026_09 SUBSET_PRODUCTS=cosmos     # resume unfinished
make rebuild-subsets POSTGRES_SCHEMA=release_2026_09 SUBSET_PRODUCTS=cosmos   # rebuild explicitly
```

## Products

| Product | What it is | Key rules |
| --- | --- | --- |
| **Exclusion metabolites** | Ions, atoms, small molecules, xenobiotics and cofactors that mass spectrometry does not measure | Curated seeds plus ChEBI derivation. Variants that share a structure key get rows of their own. A repository file adds or vetoes members. Runs first |
| **MetSigDB** | Metabolite memberships in Reactome, WikiPathways, KEGG, MACdb and ClassyFire sets | Uses published chemical entities, including native fallbacks |
| **Network views** | The MetaLinksDB, LIANA and Reactions preset definitions | Registers recipes, not materialized networks. Keeps qualified statement keys separate. LIANA roles need matching payload IDs; transport needs catalytic and compartment evidence |
| **COSMOS** | Metabolite–enzyme edges from the shared reaction tables | GeneID catalysts stay GeneID. No protein is picked from a gene's UniProt aliases. Reversible reactions get both halves; unknown direction does not |

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
> [resolution](resolution.md#gene-and-product-level).

## How a product build works

All products share one runner (`omnipath_subsets.runner`), used by both
`omnipath-postgres --finish` and `omnipath-subsets build`.

- The target schema must contain exactly one release. Each product checks the
  release ID and manifest digest before and after its work.
- A product's tables and its `subset_build_metadata` row commit together.
  A failed product rolls back only itself.
- `finish` skips products already committed for this release; an explicit
  `rebuild` replaces the selected ones.
- The same schema advisory lock as the loader prevents overlapping work.

## Where the code is

| Product | Package |
| --- | --- |
| Exclusion metabolites | `packages/subsets/omnipath_subsets/exclusion_metabolites/` |
| MetSigDB | `packages/subsets/omnipath_subsets/metsigdb/` |
| Network views | `packages/subsets/omnipath_subsets/network_views/` |
| COSMOS | `packages/subsets/omnipath_subsets/cosmos/` |
| Shared runner | `packages/subsets/omnipath_subsets/runner.py` |
