# OmniPath subsets

Build the exclusion metabolite list, MetSigDB memberships, network presets and
COSMOS reaction networks from one fixed PostgreSQL release loaded by omnipath-postgres. Products use published
identities and aliases; no resource builds or identifier resolution run.

## Build selected products

```sh
export OMNIPATH_DATABASE_URL='postgresql:///omnipath'
uv run --frozen omnipath-subsets build --schema release_2026_09
uv run --frozen omnipath-subsets build --schema release_2026_09 --products metsigdb
```

The command requires exactly one release in the target schema. It takes the
same schema advisory lock as the loader, rebuilds the selected products in one
transaction and records their counts, release ID and canonical manifest SHA-256
in `subset_build_metadata`. A failure rolls back every selected rebuild,
including changes to existing products. The Python entry point is
`omnipath_subsets.build.build_subsets`.

For long product builds, add `--checkpoint-products` (Python:
`checkpoint_products=True`). Each product's tables and metadata commit together;
a later failure or cancellation rolls back only the current product. A session
advisory lock spans all commits, and each product checks the same single release
ID and manifest digest before and after its work. Resume by selecting only the
remaining products with `--products`; existing metadata is not skipped
implicitly.

```sh
uv run --frozen omnipath-subsets build --schema release_2026_09 --checkpoint-products
uv run --frozen omnipath-subsets build --schema release_2026_09 --checkpoint-products --products network_views cosmos
```

The optional Python `on_product_committed(product, result)` callback requires
checkpoint mode and runs only after the durable commit. `result` is a
`BuildResult` snapshot containing release identity and the stats/timings of
products committed by this call. A callback failure stops further builds while
preserving that commit; `subset_build_metadata` is authoritative if observer
persistence fails.

The maintained scientific implementations live in `omnipath_subsets.exclusion_metabolites`,
`omnipath_subsets.metsigdb`, `omnipath_subsets.network_views` and
`omnipath_subsets.cosmos`. The standalone
command and PostgreSQL importer share `omnipath_subsets.runner`; they publish
`subset_build_metadata` and `parquet_phase` in the same product transaction.
Explicit selection rebuilds each product. Use `omnipath-postgres --finish` to
skip products already committed for the loaded monthly release.

## Products

- **Exclusion metabolites:** the chemical entities that mass spectrometry does not
  measure (ions, atoms, small molecules, xenobiotics, cofactors and extended cofactors), in classes.
  MetSigDB and COSMOS read it. See [Exclusion metabolites](#exclusion-metabolites).
- **MetSigDB:** metabolite memberships in Reactome, WikiPathways, KEGG, MACdb
  and ClassyFire sets, read from the normalized entity, evidence and ontology tables.
- **Network views:** the MetaLinksDB, LIANA and Reactions preset definitions and
  registry. Registering a preset does not create a materialized network; consumers
  apply the stored recipe.
- **COSMOS:** metabolite/enzyme edges projected from the shared reaction tables,
  including reversible reactions and labels drawn from published identifiers.
  GeneID-centered catalysts keep their GeneID label and `entrez` namespace,
  including records whose source reported a protein. UniProt catalogue aliases
  of a gene do not select a protein participant or create extra reaction nodes.
  Canonical UniProt product identities keep their UniProt labels. The optional
  legacy utils translation also preserves GeneID and accepts other catalyst
  mappings only when exactly one UniProt target answers. Connectors use the
  same identifier and namespace as their reaction node.

All four products use the current PostgreSQL schema. The earlier resource-record
subset implementations and their verification command have been removed.

## Exclusion metabolites

The product publishes two tables in the release schema. `exclusion_metabolite`
holds one row for each entity, class and source of the assignment.
`exclusion_metabolite_report` says what the list covers and what it does not. It
reads one loaded release. It resolves no identifier and downloads nothing. It
runs first, and you can rebuild it alone:

```sh
make rebuild-subsets POSTGRES_SCHEMA=release_2026_09 SUBSET_PRODUCTS=exclusion_metabolites
```

### The classes

| Class | Where its members come from |
|---|---|
| `ions` | The MetaProViz seeds, and every ChEBI `subclass_of` descendant of `CHEBI:24867`. |
| `atoms` | The MetaProViz seeds, and every ChEBI descendant of `CHEBI:33250`. |
| `small_molecules` | The MetaProViz KEGG seeds. Nothing derives them. |
| `cofactors` | The core list and its structure-key variants. Nothing derives them. |
| `cofactors_extended` | Chemicals with the ChEBI cofactor role (`CHEBI:23357`) and their variants. It adds to `cofactors`, so leave out both for the whole set. |
| `xenobiotics` | Valid and empty until a curator adds a member. |

A variant of a listed member shares its structure key. It can be a charge,
tautomer or stereo form. The build lists it in the member's classes as a row of
its own. Its `match_kind` is `structure_key`, and `via_entity_id` names the
member. An entity that matched in its own right gets no such row.

### The three files you edit

The files are in `omnipath_subsets/exclusion_metabolites/seed/`. Edit one,
commit it and rebuild the product. No code changes.

- `metaproviz_seed.yaml` holds the MetaProViz lists and their license. Change it
  only for a newer MetaProViz version.
- `project_cofactors.yaml` holds the reviewed core of cofactors. Add a line when
  the role misses one: `{id: "CHEBI:16761", name: ADP}`. The name is for the reader.
- `overrides.yaml` holds your changes. Give each one a reason.

```yaml
add:  [{class: ions, namespace: chebi, id: "CHEBI:29101", note: why it belongs}]
veto: [{class: cofactors, namespace: chebi, id: "CHEBI:15422", note: why it does not}]
```

The classes are `ions`, `atoms`, `small_molecules`, `xenobiotics`,
`cofactors` and `cofactors_extended`.

The namespaces are `chebi`, `kegg`, `hmdb`, `pubchem`, `lipidmaps`,
`swisslipids`, `chembl` and `drugbank`.

Write a ChEBI id with or without `CHEBI:`. Write an HMDB id in either width. A malformed file stops the build and
says why.

**A veto wins.** It removes the entity from that class, whatever put it there.
That can be a seed, the ontology, an addition or a structure-key match. A
vetoed member does not pull in its variants, and a vetoed variant stays out.

A veto leaves other classes alone. The report shows each veto with the rows it
removed. It shows a veto that removed nothing too, because that veto may be
stale.

### Read the list

```sql
-- the members of a class, with where each came from
SELECT label, chebi, kegg, hmdb, assignment_source, origin, match_kind
FROM exclusion_metabolite_members WHERE exclusion_class = 'cofactors' ORDER BY label;

-- what the build covered, and what it did not
SELECT kind, exclusion_class, namespace, resolved, fallback_only, unresolved
FROM exclusion_metabolite_report ORDER BY kind, exclusion_class;
```

`origin` is `curated` for the three seed files and `derived` for the ontology.
An entity with both has two rows. In the report, `seed_coverage` counts the seeds
that resolved and names the ones that did not. `shared_key` names the variants a
key pulled in. `veto_applied` names your vetoes.

Every row carries `list_version`, the commit of the three files. A build from
uncommitted edits gets the stamp `uncommitted:<sha256>`, and the report says so.
Commit before you build anything you intend to keep.

## Validation scope

```sh
make test-subsets
```

Tests use disposable PostgreSQL and bounded synthetic fixtures. Current product
parity tests under `packages/postgres/tests` compare complete rows, scientific
semantics, preset definitions and database objects with frozen reference code.
Runner tests cover release identity, locking, durable product checkpoints,
rollback, resume and explicit rebuild. These developer tests do not add a
verification or rebuild phase to production runs.
