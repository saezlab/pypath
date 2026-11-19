# OmniPath subsets

A published release gives four products: the exclusion metabolite list, MetSigDB
memberships, network presets and COSMOS reaction networks. Products use published
identities and aliases. They run no resource builds and no identifier resolution.

The exclusion metabolites and MetSigDB are **Parquet subsets**. They read the whole
release: the entities, identifiers, relations and ontology of many resources. The build
reads the published Parquet once, after the release. The API and PostgreSQL both take
the subsets from there. PostgreSQL builds the network views and COSMOS from a loaded
release.

## The order

```
make build                    resources                        resources/<source>/<version>/
make publish-release          the pinned release               releases/<version>.json
make build-subsets RELEASE=v  the Parquet subsets              subsets/<release>/<product>/
make release-db RELEASE=v     the serving database             serving/<version>.duckdb
make load-postgres            base tables, then the subsets    PostgreSQL schema
```

`publish-release --build-database` builds the subsets before the serving database.
Add `--skip-subsets` to serve the release without them.

## Build the Parquet subsets

```sh
uv run --frozen omnipath-subsets build --release 2026.11 --data-root data
uv run --frozen omnipath-subsets build --release 2026.11 --data-root data --products exclusion_metabolites
uv run --frozen omnipath-subsets check --release 2026.11 --data-root data
```

The build runs the products in dependency order. The exclusion metabolites come first.
MetSigDB follows and carries each metabolite's exclusion classes. The build writes each
product in a temporary directory and validates it. Then it moves the product to
`<data_root>/subsets/<release>/<product>/` with a `subset_manifest.json`. The manifest
holds:

- the SHA-256 of the release manifest and the digests of the dependencies
- the code version
- the rows and SHA-256 of each table
- the row contract and the counts the build reported

A subset is immutable. A second build for the same release and product fails. A release
that pins sample builds gives subsets marked partial.

`check` runs the release checks on the subsets and on their release. It
tests for duplicates, mandatory fields, the set type of each resource, set sizes, the
build stamp, organism and overview-map rules, ChemOnt sets, MACdb ids, and the
exclusion classes and their report. It fails on an error. It uses no frozen count.

The Python entry point is `omnipath_subsets.parquet_products.build_release_subsets`.
The layout and its rules are in `omnipath_core.subset_layout`, which the API also reads.

## Build products in PostgreSQL

The network views and COSMOS read a loaded release:

```sh
export OMNIPATH_DATABASE_URL='postgresql:///omnipath'
uv run --frozen omnipath-subsets build-postgres --schema release_2026_09
uv run --frozen omnipath-subsets build-postgres --schema release_2026_09 --products cosmos
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
uv run --frozen omnipath-subsets build-postgres --schema release_2026_09 --checkpoint-products
```

The optional Python `on_product_committed(product, result)` callback requires
checkpoint mode and runs only after the durable commit. `result` is a
`BuildResult` snapshot containing release identity and the stats/timings of
products committed by this call. A callback failure stops further builds while
preserving that commit; `subset_build_metadata` is authoritative if observer
persistence fails.

## Load the Parquet subsets into PostgreSQL

After the shared derivations, `omnipath-postgres` loads each present subset through its
staging and COPY path. It creates the tables and a `subset_build_metadata` row that
names the digest of the subset manifest. Then it creates the indexes. No product logic
runs there.

The loader reads the subsets from the data root. Use `make load-postgres` or
`omnipath-postgres --finish --data-root`. It loads a subset only for the published
release that the subset came from. The subset's entities must be entities of the loaded
release. Load the published manifest (`<data_root>/releases/<version>.json`).

A release is complete when the loader and the build have made all four products. Use `--products` to
leave a product out.

## Products

- **Exclusion metabolites:** the chemical entities that mass spectrometry does not
  measure (ions, atoms, small molecules, xenobiotics, cofactors and extended cofactors), in classes.
  MetSigDB and COSMOS read it. See [Exclusion metabolites](#exclusion-metabolites).
- **MetSigDB:** metabolite memberships in Reactome, WikiPathways, KEGG, MACdb
  and ClassyFire sets, read from the release's per-resource Parquet. See
  [MetSigDB](#metsigdb).
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

The exclusion metabolites and MetSigDB are Parquet. The network views and COSMOS are
PostgreSQL products. We removed the earlier resource-record subset implementations and
their verification command.

## MetSigDB

The product has one table, `metsigdb_membership`. One row is one metabolite in one set of
one resource. It comes from the per-resource Parquet of the release. The Parquet names
entities by their release `entity_key`. PostgreSQL gets the ids that the release gave
them. Each row carries:

- the metabolite's identifiers from every pinned resource that states them
- its structure key, and whether it has a structure at all (`has_structure`)
- its exclusion classes, as a sorted list that is empty when there are none
- the set: id, label, type, sub-type, organism and size
- the provenance of the membership

The manifest versions the row contract (`row_contract`). Version 2 prefixes MACdb set
ids (`MACDB:<n>`), names the entities by key and adds `has_structure` and
`exclusion_classes`.

Three rules decide what a set is:

- A ClassyFire set is a ChemOnt term. HMDB also links metabolites to disease terms, which
  are not chemical classes.
- A KEGG reference map (`rn*`) has no organism.
- The map number identifies a KEGG overview map, whatever the organism prefix.

The API serves the product as `/sets/metsigdb`. Its OpenAPI description gives the
identifier coverage that the build measured for the release.

## Exclusion metabolites

The product publishes two tables. `exclusion_metabolite`
holds one row for each entity, class and source of the assignment.
`exclusion_metabolite_report` says what the list covers and what it does not. It
reads the published release. It resolves no identifier and downloads nothing. It
runs first, and you can build it alone:

```sh
uv run --frozen omnipath-subsets build --release 2026.11 --data-root data --products exclusion_metabolites
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

The tables are Parquet. The loader creates `exclusion_metabolite_members` and the queries
below in PostgreSQL. In the Parquet, `entity_key` names the entities.

The API serves both tables. `/sets/exclusion-metabolites` returns the rows with each
metabolite's label and identifiers (HMDB, ChEBI, KEGG, PubChem, InChIKey), filtered by class,
source, origin, match kind or metabolite. `/sets/exclusion-metabolites/report` returns the
report. The clients call them `exclusion_metabolites()` and `exclusion_metabolite_report()`
in Python and `op_exclusion_metabolites()` and `op_exclusion_metabolite_report()` in R. The
`exclude_metabolites` filter of `/sets/metsigdb` stays as it is and reads the same table.

## Validation scope

```sh
make test-subsets
```

Tests use disposable PostgreSQL and bounded synthetic fixtures. They compare the Parquet
products with rows frozen in `tests/golden/` (see its README). The parity tests for the
PostgreSQL products are under `packages/postgres/tests`. They compare preset definitions
and database objects with frozen reference code.
Runner tests cover release identity, locking, durable product checkpoints,
rollback, resume and explicit rebuild. These developer tests do not add a
verification or rebuild phase to production runs.
