---
name: omnipath-parquet
description: Installs versioned OmniPath resources locally with the omnipath_client Python package and queries them as a DuckDB database, with SQL, entity search and relation exports. Use when finding OmniPath resources, getting OmniPath data, or answering questions about entities, relations, sources and evidence in OmniPath.
---

# Use OmniPath

Use the `omnipath_client` Python package. It downloads the resources you choose
once, writes a DuckDB database of them, and answers every later question from
your machine, with SQL or with the same search and relation export the
OmniPath website runs. The API is contacted only while installing.

```python
import omnipath_client as op

op.init(["signor"], progress=print)        # once: download and write the database

with op.connect() as db:
    db.sql("SELECT predicate, count(*) n FROM signor.relation GROUP BY 1").fetchall()
    db.search("TP53")                      # ranked entities, as the website searches
    db.relations("TP53", sources="signor") # relations, as the explorer exports them
```

Use the HTTP API directly (last section) only when Python is not available.

## Install

Python 3.11 or newer. Install the client from GitHub; it needs only DuckDB, httpx and
platformdirs:

```bash
R=git+https://github.com/saezlab/omnipath-build@parquet-migration
pip install "omnipath-client @ $R#subdirectory=packages/query/python"
# add pandas for .df()
```

The R package `omnipath.client` (`packages/query/r`) reads and writes the same data
directory with the same queries.

With uv, use `uv pip install` with the same arguments. Do **not** run a bare
`pip install omnipath-client`: that name on PyPI is a different package
(saezlab/omnipath-client) and does not have `init` or `connect`. The older
`omnipath` package is unrelated and can coexist.

## Choose resources

```bash
omnipath-client list        # name, version, size, installed state, tags
```

```python
for r in op.available():    # the catalog records
    print(r["resource_id"], r["version"], r["total_size_bytes"],
          [(t["dimension"], t["id"]) for t in r["tags"]])
```

Each record also has `display_name`, `description`, `license`, `website`,
`entity_count`, `interaction_count` and `files`.

Resources are tagged in four dimensions, which `init` takes as arguments:

- `topic`: `signaling`, `gene_regulation`, `metabolism`, `pharmacology`,
  `cell_cell_communication`, `disease_phenotype`, `subcellular_localization`,
  `cell_types`, `molecular_classification`, `nutrition`
- `content`: `physical_interactions`, `causal_regulation`,
  `enzyme_ptm_relationships`, `ligand_receptor_interactions`,
  `compound_target_interactions`, `protein_complexes`, `pathways`,
  `biochemical_reactions`, `molecular_transport`, `molecular_annotations`,
  `identifier_mappings`, `ontologies`, `trait_associations`,
  `literature_mentions`, `food_composition`
- `knowledge_origin`: `literature_curated`, `experimental`,
  `computationally_predicted`, `text_mined`
- `molecule`: currently unused

Tags can change between releases, so check `available()` when a tag is not
accepted. Install only what the task needs: every resource is about 5.6 GB;
most single resources are a few MB to a few hundred MB.

## Install data

```python
op.init(["signor", "intact"])                  # by name
op.init(topic="cell_cell_communication")       # by tag
op.init(topic=["signaling", "gene_regulation"], knowledge_origin="literature_curated")
op.init()                                      # every resource of the release
```

```bash
omnipath-client init signor intact
omnipath-client init --topic signaling --knowledge-origin literature_curated
```

- **Selection.** Several tags of one dimension are alternatives; different
  dimensions must all match.
- **Adding.** Resources already installed stay installed. Running `init` again
  with more names adds them and downloads only what is missing.
- **Releases.** The default `release="latest"` takes the newest version of each
  resource. `init(release="2026.9.6.2")` installs the versions that OmniPath
  release pins; a release and the resource versions it pins have different
  numbers. A data directory always holds one release: `init` with another
  release moves every installed resource to it.
- **Location.** The per-user data directory (`op.data_directory()`), or
  `OMNIPATH_DATA_DIR`, or `data_dir=...` (`--data-dir`). Pass the same
  `data_dir` to `connect`. Use a project directory when the work must be
  reproducible or kept apart from other work.
- **Raw evidence.** Original source records (`evidence_payloads`) are
  downloaded only with `include_evidence=True` (`--evidence`). Ask for them
  only when the task needs the source rows themselves.
- **Speed.** `materialize=True` (`--materialize`) copies the data into native
  DuckDB tables: two to five times faster queries, about twice the files' size
  in extra disk. The default views are enough for most work.
- **Interruptions.** An interrupted `init` keeps finished resources; run it
  again to continue.

`init` contacts `OMNIPATH_API_URL` when set, otherwise the hosted API
`https://omnipath-metabo-dev.schaul.click/api`.

## Query

```python
db = op.connect()                    # or op.connect(data_dir=...)
db.resources()                       # installed resources, versions, tags, licences
```

`sql`, `search` and `relations` return DuckDB relations: rows are read when you
fetch or export them.

```python
rel = db.sql("SELECT * FROM relation WHERE predicate = ?", ["affects"])
rel.fetchall()                       # list of tuples
rel.df()                             # pandas (needs pandas)
rel.to_arrow_table()
rel.write_parquet("affects.parquet")
rel.filter("sign = -1").aggregate("count(*)")   # refine a result without rewriting the SQL
```

### Prebuilt queries

These run the engine behind the hosted API and give the website's answers on
the installed resources.

```python
db.search("TP53")                                      # entity_key, label, entity_type, namespace, identifier, taxon, relation_count, ...
db.search("gluc", entity_types="chemical_entity", limit=50)
key = db.search("TP53").fetchone()[0]                  # entity_key of the best match

db.relations("TP53", sources="signor")
db.relations(["TP53", "MDM2"], predicates="affects", taxonomy_ids="9606")
db.relations("gene:entrez:7157")                       # every product of a gene

db.entity(key)                                         # dict: entity, annotations, summary, relationships
```

- `relations` takes exact names, identifiers, entity keys or group keys
  (`gene:entrez:7157`) at either endpoint. It returns one row per relation,
  with the `relation` columns plus `annotations`, `evidence` (with molecular
  forms) and `referenced_product_records` as nested lists. At most 100,000
  relations; use SQL on `relation` for more.
- Filters, for `relations` and `search`: `sources`, `predicates`, `signs`,
  `taxonomy_ids`, `entity_types`, `partner_types`, `relation_categories`,
  `interaction_types`, `is_directed`, `min_evidence_count`,
  `reference_entity_keys`, `protein_entity_keys`, `transcript_entity_keys`,
  `isoform_identifiers`, `annotation_term_ids` and the qualifier filters. One
  value stands for a list of one. An unknown filter name raises an error that
  lists the valid ones.
- `entity` takes the exact `entity_key`. It returns None when no installed
  resource has it.
- With a subset installed, counts and merged entities cover only that subset.
- `db.engine` exposes the rest of the API (facets, grouped search, ontology
  browsing, molecular context).

### Tables

| Name | Contents |
| --- | --- |
| `entity`, `relation`, ... | one table over every installed resource, with `_resource` and `_resource_version` columns |
| `signor.relation`, `chembl.entity`, ... | the same tables of one resource, in a schema named after it |
| `resources` | installed resources, versions, and whether raw evidence is installed |

In the tables over every resource, shared entities and relations are **not
deduplicated**: there is one row per resource that has them. Use
`entity_key` / `relation_key` to match them across resources, or `search` and
`relations`, which merge them.

Join child tables to their parent on `entity_id` / `relation_id` **and**
`_resource` (these ids number rows within one resource version):

```sql
SELECT e.label, r.predicate, r.object_label, r._resource
FROM entity AS e
JOIN relation AS r
  ON e.entity_key = r.subject_entity_key
 AND e._resource = r._resource
WHERE e.reference_entity_key = 'entrez:7157'   -- TP53 and its products, whatever the canonical id
```

Lookup tables `init` derives, for fast filtering without scanning:

- `entity_term` (`term`, `kind`, `entity_id`): lowercase labels and
  identifiers, to find entities by name: `WHERE term = 'tp53'`
- `relation_endpoint` (`key`, `key_kind`, `side`, `relation_id`, `predicate`,
  `category`): endpoint keys of relations, to find relations by endpoint
- `entity_group` (`group_key`, `kind`, `entity_id`): gene, reference and
  structure (`connectivity`) groups of entities
- `relation_facet` (`category`, `predicate`, `subject_type`, `object_type`,
  `taxon`, `sign`, `is_directed`, `sources`, the three qualifiers,
  `evidence_count`, `relation_count`): each distinct combination of the
  relation filter columns, with the largest evidence count and the number of
  relations, to see which values occur together without scanning `relation`

The database is read-only. Make your own tables with `CREATE TEMP TABLE`, or
write results to Parquet.

## Schema

`entity_key` and `relation_key` are 64-character keys that identify an entity
or relation across resources. `ordinal` keeps the original order of child rows.

### `entity`

| column | type | meaning |
| --- | --- | --- |
| `entity_id` | int32 | Row number within the resource version |
| `entity_key` | string | Canonical entity identity |
| `entity_type` | string | `protein`, `chemical_entity`, `macromolecular_complex`, … |
| `namespace`, `identifier` | string | Canonical identifier, e.g. `uniprot` / `P04637` |
| `taxon` | string | NCBI taxonomy id when known |
| `label` | string | Display name |
| `reference_entity_key` | string | Shared reference, normally `entrez:<GeneID>` |
| `gene_reference_keys` | list<string> | Catalogue gene links of a product |
| `group_connectivity` | string | First InChIKey block (chemical structure group) |
| `has_hierarchy`, `parent_count`, `child_count` | bool, int64 | Ontology-style hierarchy |
| `*_count` | int64 | Number of identifiers, annotations, evidence, relations |

Child tables (each with `entity_id`, `ordinal`):

- `entity_identifier`: `ns`, `id`, `is_canonical`, `source`
- `entity_annotation`: `term`, `value`, `quantity_*` (measurement fields), `source`, `dataset`
- `entity_evidence`: `source`, `dataset`, `row_id`, `upstream_id`, `annotations`, `molecular_form`

### `relation`

| column | type | meaning |
| --- | --- | --- |
| `relation_id` | int32 | Row number within the resource version |
| `relation_key` | string | Canonical relation identity |
| `statement_kind` | string | `relation` or `ontology` |
| `subject_entity_key`, `object_entity_key` | string | Endpoint `entity_key`s |
| `subject_label`, `subject_type`, `object_label`, `object_type` | string | Endpoint names and types |
| `*_reference_entity_key` | string | Endpoint references |
| `predicate` | string | Relation predicate, e.g. `affects`, `has_member` |
| `taxon` | string | Shared or primary taxon |
| `is_directed` | bool | Directed edge |
| `sign` | int32 | `-1`, `0` or `1` |
| `category` | string | `interaction`, `association`, `ontology`, … |
| `interaction_class` | string | Optional interaction class |
| `sources` | list<string> | Provenance sources |
| `evidence_count` | int64 | Number of `relation_evidence` rows |
| `*_qualifier` | list<string> | Biolink qualifiers |

Child tables (each with `relation_id`, `ordinal`):

- `relation_annotation`: `term`, `value`, `quantity_*`, `source`, `dataset`, `scope`
- `relation_evidence`: `source`, `dataset`, `row_id`, `upstream_id`, `annotations`, `subject_molecular_form`, `object_molecular_form`

### `evidence_payloads` (only with raw evidence)

| column | type | meaning |
| --- | --- | --- |
| `relation_key` | string | Parent relation, null for entity payloads |
| `entity_key` | string | Parent entity, null for relation payloads |
| `source` | string | Upstream resource |
| `row_id` | string | Source row identity |
| `payload_json` | string | Original source record as JSON |

Join a payload to its evidence row on `relation_key` (or `entity_key`) +
`source` + `row_id`, within the same resource.

## Without Python

The database works with any DuckDB client, for SQL only (`search`,
`relations` and `entity` need Python):

```bash
duckdb -readonly "$(python -c 'import omnipath_client as op; print(op.data_directory())')/omnipath.duckdb"
```

```sql
ATTACH '/path/to/omnipath.duckdb' AS omnipath (READ_ONLY);   -- from another database, R, ...
```

Without the client at all, the HTTP API serves the catalog and the files
(`https://omnipath-metabo-dev.schaul.click/api`):

```bash
API=https://omnipath-metabo-dev.schaul.click/api
curl -sS "$API/releases"                                   # published releases
curl -sS "$API/resources?shape=svelte"                     # catalog: versions, sizes, files[].url
curl -fsSL -o signor-relation.parquet "$API/resources/signor/files/relation.parquet"
curl -fsSL -o signor.zip "$API/resources/signor/download"  # every table but evidence_payloads
```

Pass `?release=<release>` to every request to pin a release. The files are the
tables above, one Parquet file each, without `_resource` columns and without
the lookup tables.

## Rules

- Use the client. Fall back to the HTTP API only without Python.
- Install the smallest selection that answers the question; add resources
  later with another `init`.
- Use `search` to turn a name into an `entity_key`, and `relations` for an
  entity's relations, before writing SQL for the same thing. Use SQL for
  counts, aggregates and anything over 100,000 relations.
- In SQL over every resource, remember there is one row per resource; join
  on `_resource` as well as the ids.
- For reproducible work, pin `release=...`, use a dedicated `data_dir`, and
  record `db.resources()` (resource versions) with the results.
- Request raw evidence only when the task needs original source records.
- Write result files as Parquet.
