# MetalinksDB from Parquet with DuckDB

2026-10-26

**Result: MetalinksDB can be rebuilt faithfully from the Parquet release with DuckDB.**

For release 2026.10.19-b, the release dev5 has loaded, the DuckDB rebuild matches PostgreSQL exactly. The checks below cover every column PostgreSQL stores. The one-time build takes about 2.5 minutes on 3 threads.

## What was rebuilt

PostgreSQL builds MetalinksDB in three layers. The first two exist only in PostgreSQL today. All three are now ported to DuckDB in `scripts/metalinksdb_from_parquet.py`:

| Layer | In PostgreSQL | Rebuilt in DuckDB |
|---|---|---|
| **1. Annotation linking** | `projection.py` | Links each evidence row to its own annotations. Also links each statement annotation whose source and dataset match the evidence row, and each row without evidence to a synthetic evidence row. |
| **2. The interaction record** (`interaction_fact_resource`) | `derived_tables.py`, the 26-minute interaction projection | Interaction class per canonical relation, with the same precedence, read from every resource of the release. Then direction, sign, PubMed IDs and DOIs, affinity with unit scaling, pChEMBL, STITCH score, curation flags, and organisms with the per-source fallback. Entity IDs are the projection's own UUIDs, so rows compare directly. |
| **3. The preset** | The `METALINKSDB` recipe and main's reference fold (`tests/fixtures/collapse.py`) | ChEMBL mechanism pairs, the curated sources, and the reaction sources' transport and ligand–receptor classes, each through the chemical, partner-type and class gate. BindingDB is excluded. Then the fold on subject, object and class. |

The script reads its constants from the PostgreSQL code: the class map `interaction_class.yaml`, and the directed predicates, annotation-class terms, affinity terms and curation terms in `derived_tables.py`. It reads the recipe from the preset definition. A change to any of them carries over automatically.

## How it was checked

1. The pins were read from dev5's `parquet_release`. They are identical to 2026.10.19-b.
2. The PostgreSQL side was exported from dev5's `interaction_fact_resource`: every row of MetalinksDB's sources, and for ChEMBL the rows flagged mechanism-of-action, with endpoint types.
3. The DuckDB record was compared on the full record key: endpoints, class, source, direction and sign. Then the folded MetalinksDB was compared on subject, object and class.

| Check | PostgreSQL | DuckDB | Result |
|---|---|---|---|
| Record rows | 1,418,198 | 1,418,198 | All matched on the full key |
| Record columns (10) | | | All identical |
| MetalinksDB pairs | 820,516 | 820,516 | All matched |
| MetalinksDB columns (15) | | | All identical |

The 10 record columns are endpoint types, organisms, affinity, pChEMBL, score, curation flags, PubMed IDs and DOIs.

The 15 MetalinksDB columns are organisms, direction, stimulation, inhibition, the sign and direction source counts, affinity, pChEMBL, score, sources, source count, curation flags, PubMed IDs and DOIs.

**Time per step**, on nicesrv with 3 threads and 8 GB:

| Step | Time |
|---|---|
| Class evidence over 46 resources | 18 s |
| Annotations of the scope | 49 s |
| Record | 46 s |
| Preset and fold | 3 s |
| Everything else | about 20 s |

PostgreSQL needs the whole 3-hour load first.

## Not covered

- **The `intercell` attribute layer.** Its implementation is in neither repo: the build stores only the layer's name, and the api-service that applies it is not here. It can't be checked until it is specified.
- **The port-only predicate gate** (`predicates` in `_definitions.py`). The PostgreSQL record has no predicate column, so PostgreSQL cannot apply it. Neither side applied it here. On Parquet it is a one-line filter on the relation's predicate.
- **Labels.** PostgreSQL's `label` attribute comes from its own labels step. Parquet carries the resolver's label. They were not compared.
- **Record rows with no interaction header**, which PostgreSQL drops. They aren't reproduced, and no difference showed up.
- **MetaChatDB.** It is not in dev5's pins, so it wasn't compared. It goes through the same code path.

## What it would take to use this

1. Decide where the output lives: a release-level product file (e.g. `releases/<version>/products/metalinksdb.parquet`), or a query the engine runs.
2. Specify the `intercell` layer and decide on the predicate gate. Neither has a PostgreSQL reference to match.
3. Turn the script into a package step with tests. Its parts map one-to-one to the PostgreSQL functions named in its docstring.
4. Rerun the comparison whenever the record logic in PostgreSQL changes, for as long as both exist.

Script: `scripts/metalinksdb_from_parquet.py`. Run output and the DuckDB working database: nicesrv `/root/data/scratch/build-20261025/ml-faithful/`.
