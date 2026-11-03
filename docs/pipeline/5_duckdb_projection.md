> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 5. DuckDB Projection

Projection streams pypath source records into DuckDB evidence tables.

## Command

```bash
make load
```

Common source-scoped variants:

```bash
make load SOURCE=bindingdb
make load SOURCES=uniprot,bindingdb,intact
make load SOURCE=bindingdb MAX_RECORDS=200000
```

## What It Does

For each selected dataset, the loader:

1. Opens the raw pypath dataset.
2. Splits records into batches using `BATCH_SIZE`, default `50000`.
3. Creates DuckDB resolver views from resolver parquet files.
4. Creates DuckDB raw evidence tables.
5. Projects source records into evidence rows.

## Evidence Tables In DuckDB

Projection writes these raw tables:

- `entity_evidence_raw`
- `entity_identifier_raw`
- `entity_annotation_raw`
- `relation_annotation_raw`
- `annotation_value`
- `relation_evidence_raw`
- `annotation_relation_evidence_raw`
- `ontology_terms_raw`

## Projection Rules

The shared projector flattens pypath silver entity trees:

- Ordinary entities become entity evidence rows.
- Identifiers become entity identifier rows.
- Non-ontology annotations become annotation rows.
- Ontology-valued annotations can become relation evidence to CV-term entities.
- Two-member interaction-like records become relation evidence.
- Other member trees produce membership relation evidence.

## Main Files

- `omnipath_build/duckdb_direct_pipeline.py`: dataset loop, batching, and
  direct load orchestration.
- `omnipath_build/evidence_projector.py`: shared projection rules.
- `omnipath_build/duckdb_load.py`: DuckDB table creation and row writers.
