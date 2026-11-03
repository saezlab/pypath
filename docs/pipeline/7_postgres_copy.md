> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 7. PostgreSQL Copy

The copy phase moves projected evidence and canonical graph rows from DuckDB
into PostgreSQL.

## Where It Runs

The copy phase is part of `make load`, after DuckDB projection and
canonicalization.

## What It Does

For each batch or staged DuckDB file, the loader:

1. Loads DuckDB's PostgreSQL extension.
2. Attaches the target PostgreSQL database.
3. Creates load views over the DuckDB tables.
4. Verifies that content tables are empty when append mode is disabled.
5. Inserts small dimension values into PostgreSQL.
6. Materializes dimension IDs used by evidence and canonical rows.
7. Drops high-volume load constraints and indexes when configured.
8. Bulk-copies evidence rows.
9. Bulk-copies canonical graph rows.
10. Resets PostgreSQL sequences after load completion.

## Append Behavior

The Makefile passes `--append` to the loader. That means normal Makefile loads
append into an existing schema after source-level skip or reload checks have
decided what should be loaded.

## Source Partitions

Evidence tables are source-partitioned in PostgreSQL. The lower-level loader
creates or attaches source partitions where needed before inserting rows.

## Main Files

- `omnipath_build/duckdb_direct_pipeline.py`: PostgreSQL preparation and copy
  orchestration.
- `omnipath_build/duckdb_load.py`: bulk-copy helpers.
- `omnipath_build/db/schema.py`: source-partitioned table definitions.

