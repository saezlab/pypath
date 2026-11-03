# Pipeline Guide

The `omnipath_build` pipeline builds PostgreSQL tables from `pypath.inputs_v2`
datasets through a DuckDB staging and canonicalization path.

This page is the short entry point. The detailed workflow is split into ordered
step files under [`docs/pipeline/`](pipeline/):

1. [Setup](pipeline/1_setup.md)
2. [Resolver Materialization](pipeline/2_resolver_materialization.md)
3. [Database Schema](pipeline/3_database_schema.md)
4. [Source Discovery](pipeline/4_source_discovery.md)
5. [DuckDB Projection](pipeline/5_duckdb_projection.md)
6. [Canonicalization](pipeline/6_canonicalization.md)
7. [PostgreSQL Copy](pipeline/7_postgres_copy.md)
8. [Refresh And Reload](pipeline/8_refresh_reload.md)
9. [Derivation](pipeline/9_derivation.md)

## Pipeline Shape

At a high level, the pipeline is:

```text
setup
  -> resolver parquet files
  -> PostgreSQL schema
  -> pypath source discovery
  -> DuckDB evidence projection
  -> DuckDB canonicalization
  -> PostgreSQL COPY
  -> derived query tables and bitmaps
```

The main Makefile targets are:

```bash
make setup
make resolver
make db-setup
make load
make derive
```

For a scratch build:

```bash
make all DROP_EXISTING=1
```

For a bounded source-level smoke test:

```bash
make load SOURCE=bindingdb MAX_RECORDS=200000 SCHEMA=omnipath_test
```

## Main Code Entry Points

- `Makefile`: user-facing commands and defaults.
- `omnipath_build/cli.py`: administrative commands such as schema setup,
  resolver materialization, source deletion, and derivation.
- `omnipath_build/duckdb_direct_pipeline.py`: active load orchestration.
- `omnipath_build/resources.py`: dynamic `pypath.inputs_v2` source discovery.
- `omnipath_build/evidence_projector.py`: shared source-record projection rules.
- `omnipath_build/duckdb_load.py`: DuckDB resolver views, canonicalization SQL,
  and PostgreSQL COPY helpers.
- `omnipath_build/db/`: PostgreSQL schema, refresh, indexes, derived tables,
  bitmaps, and resource metadata sync.
