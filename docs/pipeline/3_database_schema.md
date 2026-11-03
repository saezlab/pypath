> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 3. Database Schema

Schema setup creates the PostgreSQL objects needed for evidence, canonical graph
rows, refresh operations, and derived query artifacts.

## Command

```bash
make db-setup
```

For a clean schema:

```bash
make db-setup DROP_EXISTING=1
```

To drop and recreate the target schema without running the resolver:

```bash
make db-reset
```

## What It Creates

The schema contains:

- Dimension and vocabulary tables such as `data_source`, `dataset`, identifier
  types, entity types, predicates, roles, and resolution statuses.
- Source-partitioned evidence tables such as `entity_evidence`,
  `entity_evidence_identifier`, `relation_evidence`, evidence annotations, and
  evidence-to-canonical mapping tables.
- Canonical graph tables such as `entity` and `relation`.
- Derived tables and bitmap tables, created or refreshed later by `make derive`,
  including `entity_identifier_lookup` for indexed app search.

## Deferred Indexes

When `DROP_EXISTING=1` is used, the Makefile calls `init-db --drop-existing
--no-indexes`. This defers expensive secondary indexes until after bulk loading,
which avoids maintaining those indexes row by row during ingest.

## Main Files

- `Makefile`: `db-setup` and `db-reset` targets.
- `omnipath_build/cli.py`: `init-db` command.
- `omnipath_build/db/schema.py`: PostgreSQL DDL.
- `omnipath_build/db/indexes.py`: deferred and secondary indexes.

## Phase Boundary

Schema setup prepares empty structures. It does not discover pypath sources,
load source evidence, canonicalize records, or build query-ready derived data.
