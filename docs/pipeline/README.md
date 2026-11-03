> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# Pipeline Steps

This folder records the previous `omnipath_build` pipeline in execution order.
Each numbered file focuses on one phase and points to the code that owns it.

## Ordered Steps

1. [Setup](1_setup.md)
2. [Resolver Materialization](2_resolver_materialization.md)
3. [Database Schema](3_database_schema.md)
4. [Source Discovery](4_source_discovery.md)
5. [DuckDB Projection](5_duckdb_projection.md)
6. [Canonicalization](6_canonicalization.md)
7. [PostgreSQL Copy](7_postgres_copy.md)
8. [Refresh And Reload](8_refresh_reload.md)
9. [Derivation](9_derivation.md)

## Data Layers

The build produces three main layers:

- Source evidence tables: source rows, entity occurrences, identifiers,
  relation evidence, and annotations.
- Canonical graph tables: deduplicated `entity` and `relation` rows plus the
  evidence-to-canonical mappings.
- Derived query artifacts: count tables, ontology-term summaries, roaring bitmap
  indexes, and resource metadata.

## Operational Defaults

- Default database URL:
  `postgresql://omnipath:omnipath@localhost:55432/omnipath`
- Default schema: `public`
- Default inputs package: `pypath.inputs_v2`
- Default resolver/output root: `data`
- Default batch size: `50000`
- Default DuckDB thread count: `4`

## Common Workflows

Full scratch build:

```bash
make all DROP_EXISTING=1 DERIVE=1
```

Scratch build with existing resolver files:

```bash
make db-setup DROP_EXISTING=1
make load
make derive
```

Load one source:

```bash
make load SOURCE=bindingdb
make derive
```

Reload one source:

```bash
make reload SOURCE=bindingdb
make derive
```

Smoke test in another schema:

```bash
make db-setup SCHEMA=omnipath_test DROP_EXISTING=1
make load SOURCE=bindingdb MAX_RECORDS=200000 SCHEMA=omnipath_test
make derive SCHEMA=omnipath_test
```
