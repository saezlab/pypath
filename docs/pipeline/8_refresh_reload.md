> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 8. Refresh And Reload

Refresh behavior is source-scoped. A normal load is additive; reload deletes one
or more sources before loading current parser output.

## Additive Load

```bash
make load SOURCE=bindingdb
```

Before loading a selected source, the pipeline checks whether that source
already has source-scoped content in the target schema. If content exists, the
source is skipped.

## Reload

```bash
make reload SOURCE=bindingdb
```

`make reload` sets `RELOAD_EXISTING=1` and then runs the same load path. Before
streaming parser output, existing source content is deleted.

## Delete Strategy

Reload deletion prefers dropping source partitions when all expected partitions
exist. If partitions are missing or the database uses older default-partition
rows, it falls back to row deletion.

After source evidence is removed, the refresh code garbage-collects canonical
relations, entities, identifiers, and annotations that are no longer supported
by any source. Shared graph rows backed by other sources are preserved.

## Manual Source Drop

```bash
make drop-source SOURCE=bindingdb
```

This performs the source deletion and garbage-collection part without reloading
the source.

## Main Files

- `Makefile`: `load`, `reload`, and `drop-source` targets.
- `omnipath_build/duckdb_direct_pipeline.py`: skip/reload decision.
- `omnipath_build/db/refresh.py`: source deletion and garbage collection.

## Phase Boundary

Reload refreshes source evidence and canonical graph rows. It does not rebuild
query-ready derived tables unless `make derive` is run afterward.

