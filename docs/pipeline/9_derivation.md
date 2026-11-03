> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 9. Derivation

Derivation rebuilds query-facing artifacts after selected sources have been
loaded and canonicalized.

## Command

```bash
make derive
```

## What It Does

The `derive` command:

1. Ensures content primary keys.
2. Ensures the base schema exists.
3. Creates deferred and secondary indexes when enabled.
4. Rebuilds derived search/count tables.
5. Rebuilds the entity identifier lookup used by app search.
6. Rebuilds roaring bitmap tables.
7. Discovers pypath resources again.
8. Syncs the `resources` summary table.

## Derived Tables

Derived tables summarize canonical graph content for search and counts. Core
derived tables include `entity_identifier_lookup`, which joins canonical
entities to deduplicated rows in `identifier_evidence`, and
`entity_relation_counts`; ontology-term summary data is also maintained in this
phase.

## Bitmap Tables

Bitmap tables encode ontology annotations, source facets, entity types, and
relation categories as compressed sets of canonical entity or relation bitmap
IDs.

## Main Files

- `omnipath_build/cli.py`: `derive` command.
- `omnipath_build/db/derived_tables.py`: derived table rebuild.
- `omnipath_build/db/bitmaps.py`: roaring bitmap rebuild.
- `omnipath_build/db/resources.py`: resource summary sync.
- `omnipath_build/db/indexes.py`: deferred and secondary indexes.

## Phase Boundary

`make load` and `make reload` do not make the database fully query-ready on
their own. Run `make derive` once after a source batch is complete.
