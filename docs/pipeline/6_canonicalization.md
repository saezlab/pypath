> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 6. Canonicalization

Canonicalization turns source evidence occurrences into deduplicated canonical
entities and relations.

## Where It Runs

Canonicalization is part of `make load`; it runs after projection for each
DuckDB batch or staged DuckDB file.

## Entity Canonicalization

The DuckDB canonicalization step:

1. Builds evidence identifier keys from projected identifiers.
2. Reads resolver parquet files through DuckDB views.
3. Matches evidence identifiers to resolver candidates by entity type and
   taxonomy where applicable.
4. Resolves CV-term entities directly from their controlled-vocabulary
   accession.
5. Chooses canonical identifiers for resolved entities.
6. Produces unresolved fallback identifiers when no resolver candidate is
   accepted.
7. Builds special signatures for complex and reaction member structures.

Canonical entity IDs are deterministic UUIDs based on canonical entity keys.
This keeps repeated ingests of the same source snapshot stable.

## Relation Canonicalization

After entities are resolved, relation canonicalization replaces relation
evidence endpoints with canonical entity IDs and deduplicates graph relations by
subject, predicate, and object.

## Main Files

- `omnipath_build/duckdb_direct_pipeline.py`: calls canonicalization.
- `omnipath_build/duckdb_load.py`: `_canonicalize_loaded_duckdb` and supporting
  SQL.
- `omnipath_build/relation_rules.py`: relation predicate/category rules.
- `omnipath_build/cv_terms.py`: controlled vocabulary constants.

## Phase Boundary

Canonicalization still happens in DuckDB. PostgreSQL is not populated until the
copy phase.

