# Direct PostgreSQL Schema Decisions

This document records the current direction for simplifying the `omnipath_build`
PostgreSQL schema. It supersedes older assumptions in `PLAN.md` where they
conflict.

## Evidence Layer

- Keep `entity_evidence` and `relation_evidence` as source-scoped occurrence
  tables. They are still needed for refresh, canonicalization, and source
  provenance.
- Make evidence tables leaner. Remove legacy or redundant fields:
  - remove `snapshot_id`
  - remove `occurrence_id` and `relation_occurrence_id`
  - keep provenance as only `source_id`, `dataset_id`, and `row_id`
- Store taxonomy IDs as integers across the database.
- Normalize repeated dimensions where useful:
  - source, via `data_source`
  - dataset, via `dataset`
  - entity type, via `vocab_entity_type`
  - entity role, via `vocab_entity_role`, if still needed
  - predicate, via `vocab_relation_predicate`
  - relation category, via `vocab_relation_category`
  - annotation scope, via `vocab_annotation_scope`

## Entity Evidence

`entity_evidence` should remain the bridge between raw source rows and
canonical entities.

Its purpose is to hold the evidence participant, its source row provenance,
type, taxonomy, and links to identifiers and annotations. Identifier values
belong in `identifier_evidence` and `entity_evidence_identifier`, not repeated
on the entity evidence row.

## Relation Evidence

`relation_evidence` should remain the bridge between raw source rows and
canonical relations.

Its provenance only needs `source_id`, `dataset_id`, and `row_id`. Endpoints should
continue to point either to source entity evidence or directly to already
canonical entities when the object is known at ingest time, such as ontology
term relations.

Predicate and relation category should be normalized instead of stored as
repeated text labels.

## Annotations

The primary annotation model should be evidence-first:

- `annotation` stores the deduplicated annotation value.
- evidence annotation link tables attach annotation values to entity or
  relation evidence.
- canonical `entity_annotation` and `relation_annotation` should be removed
  from the primary schema for now.
- If canonical annotation tables are needed later for query performance, they
  should be rebuilt as derived tables.

Relation annotation scope should be normalized and limited to:

- `relation`
- `subject`
- `object`

This is enough to distinguish attributes of the assertion itself from
participant-specific attributes.

## Source Refresh And Partitioning

Partitioning by source is a good fit for the evidence layer. Partition these
tables, or their normalized equivalents, by `source_id`:

- `entity_evidence`
- `entity_evidence_identifier`
- `relation_evidence`
- evidence annotation link tables

Canonical `entity` and `relation` tables should stay global and deduplicated.
Refreshing a source should remove that source's evidence partitions, recompute
affected evidence-to-canonical mappings, and garbage-collect canonical graph
rows that no longer have supporting evidence.

## Implementation Notes

- Source partitioning should be native PostgreSQL partitioning
- `data_source` and `dataset` are provenance dimensions, not controlled
  vocabularies. Keep them separate from the `vocab_` naming convention.
- Use separate `vocab_` tables for controlled vocabularies rather than one
  generic vocabulary table, so foreign keys remain precise.
- Update canonicalization so it does not derive canonical annotation tables.
