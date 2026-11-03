# DuckDB Canonical Load Plan

## Current Situation

- The current DuckDB direct path lives in `omnipath_build/duckdb_direct_pipeline.py`
  and `omnipath_build/duckdb_load.py`.
  - Projects source records into DuckDB evidence tables.
  - Canonicalizes entities and annotation-derived relations in DuckDB.
  - Loads Postgres through either DuckDB's Postgres extension or experimental `COPY` payloads.
  - Contains the current source-partition staging + `ATTACH PARTITION` experiment.

- `omnipath_build/evidence_projector.py` contains the shared source-to-evidence
  flattening logic.
  - `DuckDBEvidenceProjector` writes those projected rows directly into DuckDB
    evidence tables.

- `omnipath_build/db/schema.py` defines the Postgres schema.
  - Canonical `entity` and `relation` IDs have been changed experimentally to deterministic UUIDs.
  - Source-specific evidence tables are partitioned by `source_id`.
  - The experimental fast path drops constraints/indexes before bulk loading and attaches source staging partitions.

- `omnipath_build/resolver/sources/chemicals.py` builds chemical resolver parquet files.
  - HMDB and ChEBI now populate `data/chemicals/chemical_identifier_lookup.parquet`.
  - DuckDB resolver views now union protein and chemical resolver rows when the chemical cache exists.

- `todo.md` is a short running list of immediate loose ends.

## Current Gaps

- ChEMBL activities currently project into `relation_evidence`, but DuckDB canonicalization does not yet create canonical `relation` rows from member-style `relation_evidence_raw`.
- The old build pipeline handles membership-style interactions; the DuckDB path needs parity for:
  - resolving subject/object members from relation evidence
  - creating canonical relation keys from resolved member entity IDs
  - populating `relation`
  - populating `relation_evidence_relation`
- The cache-first entity/relation anti-join flow is not implemented yet.
- The 50k ChEMBL batch benchmark is currently ad hoc; it should become a real runner.
- Constraints/index rebuild or validation after source-partition attach is still a separate phase to design.

## Target Direction

- DuckDB is the ingest engine for projection, canonicalization, cache checks, and load-ready payload construction.
- Postgres only receives typed `COPY` payloads and source-partition attach operations.
- The DuckDB/Parquet caches are part of the canonical build state. If that state is discarded, the database is rebuilt from scratch.

## Persistent DuckDB/Parquet Caches

- `entity_key_cache`
  - `entity_id`
  - `entity_key`
  - `entity_type`
  - `taxonomy_id`
  - `canonical_identifier_type`
  - `canonical_identifier`
  - `identifiers_json`
  - `sources`
  - `first_seen_at`
  - `last_seen_at`

- `relation_key_cache`
  - `relation_id`
  - `relation_key`
  - `subject_entity_id`
  - `predicate`
  - `object_entity_id`
  - `sources`
  - `first_seen_at`
  - `last_seen_at`

- dimension caches
  - `data_source`
  - `dataset`
  - identifier types
  - entity types
  - predicates
  - relation categories
  - resolution statuses
  - annotation scopes

## Entity And Relation Keys

- Resolved entity key:
  - `entity_type | taxonomy_id | canonical_identifier_type | canonical_identifier`

- Unresolved entity key:
  - hash of sorted identifier set for the evidence group

- Relation key:
  - `subject_entity_id | predicate | object_entity_id`

- IDs are deterministic UUIDs derived from these keys.

## Batch Flow

1. Project incoming source rows into DuckDB evidence tables.
   - `entity_evidence_raw`
   - `entity_identifier_raw`
   - `entity_annotation_raw`
   - `annotation_value`
   - `relation_evidence_raw`
   - `relation_annotation_raw`
   - `annotation_relation_evidence_raw`

2. Prefilter resolver cache to batch identifiers.
   - Build `evidence_identifier_key`.
   - Build `needed_resolver_lookup` by joining resolver rows only to keys present in the batch.

3. Resolve entities in DuckDB.
   - Produce `entity_resolution`.
   - Produce batch canonical entity candidates.
   - Assign deterministic `entity_id` from entity key.

4. Compare entity candidates against `entity_key_cache`.
   - `new_entities = batch_entity_candidates ANTI JOIN entity_key_cache`
   - `existing_entities = batch_entity_candidates JOIN entity_key_cache`
   - Build `entity` COPY payload only from `new_entities`.
   - Update `entity_key_cache.sources` to include the current source for all batch entity keys.

5. Build entity evidence resolution payload.
   - Every entity evidence row maps to an `entity_id`.
   - Build source-partition COPY payload for `entity_evidence_resolution`.

6. Build relation candidates in DuckDB.
   - Use resolved subject/object `entity_id`.
   - Build deterministic relation key and `relation_id`.
   - Include member-style `relation_evidence_raw` rows and annotation-derived relation rows.

7. Compare relation candidates against `relation_key_cache`.
   - `new_relations = batch_relation_candidates ANTI JOIN relation_key_cache`
   - `existing_relations = batch_relation_candidates JOIN relation_key_cache`
   - Build `relation` COPY payload only from `new_relations`.
   - Update `relation_key_cache.sources` to include the current source for all batch relation keys.

8. Build relation evidence link payload.
   - Every source relation evidence row maps to a `relation_id`.
   - Build source-partition COPY payload for `relation_evidence_relation`.

9. Build evidence and annotation payloads.
   - Source-partition payloads:
     - `entity_evidence`
     - `entity_evidence_identifier`
     - `entity_evidence_annotation`
     - `relation_evidence`
     - `relation_evidence_annotation`
   - Global deduplicated payloads:
     - `identifier_evidence`
     - `annotation`

10. Load into Postgres.
    - Global tables receive `COPY` payloads for new canonical/deduplicated rows only.
    - Source-specific tables load into staging tables created with `LIKE`.
    - Attach source staging tables with `ALTER TABLE ... ATTACH PARTITION ... FOR VALUES IN (<source_id>)`.

## Source Update Flow

1. Process the refreshed source into new DuckDB batch state.
2. Build canonical entity and relation keys against the local caches.
3. COPY only new canonical entities, relations, identifiers, and annotations.
4. Build complete replacement source partitions for evidence and evidence-link tables.
5. Swap source partitions in Postgres.
   - detach/drop old source partition
   - attach new staging partition
6. Update cache source membership.
   - add current source to keys present in the refreshed source
   - remove current source from keys no longer referenced by that source
7. Optionally garbage-collect canonical keys whose `sources` set is empty.

## Immediate Implementation Steps

- Add persistent `entity_key_cache` and `relation_key_cache` tables/files.
- Add batch entity key construction and anti-join against `entity_key_cache`.
- Add batch relation key construction for both member-style and annotation-derived relations.
- Add relation anti-join against `relation_key_cache`.
- Change canonical `entity` and `relation` COPY payloads to emit only cache misses.
- Keep source-partition staging and attach for all source-specific evidence tables.
- Add a 50k-batch ChEMBL runner using this cache-first flow.
