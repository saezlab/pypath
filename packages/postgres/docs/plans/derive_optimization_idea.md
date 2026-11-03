Yes. The best next work is to add explicit “dirty sets” from refresh/canonicalize, then have derive consume only those. Some derive work is safely skippable now; some is global unless we track what changed.

Highest value candidates:

1. `entity_relation_counts`
Currently truncates and rebuilds from all `relation` endpoints. This can be scoped cleanly.
Track dirty `relation_id`s and dirty endpoint `entity_id`s during canonicalize/refresh, then recompute counts only for affected entities:
- endpoints of inserted/deleted/changed relations
- entities whose evidence resolution changed
- entities removed from canonical graph

2. `resources`
This loops all discovered resources and recounts every present source. Counts are source-local, so derive can accept `--source` or consume `derive_dirty_source`.
For a one-source refresh, update only that source’s resource row. Static metadata sync can stay global but counts do not need to.

3. Bitmap source facets
These are very scopeable:
- `facet_entity_bitmap` where `facet_name = 'source'`
- `facet_relation_bitmap` where `facet_name = 'source'`

If one source was refreshed, delete and rebuild only that source’s bitmap rows.

4. Other bitmap facets
These need dirty entity/relation sets, but are still feasible:
- entity `entity_type` / `taxonomy_id`: recompute only facet values touched by dirty entities
- relation `predicate` / `participant_type` / `taxonomy_id`: recompute only facet buckets touched by dirty relations/endpoints
- annotation-term bitmaps: hardest, because changed relations can affect term-to-relation bitmaps through neighboring annotated entities

5. `ontology_terms`
Can be scoped to dirty CV-term entities. Recompute only terms whose evidence annotations changed or whose term entity was newly created/removed.

6. Index creation
`CREATE INDEX IF NOT EXISTS` is mostly metadata work once indexes exist. Not urgent, but we can avoid even issuing these statements by checking `pg_class` first or making `--indexes` default smarter.

The right architecture is probably:
- add small tables like `derive_dirty_source`, `derive_dirty_entity`, `derive_dirty_relation`, `derive_dirty_term`
- populate them from `delete_source_content` before deletes and from canonicalize after inserts/upserts
- derive consumes and clears them transactionally
- keep a full rebuild path as fallback

One caveat on the identifier guard we just added: if resolver tables are reloaded, resolved entities with existing identifier arrays may become stale. That is fine for “fill missing once”, but if resolver data changes we need either a `--refresh-entity-identifiers` flag or a resolver-version invalidation that marks resolved entities dirty.