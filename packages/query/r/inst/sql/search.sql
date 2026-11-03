-- The explorer's entity search as a DuckDB macro over the installed database.
--
-- Mirrors omnipath_api.queries.services.entities.search_entities_api and the columns
-- omnipath_client.Client.search returns. Needs relations.sql (op_unset, op_biolink,
-- op_entity_keys).
--
-- Phases, as in the API's queries: an empty query lists entities by key; a 64-character key
-- finds that entity; otherwise label and identifier prefix matches (``entity_term``),
-- and only when there are none, labels containing the query.
--
-- One difference by design: the API's search reads prefix hits in growing batches and stops
-- once the page is certain; it counts an entity's relations over the hits read so far.
-- Here every hit is ranked at once, so ties between equally ranked entities can come
-- out in a different order when a query has very many hits.

-- An entity row passes the search filters; NULL or empty lists filter nothing, and
-- ``keys`` is NULL when no entity filter is given.
CREATE OR REPLACE MACRO op_entity_ok(e, scope, type_names, taxa, references_, keys) AS (
    (op_unset(scope) OR list_contains(scope, e._resource))
    AND (op_unset(type_names) OR list_contains(type_names, e.entity_type))
    AND (op_unset(taxa) OR list_contains(taxa, e.taxon))
    AND (op_unset(references_) OR list_contains(references_, e.reference_entity_key))
    AND (keys IS NULL OR list_contains(keys, e.entity_key))
);

-- The query's phase and parts, from the macro argument itself: as constants, the
-- conditions on them are folded while planning, and the term range reaches the scan.
CREATE OR REPLACE MACRO op_search_mode(query) AS (CASE
    WHEN trim(coalesce(query, '')) = '' THEN 'all'
    WHEN regexp_full_match(trim(query), '[0-9a-fA-F]{64}') THEN 'key'
    ELSE 'text' END);
CREATE OR REPLACE MACRO op_search_text(query) AS (lower(trim(coalesce(query, ''))));
CREATE OR REPLACE MACRO op_prefix_upper(ql) AS (CASE WHEN ql = '' THEN NULL
    ELSE left(ql, length(ql) - 1) || chr(unicode(right(ql, 1)) + 1) END);
CREATE OR REPLACE MACRO op_curie_ns(query) AS (
    nullif(lower(regexp_extract(trim(coalesce(query, '')), '^([A-Za-z][\w.-]*)[:|](.+)$', 1)), ''));
CREATE OR REPLACE MACRO op_curie_id(query) AS (
    nullif(regexp_extract(trim(coalesce(query, '')), '^([A-Za-z][\w.-]*)[:|](.+)$', 2), ''));

CREATE OR REPLACE MACRO op_search(
    query,
    sources := NULL,
    entity_types := NULL,
    taxonomy_ids := NULL,
    reference_entity_keys := NULL,
    entity_ids := NULL,
    max_results := 20
) AS TABLE
WITH s AS MATERIALIZED (
    SELECT
        trim(coalesce(query, '')) AS q,
        lower(trim(coalesce(query, ''))) AS ql,
        CASE
            WHEN trim(coalesce(query, '')) = '' THEN 'all'
            WHEN regexp_full_match(trim(query), '[0-9a-fA-F]{64}') THEN 'key'
            ELSE 'text'
        END AS mode,
        nullif(lower(regexp_extract(trim(coalesce(query, '')), '^([A-Za-z][\w.-]*)[:|](.+)$', 1)), '') AS curie_ns,
        nullif(regexp_extract(trim(coalesce(query, '')), '^([A-Za-z][\w.-]*)[:|](.+)$', 2), '') AS curie_id,
        sources::VARCHAR[] AS scope,
        (SELECT list(name) FROM op_biolink('class', entity_types)) AS type_names,
        taxonomy_ids::VARCHAR[] AS taxa,
        reference_entity_keys::VARCHAR[] AS references_,
        CASE WHEN op_unset(entity_ids::VARCHAR[]) THEN NULL ELSE (
            SELECT coalesce(list(entity_key), []) FROM op_entity_keys(entity_ids, scope := sources::VARCHAR[])
        ) END AS keys,
        max_results::BIGINT AS wanted
), bounds AS MATERIALIZED (
    -- [ql, upper): the terms starting with the query, as a range on the sorted terms
    SELECT *, CASE WHEN ql = '' THEN NULL
        ELSE left(ql, length(ql) - 1) || chr(unicode(right(ql, 1)) + 1) END AS upper
    FROM s
), entities AS NOT MATERIALIZED (
    -- every entity passing the filters: read only by the substring and listing phases
    SELECT e.*, e._resource || '/' || e._resource_version AS resource
    FROM entity e, s
    WHERE (op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], e._resource))
        AND op_entity_ok(e, s.scope, s.type_names, s.taxa, s.references_, s.keys)
), terms AS NOT MATERIALIZED (
    SELECT t.* FROM entity_term t
    WHERE op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], t._resource)
), hits AS (
    -- (resource, entity_id, rank, term, namespace the entity must have or NULL)
    SELECT t._resource, t.entity_id,
        CASE WHEN t.term = b.ql THEN (CASE WHEN t.kind = 'label' THEN 1 ELSE 2 END)
            ELSE (CASE WHEN t.kind = 'label' THEN 3 ELSE 4 END) END AS match_rank,
        t.term, NULL::VARCHAR AS namespace
    FROM terms t, bounds b
    WHERE op_search_mode(query) = 'text' AND t.term >= op_search_text(query)
        AND t.term < op_prefix_upper(op_search_text(query))
    UNION ALL
    SELECT t._resource, t.entity_id,
        CASE WHEN t.kind = 'label' THEN 1 ELSE 2 END, t.term, b.curie_ns
    FROM terms t, bounds b
    WHERE op_search_mode(query) = 'text' AND op_curie_ns(query) IS NOT NULL
        AND t.term = lower(op_curie_id(query))
    UNION ALL
    SELECT g._resource, g.entity_id, 2, b.ql, NULL
    FROM entity_group g, bounds b
    WHERE op_search_mode(query) = 'text' AND op_curie_ns(query) IS NOT NULL AND g.kind = 'reference'
        AND g.group_key = op_curie_ns(query) || ':' || op_curie_id(query)
        AND (op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], g._resource))
), hit_rows AS MATERIALIZED (
    SELECT _resource, entity_id,
        min(match_rank) AS match_rank, min(length(term)) AS term_length, min(term) AS term,
        list(DISTINCT namespace) FILTER (WHERE namespace IS NOT NULL) AS namespaces,
        bool_or(namespace IS NULL) AS any_namespace
    FROM hits GROUP BY ALL
), hit_entities AS (
    SELECT e.* FROM op_entity_at((SELECT list(struct_pack(r := _resource, v := entity_id)) FROM hit_rows)) e, s
    WHERE op_entity_ok(e, s.scope, s.type_names, s.taxa, s.references_, s.keys)
), prefix_matched AS (
    SELECT e.entity_key, e.relation_count,
        struct_pack(r := h.match_rank, l := h.term_length, t := h.term) AS hit
    FROM hit_rows h JOIN hit_entities e USING (_resource, entity_id)
    WHERE h.any_namespace OR list_contains(h.namespaces, lower(coalesce(e.namespace, '')))
), prefix_page AS MATERIALIZED (
    SELECT entity_key,
        row_number() OVER (ORDER BY min(hit).r, min(hit).l,
            sum(coalesce(relation_count, 0)) DESC, min(hit).t, entity_key) AS position
    FROM prefix_matched GROUP BY entity_key
    QUALIFY position <= (SELECT wanted FROM s)
), prefix_rows AS (
    -- every resource's copy of the page's entities that passes the filters
    SELECT e.*, e._resource || '/' || e._resource_version AS resource, p.position
    FROM prefix_page p
    JOIN op_entity_with((SELECT list(entity_key) FROM prefix_page)) e USING (entity_key), s
    WHERE op_entity_ok(e, s.scope, s.type_names, s.taxa, s.references_, s.keys)
), contains_gate AS (
    -- one row only when the prefix phase found nothing: joined first, an empty gate
    -- ends the scan of every entity's label
    SELECT 1 AS open WHERE op_search_mode(query) = 'text' AND NOT EXISTS (SELECT 1 FROM prefix_page)
), contains_matched AS (
    SELECT e.* FROM contains_gate g, entities e
    WHERE contains(lower(e.label), op_search_text(query))
), contains_page AS (
    SELECT entity_key,
        row_number() OVER (ORDER BY length(coalesce(min(label), '')),
            sum(relation_count) DESC, coalesce(min(label), ''), entity_key) AS position
    FROM contains_matched GROUP BY entity_key
    QUALIFY position <= (SELECT wanted FROM s)
), contains_rows AS (
    SELECT m.*, p.position FROM contains_page p JOIN contains_matched m USING (entity_key)
), listed AS (
    SELECT e.* FROM entities e
    WHERE op_search_mode(query) = 'all'
        OR (op_search_mode(query) = 'key' AND e.entity_key = trim(query))
), listed_page AS (
    SELECT entity_key, row_number() OVER (ORDER BY entity_key) AS position
    FROM listed GROUP BY entity_key
    QUALIFY position <= (SELECT wanted FROM s)
), listed_rows AS (
    SELECT l.*, p.position FROM listed_page p JOIN listed l USING (entity_key)
), copies AS (
    SELECT * FROM prefix_rows
    UNION ALL BY NAME SELECT * FROM contains_rows
    UNION ALL BY NAME SELECT * FROM listed_rows
), merged AS (
    -- One entity from its resource copies: the richest copy's columns, the relation
    -- counts summed, gene reference keys united and one reference key chosen.
    SELECT entity_key, min(position) AS position,
        arg_max(struct_pack(
            namespace := namespace, identifier := identifier, label := label,
            entity_type := entity_type, taxon := taxon), struct_pack(
            a := coalesce(annotation_count, 0),
            i := coalesce(identifier_count, 0),
            g := regexp_full_match(coalesce(label, ''), '[A-Z][A-Z0-9-]{1,14}')::INTEGER,
            n := (coalesce(label, '') <> '' AND coalesce(label, '') <> coalesce(identifier, '')
                AND NOT regexp_full_match(coalesce(label, ''), '[0-9]+'))::INTEGER,
            f := lower(coalesce(label, '')),
            l := coalesce(label, ''),
            k := entity_key,
            r := resource,
            h := coalesce(has_hierarchy, false),
            p := coalesce(parent_count, 0),
            c := coalesce(child_count, 0))) AS best,
        sum(coalesce(relation_count, 0))::BIGINT AS relation_count,
        list_sort(list_distinct(flatten(list(coalesce(gene_reference_keys, []))))) AS gene_reference_keys,
        list(DISTINCT reference_entity_key)
            FILTER (WHERE coalesce(reference_entity_key, '') <> '') AS reference_keys
    FROM copies GROUP BY entity_key
), shaped AS (
    SELECT *,
        trim(coalesce(nullif(best.identifier, ''), entity_key)) AS shown_identifier,
        coalesce(best.namespace, 'None') || ':' || coalesce(best.identifier, 'None') AS native
    FROM merged
)
SELECT
    entity_key,
    trim(coalesce(nullif(best.label, ''), shown_identifier)) AS label,
    lower(trim(coalesce(nullif(best.entity_type, ''), 'unknown'))) AS entity_type,
    lower(trim(coalesce(nullif(best.namespace, ''), 'unknown'))) AS namespace,
    shown_identifier AS identifier,
    nullif(coalesce(best.taxon, ''), '') AS taxon,
    relation_count,
    CASE
        WHEN list_contains(coalesce(reference_keys, []), native) OR len(coalesce(reference_keys, [])) > 1 THEN native
        WHEN len(coalesce(reference_keys, [])) = 1 THEN reference_keys[1]
    END AS reference_entity_key,
    gene_reference_keys
FROM shaped
ORDER BY position;
