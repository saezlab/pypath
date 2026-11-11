-- query: ontology_search
-- Ontology terms with the entities and relations they annotate, most annotated first;
-- scoped to entities (entity_pks, term_ids, or a selection scope) when given. query
-- filters by label or term id, ontology_ids and prefixes by ontology; a page of limit
-- from offset, or (counts) the number of terms per ontology.
-- param: query VARCHAR = ""
-- param: limit INTEGER = 24
-- param: offset INTEGER = 0
-- param: ontology_ids VARCHAR[] = NULL
-- param: prefixes VARCHAR[] = NULL
-- param: counts BOOLEAN = false
-- param: scoped BOOLEAN = false
-- param: entity_pks VARCHAR[] = NULL
-- param: term_ids VARCHAR[] = NULL
-- param: has_selection BOOLEAN = false
-- param: selection_entity_pks VARCHAR[] = NULL
-- param: selection_include_associated BOOLEAN = true
-- param: selection_term_ids VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"query": "apoptosis", "limit": 5}

-- step: scope_entity_pks AS value
SELECT CASE WHEN {{has_selection}} THEN {{selection_entity_pks}} END

-- step: scope_include_associated AS value
SELECT {{selection_include_associated}}

-- steps: selection_scope

-- step: keys AS value
SELECT CASE WHEN {{has_selection}} THEN coalesce({{scope_keys}}, []) ELSE coalesce({{entity_pks}}, []) END

-- step: terms AS value
SELECT CASE WHEN {{has_selection}} THEN coalesce({{selection_term_ids}}, []) ELSE coalesce({{term_ids}}, []) END

-- step: is_scoped AS value
SELECT {{has_selection}} OR {{scoped}} OR len({{keys}}) > 0 OR len({{terms}}) > 0

-- step: scope_rows AS table
-- the scope's entities as a table: a selection with its neighbours has thousands, and
-- checking each of the ten million term_entity rows against a list of them took 75 s
-- (a gene's selection); joining them takes 2-3 s
SELECT e.entity_n, e.entity_key FROM entity e
JOIN (SELECT DISTINCT k FROM (SELECT unnest({{keys}}) AS k) WHERE k IS NOT NULL) s ON s.k = e.entity_key

-- step: term_rows AS table
SELECT DISTINCT t AS term FROM (SELECT unnest({{terms}}) AS t) WHERE t IS NOT NULL

-- step: items AS table cache
-- every term with its counts (or the scope's): kept, so that searches only filter it
WITH terms AS (
    SELECT entity_n, entity_key, identifier, label, namespace
    FROM entity
    WHERE entity_type = 'ontology_class' AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))
        AND NOT ({{is_scoped}} AND len({{keys}}) = 0 AND len({{terms}}) = 0)
), edges AS (
    SELECT term_n, entity_n, relation_n FROM term_entity
    WHERE (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))
        AND (NOT {{is_scoped}} OR entity_n IN (SELECT entity_n FROM __scope_rows))
), counts AS (
    SELECT e.term_n, count(DISTINCT e.entity_n) AS entities, count(DISTINCT e.relation_n) AS relations
    FROM edges e JOIN terms t ON t.entity_n = e.term_n GROUP BY e.term_n
)
SELECT t.* EXCLUDE (entity_n), coalesce(a.entities, 0) AS entity_count, coalesce(a.relations, 0) AS relation_count,
    row_number() OVER (ORDER BY coalesce(a.entities, 0) DESC, t.identifier, t.entity_key) AS position
FROM terms t LEFT JOIN counts a ON a.term_n = t.entity_n
WHERE NOT {{is_scoped}} OR a.term_n IS NOT NULL
    OR t.entity_key IN (SELECT entity_key FROM __scope_rows)
    OR t.entity_key IN (SELECT term FROM __term_rows) OR t.identifier IN (SELECT term FROM __term_rows)

-- step: matching AS table
SELECT *, coalesce(identifier, entity_key) AS term_id,
    CASE WHEN contains(coalesce(identifier, ''), ':') THEN split_part(identifier, ':', 1) ELSE namespace END AS prefix,
    coalesce(namespace, 'unknown') AS ontology
FROM __items
WHERE trim(coalesce({{query}}, '')) = ''
    OR contains(lower(coalesce(label, '')), lower(trim({{query}})))
    OR contains(lower(coalesce(identifier, entity_key)), lower(trim({{query}})))

-- step: response AS value
SELECT CASE WHEN {{counts}} THEN (
    SELECT to_json(coalesce(list(struct_pack(ontologyId := ontology, scopedCount := n) ORDER BY n DESC, ontology), []))
    FROM (SELECT ontology, count(*) AS n FROM __matching GROUP BY ontology))
ELSE (
    SELECT to_json(coalesce(list(struct_pack(
            entityPk := entity_key, termId := term_id, ontologyPrefix := prefix, ontologyId := ontology,
            label := label, definition := NULL::VARCHAR, synonyms := []::VARCHAR[], sources := []::VARCHAR[],
            annotatedEntityCount := entity_count, annotatedRelationCount := relation_count,
            annotatedItemCount := entity_count + relation_count) ORDER BY position), []))
    FROM (SELECT * FROM __matching
          WHERE (op_unset({{ontology_ids}}) OR list_contains({{ontology_ids}}, ontology))
              AND (op_unset({{prefixes}}) OR list_contains({{prefixes}}, prefix))
          ORDER BY position
          LIMIT greatest(1, least(coalesce({{limit}}, 24), 100)) OFFSET greatest(0, coalesce({{offset}}, 0))))
END
