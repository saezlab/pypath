-- query: entity_facets
-- The entity filter options a search leaves: entity types and taxa of the matching
-- entities (most first), and every selected resource's count of them, which ignores
-- the source filter so that other sources stay selectable.
-- param: query VARCHAR = ""
-- param: sources VARCHAR[] = NULL
-- param: entity_types VARCHAR[] = NULL
-- param: taxonomy_ids VARCHAR[] = NULL
-- param: reference_entity_keys VARCHAR[] = NULL
-- param: entity_ids VARCHAR[] = NULL
-- param: facet_limit INTEGER = 50
-- param: selection VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"query": "TP53"}

-- step: q AS value
SELECT trim(coalesce({{query}}, ''))

-- step: ql AS value
SELECT lower({{q}})

-- step: q_upper AS value
SELECT CASE WHEN {{ql}} = '' THEN NULL
    ELSE left({{ql}}, length({{ql}}) - 1) || chr(unicode(right({{ql}}, 1)) + 1) END

-- step: curie_ns AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN lower(regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 1)) END

-- step: curie_id AS value
SELECT CASE WHEN regexp_full_match({{q}}, '([A-Za-z][\w.-]*)[:|](.+)')
    THEN regexp_extract({{q}}, '^([A-Za-z][\w.-]*)[:|](.+)$', 2) END

-- step: phases AS value
SELECT CASE WHEN {{q}} = '' THEN ['all'] WHEN regexp_full_match({{q}}, '[0-9a-fA-F]{64}') THEN ['key']
    ELSE ['prefix', 'contains'] END

-- step: type_names AS value
SELECT (SELECT list(name) FROM op_biolink('class', {{entity_types}}))

-- step: selected AS value
-- the resources a request names (selection), of those the query policy leaves
SELECT CASE WHEN op_unset({{selection}}) THEN {{resources}}
    WHEN op_unset({{resources}}) THEN {{selection}}
    ELSE coalesce(nullif(list_intersect({{selection}}, {{resources}}), []), ['']) END

-- step: typed_scope AS value
-- the source filter, else the selected resources
SELECT CASE WHEN op_unset({{sources}}) THEN {{selected}} ELSE {{sources}} END

-- step: key_tokens AS list
SELECT unnest(coalesce({{entity_ids}}, []))

-- step: key_scope AS value
SELECT {{typed_scope}}

-- steps: entity_keys

-- step: match_keys AS value
SELECT CASE WHEN op_unset({{entity_ids}}) THEN NULL ELSE coalesce({{entity_keys}}, []) END

-- step: match_scope AS value
SELECT {{typed_scope}}

-- steps: entity_matches

-- step: typed_counts AS table
-- the matching entities by type, taxon and resource: without a source filter the same
-- entities give the source counts
WITH m AS {{include:entity_match_rows m}}
SELECT 'entity_type' AS facet, entity_type::VARCHAR AS value, count(*) AS n FROM m GROUP BY entity_type
UNION ALL
SELECT 'taxonomy_id', taxon::VARCHAR, count(*) FROM m GROUP BY taxon
UNION ALL
SELECT 'source', r, count(*) FROM (SELECT unnest(resources) AS r FROM m) GROUP BY r

-- step: key_scope AS value
SELECT {{selected}}

-- steps: entity_keys

-- step: match_keys AS value
SELECT CASE WHEN op_unset({{entity_ids}}) THEN NULL ELSE coalesce({{entity_keys}}, []) END

-- step: match_scope AS value
SELECT {{selected}}

-- steps: entity_matches

-- step: filtered_source_counts AS table when NOT op_unset({{sources}})
SELECT r AS value, count(*) AS n
FROM (SELECT unnest(resources) AS r FROM {{include:entity_match_rows m}} m) GROUP BY r

-- step: source_counts AS table
SELECT r.resource_id AS value, coalesce(c.n, 0) AS n
FROM resources r LEFT JOIN (
    SELECT value, n FROM __filtered_source_counts WHERE NOT op_unset({{sources}})
    UNION ALL
    SELECT value, n FROM __typed_counts WHERE facet = 'source' AND op_unset({{sources}})) c
    ON c.value = r.resource_id
WHERE op_unset({{resources}}) OR list_contains({{resources}}, r.resource_id)

-- step: typed_any AS value
-- a selection without any resource of the release has no options
SELECT op_unset({{typed_scope}}) OR (SELECT count(*) FROM resources WHERE list_contains({{typed_scope}}, resource_id)) > 0

-- step: response AS value when {{typed_any}}
WITH typed AS (
    SELECT facet, value, n, row_number() OVER (ORDER BY n DESC, facet, value) AS position
    FROM __typed_counts WHERE facet <> 'source' AND value IS NOT NULL AND value <> ''
), shown AS (
    SELECT * FROM typed WHERE position <= greatest(coalesce({{facet_limit}}, 50) * 2, 20)
), items AS (
    SELECT position, CASE WHEN facet = 'taxonomy_id'
        THEN to_json(struct_pack(facetName := facet, facetValue := value, scopedCount := n,
            facetLabel := (SELECT name FROM op_taxonomy WHERE taxon = value)))
        ELSE to_json(struct_pack(facetName := facet, facetValue := value, scopedCount := n)) END AS item
    FROM shown
    UNION ALL
    SELECT 1000000 + row_number() OVER (ORDER BY n DESC, value),
        to_json(struct_pack(facetName := 'source', facetValue := value, scopedCount := n))
    FROM __source_counts
)
SELECT to_json(coalesce(list(item ORDER BY position), []::JSON[])) FROM items
