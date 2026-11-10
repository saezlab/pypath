-- query: example_entity
-- The entity an example identifier names: of those resolved from it, the most
-- connected (a ChEBI id names the chemical and the ontology term it also is), then the
-- one with that identifier; a source's generic participant record comes last. With the
-- key of its group as the grouped view forms it. NULL when the release has none.
-- param: token VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"token": "entrez:7157"}

-- step: key_tokens AS list
SELECT {{token}}

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: copies AS table
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource, 0 AS position,
    NULL::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
    NULL::STRUCT(term VARCHAR, "value" VARCHAR, quantity STRUCT(has_numeric_value DOUBLE,
        has_unit VARCHAR, has_unit_prefix VARCHAR, has_binary_relation VARCHAR,
        source_field VARCHAR, comparator VARCHAR), "source" VARCHAR, dataset VARCHAR)[] AS annotations
FROM entity_by_key
WHERE entity_key IN ({{entity_keys*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: chosen AS table
WITH merged AS (SELECT * FROM {{include:merged_entity __copies}})
SELECT * FROM merged
ORDER BY coalesce(raw_type IN ('physical_entity', 'named_thing', 'entity'), false),
    -coalesce(relation_count, 0),
    lower(coalesce(raw_identifier, '')) NOT IN (lower({{token}}), lower(substr({{token}}, strpos({{token}}, ':') + 1))),
    pk
LIMIT 1

-- step: response AS value when (SELECT count(*) FROM __chosen) > 0
WITH shaped AS (SELECT * FROM {{include:entity_summary __chosen}})
SELECT to_json(struct_pack(
    entity := (SELECT summary FROM shaped),
    group_key := (SELECT CASE
        WHEN list_contains((SELECT list(entity_type) FROM op_chemical_types), raw_type) AND group_connectivity IS NOT NULL
            THEN 'connectivity:' || group_connectivity
        WHEN NOT list_contains((SELECT list(entity_type) FROM op_chemical_types), coalesce(raw_type, ''))
            AND starts_with(coalesce(reference_entity_key, ''), 'entrez:') THEN 'gene:' || reference_entity_key
        ELSE 'entity:' || pk END FROM __chosen)))
