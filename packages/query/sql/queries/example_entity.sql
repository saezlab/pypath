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

-- step: chosen AS table
SELECT * FROM entity_merged
WHERE entity_key IN ({{entity_keys*}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))
ORDER BY coalesce(entity_type IN ('physical_entity', 'named_thing', 'entity'), false),
    -coalesce(relation_count, 0),
    lower(coalesce(identifier, '')) NOT IN (lower({{token}}), lower(substr({{token}}, strpos({{token}}, ':') + 1))),
    entity_key
LIMIT 1

-- step: response AS value when (SELECT count(*) FROM __chosen) > 0
SELECT to_json(struct_pack(
    entity := (SELECT op_summary(summary) FROM __chosen),
    group_key := (SELECT CASE
        WHEN list_contains((SELECT list(entity_type) FROM op_chemical_types), entity_type) AND group_connectivity IS NOT NULL
            THEN 'connectivity:' || group_connectivity
        WHEN NOT list_contains((SELECT list(entity_type) FROM op_chemical_types), coalesce(entity_type, ''))
            AND starts_with(coalesce(reference_entity_key, ''), 'entrez:') THEN 'gene:' || reference_entity_key
        ELSE 'entity:' || entity_key END FROM __chosen)))
