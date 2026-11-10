-- query: entity_filter_options
-- The values the entity filters offer: entity types, resources and taxa.
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {}

-- step: types AS list
SELECT DISTINCT entity_type FROM entity
WHERE entity_type IS NOT NULL AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY 1

-- step: taxa AS list
SELECT DISTINCT taxon::VARCHAR FROM entity
WHERE taxon IS NOT NULL AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY 1

-- step: sources AS list
SELECT resource_id FROM resources
WHERE op_unset({{resources}}) OR list_contains({{resources}}, resource_id)
ORDER BY 1

-- step: response AS value
SELECT to_json(struct_pack(entity_types := coalesce({{types}}, []::VARCHAR[]),
    sources := coalesce({{sources}}, []::VARCHAR[]), taxonomy_ids := coalesce({{taxa}}, []::VARCHAR[])))
