-- query: stats_entity_types
-- Entity rows per entity type, most first.
-- param: resources VARCHAR[] = NULL
-- output: json response
-- output: table counts
-- example: {}

-- step: counts AS table
SELECT entity_type, count(*) AS n FROM entity
WHERE entity_type IS NOT NULL AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
GROUP BY 1 ORDER BY n DESC, entity_type

-- step: response AS value
SELECT to_json(coalesce(list(struct_pack(entityType := entity_type, count := n) ORDER BY n DESC, entity_type), []))
FROM __counts
