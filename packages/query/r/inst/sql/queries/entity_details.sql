-- query: entity_details
-- The explorer's entity details: the entity, its first
-- 50 structural relationships with both endpoints, its relation count and up to 100
-- association relations, as the API's JSON. NULL when no selected resource has it.
-- param: entity_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- steps: entity_main

-- step: endpoint_entries AS table
SELECT DISTINCT relation_n, predicate, category, page_rank FROM relation_endpoint
WHERE entity_n IN ({{main_ids*}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: association_ids AS list
SELECT relation_n FROM (SELECT * FROM __endpoint_entries WHERE category = 'association'
                        ORDER BY predicate, page_rank LIMIT 100)
ORDER BY relation_n

-- step: association_rows AS table
SELECT relation_key, predicate FROM relation WHERE relation_n IN ({{association_ids*}})

-- step: rel_keys AS list
SELECT {{entity_key}}

-- step: rel_limit AS value
SELECT 50

-- step: rel_offset AS value
SELECT 0

-- steps: structural_relationships

-- step: response AS value when {{entity_json}} IS NOT NULL
SELECT to_json(struct_pack(
    entity := {{entity_json}}::JSON,
    relationships := coalesce((SELECT list(item ORDER BY position) FROM __rel_items), []),
    relationshipsTotal := coalesce((SELECT any_value(total) FROM __rel_page), 0),
    summary := struct_pack(interactionCount := (SELECT count(DISTINCT relation_n) FROM __endpoint_entries)),
    annotations := coalesce((SELECT list(struct_pack(relationPk := relation_key, predicate := predicate)
            ORDER BY relation_key) FROM __association_rows), [])
))
