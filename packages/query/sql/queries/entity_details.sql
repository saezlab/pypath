-- query: entity_details
-- The explorer's entity details: the entity merged from its resource copies, its first
-- 50 structural relationships with both endpoints, its relation count and up to 100
-- association relations, as the API's JSON. NULL when no selected resource has it.
-- param: entity_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- steps: entity_main

-- step: endpoint_entries AS table
SELECT DISTINCT _rid, resource, relation_id, predicate, category
FROM relation_endpoint
WHERE key = {{entity_key}} AND key_kind = 'entity'
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: association_page AS table
SELECT _rid FROM (SELECT DISTINCT _rid, resource, relation_id, predicate FROM __endpoint_entries
                  WHERE category = 'association')
ORDER BY predicate, resource, relation_id LIMIT 100

-- step: association_ids AS list
SELECT _rid FROM __association_page ORDER BY _rid

-- step: association_rows AS table
SELECT relation_key, predicate, _resource, _resource_version FROM relation WHERE _rid IN ({{association_ids*}})

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
    summary := struct_pack(interactionCount := (SELECT count(*) FROM __endpoint_entries)),
    annotations := coalesce((SELECT list(struct_pack(relationPk := relation_key, predicate := predicate)
            ORDER BY relation_key, _resource || '/' || _resource_version)
        FROM __association_rows), [])
))
