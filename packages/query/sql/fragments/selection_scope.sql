-- Steps: an entity selection (the values scope_entity_pks and scope_include_associated):
-- the seeds, their keys (selected_keys) and the relations to their neighbours
-- (__neighbour_rows); the scope_keys value lists the keys, then the neighbours.

-- step: seeds AS list
SELECT t FROM (SELECT trim(unnest(coalesce({{scope_entity_pks}}, []))) AS t, generate_subscripts(coalesce({{scope_entity_pks}}, []), 1) AS i)
WHERE t <> '' ORDER BY i

-- step: key_tokens AS list
SELECT unnest({{seeds}})

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: selected_keys AS value
SELECT CASE WHEN len(coalesce({{entity_keys}}, [])) > 0 THEN {{entity_keys}} ELSE {{seeds}} END

-- step: selected_ns AS list when {{scope_include_associated}} AND len(coalesce({{selected_keys}}, [])) > 0
SELECT entity_n FROM entity WHERE entity_key IN ({{selected_keys*}}) ORDER BY entity_n

-- step: neighbour_endpoints AS table
SELECT relation_n FROM relation_endpoint
WHERE entity_n IN ({{selected_ns+}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: neighbour_ids AS list
SELECT DISTINCT relation_n FROM __neighbour_endpoints ORDER BY 1

-- step: neighbour_rows AS table
SELECT relation_n, subject_n, object_n FROM relation WHERE relation_n IN ({{neighbour_ids*}})

-- step: neighbour_ns AS list
-- the other ends, in the order of their first relation (by number), subject first
SELECT n FROM (
    SELECT n, min(struct_pack(r := relation_n, s := side)) AS first
    FROM (SELECT relation_n, 1 AS side, subject_n AS n FROM __neighbour_rows
          UNION ALL SELECT relation_n, 2, object_n FROM __neighbour_rows)
    WHERE NOT list_contains({{selected_ns}}, n)
    GROUP BY n)
ORDER BY first

-- step: neighbour_keys AS table
SELECT entity_n, entity_key FROM entity WHERE entity_n IN ({{neighbour_ns*}})

-- step: scope_keys AS value
SELECT coalesce({{selected_keys}}, []) || coalesce((
    SELECT list(k.entity_key ORDER BY o.i) FROM __neighbour_keys k
    JOIN (SELECT unnest({{neighbour_ns}}::INTEGER[]) AS entity_n, generate_subscripts({{neighbour_ns}}::INTEGER[], 1) AS i) o
        USING (entity_n)), [])
