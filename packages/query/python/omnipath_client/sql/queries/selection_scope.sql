-- query: selection_scope
-- An entity selection: the seed entities' keys (resolved from names, identifiers and
-- keys) and, unless include_associated is false, the entities at the other end of
-- their relations, in the order the relations list them.
-- param: entity_pks VARCHAR[] = NULL
-- param: include_associated BOOLEAN = true
-- param: annotation_term_ids VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_pks": ["TP53"]}

-- step: seeds AS list
SELECT t FROM (SELECT trim(unnest(coalesce({{entity_pks}}, []))) AS t, generate_subscripts(coalesce({{entity_pks}}, []), 1) AS i)
WHERE t <> '' ORDER BY i

-- step: key_tokens AS list
SELECT unnest({{seeds}})

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: selected_keys AS value
SELECT CASE WHEN len(coalesce({{entity_keys}}, [])) > 0 THEN {{entity_keys}} ELSE {{seeds}} END

-- step: neighbour_ids AS list when {{include_associated}} AND len(coalesce({{selected_keys}}, [])) > 0
SELECT DISTINCT _rid FROM relation_endpoint
WHERE key IN ({{selected_keys*}}) AND key_kind = 'entity'
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY _rid

-- step: neighbour_rows AS table
SELECT _rid, _resource || '/' || _resource_version AS resource, relation_id, subject_entity_key, object_entity_key
FROM relation WHERE _rid IN ({{neighbour_ids*}})

-- step: response AS value
WITH ends AS (
    SELECT k, min(struct_pack(r := resource, i := relation_id, s := side)) AS first
    FROM (SELECT resource, relation_id, 1 AS side, subject_entity_key AS k FROM __neighbour_rows
          UNION ALL SELECT resource, relation_id, 2, object_entity_key FROM __neighbour_rows)
    WHERE k IS NOT NULL AND NOT list_contains(coalesce({{selected_keys}}, []), k)
    GROUP BY k
), keys AS (
    SELECT coalesce({{selected_keys}}, []) || coalesce((SELECT list(k ORDER BY first) FROM ends), []) AS all_keys
)
SELECT to_json(struct_pack(
    entityPks := (SELECT all_keys FROM keys),
    seedEntityPks := coalesce({{seeds}}, []::VARCHAR[]),
    termEntityPks := []::VARCHAR[],
    ontologyTermIds := coalesce({{annotation_term_ids}}, []::VARCHAR[]),
    criteriaCount := len(coalesce({{seeds}}, [])),
    expandedEntityCount := greatest(0, len((SELECT all_keys FROM keys)) - len(coalesce({{seeds}}, [])))))
