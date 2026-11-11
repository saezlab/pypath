-- query: entities_by_keys
-- Entities by names, identifiers or keys (resolved as resolve_keys does) as entity
-- summaries, in the resolved order; slim
-- leaves out identifiers and annotations.
-- param: tokens VARCHAR[] = NULL
-- param: slim BOOLEAN = false
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"tokens": ["TP53"]}

-- step: key_tokens AS list
SELECT t FROM (SELECT trim(unnest(coalesce({{tokens}}, []))) AS t, generate_subscripts(coalesce({{tokens}}, []), 1) AS i)
WHERE t <> '' ORDER BY i

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: merged AS table
-- the entities, in the order of the keys
SELECT e.*, o.position
FROM (SELECT * FROM entity WHERE entity_key IN ({{entity_keys*}})
      AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))) e
JOIN (SELECT unnest({{entity_keys}}::VARCHAR[]) AS entity_key,
        generate_subscripts({{entity_keys}}::VARCHAR[], 1) AS position) o USING (entity_key)

-- step: list_ids AS list when NOT {{slim}}
SELECT entity_n FROM __merged ORDER BY entity_n

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE entity_n IN ({{list_ids+}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE entity_n IN ({{list_ids*}})

-- step: summaries AS table
-- slim: the stored summary (no identifier or annotation lists); else built with them
WITH full_rows AS (SELECT * FROM __merged WHERE NOT {{slim}}),
listed AS (SELECT * FROM {{include:merged_rows full_rows}})
SELECT position, op_summary(summary) AS summary FROM __merged WHERE {{slim}}
UNION ALL
SELECT position, summary FROM {{include:entity_summary listed}}

-- step: response AS value
SELECT to_json(struct_pack(entities := coalesce((SELECT list(summary ORDER BY position) FROM __summaries), [])))
