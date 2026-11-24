-- query: entities_by_keys
-- Entities by names, identifiers or keys (resolved as resolve_keys does) as entity
-- summaries, in the resolved order; slim
-- leaves out identifiers and annotations. matches adds, for each token given (in order,
-- repeats kept), the keys of the entities it resolved to.
-- param: tokens VARCHAR[] = NULL
-- param: slim BOOLEAN = false
-- param: matches BOOLEAN = false
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

-- step: token_matches AS table when {{matches}}
-- each token and the entities read that it found: by a name, identifier or group key
-- (as entity_keys keeps them), or as an entity key itself
WITH found AS (
    SELECT DISTINCT h.token, e.entity_key
    FROM __key_hits h JOIN __key_hit_entities e USING (entity_n)
    WHERE h.member OR h.namespace IS NULL OR h.namespace = lower(coalesce(e.namespace, ''))
    UNION
    SELECT c.token, m.entity_key FROM __key_classes c
    JOIN __merged m ON m.entity_key = CASE WHEN c.kind = 'key' THEN c.token ELSE substr(c.token, 8) END
    WHERE c.kind = 'key' OR (c.kind = 'group' AND starts_with(c.token, 'entity:'))
)
SELECT c.position, c.token,
    coalesce(list(f.entity_key ORDER BY f.entity_key) FILTER (WHERE f.entity_key IS NOT NULL), []) AS entity_keys
FROM __key_classes c LEFT JOIN found f USING (token)
GROUP BY c.position, c.token

-- step: response AS value
SELECT CASE WHEN {{matches}} THEN to_json(struct_pack(
        entities := coalesce((SELECT list(summary ORDER BY position) FROM __summaries), []),
        matches := coalesce((SELECT list(struct_pack(identifier := token, entityKeys := entity_keys) ORDER BY position)
                             FROM __token_matches), [])))
    ELSE to_json(struct_pack(entities := coalesce((SELECT list(summary ORDER BY position) FROM __summaries), [])))
END
