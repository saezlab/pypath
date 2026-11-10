-- query: entities_by_keys
-- Entities by names, identifiers or keys (resolved as resolve_keys does), each merged
-- from its resource copies as an entity summary, in the resolved order; slim leaves
-- out identifiers and annotations.
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

-- step: copies AS table
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key IN ({{entity_keys*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: copy_ids AS list when NOT {{slim}}
SELECT _eid FROM __copies ORDER BY _eid

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE _eid IN ({{copy_ids*}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE _eid IN ({{copy_ids*}})

-- step: summaries AS table
WITH copy_identifiers AS (
    SELECT _eid, list(struct_pack(ns := ns, id := id, is_canonical := is_canonical, source := source)
        ORDER BY ordinal) AS identifiers
    FROM __identifier_rows GROUP BY ALL
), copy_annotations AS (
    SELECT _eid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset) ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY ALL
), ordered AS (
    SELECT unnest({{entity_keys}}::VARCHAR[]) AS entity_key, generate_subscripts({{entity_keys}}::VARCHAR[], 1) AS position
), nested AS (
    SELECT c.*, o.position,
        CASE WHEN {{slim}} THEN NULL ELSE coalesce(i.identifiers, []) END AS identifiers,
        CASE WHEN {{slim}} THEN NULL ELSE coalesce(a.annotations, []) END AS annotations
    FROM __copies c JOIN ordered o USING (entity_key)
    LEFT JOIN copy_identifiers i USING (_eid) LEFT JOIN copy_annotations a USING (_eid)
), merged AS (
    SELECT * FROM {{include:merged_entity nested}}
)
SELECT * FROM {{include:entity_summary merged}}

-- step: response AS value
SELECT to_json(struct_pack(entities := coalesce((SELECT list(summary ORDER BY position) FROM __summaries), [])))
