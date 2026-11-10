-- Steps: the entity {{entity_key}} merged from its copies in the resources {{resources}}
-- with their identifiers and annotations: the table __main (merged_entity), and the
-- value entity_json, its summary with molecularEvidenceTotal (NULL without copies).

-- step: copies AS table
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: copy_ids AS list
SELECT _eid FROM __copies ORDER BY _eid

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE _eid IN ({{copy_ids*}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE _eid IN ({{copy_ids*}})

-- step: main AS table
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
), nested AS (
    SELECT c.*, 0 AS position, coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations
    FROM __copies c LEFT JOIN copy_identifiers i USING (_eid) LEFT JOIN copy_annotations a USING (_eid)
)
SELECT * FROM {{include:merged_entity nested}}

-- step: main_summary AS table
SELECT * FROM {{include:entity_summary __main}}

-- step: entity_json AS value when (SELECT count(*) FROM __copies) > 0
-- the summary with the evidence count, merged as JSON: DuckDB 2.0 drops the name of
-- struct_insert's field inside a macro
SELECT json_merge_patch(to_json(s.summary),
    to_json(struct_pack(molecularEvidenceTotal := (SELECT evidence_total FROM __main))))
FROM __main_summary s
