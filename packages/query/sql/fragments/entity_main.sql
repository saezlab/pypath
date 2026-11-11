-- Steps: the entity {{entity_key}} when a resource of {{resources}} has it, with its
-- identifiers and annotations: the table __main (the columns entity_summary reads),
-- __main_summary, the value entity_json (its summary with molecularEvidenceTotal; NULL
-- without the entity) and the list main_ids (its number, for its evidence).

-- step: merged AS table
SELECT * FROM entity
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: main_ids AS list
SELECT entity_n FROM __merged

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE entity_n IN ({{main_ids+}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE entity_n IN ({{main_ids+}})

-- step: main AS table
WITH placed AS (SELECT *, 0 AS position FROM __merged)
SELECT * FROM {{include:merged_rows placed}}

-- step: main_summary AS table
SELECT * FROM {{include:entity_summary __main}}

-- step: entity_json AS value when (SELECT count(*) FROM __merged) > 0
-- the summary with the evidence count, merged as JSON: DuckDB 2.0 drops the name of
-- struct_insert's field inside a macro
SELECT op_json_patch(to_json(s.summary),
    to_json(struct_pack(molecularEvidenceTotal := (SELECT evidence_total FROM __main))))
FROM __main_summary s
