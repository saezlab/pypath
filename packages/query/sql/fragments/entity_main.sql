-- Steps: the entity {{entity_key}} (merged over its copies at build) when one of its
-- copies is in the resources {{resources}}, with its identifiers and annotations: the
-- table __main (the columns entity_summary reads), __main_summary, the value
-- entity_json (its summary with molecularEvidenceTotal; NULL without the entity) and
-- __copies (its copies in those resources, for their evidence).

-- step: merged AS table
SELECT * FROM entity_merged
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: copies AS table
SELECT _eid, entity_id, _resource, _resource || '/' || _resource_version AS resource FROM entity_by_key
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

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
