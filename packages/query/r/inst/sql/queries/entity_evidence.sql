-- query: entity_evidence
-- A page of an entity's evidence items (the selected resources'), by resource and
-- position: the total and the next page's offset; NULL when no selected resource has
-- the entity.
-- param: entity_key VARCHAR
-- param: limit INTEGER = 20
-- param: offset INTEGER = 0
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- step: found AS table
SELECT entity_n FROM entity
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: found_ids AS list
SELECT entity_n FROM __found

-- step: page AS table
-- the selected resources' evidence records
SELECT * FROM entity_evidence WHERE entity_n IN ({{found_ids*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))
ORDER BY resource, resource_version, entity_id, ordinal
LIMIT {{limit}} OFFSET {{offset}}

-- step: total AS value
SELECT count(*) FROM entity_evidence WHERE entity_n IN ({{found_ids*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))

-- step: response AS value when (SELECT count(*) FROM __found) > 0
WITH total AS (SELECT {{total}} AS t)
SELECT to_json(struct_pack(
    evidence := coalesce((SELECT list(struct_pack(source := source, dataset := dataset, row_id := row_id,
        upstream_id := upstream_id, annotations := annotations, molecular_form := molecular_form)
        ORDER BY resource, resource_version, entity_id, ordinal) FROM __page), []),
    evidenceTotal := (SELECT t FROM total),
    nextCursor := CASE WHEN (SELECT count(*) FROM __page) > 0
        AND {{offset}} + (SELECT count(*) FROM __page) < (SELECT t FROM total)
        THEN ({{offset}} + (SELECT count(*) FROM __page))::VARCHAR END))
