-- query: entity_evidence
-- A page of an entity's evidence items, by resource and position: the total and the
-- next page's offset; NULL when no selected resource has the entity.
-- param: entity_key VARCHAR
-- param: limit INTEGER = 20
-- param: offset INTEGER = 0
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- step: copies AS table
SELECT _eid, _resource || '/' || _resource_version AS resource, entity_id, coalesce(evidence_count, 0) AS n
FROM entity_by_key
WHERE entity_key = {{entity_key}} AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: ranges AS table
SELECT *, coalesce(sum(n) OVER (ORDER BY resource, entity_id ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), 0) AS start
FROM __copies

-- step: page_ids AS list
-- the copies the page overlaps
SELECT _eid FROM __ranges WHERE start < {{offset}} + {{limit}} AND start + n > {{offset}} ORDER BY _eid

-- step: evidence_rows AS table
SELECT * FROM entity_evidence WHERE _eid IN ({{page_ids*}})

-- step: response AS value when (SELECT count(*) FROM __copies) > 0
WITH numbered AS (
    SELECT e.* EXCLUDE (_eid, _resource, _resource_version, entity_id, ordinal),
        r.resource, r.entity_id, e.ordinal,
        r.start + row_number() OVER (PARTITION BY e._eid ORDER BY e.ordinal) - 1 AS position
    FROM __evidence_rows e JOIN __ranges r USING (_eid)
), page AS (
    SELECT * FROM numbered WHERE position >= {{offset}} AND position < {{offset}} + {{limit}}
), total AS (SELECT sum(n) AS t FROM __copies)
SELECT to_json(struct_pack(
    evidence := coalesce((SELECT list(struct_pack(source := source, dataset := dataset, row_id := row_id,
        upstream_id := upstream_id, annotations := annotations, molecular_form := molecular_form)
        ORDER BY position) FROM page), []),
    evidenceTotal := (SELECT t FROM total),
    nextCursor := CASE WHEN (SELECT count(*) FROM page) > 0
        AND {{offset}} + (SELECT count(*) FROM page) < (SELECT t FROM total)
        THEN ({{offset}} + (SELECT count(*) FROM page))::VARCHAR END))
