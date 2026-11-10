-- query: ontology_terms
-- Ontology term metadata by identifier or entity key: id, label, namespace and the
-- definition among its annotations; keyed by identifier and by label, NULL for ids
-- without a term.
-- param: term_ids VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"term_ids": ["GO:0006915"]}

-- step: rows AS table
SELECT _eid, _resource || '/' || _resource_version AS resource, entity_id, entity_key, label, identifier, namespace
FROM entity
WHERE (identifier IN ({{term_ids*}}) OR entity_key IN ({{term_ids*}}))
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: ids AS list
SELECT _eid FROM __rows ORDER BY _eid

-- step: annotation_rows AS table
SELECT _eid, term, value, ordinal FROM entity_annotation WHERE _eid IN ({{ids*}})

-- step: response AS value
WITH defined AS (
    SELECT r.*, (SELECT value FROM __annotation_rows a WHERE a._eid = r._eid AND coalesce(a.value, '') <> ''
            AND (a.term IN ('description', 'biolink:description', 'IAO:0000115') OR contains(a.term, '0801')
                OR contains(lower(a.term), 'def') OR contains(lower(a.term), 'comment'))
            ORDER BY a.ordinal LIMIT 1) AS definition,
        row_number() OVER (ORDER BY resource, entity_id) AS n
    FROM __rows r
), items AS (
    SELECT n, coalesce(identifier, '') AS ident, coalesce(nullif(label, ''), coalesce(identifier, '')) AS lbl,
        to_json(struct_pack(id := coalesce(nullif(identifier, ''), coalesce(nullif(label, ''), coalesce(identifier, ''))),
            termId := coalesce(nullif(identifier, ''), coalesce(nullif(label, ''), coalesce(identifier, ''))),
            label := coalesce(nullif(label, ''), coalesce(identifier, '')),
            name := coalesce(nullif(label, ''), coalesce(identifier, '')),
            namespace := namespace, definition := definition)) AS item
    FROM defined
), keyed AS (
    -- by identifier and by label; a later row wins
    SELECT k, arg_max(item, struct_pack(n := n, o := o)) AS item FROM (
        SELECT ident AS k, item, n, 1 AS o FROM items WHERE ident <> ''
        UNION ALL SELECT lbl, item, n, 2 FROM items WHERE lbl <> '')
    GROUP BY k
), wanted AS (SELECT DISTINCT unnest(coalesce({{term_ids}}, [])) AS k)
SELECT to_json(struct_pack(terms := coalesce(map_from_entries(list(struct_pack(k := k, v := item) ORDER BY k)), map([]::VARCHAR[], []::JSON[]))))
FROM (SELECT k, item FROM keyed UNION ALL SELECT k, NULL::JSON FROM wanted WHERE k NOT IN (SELECT k FROM keyed))
