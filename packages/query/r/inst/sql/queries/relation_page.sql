-- query: relation_page
-- A page of the relations matching the explorer's filters, in page order (category,
-- predicate, labels, resource, id), with both endpoints, the total and the next page's
-- offset. A plain entity or reference key selection is one range of the endpoint
-- table, in page order already; other selections read the candidates' relations.
-- param: entities VARCHAR[] = NULL
-- param: gene_mode VARCHAR = "any"
-- param: scope_entity_ids VARCHAR[] = NULL
-- param: scope_endpoint_mode VARCHAR = "any"
-- param: reference_entity_keys VARCHAR[] = NULL
-- param: sources VARCHAR[] = NULL
-- param: predicates VARCHAR[] = NULL
-- param: relation_categories VARCHAR[] = NULL
-- param: taxonomy_ids VARCHAR[] = NULL
-- param: interaction_types VARCHAR[] = NULL
-- param: entity_types VARCHAR[] = NULL
-- param: partner_types VARCHAR[] = NULL
-- param: is_directed BOOLEAN = NULL
-- param: sign INTEGER = NULL
-- param: signs INTEGER[] = NULL
-- param: min_evidence_count INTEGER = NULL
-- param: object_aspect_qualifier VARCHAR[] = NULL
-- param: object_direction_qualifier VARCHAR[] = NULL
-- param: causal_mechanism_qualifier VARCHAR[] = NULL
-- param: protein_entity_keys VARCHAR[] = NULL
-- param: transcript_entity_keys VARCHAR[] = NULL
-- param: isoform_identifiers VARCHAR[] = NULL
-- param: molecular_endpoint_mode VARCHAR = "any"
-- param: resources VARCHAR[] = NULL
-- param: limit INTEGER = 50
-- param: offset INTEGER = 0
-- output: json response
-- output: table items
-- example: {"entities": ["TP53"], "limit": 20}

-- step: fast AS value
-- entity keys or reference keys alone, without a scope or molecular forms: one range of
-- the endpoint table, in page order
SELECT (NOT op_unset({{entities}})) <> (NOT op_unset({{reference_entity_keys}}))
    AND op_unset({{scope_entity_ids}}) AND op_unset({{protein_entity_keys}})
    AND op_unset({{transcript_entity_keys}}) AND op_unset({{isoform_identifiers}})

-- step: read_candidates AS value
SELECT NOT {{fast}}

-- steps: relation_candidates

-- step: range_keys AS list when {{fast}}
SELECT unnest(CASE WHEN op_unset({{reference_entity_keys}}) THEN coalesce({{entity_key_list}}, [])
    ELSE {{reference_entity_keys}} END) AS k ORDER BY k

-- step: range_rows AS table when {{fast}}
SELECT * FROM relation_endpoint
WHERE key IN ({{range_keys*}})
    AND key_kind = CASE WHEN op_unset({{reference_entity_keys}}) THEN 'entity' ELSE 'reference' END

-- step: fast_order AS table when {{fast}}
-- a relation with both ends among the keys is one relation; a reference key selects
-- relations at either end
WITH mode AS (SELECT CASE WHEN op_unset({{reference_entity_keys}}) THEN {{gene_mode}} ELSE 'any' END AS m),
sides AS (
    SELECT _rid FROM __range_rows, mode GROUP BY _rid
    HAVING any_value(m) <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
)
SELECT DISTINCT ON (r._rid) r._rid, r.page_rank FROM __range_rows r, mode
WHERE (op_unset({{scope}}) OR list_contains({{scope}}, r._resource))
    AND (m NOT IN ('source', 'source_only') OR r.side = 'subject')
    AND (m NOT IN ('target', 'target_only') OR r.side = 'object')
    AND r._rid IN (SELECT _rid FROM sides)
    AND {{include:relation_conditions r}}

-- step: general_order AS table when NOT {{fast}} AND {{keyed}}
SELECT r._rid, row_number() OVER (ORDER BY r.category, r.predicate, r.subject_label, r.object_label,
    r._resource || '/' || r._resource_version, r.relation_id) AS page_rank
FROM __candidate_rows r
WHERE {{include:relation_conditions r}}

-- step: page_ids AS list
-- the page: of the ordered selection, or straight from a scan without endpoint filters
SELECT _rid FROM (
    SELECT _rid, page_rank AS o1, NULL::VARCHAR AS o2, NULL::VARCHAR AS o3, NULL::VARCHAR AS o4,
        NULL::VARCHAR AS o5, NULL::VARCHAR AS o6, NULL::INTEGER AS o7
    FROM __fast_order WHERE {{fast}}
    UNION ALL
    SELECT _rid, page_rank, NULL, NULL, NULL, NULL, NULL, NULL FROM __general_order WHERE NOT {{fast}} AND {{keyed}}
    UNION ALL
    SELECT r._rid, NULL, r.category, r.predicate, r.subject_label, r.object_label,
        r._resource || '/' || r._resource_version, r.relation_id
    FROM relation r
    WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, r._resource))
        AND {{include:relation_conditions r}}
)
ORDER BY o1, o2, o3, o4, o5, o6, o7
LIMIT {{limit}} OFFSET {{offset}}

-- step: total AS value
SELECT CASE
    WHEN {{fast}} THEN (SELECT count(*) FROM __fast_order)
    WHEN {{keyed}} THEN (SELECT count(*) FROM __general_order)
    WHEN {{#page_ids}} = {{limit}} OR {{offset}} > 0 THEN (
        SELECT count(*) FROM relation r
        WHERE (op_unset({{scope}}) OR list_contains({{scope}}, r._resource))
            AND {{include:relation_conditions r}})
    ELSE {{offset}} + {{#page_ids}}
END

-- step: endpoint_relations AS table
SELECT * FROM relation WHERE _rid IN ({{page_ids*}})

-- step: page_annotation_rows AS table when {{forms}}
SELECT * FROM relation_annotation WHERE _rid IN ({{page_ids*}})

-- step: page_evidence_rows AS table when {{forms}}
SELECT * FROM relation_evidence WHERE _rid IN ({{page_ids*}})

-- steps: endpoint_rows

-- step: items AS table
WITH annotations AS (
    SELECT _rid, list(struct_pack(term := term, value := value, scope := scope) ORDER BY ordinal) AS annotations
    FROM __page_annotation_rows GROUP BY ALL
), evidence AS (
    -- under molecular form filters, the matching evidence rows count
    SELECT _rid, count(*) AS n FROM __page_evidence_rows ev
    WHERE op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
        {{isoform_identifiers}}, {{molecular_endpoint_mode}})
    GROUP BY ALL
), ordered AS (
    SELECT unnest({{page_ids}}::BIGINT[]) AS _rid, generate_subscripts({{page_ids}}::BIGINT[], 1) AS position
), relations AS (
    SELECT r.* REPLACE (CASE WHEN {{forms}} THEN coalesce(e.n, 0) ELSE r.evidence_count END AS evidence_count),
        o.position,
        CASE WHEN {{forms}} THEN coalesce(a.annotations, []) ELSE op_qualifier_annotations(r) END AS annotations
    FROM __endpoint_relations r JOIN ordered o USING (_rid)
    LEFT JOIN annotations a USING (_rid) LEFT JOIN evidence e USING (_rid)
)
SELECT * FROM {{include:relation_items relations}}

-- step: response AS value
SELECT to_json(struct_pack(
    relations := coalesce((SELECT list(item.relation ORDER BY position) FROM __items), []),
    rows := coalesce((SELECT list(item ORDER BY position) FROM __items), []),
    total := {{total}},
    nextCursor := CASE WHEN {{offset}} + {{#page_ids}} < {{total}} THEN ({{offset}} + {{#page_ids}})::VARCHAR END))
