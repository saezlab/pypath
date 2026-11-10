-- Steps: a page of the relations matching the relation filters (the parameters of the
-- relations query, and limit and offset), in page order: the list page_ids, the value
-- total, the table __items (position, item: the relation with both endpoints) and
-- __page_evidence (each relation's evidence, when details are read). The includer's
-- value details (true) reads the page's annotations and evidence also without
-- molecular form filters.

-- step: selection_sources AS value
SELECT {{sources}}

-- steps: relation_fast

-- step: read_candidates AS value
SELECT NOT {{fast}}

-- steps: relation_candidates

-- steps: relation_range

-- step: include_details AS value
-- the page's annotations and evidence: under molecular form filters, or when asked
SELECT {{forms}} OR {{details}}

-- step: fast_order AS table when {{fast}}
SELECT r._rid, r.page_rank FROM __range_selected r WHERE {{include:relation_conditions r}}

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

-- step: page_annotation_rows AS table when {{include_details}}
SELECT * FROM relation_annotation WHERE _rid IN ({{page_ids*}})

-- step: page_evidence_rows AS table when {{include_details}}
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
        CASE WHEN {{include_details}} THEN coalesce(a.annotations, []) ELSE op_qualifier_annotations(r) END AS annotations
    FROM __endpoint_relations r JOIN ordered o USING (_rid)
    LEFT JOIN annotations a USING (_rid) LEFT JOIN evidence e USING (_rid)
)
SELECT * FROM {{include:relation_items relations}}

-- step: page_evidence AS table
-- each page relation's evidence (the matching rows, under molecular form filters)
SELECT _rid, list(struct_pack(source := source, dataset := dataset, row_id := row_id, upstream_id := upstream_id,
        annotations := annotations, subject_molecular_form := subject_molecular_form,
        object_molecular_form := object_molecular_form) ORDER BY ordinal) AS evidence
FROM __page_evidence_rows ev
WHERE NOT {{forms}} OR op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
GROUP BY ALL
