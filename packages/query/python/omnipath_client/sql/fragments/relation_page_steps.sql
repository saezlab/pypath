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
SELECT r.relation_n, r.page_rank FROM __range_selected r WHERE {{include:relation_conditions r}}

-- step: general_order AS table when NOT {{fast}} AND {{keyed}}
SELECT r.relation_n, r.page_rank FROM __candidate_rows r WHERE {{include:relation_conditions r}}

-- step: page_ids AS list
-- the page: of the ordered selection, or straight from the relations without endpoint
-- filters
SELECT relation_n FROM (
    SELECT relation_n, page_rank FROM __fast_order WHERE {{fast}}
    UNION ALL
    SELECT relation_n, page_rank FROM __general_order WHERE NOT {{fast}} AND {{keyed}}
    UNION ALL
    SELECT r.relation_n, r.page_rank FROM relation r
    WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_has_any(r.resources, {{scope}}))
        AND {{include:relation_conditions r}}
)
ORDER BY page_rank
LIMIT {{limit}} OFFSET {{offset}}

-- step: total AS value
SELECT CASE
    WHEN {{fast}} THEN (SELECT count(*) FROM __fast_order)
    WHEN {{keyed}} THEN (SELECT count(*) FROM __general_order)
    WHEN {{#page_ids}} = {{limit}} OR {{offset}} > 0 THEN (
        SELECT count(*) FROM relation r
        WHERE (op_unset({{scope}}) OR list_has_any(r.resources, {{scope}}))
            AND {{include:relation_conditions r}})
    ELSE {{offset}} + {{#page_ids}}
END

-- step: endpoint_relations AS table
SELECT * FROM relation WHERE relation_n IN ({{page_ids*}})

-- step: page_annotation_rows AS table when {{include_details}}
SELECT * FROM relation_annotation WHERE relation_n IN ({{page_ids+}})

-- step: page_evidence_rows AS table when {{include_details}}
SELECT * FROM relation_evidence WHERE relation_n IN ({{page_ids*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))

-- steps: endpoint_rows

-- step: items AS table
WITH annotations AS (
    SELECT relation_n, list(struct_pack(term := term, value := value, scope := scope) ORDER BY ordinal) AS annotations
    FROM __page_annotation_rows GROUP BY ALL
), evidence AS (
    -- under molecular form filters, the matching evidence records count
    SELECT relation_n, count(*) AS n FROM __page_evidence_rows ev
    WHERE op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
        {{isoform_identifiers}}, {{molecular_endpoint_mode}})
    GROUP BY ALL
), ordered AS (
    SELECT unnest({{page_ids}}::INTEGER[]) AS relation_n, generate_subscripts({{page_ids}}::INTEGER[], 1) AS position
), relations AS (
    SELECT r.* REPLACE (CASE WHEN {{forms}} THEN coalesce(e.n, 0) ELSE r.evidence_count END AS evidence_count),
        o.position,
        CASE WHEN {{include_details}} THEN coalesce(a.annotations, []) ELSE op_qualifier_annotations(r) END AS annotations
    FROM __endpoint_relations r JOIN ordered o USING (relation_n)
    LEFT JOIN annotations a USING (relation_n) LEFT JOIN evidence e USING (relation_n)
)
SELECT * FROM {{include:relation_items relations}}

-- step: page_evidence AS table
-- each page relation's evidence (the matching records, under molecular form filters)
SELECT relation_n, list(struct_pack(source := source, dataset := dataset, row_id := row_id, upstream_id := upstream_id,
        annotations := annotations, subject_molecular_form := subject_molecular_form,
        object_molecular_form := object_molecular_form) ORDER BY resource, resource_version, relation_id, ordinal) AS evidence
FROM __page_evidence_rows ev
WHERE NOT {{forms}} OR op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
GROUP BY ALL
