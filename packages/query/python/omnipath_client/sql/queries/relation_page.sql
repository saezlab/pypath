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
-- param: details BOOLEAN = false
-- param: records BOOLEAN = false
-- output: json response
-- output: table items
-- example: {"entities": ["TP53"], "limit": 20}

-- steps: relation_page_steps

-- step: response AS value
-- records: also the page's relation rows as stored (with their annotations and evidence
-- under details)
WITH annotations AS (
    SELECT relation_n, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __page_annotation_rows GROUP BY ALL
), ordered AS (
    SELECT unnest({{page_ids}}::INTEGER[]) AS relation_n, generate_subscripts({{page_ids}}::INTEGER[], 1) AS position
), stored AS (
    SELECT o.position AS position_, a.annotations AS annotations_, e.evidence AS evidence_,
        r.* EXCLUDE (relation_n, subject_n, object_n, page_rank, evidence_count),
        s.entity_key AS subject_entity_key, s.label AS subject_label,
        b.entity_key AS object_entity_key, b.label AS object_label,
        CASE WHEN {{forms}} THEN len(coalesce(e.evidence, [])) ELSE r.evidence_count END AS evidence_count
    FROM __endpoint_relations r JOIN ordered o USING (relation_n)
    JOIN __endpoint_rows s ON s.entity_n = r.subject_n JOIN __endpoint_rows b ON b.entity_n = r.object_n
    LEFT JOIN annotations a USING (relation_n) LEFT JOIN __page_evidence e USING (relation_n)
), records AS (
    -- each row's own columns (the helper fields dropped)
    SELECT position_ AS position, annotations_ AS annotations, evidence_ AS evidence,
        json_merge_patch(to_json(t), '{"position_": null, "annotations_": null, "evidence_": null}') AS row_
    FROM stored t
)
-- (json_merge_patch: no records field unless asked)
SELECT json_merge_patch(to_json(struct_pack(
    relations := coalesce((SELECT list(item.relation ORDER BY position) FROM __items), []),
    rows := coalesce((SELECT list(item ORDER BY position) FROM __items), []),
    total := {{total}},
    nextCursor := CASE WHEN {{offset}} + {{#page_ids}} < {{total}} THEN ({{offset}} + {{#page_ids}})::VARCHAR END)),
    to_json(struct_pack(records := CASE WHEN {{records}} THEN coalesce((SELECT list(CASE WHEN {{details}}
            THEN op_json_patch(row_, to_json(struct_pack(annotations := coalesce(annotations, []), evidence := coalesce(evidence, []))))
            ELSE row_ END ORDER BY position) FROM records), []) END)))
