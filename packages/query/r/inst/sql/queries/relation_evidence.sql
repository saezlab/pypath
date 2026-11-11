-- query: relation_evidence
-- The evidence of a relation, from the selected resources stating it, as the API lists it: each record with its molecular forms and its annotations by scope (a record
-- without its own has its relation's), and the relation-scope annotations. Molecular
-- form filters keep the matching records.
-- param: relation_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- param: protein_entity_keys VARCHAR[] = NULL
-- param: transcript_entity_keys VARCHAR[] = NULL
-- param: isoform_identifiers VARCHAR[] = NULL
-- param: molecular_endpoint_mode VARCHAR = "any"
-- output: json response
-- example: {"relation_key": "missing"}

-- step: relation_keys AS list
SELECT {{relation_key}}

-- steps: relation_by_keys

-- step: ids AS list
SELECT relation_n FROM __found ORDER BY relation_n

-- step: annotation_rows AS table
SELECT * FROM relation_annotation WHERE relation_n IN ({{ids+}})

-- step: evidence_rows AS table
SELECT * FROM relation_evidence WHERE relation_n IN ({{ids*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))

-- step: response AS value
WITH relation_annotations AS (
    SELECT relation_n, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY ALL
), records AS (
    SELECT e.* EXCLUDE (resource), e.resource || '/' || e.resource_version AS resource,
        row_number() OVER (PARTITION BY e.relation_n, e.resource, e.resource_version, e.relation_id ORDER BY e.ordinal) - 1 AS index_,
        coalesce(a.annotations, []) AS relation_anns
    FROM __evidence_rows e LEFT JOIN relation_annotations a USING (relation_n)
), kept AS (
    SELECT * FROM records
    WHERE op_evidence_ok(records, {{protein_entity_keys}}, {{transcript_entity_keys}},
        {{isoform_identifiers}}, {{molecular_endpoint_mode}})
        OR (op_unset({{protein_entity_keys}}) AND op_unset({{transcript_entity_keys}})
            AND op_unset({{isoform_identifiers}}))
), items AS (
    SELECT resource, relation_id, ordinal, struct_pack(
        relationPk := {{relation_key}},
        relationEvidencePk := {{relation_key}} || ':' || sha256(
            '[' || to_json(resource) || ', ' || relation_id || ', ' || coalesce(to_json(source), 'null')
            || ', ' || coalesce(to_json(dataset), 'null') || ', ' || coalesce(to_json(row_id), 'null')
            || ', ' || index_ || ']'),
        source := coalesce(nullif(source, ''), 'unknown'),
        dataset := dataset,
        rowId := row_id,
        upstreamId := upstream_id,
        subjectMolecularForm := subject_molecular_form,
        objectMolecularForm := object_molecular_form,
        subjectAttributes := op_scoped_annotations(coalesce(annotations, relation_anns), 'subject'),
        objectAttributes := op_scoped_annotations(coalesce(annotations, relation_anns), 'object'),
        recordAttributes := op_scoped_annotations(coalesce(annotations, relation_anns), 'relation'),
        evidence := op_scoped_annotations(coalesce(annotations, relation_anns), 'relation')) AS item
    FROM kept
)
SELECT to_json(struct_pack(
    evidence := coalesce((SELECT list(item ORDER BY resource, relation_id, ordinal) FROM items), []),
    annotations := coalesce((SELECT flatten(list(op_scoped_annotations(a.annotations, 'relation')))
        FROM relation_annotations a), [])))
