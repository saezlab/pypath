-- query: relations
-- Relations as the explorer exports them: one row per relation with its annotations,
-- evidence and the product records its evidence names. Endpoint filters (entities,
-- scope entities, reference keys) read the endpoint table's ranges; the other filters
-- are conditions on the relations. At most max_relations relations, in global id order.
-- A relation lists only its own resource among its sources, so a source filter reads
-- only those resources.
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
-- param: max_relations INTEGER = 100000
-- param: resources VARCHAR[] = NULL
-- output: table rows
-- example: {"entities": ["TP53"]}

-- step: read_candidates AS value
SELECT true

-- steps: relation_candidates

-- step: selected AS table
SELECT r.* FROM (
    SELECT * FROM __candidate_rows
    UNION ALL BY NAME
    SELECT * FROM relation WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))
) r
WHERE {{include:relation_conditions r}}
ORDER BY r._rid
LIMIT {{max_relations}}

-- step: selected_ids AS list
SELECT _rid FROM __selected ORDER BY _rid

-- step: annotation_rows AS table
SELECT * FROM relation_annotation WHERE _rid IN ({{selected_ids*}})

-- step: evidence_rows AS table
SELECT * FROM relation_evidence WHERE _rid IN ({{selected_ids*}})

-- step: evidence AS table
SELECT _rid, list(struct_pack(
    source := source, dataset := dataset, row_id := row_id,
    upstream_id := upstream_id, annotations := annotations,
    subject_molecular_form := subject_molecular_form,
    object_molecular_form := object_molecular_form) ORDER BY ordinal) AS evidence
FROM __evidence_rows ev
WHERE NOT {{forms}} OR op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
GROUP BY ALL

-- step: product_keys AS list
SELECT DISTINCT k FROM (
    SELECT unnest([f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key]) AS k
    FROM __evidence, unnest(__evidence.evidence) AS occurrences(f))
WHERE k IS NOT NULL ORDER BY k

-- step: product_rows AS table
SELECT * EXCLUDE (key_rank, _eid, _resource, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key IN ({{product_keys*}}) AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))

-- step: rows AS table
WITH annotations AS (
    SELECT _rid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY ALL
), nested AS (
    SELECT s.* EXCLUDE (evidence_count),
        coalesce(a.annotations, []) AS annotations,
        coalesce(e.evidence, []) AS evidence,
        CASE WHEN {{forms}} THEN len(coalesce(e.evidence, [])) ELSE s.evidence_count END AS evidence_count
    FROM __selected s LEFT JOIN annotations a USING (_rid) LEFT JOIN __evidence e USING (_rid)
), product_keys AS (
    SELECT DISTINCT n._rid, unnest([
        f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key
    ]) AS entity_key
    FROM nested n, unnest(n.evidence) AS occurrences(f)
), products AS (
    SELECT k._rid, list(pe ORDER BY pe.entity_key, pe.resource) AS referenced_product_records
    FROM product_keys k JOIN __product_rows pe USING (entity_key)
    GROUP BY ALL
)
SELECT nested.* EXCLUDE (_rid, _resource, _resource_version, relation_id),
    coalesce(products.referenced_product_records, []) AS referenced_product_records
FROM nested LEFT JOIN products USING (_rid)
ORDER BY nested._rid
