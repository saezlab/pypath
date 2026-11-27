-- query: relations
-- Relations as the explorer exports them: one row per relation with its annotations,
-- evidence and the product records its evidence names. Endpoint filters (entities,
-- scope entities, reference keys) read the endpoint table's ranges; the other filters
-- are conditions on the relations. At most max_relations relations, by key.
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
-- param: metabolite_classes VARCHAR[] = NULL
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

-- step: selection_sources AS value
SELECT {{sources}}

-- step: read_candidates AS value
SELECT true

-- steps: relation_candidates

-- step: fold_ids AS list when NOT op_unset({{resources}}) AND coalesce({{min_evidence_count}}, 0) > 1 AND {{keyed}}
-- (the relation_folds fragment's, of the candidates)
SELECT DISTINCT relation_n FROM __candidate_rows WHERE NOT list_has_all({{resources}}, resources) ORDER BY relation_n

-- step: selection AS table
-- the relations to return, by key: the filters read the candidates' filter columns, so
-- that only the returned relations' rows are read below (not every candidate's), and
-- the steps below read their tables by joining this one (a literal list of this many
-- keys is read much more slowly)
SELECT relation_n FROM (
    SELECT r.relation_n FROM __candidate_rows r
    WHERE {{keyed}} AND {{include:relation_conditions r}}
    UNION ALL
    SELECT r.relation_n FROM relation r
    WHERE NOT {{keyed}} AND (op_unset({{scope}})
        OR r.relation_n IN (SELECT relation_n FROM relation_resource WHERE resource IN ({{scope*}})))
        AND {{include:relation_conditions r}})
ORDER BY relation_n
LIMIT {{max_relations}}

-- step: selected AS table
SELECT r.* FROM relation r JOIN __selection USING (relation_n)

-- step: endpoints AS table
SELECT entity_n, entity_key, label FROM entity
WHERE entity_n IN (SELECT subject_n FROM __selected UNION SELECT object_n FROM __selected)

-- step: annotation_rows AS table
SELECT * FROM relation_annotation WHERE relation_n IN (SELECT relation_n FROM __selection)
    AND op_annotation_in_scope(source, {{resources}})

-- step: evidence AS table
-- (each relation's evidence in its order: sorted by the order key, which is much faster
-- than an ORDER BY inside the aggregate)
SELECT relation_n, list_transform(list_sort(list(struct_pack(
    key := struct_pack(resource := resource, resource_version := resource_version,
        relation_id := relation_id, ordinal := ordinal),
    value := struct_pack(source := source, dataset := dataset, row_id := row_id,
        upstream_id := upstream_id, annotations := annotations,
        subject_molecular_form := subject_molecular_form,
        object_molecular_form := object_molecular_form)))), x -> x.value) AS evidence
FROM relation_evidence ev
WHERE relation_n IN (SELECT relation_n FROM __selection)
    AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))
    AND (NOT {{forms}} OR op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
        {{isoform_identifiers}}, {{molecular_endpoint_mode}}))
GROUP BY ALL

-- step: product_keys AS list
SELECT DISTINCT k FROM (
    SELECT unnest([f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key]) AS k
    FROM __evidence, unnest(__evidence.evidence) AS occurrences(f))
WHERE k IS NOT NULL ORDER BY k

-- step: product_rows AS table
SELECT * EXCLUDE (entity_n, summary, chemical, structure_group, gene_group) FROM entity
WHERE entity_key IN ({{product_keys*}}) AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))

-- step: rows AS table
WITH annotations AS (
    SELECT relation_n, list_transform(list_sort(list(struct_pack(key := ordinal,
        value := struct_pack(term := term, value := value,
            quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
                quantity_has_binary_relation, quantity_source_field, quantity_comparator),
            source := source, dataset := dataset, scope := scope)))), x -> x.value) AS annotations
    FROM __annotation_rows GROUP BY ALL
), nested AS (
    SELECT s.* EXCLUDE (evidence_count) REPLACE (op_scoped_sources(s.sources, {{resources}}) AS sources),
        sb.entity_key AS subject_entity_key, sb.label AS subject_label,
        ob.entity_key AS object_entity_key, ob.label AS object_label,
        coalesce(a.annotations, []) AS annotations,
        coalesce(e.evidence, []) AS evidence,
        -- under a resource selection the evidence is the selected resources' records
        CASE WHEN {{forms}} OR NOT op_unset({{resources}}) THEN len(coalesce(e.evidence, [])) ELSE s.evidence_count END AS evidence_count
    FROM __selected s
    JOIN __endpoints sb ON sb.entity_n = s.subject_n JOIN __endpoints ob ON ob.entity_n = s.object_n
    LEFT JOIN annotations a USING (relation_n) LEFT JOIN __evidence e USING (relation_n)
), product_keys AS (
    SELECT DISTINCT n.relation_n, unnest([
        f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key
    ]) AS entity_key
    FROM nested n, unnest(n.evidence) AS occurrences(f)
), products AS (
    SELECT k.relation_n, list(pe ORDER BY pe.entity_key) AS referenced_product_records
    FROM product_keys k JOIN __product_rows pe USING (entity_key)
    GROUP BY ALL
)
SELECT nested.* EXCLUDE (relation_n, subject_n, object_n, page_rank),
    coalesce(products.referenced_product_records, []) AS referenced_product_records
FROM nested LEFT JOIN products USING (relation_n)
ORDER BY nested.relation_n
