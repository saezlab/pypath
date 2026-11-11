-- query: molecular_context
-- An entity's (or gene group's) molecular context: the relations of its gene reference
-- (or of the entity, or, in product view, of a protein or transcript form) with their
-- evidence, the molecular forms observed, the gene's catalogue products, the entity's
-- own evidence (standalone) and the other products the evidence names. NULL for an
-- unknown entity.
-- param: entity_id VARCHAR
-- param: view VARCHAR = "reference"
-- param: isoform_identifier VARCHAR = NULL
-- param: limit INTEGER = 20
-- param: offset INTEGER = 0
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_id": "missing"}

-- step: checked AS value
SELECT CASE
    WHEN {{view}} NOT IN ('reference', 'product') THEN error('invalid: Context view must be reference or product.')
    WHEN {{isoform_identifier}} IS NOT NULL AND {{view}} <> 'product' THEN error('invalid: Isoform selection requires product view.')
    ELSE true END

-- step: gene_group AS value
SELECT starts_with({{entity_id}}, 'gene:')

-- step: entity_key AS value
SELECT CASE WHEN {{gene_group}} THEN NULL ELSE {{entity_id}} END

-- steps: entity_main

-- step: entity_type AS value
SELECT lower(coalesce(summary.entityType, '')) FROM __main_summary

-- step: reference AS value
SELECT CASE WHEN {{gene_group}} THEN substr({{entity_id}}, 6) ELSE (SELECT summary.referenceEntityKey FROM __main_summary) END

-- step: found AS value
SELECT {{gene_group}} OR {{entity_json}} IS NOT NULL

-- step: product_kind AS value
SELECT CASE WHEN {{entity_type}} IN ('protein', 'transcript') THEN {{entity_type}} END

-- step: is_product AS value
SELECT CASE WHEN {{view}} = 'product' AND {{found}} AND {{product_kind}} IS NULL
    THEN error('invalid: Product view requires a protein or transcript record.') ELSE {{view}} = 'product' END

-- step: gene_reference AS value
SELECT coalesce(starts_with({{reference}}, 'entrez:'), false)

-- step: entities AS value
SELECT CASE WHEN NOT {{is_product}} AND NOT {{gene_reference}} THEN [{{entity_id}}] END

-- step: reference_entity_keys AS value
SELECT CASE WHEN NOT {{is_product}} AND {{gene_reference}} THEN [{{reference}}] END

-- step: protein_entity_keys AS value
SELECT CASE WHEN {{is_product}} AND {{product_kind}} = 'protein' THEN [{{entity_id}}] END

-- step: transcript_entity_keys AS value
SELECT CASE WHEN {{is_product}} AND {{product_kind}} = 'transcript' THEN [{{entity_id}}] END

-- step: isoform_identifiers AS value
SELECT CASE WHEN {{isoform_identifier}} IS NOT NULL THEN [{{isoform_identifier}}] END

-- step: gene_mode AS value
SELECT 'any'

-- step: scope_entity_ids AS value
SELECT NULL::VARCHAR[]

-- step: scope_endpoint_mode AS value
SELECT 'any'

-- step: molecular_endpoint_mode AS value
SELECT 'any'

-- step: sources AS value
SELECT NULL::VARCHAR[]

-- step: predicates AS value
SELECT NULL::VARCHAR[]

-- step: relation_categories AS value
SELECT NULL::VARCHAR[]

-- step: taxonomy_ids AS value
SELECT NULL::VARCHAR[]

-- step: interaction_types AS value
SELECT NULL::VARCHAR[]

-- step: entity_types AS value
SELECT NULL::VARCHAR[]

-- step: partner_types AS value
SELECT NULL::VARCHAR[]

-- step: is_directed AS value
SELECT NULL::BOOLEAN

-- step: sign AS value
SELECT NULL::INTEGER

-- step: signs AS value
SELECT NULL::INTEGER[]

-- step: min_evidence_count AS value
SELECT NULL::INTEGER

-- step: object_aspect_qualifier AS value
SELECT NULL::VARCHAR[]

-- step: object_direction_qualifier AS value
SELECT NULL::VARCHAR[]

-- step: causal_mechanism_qualifier AS value
SELECT NULL::VARCHAR[]

-- step: details AS value
SELECT true

-- steps: relation_page_steps

-- step: observed_relation_forms AS table
-- the forms at the selected end of each relation's evidence (in product view, the
-- forms of the selected product), first seen first
SELECT f.form, min(struct_pack(p := i.position, e := f.e, s := f.s)) AS first FROM (
    SELECT r.relation_n, ev.e, side.s,
        CASE side.s WHEN 1 THEN ev.item.subject_molecular_form ELSE ev.item.object_molecular_form END AS form,
        CASE side.s WHEN 1 THEN rel.subject_reference_entity_key ELSE rel.object_reference_entity_key END AS side_reference,
        CASE side.s WHEN 1 THEN rel.subject_n ELSE rel.object_n END AS side_n
    FROM __page_evidence r
    JOIN __endpoint_relations rel USING (relation_n),
    LATERAL (SELECT unnest(r.evidence) AS item, generate_subscripts(r.evidence, 1) AS e) ev,
    (VALUES (1), (2)) AS side(s)
) f JOIN (SELECT unnest({{page_ids}}::INTEGER[]) AS relation_n, generate_subscripts({{page_ids}}::INTEGER[], 1) AS position) i USING (relation_n)
WHERE f.form IS NOT NULL
    AND CASE WHEN {{is_product}} THEN op_form_ok(f.form, {{protein_entity_keys}}, {{transcript_entity_keys}}, {{isoform_identifiers}})
        WHEN {{gene_reference}} THEN f.side_reference = {{reference}}
        ELSE list_contains({{main_ids}}, f.side_n) END
GROUP BY f.form

-- step: product_ids AS list when {{reference}} IS NOT NULL
SELECT DISTINCT entity_n FROM entity_group WHERE group_key = {{reference}} AND kind = 'gene' ORDER BY 1

-- step: product_rows AS table
-- the gene's products (not its Entrez records) that a resource read has
SELECT * FROM entity WHERE entity_n IN ({{product_ids*}}) AND namespace <> 'entrez'
    AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: product_list_ids AS list
SELECT entity_n FROM __product_rows ORDER BY entity_n

-- step: identifier_rows AS table
SELECT * FROM entity_identifier WHERE entity_n IN ({{product_list_ids+}})

-- step: annotation_rows AS table
SELECT * FROM entity_annotation WHERE entity_n IN ({{product_list_ids*}})

-- step: products AS table
-- the products, most connected first
WITH placed AS (SELECT *, row_number() OVER (ORDER BY relation_count DESC, entity_key) AS position FROM __product_rows),
listed AS (SELECT * FROM {{include:merged_rows placed}})
SELECT * FROM {{include:entity_summary listed}}

-- step: standalone_ns AS list
-- the entities whose evidence is standalone: the gene's reference members, or the
-- entity (product view reads every evidence record of the form)
SELECT DISTINCT n FROM (
    SELECT entity_n AS n FROM entity_group
    WHERE NOT {{is_product}} AND {{gene_reference}} AND group_key = {{reference}} AND kind = 'reference'
    UNION ALL
    SELECT unnest({{main_ids}}::INTEGER[]) WHERE NOT {{is_product}} AND NOT {{gene_reference}})
ORDER BY 1

-- step: standalone_rows AS table
SELECT * FROM (
    SELECT * FROM entity_evidence WHERE entity_n IN ({{standalone_ns*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))
    UNION ALL BY NAME
    SELECT * FROM entity_evidence
    WHERE {{is_product}} AND (op_unset({{resources}}) OR list_contains({{resources}}, resource))
        AND CASE {{product_kind}} WHEN 'protein' THEN molecular_form.protein_entity_key = {{entity_id}}
            WHEN 'transcript' THEN molecular_form.transcript_entity_key = {{entity_id}} ELSE false END)
WHERE {{isoform_identifier}} IS NULL
    OR molecular_form.isoform_identifier.ns || ':' || molecular_form.isoform_identifier.id = {{isoform_identifier}}
ORDER BY resource, resource_version, entity_id, ordinal
LIMIT {{limit}} + 1 OFFSET {{offset}}

-- step: standalone_row_ns AS list
SELECT DISTINCT entity_n FROM __standalone_rows ORDER BY 1

-- step: standalone_keys AS table
SELECT entity_n, entity_key FROM entity WHERE entity_n IN ({{standalone_row_ns*}})

-- step: standalone AS table
SELECT row_number() OVER (ORDER BY s.resource, s.resource_version, s.entity_id, s.ordinal) AS position,
    k.entity_key, s.molecular_form,
    struct_pack(source := s.source, dataset := s.dataset, row_id := s.row_id, upstream_id := s.upstream_id,
        annotations := s.annotations, molecular_form := s.molecular_form) AS occurrence
FROM __standalone_rows s JOIN __standalone_keys k USING (entity_n)

-- step: known_keys AS value
SELECT coalesce(list(DISTINCT k), []::VARCHAR[]) FROM (
    SELECT item.subjectEntity.entityPk AS k FROM __items UNION ALL SELECT item.objectEntity.entityPk FROM __items
    UNION ALL SELECT summary.entityPk FROM __products)

-- step: referenced_keys AS list
-- at most 100 product keys the forms name and the page does not show: the standalone
-- forms' first, then each relation's evidence, subject then object
SELECT k FROM (
    SELECT k, min(o) AS o FROM (
        SELECT unnest([molecular_form.protein_entity_key, molecular_form.transcript_entity_key]) AS k,
            struct_pack(a := 0, p := position, e := 0, s := 0, f := generate_subscripts([1, 2], 1)) AS o
        FROM __standalone WHERE position <= {{limit}}
        UNION ALL
        SELECT unnest([form.protein_entity_key, form.transcript_entity_key]),
            struct_pack(a := 1, p := position, e := e, s := s, f := generate_subscripts([1, 2], 1))
        FROM (
            SELECT i.position, ev.e, side.s,
                CASE side.s WHEN 1 THEN ev.item.subject_molecular_form ELSE ev.item.object_molecular_form END AS form
            FROM __page_evidence r
            JOIN (SELECT unnest({{page_ids}}::INTEGER[]) AS relation_n, generate_subscripts({{page_ids}}::INTEGER[], 1) AS position) i USING (relation_n),
            LATERAL (SELECT unnest(r.evidence) AS item, generate_subscripts(r.evidence, 1) AS e) ev,
            (VALUES (1), (2)) AS side(s)))
    WHERE k IS NOT NULL AND NOT list_contains({{known_keys}}, k)
    GROUP BY k)
ORDER BY o LIMIT 100

-- step: referenced AS table
-- the referenced products' stored summaries, in the order of the keys
SELECT o.position, op_summary(m.summary) AS summary
FROM (SELECT entity_key, summary FROM entity
      WHERE entity_key IN ({{referenced_keys*}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))) m
JOIN (SELECT unnest({{referenced_keys}}::VARCHAR[]) AS entity_key,
        generate_subscripts({{referenced_keys}}::VARCHAR[], 1) AS position) o USING (entity_key)

-- step: response AS value when {{found}}
WITH forms AS (
    SELECT form, struct_pack(a := 0, p := first.p, e := first.e, s := first.s) AS o FROM __observed_relation_forms
    UNION ALL
    SELECT molecular_form, struct_pack(a := 1, p := position, e := 0, s := 0) FROM __standalone
    WHERE position <= {{limit}} AND molecular_form IS NOT NULL
), distinct_forms AS (
    SELECT form, min(o) AS o FROM forms GROUP BY form
)
SELECT to_json(struct_pack(
    referenceEntityKey := {{reference}},
    standaloneEvidence := coalesce((SELECT list(struct_pack(entityPk := entity_key, occurrence := occurrence) ORDER BY position)
        FROM __standalone WHERE position <= {{limit}}), []),
    catalogueProducts := coalesce((SELECT list(summary ORDER BY position) FROM __products), []),
    referencedProducts := coalesce((SELECT list(summary ORDER BY position) FROM __referenced), []),
    observedForms := coalesce((SELECT list(form ORDER BY o) FROM distinct_forms), []),
    relations := coalesce((SELECT list(op_json_patch(to_json(i.item), to_json(struct_pack(evidence := coalesce(e.evidence, []))))
            ORDER BY i.position)
        FROM __items i LEFT JOIN __page_evidence e ON e.relation_n = list_extract({{page_ids}}::INTEGER[], i.position)), []),
    relationsTotal := {{total}},
    nextCursor := CASE WHEN {{offset}} + {{#page_ids}} < {{total}} OR (SELECT count(*) FROM __standalone) > {{limit}}
        THEN ({{offset}} + {{limit}})::VARCHAR END,
    filters := CASE
        WHEN {{is_product}} THEN to_json(map_from_entries(list_filter([
            struct_pack(k := {{product_kind}} || '_entity_keys', v := to_json([{{entity_id}}])),
            struct_pack(k := 'isoform_identifiers', v := to_json({{isoform_identifiers}}))], lambda e: e.v::VARCHAR <> 'null')))
        WHEN {{gene_reference}} THEN to_json(struct_pack(reference_entity_keys := [{{reference}}]))
        ELSE to_json(map_from_entries(list_filter([
            struct_pack(k := 'entity_pks', v := to_json([{{entity_id}}])),
            struct_pack(k := 'isoform_identifiers', v := to_json({{isoform_identifiers}}))], lambda e: e.v::VARCHAR <> 'null')))
    END))
