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
    SELECT r._rid, ev.e, side.s,
        CASE side.s WHEN 1 THEN ev.item.subject_molecular_form ELSE ev.item.object_molecular_form END AS form,
        CASE side.s WHEN 1 THEN rel.subject_reference_entity_key ELSE rel.object_reference_entity_key END AS side_reference,
        CASE side.s WHEN 1 THEN rel.subject_entity_key ELSE rel.object_entity_key END AS side_key
    FROM __page_evidence r
    JOIN __endpoint_relations rel USING (_rid),
    LATERAL (SELECT unnest(r.evidence) AS item, generate_subscripts(r.evidence, 1) AS e) ev,
    (VALUES (1), (2)) AS side(s)
) f JOIN (SELECT unnest({{page_ids}}::BIGINT[]) AS _rid, generate_subscripts({{page_ids}}::BIGINT[], 1) AS position) i USING (_rid)
WHERE f.form IS NOT NULL
    AND CASE WHEN {{is_product}} THEN op_form_ok(f.form, {{protein_entity_keys}}, {{transcript_entity_keys}}, {{isoform_identifiers}})
        WHEN {{gene_reference}} THEN f.side_reference = {{reference}}
        ELSE f.side_key = {{entity_id}} END
GROUP BY f.form

-- step: product_ids AS list when {{reference}} IS NOT NULL
SELECT DISTINCT _eid FROM entity_group
WHERE group_key = {{reference}} AND kind = 'gene' AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY 1

-- step: product_copies AS table
SELECT * EXCLUDE (_resource_version), _resource || '/' || _resource_version AS resource FROM entity
WHERE _eid IN ({{product_ids*}})

-- step: product_copy_ids AS list
SELECT _eid FROM __product_copies WHERE namespace <> 'entrez' ORDER BY 1

-- step: product_identifier_rows AS table
SELECT * FROM entity_identifier WHERE _eid IN ({{product_copy_ids*}})

-- step: product_annotation_rows AS table
SELECT * FROM entity_annotation WHERE _eid IN ({{product_copy_ids*}})

-- step: products AS table
WITH copy_identifiers AS (
    SELECT _eid, list(struct_pack(ns := ns, id := id, is_canonical := is_canonical, source := source)
        ORDER BY ordinal) AS identifiers
    FROM __product_identifier_rows GROUP BY ALL
), copy_annotations AS (
    SELECT _eid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset) ORDER BY ordinal) AS annotations
    FROM __product_annotation_rows GROUP BY ALL
), nested AS (
    -- a key's place: its first copy's, by resource and id
    SELECT c.*, min(struct_pack(r := c.resource, i := c.entity_id)) OVER (PARTITION BY c.entity_key) AS first_copy,
        coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations
    FROM __product_copies c LEFT JOIN copy_identifiers i USING (_eid) LEFT JOIN copy_annotations a USING (_eid)
    WHERE c.namespace <> 'entrez'
), placed AS (
    SELECT * EXCLUDE (first_copy), dense_rank() OVER (ORDER BY first_copy) AS position FROM nested
), merged AS (
    SELECT * FROM {{include:merged_entity placed}}
)
SELECT * FROM {{include:entity_summary merged}}

-- step: standalone_ids AS list
-- the entity rows whose evidence is standalone: the gene's reference members, or the
-- entity's copies (product view reads every evidence row of the form)
SELECT _eid FROM (
    SELECT DISTINCT _eid FROM entity_group
    WHERE NOT {{is_product}} AND {{gene_reference}} AND group_key = {{reference}} AND kind = 'reference'
        AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
    UNION
    SELECT _eid FROM __copies WHERE NOT {{is_product}} AND NOT {{gene_reference}}
) ORDER BY 1

-- step: standalone_rows AS table
SELECT * FROM (
    SELECT * FROM entity_evidence WHERE _eid IN ({{standalone_ids*}})
    UNION ALL BY NAME
    SELECT * FROM entity_evidence
    WHERE {{is_product}} AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
        AND CASE {{product_kind}} WHEN 'protein' THEN molecular_form.protein_entity_key = {{entity_id}}
            WHEN 'transcript' THEN molecular_form.transcript_entity_key = {{entity_id}} ELSE false END)
WHERE {{isoform_identifier}} IS NULL
    OR molecular_form.isoform_identifier.ns || ':' || molecular_form.isoform_identifier.id = {{isoform_identifier}}
ORDER BY _resource || '/' || _resource_version, entity_id, ordinal
LIMIT {{limit}} + 1 OFFSET {{offset}}

-- step: standalone_keys AS table
SELECT _eid, entity_key FROM entity WHERE _eid IN (SELECT DISTINCT _eid FROM __standalone_rows)

-- step: standalone AS table
SELECT row_number() OVER (ORDER BY s._resource || '/' || s._resource_version, s.entity_id, s.ordinal) AS position,
    k.entity_key, s.molecular_form,
    struct_pack(source := s.source, dataset := s.dataset, row_id := s.row_id, upstream_id := s.upstream_id,
        annotations := s.annotations, molecular_form := s.molecular_form) AS occurrence
FROM __standalone_rows s JOIN __standalone_keys k USING (_eid)

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
            JOIN (SELECT unnest({{page_ids}}::BIGINT[]) AS _rid, generate_subscripts({{page_ids}}::BIGINT[], 1) AS position) i USING (_rid),
            LATERAL (SELECT unnest(r.evidence) AS item, generate_subscripts(r.evidence, 1) AS e) ev,
            (VALUES (1), (2)) AS side(s)))
    WHERE k IS NOT NULL AND NOT list_contains({{known_keys}}, k)
    GROUP BY k)
ORDER BY o LIMIT 100

-- step: referenced_copies AS table
SELECT * EXCLUDE (key_rank, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key IN ({{referenced_keys*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: referenced AS table
WITH ordered AS (
    SELECT unnest({{referenced_keys}}::VARCHAR[]) AS entity_key, generate_subscripts({{referenced_keys}}::VARCHAR[], 1) AS position
), nested AS (
    SELECT c.*, o.position,
        NULL::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
        NULL::STRUCT(term VARCHAR, "value" VARCHAR, quantity STRUCT(has_numeric_value DOUBLE,
            has_unit VARCHAR, has_unit_prefix VARCHAR, has_binary_relation VARCHAR,
            source_field VARCHAR, comparator VARCHAR), "source" VARCHAR, dataset VARCHAR)[] AS annotations
    FROM __referenced_copies c JOIN ordered o USING (entity_key)
), merged AS (
    SELECT * FROM {{include:merged_entity nested}}
)
SELECT * FROM {{include:entity_summary merged}}

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
        FROM __items i LEFT JOIN __page_evidence e ON e._rid = list_extract({{page_ids}}::BIGINT[], i.position)), []),
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
