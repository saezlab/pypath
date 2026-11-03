-- The explorer's relation export as DuckDB macros over the installed database.
--
-- Mirrors omnipath_api.queries: resolve_entity_keys, _relation_selection, export_slice.
-- Reads the views over every installed resource (relation, relation_evidence, ...)
-- only for scans the API's queries also make; rows by id or key come through the lookup
-- macros of the installation (op_<table>_at / op_<table>_with), which the clients
-- write with the database (omnipath_client/database.py, the R package's database.R).
-- Needs the table op_names (Biolink names, from contract.json).
--
-- Conditions on whether a filter is given test the macro's parameters directly, so
-- DuckDB folds them while planning and drops the branches that do not apply.

-- True when a list filter is not given: NULL or empty, as the API's "if values:".
CREATE OR REPLACE MACRO op_unset(x) AS (x IS NULL OR len(x) = 0);

-- One endpoint's molecular form matches every given identity constraint.
CREATE OR REPLACE MACRO op_form_ok(form, proteins, transcripts, isoforms) AS (
    (op_unset(proteins) OR list_contains(proteins, form.protein_entity_key))
    AND (op_unset(transcripts) OR list_contains(transcripts, form.transcript_entity_key))
    AND (op_unset(isoforms) OR list_contains(
        isoforms, form.isoform_identifier.ns || ':' || form.isoform_identifier.id))
);

-- An evidence row matches the molecular form filters at the endpoints the mode names.
CREATE OR REPLACE MACRO op_evidence_ok(ev, proteins, transcripts, isoforms, mode) AS (
    CASE mode
        WHEN 'source' THEN op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
        WHEN 'target' THEN op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
        WHEN 'both' THEN op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
            AND op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
        ELSE op_form_ok(ev.subject_molecular_form, proteins, transcripts, isoforms)
            OR op_form_ok(ev.object_molecular_form, proteins, transcripts, isoforms)
    END
);

-- Biolink names as stored: predicates and entity types by name, storage name or CURIE.
-- An unknown name is an error naming it.
CREATE OR REPLACE MACRO op_biolink(kind_, names) AS TABLE
SELECT CASE WHEN n.name IS NULL
    THEN error('Unknown Biolink ' || kind_ || ': ' || u.v) ELSE n.name END AS name
FROM unnest(names::VARCHAR[]) AS u(v)
LEFT JOIN op_names n ON n.kind = kind_ AND n.alias = lower(trim(u.v));

-- Entity keys of names, identifiers, ns:id / ns|id CURIEs, reference keys, entity keys
-- and group keys (gene:<reference>, connectivity:<block>, entity:<key>), by exact
-- case-insensitive match, in the resources of ``scope`` (NULL or empty: all).
CREATE OR REPLACE MACRO op_entity_keys(tokens, scope := NULL) AS TABLE
WITH classified AS MATERIALIZED (
    SELECT
        token,
        CASE
            WHEN starts_with(token, 'gene:') OR starts_with(token, 'connectivity:')
                OR starts_with(token, 'entity:') THEN 'group'
            WHEN regexp_full_match(token, '[0-9a-fA-F]{64}') THEN 'key'
            ELSE 'other'
        END AS kind,
        nullif(lower(regexp_extract(token, '^([A-Za-z][\w.-]*)[:|](.+)$', 1)), '') AS curie_ns,
        nullif(regexp_extract(token, '^([A-Za-z][\w.-]*)[:|](.+)$', 2), '') AS curie_id
    FROM (SELECT DISTINCT trim(t) AS token FROM unnest(tokens::VARCHAR[]) AS u(t) WHERE trim(t) <> '')
), curies AS MATERIALIZED (
    SELECT * FROM classified WHERE kind = 'other' AND curie_ns IS NOT NULL
), wanted_groups AS (
    SELECT CASE WHEN starts_with(token, 'gene:') THEN 'reference' ELSE 'connectivity' END AS group_kind,
        substr(token, strpos(token, ':') + 1) AS group_key
    FROM classified WHERE kind = 'group' AND NOT starts_with(token, 'entity:')
), term_rows AS MATERIALIZED (
    SELECT * FROM op_entity_term_with((SELECT list(t) FROM (
        SELECT lower(token) AS t FROM classified WHERE kind = 'other'
        UNION SELECT lower(curie_id) FROM curies)))
    WHERE op_unset(scope) OR list_contains(scope, _resource)
), group_rows AS MATERIALIZED (
    SELECT * FROM op_entity_group_with((SELECT list(k) FROM (
        SELECT curie_ns || ':' || curie_id AS k FROM curies
        UNION SELECT group_key FROM wanted_groups)))
    WHERE op_unset(scope) OR list_contains(scope, _resource)
), hits AS MATERIALIZED (
    -- (resource, entity_id, namespace the entity must have, NULL for any)
    SELECT _resource, entity_id, NULL::VARCHAR AS namespace
    FROM op_entity_with((SELECT list(token) FROM classified WHERE kind = 'other'))
    WHERE op_unset(scope) OR list_contains(scope, _resource)
    UNION ALL
    SELECT t._resource, t.entity_id, NULL FROM term_rows t
    WHERE t.term IN (SELECT lower(token) FROM classified WHERE kind = 'other')
    UNION ALL
    SELECT t._resource, t.entity_id, c.curie_ns FROM term_rows t
    JOIN curies c ON t.term = lower(c.curie_id)
    UNION ALL
    SELECT g._resource, g.entity_id, NULL FROM group_rows g
    JOIN curies c ON g.kind = 'reference' AND g.group_key = c.curie_ns || ':' || c.curie_id
    UNION ALL
    SELECT g._resource, g.entity_id, NULL FROM group_rows g
    JOIN wanted_groups w ON g.kind = w.group_kind AND g.group_key = w.group_key
), hit_entities AS (
    SELECT _resource, entity_id, entity_key, namespace
    FROM op_entity_at((SELECT list(struct_pack(r := _resource, v := entity_id)) FROM hits))
)
SELECT DISTINCT entity_key FROM (
    SELECT token AS entity_key FROM classified WHERE kind = 'key'
    UNION ALL
    SELECT substr(token, 8) FROM classified WHERE kind = 'group' AND starts_with(token, 'entity:')
    UNION ALL
    SELECT e.entity_key FROM hits h JOIN hit_entities e USING (_resource, entity_id)
    WHERE h.namespace IS NULL OR h.namespace = lower(coalesce(e.namespace, ''))
);

-- Relations as the explorer exports them: one row per relation with its annotations,
-- evidence and the product records its evidence names. Filters are the explorer's, as
-- lists (except the modes, is_directed, sign and min_evidence_count). At most
-- ``max_relations`` relations. A relation lists only its own resource among its
-- sources, so a source filter reads only those resources.
CREATE OR REPLACE MACRO op_relations(
    entities := NULL,
    gene_mode := 'any',
    reference_entity_keys := NULL,
    sources := NULL,
    predicates := NULL,
    relation_categories := NULL,
    taxonomy_ids := NULL,
    interaction_types := NULL,
    entity_types := NULL,
    partner_types := NULL,
    is_directed := NULL,
    sign := NULL,
    signs := NULL,
    min_evidence_count := NULL,
    object_aspect_qualifier := NULL,
    object_direction_qualifier := NULL,
    causal_mechanism_qualifier := NULL,
    protein_entity_keys := NULL,
    transcript_entity_keys := NULL,
    isoform_identifiers := NULL,
    molecular_endpoint_mode := 'any',
    max_relations := 100000
) AS TABLE
WITH settings AS (
    SELECT
        (SELECT list(name) FROM op_biolink('predicate', predicates)) AS predicate_names,
        (SELECT list(name) FROM op_biolink('class',
            CASE WHEN op_unset(interaction_types::VARCHAR[]) THEN entity_types ELSE interaction_types END
        )) AS type_names,
        (SELECT list(name) FROM op_biolink('class', partner_types)) AS partner_names,
        CASE WHEN op_unset(signs::INTEGER[]) AND sign IS NOT NULL THEN [sign::INTEGER]
            ELSE signs::INTEGER[] END AS sign_values
), keys AS MATERIALIZED (
    SELECT entity_key FROM op_entity_keys(entities, scope := sources::VARCHAR[])
), endpoints AS MATERIALIZED (
    SELECT * FROM op_relation_endpoint_with((SELECT list(k) FROM (
        SELECT entity_key AS k FROM keys
        UNION ALL SELECT unnest(reference_entity_keys::VARCHAR[]))))
    WHERE op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], _resource)
), by_entity AS (
    SELECT p._resource, p.relation_id FROM endpoints p
    WHERE p.key IN (SELECT entity_key FROM keys) AND p.key_kind = 'entity'
        AND (gene_mode NOT IN ('source', 'source_only') OR p.side = 'subject')
        AND (gene_mode NOT IN ('target', 'target_only') OR p.side = 'object')
    GROUP BY ALL
    HAVING gene_mode <> 'both' OR (bool_or(p.side = 'subject') AND bool_or(p.side = 'object'))
), by_reference AS (
    SELECT DISTINCT p._resource, p.relation_id FROM endpoints p
    WHERE p.key IN (SELECT unnest(reference_entity_keys::VARCHAR[])) AND p.key_kind = 'reference'
), candidates AS MATERIALIZED (
    SELECT * FROM by_entity
    WHERE op_unset(reference_entity_keys::VARCHAR[])
        OR (_resource, relation_id) IN (SELECT _resource, relation_id FROM by_reference)
    UNION ALL
    SELECT * FROM by_reference WHERE op_unset(entities::VARCHAR[])
), form_rows AS (
    -- relations with an evidence row matching the molecular form filters
    SELECT DISTINCT ev._resource, ev.relation_id
    FROM op_relation_evidence_at((SELECT list(struct_pack(r := _resource, v := relation_id)) FROM candidates)) ev
    WHERE NOT (op_unset(entities::VARCHAR[]) AND op_unset(reference_entity_keys::VARCHAR[]))
        AND op_evidence_ok(ev, protein_entity_keys::VARCHAR[], transcript_entity_keys::VARCHAR[],
            isoform_identifiers::VARCHAR[], molecular_endpoint_mode)
    UNION ALL
    SELECT DISTINCT ev._resource, ev.relation_id FROM relation_evidence ev
    WHERE op_unset(entities::VARCHAR[]) AND op_unset(reference_entity_keys::VARCHAR[])
        AND (op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], ev._resource))
        AND op_evidence_ok(ev, protein_entity_keys::VARCHAR[], transcript_entity_keys::VARCHAR[],
            isoform_identifiers::VARCHAR[], molecular_endpoint_mode)
), selected AS MATERIALIZED (
    SELECT r.* FROM (
        SELECT * FROM op_relation_at((SELECT list(struct_pack(r := _resource, v := relation_id)) FROM candidates))
        WHERE NOT (op_unset(entities::VARCHAR[]) AND op_unset(reference_entity_keys::VARCHAR[]))
        UNION ALL BY NAME
        -- native tables also hold the global ids of their lookups; not part of the answer
        SELECT COLUMNS(lambda c: c NOT IN ('_rid', '_eid')) FROM relation
        WHERE op_unset(entities::VARCHAR[]) AND op_unset(reference_entity_keys::VARCHAR[])
            AND (op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], _resource))
    ) r, settings s
    WHERE (op_unset(protein_entity_keys::VARCHAR[]) AND op_unset(transcript_entity_keys::VARCHAR[])
            AND op_unset(isoform_identifiers::VARCHAR[])
            OR (r._resource, r.relation_id) IN (SELECT _resource, relation_id FROM form_rows))
        AND (op_unset(relation_categories::VARCHAR[]) OR list_contains(relation_categories::VARCHAR[], r.category))
        AND (op_unset(s.predicate_names) OR list_contains(s.predicate_names, r.predicate))
        AND (op_unset(object_aspect_qualifier::VARCHAR[])
            OR list_has_any(r.object_aspect_qualifier, object_aspect_qualifier::VARCHAR[]))
        AND (op_unset(object_direction_qualifier::VARCHAR[])
            OR list_has_any(r.object_direction_qualifier, object_direction_qualifier::VARCHAR[]))
        AND (op_unset(causal_mechanism_qualifier::VARCHAR[])
            OR list_has_any(r.causal_mechanism_qualifier, causal_mechanism_qualifier::VARCHAR[]))
        AND (is_directed IS NULL OR r.is_directed = is_directed::BOOLEAN)
        AND (op_unset(s.sign_values) OR list_contains(s.sign_values, r.sign))
        AND (coalesce(min_evidence_count::BIGINT, 0) <= 1 OR r.evidence_count >= min_evidence_count::BIGINT)
        AND (op_unset(s.type_names)
            OR (op_unset(s.partner_names)
                AND (list_contains(s.type_names, r.subject_type) OR list_contains(s.type_names, r.object_type)))
            OR (NOT op_unset(s.partner_names)
                AND ((list_contains(s.type_names, r.subject_type) AND list_contains(s.partner_names, r.object_type))
                  OR (list_contains(s.type_names, r.object_type) AND list_contains(s.partner_names, r.subject_type)))))
        AND (op_unset(taxonomy_ids::VARCHAR[]) OR list_contains(taxonomy_ids::VARCHAR[], r.taxon))
        AND (op_unset(sources::VARCHAR[]) OR list_has_any(r.sources, sources::VARCHAR[]))
    LIMIT max_relations
), selected_pairs AS MATERIALIZED (
    SELECT list(struct_pack(r := _resource, v := relation_id)) AS pairs FROM selected
), annotations AS (
    -- Lists in ordinal order: collected unordered and sorted by their ordinals after,
    -- since list(... ORDER BY ...) over these structs is ten times slower.
    SELECT a._resource, a.relation_id, list_select(list(struct_pack(
        term := a.term, value := a.value,
        quantity := CASE WHEN a.quantity_has_numeric_value IS NULL AND a.quantity_has_unit IS NULL
            AND a.quantity_has_unit_prefix IS NULL AND a.quantity_has_binary_relation IS NULL
            AND a.quantity_source_field IS NULL AND a.quantity_comparator IS NULL THEN NULL
            ELSE struct_pack(
                has_numeric_value := a.quantity_has_numeric_value,
                has_unit := a.quantity_has_unit,
                has_unit_prefix := a.quantity_has_unit_prefix,
                has_binary_relation := a.quantity_has_binary_relation,
                source_field := a.quantity_source_field,
                comparator := a.quantity_comparator) END,
        source := a.source, dataset := a.dataset, scope := a.scope)), list_grade_up(list(a.ordinal))) AS annotations
    FROM op_relation_annotation_at((SELECT pairs FROM selected_pairs)) a
    GROUP BY ALL
), evidence AS (
    SELECT ev._resource, ev.relation_id, list_select(list(struct_pack(
        source := ev.source, dataset := ev.dataset, row_id := ev.row_id,
        upstream_id := ev.upstream_id, annotations := ev.annotations,
        subject_molecular_form := ev.subject_molecular_form,
        object_molecular_form := ev.object_molecular_form)), list_grade_up(list(ev.ordinal))) AS evidence
    FROM op_relation_evidence_at((SELECT pairs FROM selected_pairs)) ev
    WHERE op_unset(protein_entity_keys::VARCHAR[]) AND op_unset(transcript_entity_keys::VARCHAR[])
            AND op_unset(isoform_identifiers::VARCHAR[])
        OR op_evidence_ok(ev, protein_entity_keys::VARCHAR[], transcript_entity_keys::VARCHAR[],
            isoform_identifiers::VARCHAR[], molecular_endpoint_mode)
    GROUP BY ALL
), rows AS MATERIALIZED (
    SELECT sel.* EXCLUDE (evidence_count),
        coalesce(a.annotations, []) AS annotations,
        coalesce(e.evidence, []) AS evidence,
        CASE WHEN op_unset(protein_entity_keys::VARCHAR[]) AND op_unset(transcript_entity_keys::VARCHAR[])
                AND op_unset(isoform_identifiers::VARCHAR[])
            THEN sel.evidence_count ELSE len(coalesce(e.evidence, [])) END AS evidence_count
    FROM selected sel
    LEFT JOIN annotations a USING (_resource, relation_id)
    LEFT JOIN evidence e USING (_resource, relation_id)
), product_keys AS MATERIALIZED (
    SELECT DISTINCT rows._resource, rows.relation_id, unnest([
        f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key
    ]) AS entity_key
    FROM rows, unnest(rows.evidence) AS occurrences(f)
), product_entities AS (
    SELECT * EXCLUDE (_resource, _resource_version),
        _resource || '/' || _resource_version AS resource
    FROM op_entity_with((SELECT list(DISTINCT entity_key) FROM product_keys WHERE entity_key IS NOT NULL))
    WHERE op_unset(sources::VARCHAR[]) OR list_contains(sources::VARCHAR[], _resource)
), products AS (
    SELECT k._resource, k.relation_id,
        list_select(list(pe), list_grade_up(list(struct_pack(k := pe.entity_key, r := pe.resource))))
            AS referenced_product_records
    FROM product_keys k JOIN product_entities pe USING (entity_key)
    GROUP BY ALL
)
SELECT rows.* EXCLUDE (_resource, _resource_version, relation_id),
    coalesce(products.referenced_product_records, []) AS referenced_product_records
FROM rows LEFT JOIN products USING (_resource, relation_id);
