-- The explorer's entity details as a DuckDB macro over the installed database.
--
-- Mirrors omnipath_api.queries.entity_details.details, the API's entity details:
-- the entity merged from its resource copies, its first 50 structural
-- relationships with both endpoints, its relation count and up to 100 association
-- relations, shaped as the API's JSON (shape/entity.py, shape/relation.py,
-- shape/display_name.py). One row with one JSON column, ``details``; no row when no
-- installed resource has the entity. Needs relations.sql (lookups through op_*_at,
-- op_*_with). Taxon names are null, as in the client: the taxonomy reference is not
-- installed.

-- Values in first-seen order without repeats (positions are 1-based).
CREATE OR REPLACE MACRO op_ordered_distinct(xs) AS (
    list_transform(
        list_filter(list_transform(xs, lambda x, i: struct_pack(x := x, i := i)),
            lambda s: list_position(xs, s.x) = s.i),
        lambda s: s.x)
);

-- The last value of a relation-scope annotation term, as a dict built from the
-- annotations in order keeps it (interaction_label).
CREATE OR REPLACE MACRO op_last_value(anns, term_) AS (
    list_reverse(list_filter(anns, lambda a: a.scope = 'relation' AND a.term = term_))[1].value
);

-- One relation as the API's relation summary (shape_relation_summary).
CREATE OR REPLACE MACRO op_relation_json(r, anns) AS (
    struct_pack(
        relationPk := coalesce(r.relation_key, ''),
        subjectEntityPk := coalesce(r.subject_entity_key, ''),
        subjectReferenceEntityKey := r.subject_reference_entity_key,
        objectReferenceEntityKey := r.object_reference_entity_key,
        predicate := coalesce(r.predicate, ''),
        displayLabel := CASE
            WHEN r.predicate = 'interacts_with'
                AND op_last_value(anns, 'causal_mechanism_qualifier') = 'binding' THEN 'Binds'
            WHEN r.predicate = 'affects'
                AND op_last_value(anns, 'object_aspect_qualifier') = 'transport'
                AND op_last_value(anns, 'causal_mechanism_qualifier') = 'relocalization'
                AND coalesce(op_last_value(anns, 'object_direction_qualifier'), '') = ''
                THEN 'Transports'
        END,
        -- effect qualifiers, the first value of each
        qualifiers := to_json(map_from_entries(coalesce(list_transform(
            list_filter(list_transform(
                list_filter(anns, lambda a: a.scope = 'relation' AND coalesce(a.value, '') <> ''
                    AND a.term IN ('object_direction_qualifier', 'object_aspect_qualifier',
                        'causal_mechanism_qualifier')),
                lambda a, i: struct_pack(a := a, i := i)),
                lambda s: list_position(list_transform(
                    list_filter(anns, lambda a: a.scope = 'relation' AND coalesce(a.value, '') <> ''
                        AND a.term IN ('object_direction_qualifier', 'object_aspect_qualifier',
                            'causal_mechanism_qualifier')),
                    lambda a: a.term), s.a.term) = s.i),
            lambda s: struct_pack(k := s.a.term, v := s.a.value)), []))),
        objectEntityPk := coalesce(r.object_entity_key, ''),
        relationCategory := r.category,
        interactionClass := r.interaction_class,
        isDirected := coalesce(r.is_directed, false),
        sign := coalesce(r.sign, 0),
        participantTypes := op_ordered_distinct(list_filter(
            [r.subject_type, r.object_type], lambda t: coalesce(t, '') <> '')),
        evidenceCount := coalesce(r.evidence_count, 0),
        sources := list_filter(coalesce(r.sources, []), lambda x: coalesce(x, '') <> '')
    )
);

CREATE OR REPLACE MACRO op_entity(entity_key_) AS TABLE
WITH copies AS MATERIALIZED (
    -- the entity in every installed resource
    SELECT * EXCLUDE (_resource_version), _resource || '/' || _resource_version AS resource
    FROM op_entity_with([entity_key_])
), copy_pairs AS MATERIALIZED (
    SELECT list(struct_pack(r := _resource, v := entity_id)) AS pairs FROM copies
), copy_identifiers AS (
    SELECT _resource, entity_id, list(struct_pack(ns := ns, id := id,
        is_canonical := is_canonical, source := source) ORDER BY ordinal) AS identifiers
    FROM op_entity_identifier_at((SELECT pairs FROM copy_pairs)) GROUP BY ALL
), copy_annotations AS (
    SELECT _resource, entity_id, list(struct_pack(term := term, value := value,
        quantity := CASE WHEN quantity_has_numeric_value IS NULL AND quantity_has_unit IS NULL
            AND quantity_has_unit_prefix IS NULL AND quantity_has_binary_relation IS NULL
            AND quantity_source_field IS NULL AND quantity_comparator IS NULL THEN NULL
            ELSE struct_pack(
                has_numeric_value := quantity_has_numeric_value,
                has_unit := quantity_has_unit,
                has_unit_prefix := quantity_has_unit_prefix,
                has_binary_relation := quantity_has_binary_relation,
                source_field := quantity_source_field,
                comparator := quantity_comparator) END,
        source := source, dataset := dataset) ORDER BY ordinal) AS annotations
    FROM op_entity_annotation_at((SELECT pairs FROM copy_pairs)) GROUP BY ALL
), nested AS MATERIALIZED (
    SELECT c.*, coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations
    FROM copies c
    LEFT JOIN copy_identifiers i USING (_resource, entity_id)
    LEFT JOIN copy_annotations a USING (_resource, entity_id)
), merged AS MATERIALIZED (
    -- The richest copy's columns, relation counts summed, identifiers united by (ns, id)
    -- and annotations by value when there are several copies.
    SELECT
        any_value(entity_key) AS entity_key,
        arg_max(struct_pack(namespace := namespace, identifier := identifier, label := label,
            entity_type := entity_type, taxon := taxon, has_hierarchy := has_hierarchy,
            parent_count := parent_count, child_count := child_count), struct_pack(
            a := len(annotations),
            i := len(identifiers),
            g := regexp_full_match(coalesce(label, ''), '[A-Z][A-Z0-9-]{1,14}')::INTEGER,
            n := (coalesce(label, '') <> '' AND coalesce(label, '') <> coalesce(identifier, '')
                AND NOT regexp_full_match(coalesce(label, ''), '[0-9]+'))::INTEGER,
            f := lower(coalesce(label, '')),
            l := coalesce(label, ''),
            k := entity_key,
            r := resource,
            h := coalesce(has_hierarchy, false),
            p := coalesce(parent_count, 0),
            c := coalesce(child_count, 0))) AS best,
        sum(coalesce(relation_count, 0))::BIGINT AS relation_count,
        sum(coalesce(evidence_count, 0))::BIGINT AS evidence_total,
        list_sort(list_distinct(flatten(list(coalesce(gene_reference_keys, []))))) AS gene_reference_keys,
        list(DISTINCT reference_entity_key)
            FILTER (WHERE coalesce(reference_entity_key, '') <> '') AS reference_keys,
        flatten(list(identifiers ORDER BY resource, entity_id)) AS all_identifiers,
        flatten(list(annotations ORDER BY resource, entity_id)) AS all_annotations,
        count(*) AS copies
    FROM nested
    HAVING count(*) > 0
), main_entity AS (
    SELECT
        'main' AS item, 0 AS position,
        entity_key AS pk, best.namespace AS raw_namespace, best.identifier AS raw_identifier,
        best.label AS raw_label, best.entity_type AS raw_type, best.taxon AS taxon,
        best.has_hierarchy AS has_hierarchy, best.parent_count AS parent_count,
        best.child_count AS child_count, relation_count, gene_reference_keys,
        CASE
            WHEN list_contains(coalesce(reference_keys, []),
                    coalesce(best.namespace, 'None') || ':' || coalesce(best.identifier, 'None'))
                OR len(coalesce(reference_keys, [])) > 1
            THEN coalesce(best.namespace, 'None') || ':' || coalesce(best.identifier, 'None')
            WHEN len(coalesce(reference_keys, [])) = 1 THEN reference_keys[1]
        END AS reference_entity_key,
        list_transform(list_filter(
            list_transform(all_identifiers, lambda x, i: struct_pack(x := x, i := i)),
            lambda s: list_position(list_transform(all_identifiers,
                lambda x: struct_pack(n := coalesce(x.ns, chr(1)), d := coalesce(x.id, chr(1)))),
                struct_pack(n := coalesce(s.x.ns, chr(1)), d := coalesce(s.x.id, chr(1)))) = s.i),
            lambda s: s.x) AS identifiers,
        CASE WHEN copies > 1 THEN list_transform(list_filter(
                list_transform(all_annotations, lambda x, i: struct_pack(x := x, i := i)),
                lambda s: list_position(all_annotations, s.x) = s.i), lambda s: s.x)
            ELSE all_annotations END AS annotations,
        [entity_key] AS source_keys
    FROM merged
), endpoint_entries AS MATERIALIZED (
    -- relations with the entity at an endpoint: (resource, relation, predicate, category)
    SELECT DISTINCT _resource, _resource || '/' || _resource_version AS resource,
        relation_id, predicate, category
    FROM op_relation_endpoint_with([entity_key_]) WHERE key_kind = 'entity'
), counts AS (
    SELECT count(*) AS interactions, count(*) FILTER (WHERE category = 'association') AS associations
    FROM endpoint_entries
), structural_page AS MATERIALIZED (
    SELECT _resource, resource, relation_id, count(*) OVER () AS total,
        row_number() OVER (ORDER BY predicate, resource, relation_id) AS position
    FROM (SELECT DISTINCT _resource, resource, relation_id, predicate FROM endpoint_entries
          WHERE predicate IN ('has_input', 'has_output', 'enabled_by', 'has_member', 'has_part'))
    QUALIFY position <= 50
), association_page AS MATERIALIZED (
    SELECT _resource, resource, relation_id
    FROM (SELECT DISTINCT _resource, resource, relation_id, predicate FROM endpoint_entries
          WHERE category = 'association')
    ORDER BY predicate, resource, relation_id LIMIT 100
), page_relations AS MATERIALIZED (
    SELECT * FROM op_relation_at((SELECT list(struct_pack(r := _resource, v := relation_id)) FROM (
        SELECT _resource, relation_id FROM structural_page
        UNION SELECT _resource, relation_id FROM association_page)))
), structural AS MATERIALIZED (
    SELECT p.position, r.* FROM structural_page p JOIN page_relations r USING (_resource, relation_id)
), structural_annotations AS (
    SELECT _resource, relation_id, list(struct_pack(term := term, value := value,
        quantity := CASE WHEN quantity_has_numeric_value IS NULL AND quantity_has_unit IS NULL
            AND quantity_has_unit_prefix IS NULL AND quantity_has_binary_relation IS NULL
            AND quantity_source_field IS NULL AND quantity_comparator IS NULL THEN NULL
            ELSE struct_pack(
                has_numeric_value := quantity_has_numeric_value,
                has_unit := quantity_has_unit,
                has_unit_prefix := quantity_has_unit_prefix,
                has_binary_relation := quantity_has_binary_relation,
                source_field := quantity_source_field,
                comparator := quantity_comparator) END,
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM op_relation_annotation_at((SELECT list(struct_pack(r := _resource, v := relation_id)) FROM structural))
    GROUP BY ALL
), endpoint_rows AS MATERIALIZED (
    -- each structural relation's endpoint entities, read in the relation's own resource
    SELECT * FROM op_entity_with((SELECT list(k) FROM (
        SELECT subject_entity_key AS k FROM structural WHERE subject_entity_key IS NOT NULL
        UNION SELECT object_entity_key FROM structural WHERE object_entity_key IS NOT NULL)))
), endpoints AS (
    SELECT side AS item, s.position,
        CASE side WHEN 'subject' THEN s.subject_entity_key ELSE s.object_entity_key END AS key,
        CASE side WHEN 'subject' THEN s.subject_type ELSE s.object_type END AS side_type,
        CASE side WHEN 'subject' THEN s.subject_label ELSE s.object_label END AS side_label,
        CASE side WHEN 'subject' THEN s.subject_reference_entity_key
            ELSE s.object_reference_entity_key END AS side_reference,
        s._resource, s.taxon AS relation_taxon
    FROM structural s, (VALUES ('subject'), ('object')) AS sides(side)
), endpoint_entities AS (
    -- the endpoint's entity row, or the relation's own columns when its resource lacks it
    SELECT p.item, p.position,
        coalesce(e.entity_key, p.key) AS pk,
        CASE WHEN e.entity_key IS NULL THEN coalesce(nullif(p.side_type, ''), 'unknown') ELSE e.namespace END AS raw_namespace,
        CASE WHEN e.entity_key IS NULL THEN coalesce(nullif(p.side_label, ''), p.key) ELSE e.identifier END AS raw_identifier,
        CASE WHEN e.entity_key IS NULL THEN p.side_label ELSE e.label END AS raw_label,
        CASE WHEN e.entity_key IS NULL THEN p.side_type ELSE e.entity_type END AS raw_type,
        CASE WHEN e.entity_key IS NULL THEN p.relation_taxon ELSE e.taxon END AS taxon,
        e.has_hierarchy, e.parent_count, e.child_count, e.relation_count,
        coalesce(e.gene_reference_keys, []) AS gene_reference_keys,
        CASE WHEN e.entity_key IS NULL THEN p.side_reference ELSE e.reference_entity_key END AS reference_entity_key,
        []::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
        NULL AS annotations,
        []::VARCHAR[] AS source_keys
    FROM endpoints p
    LEFT JOIN endpoint_rows e ON e._resource = p._resource AND e.entity_key = p.key
), to_shape AS (
    SELECT * FROM main_entity
    UNION ALL BY NAME SELECT * FROM endpoint_entities
), shape_1 AS (
    SELECT *,
        coalesce(pk, '') AS pk_s,
        trim(coalesce(nullif(raw_namespace, ''), 'unknown')) AS namespace_s,
        trim(coalesce(nullif(raw_identifier, ''), coalesce(pk, ''))) AS identifier_s,
        list_transform(coalesce(annotations, []), lambda a: struct_pack(
            term := trim(coalesce(a.term, '')), value := trim(coalesce(a.value, '')),
            quantity := a.quantity, source := trim(coalesce(a.source, '')),
            dataset := trim(coalesce(a.dataset, '')))) AS annotations_s
    FROM to_shape
), shape_2 AS (
    SELECT *,
        lower(namespace_s) AS slug,
        trim(coalesce(nullif(raw_label, ''), identifier_s)) AS label_s
    FROM shape_1
), shape_3 AS (
    -- identifier candidates: the canonical one, the label, then the nested ones
    SELECT *, [struct_pack(type := slug, id := identifier_s, source := '', canon := true)]
        || CASE WHEN label_s <> '' AND label_s <> identifier_s THEN [struct_pack(
            type := CASE WHEN slug IN ('uniprot', 'entrez', 'ensg') THEN 'genesymbol' ELSE 'name' END,
            id := label_s, source := '', canon := false)] ELSE [] END
        || list_transform(coalesce(identifiers, []), lambda x: struct_pack(
            type := lower(trim(coalesce(x.ns, ''))), id := trim(coalesce(x.id, '')),
            source := coalesce(x.source, ''), canon := coalesce(x.is_canonical, false))) AS candidates
    FROM shape_2
), shape_4 AS (
    SELECT *, list_transform(list_filter(
            list_transform(candidates, lambda c, i: struct_pack(c := c, i := i)),
            lambda s: s.c.type <> '' AND s.c.id <> ''
                AND list_position(list_transform(candidates, lambda c: struct_pack(t := c.type, d := c.id)),
                    struct_pack(t := s.c.type, d := s.c.id)) = s.i),
        lambda s: s.c) AS kept
    FROM shape_3
), shape_5 AS (
    SELECT *,
        op_ordered_distinct(
            list_transform(list_filter(kept, lambda k: k.source <> ''), lambda k: k.source)
            || list_transform(list_filter(annotations_s, lambda a: a.source <> ''), lambda a: a.source)) AS sources,
        list_transform(annotations_s, lambda a: CASE WHEN a.quantity IS NULL
            THEN to_json(struct_pack(term := a.term, value := a.value, source := a.source, dataset := a.dataset))
            ELSE to_json(struct_pack(term := a.term, value := a.value, quantity := a.quantity,
                source := a.source, dataset := a.dataset)) END) AS annotations_json,
        lower(regexp_replace(coalesce(raw_type, ''), '[\s_]', '', 'g')) AS kind,
        flatten([list_transform(list_filter(kept, lambda k: k.type = ns_), lambda k: k.id)
            FOR ns_ IN ['chebi', 'chembl', 'hmdb', 'pubchem']]) AS chemical_names,
        coalesce(has_hierarchy, false) AS hierarchy,
        coalesce(parent_count, 0) AS parents,
        coalesce(child_count, 0) AS children_
    FROM shape_4
), shaped AS (
    SELECT item, position, struct_pack(
        entityPk := pk_s,
        referenceEntityKey := reference_entity_key,
        geneReferenceKeys := coalesce(gene_reference_keys, []),
        molecularEvidence := []::JSON[],
        resolutionStatus := 'resolved',
        sourceEntityPks := op_ordered_distinct(
            list_filter(coalesce(source_keys, []), lambda k: trim(coalesce(k, '')) <> '') || [pk_s]),
        entityFacetHints := struct_pack(chemicalClasses := []::JSON[],
            metabolicDomains := []::JSON[], structuralSpecificities := []::JSON[]),
        entityType := lower(trim(coalesce(nullif(raw_type, ''), 'unknown'))),
        primaryNamespace := slug,
        canonicalIdentifier := identifier_s,
        canonicalIdentifierType := slug,
        primaryIdentifier := identifier_s,
        label := label_s,
        taxon := taxon,
        taxonomyId := nullif(coalesce(taxon, ''), ''),
        hasHierarchy := hierarchy,
        parentCount := parents,
        childCount := children_,
        relationCount := coalesce(relation_count, 0),
        sources := sources,
        sourceCount := len(sources),
        resources := sources,
        entityAttributes := CASE WHEN len(annotations_json) > 0 THEN annotations_json END,
        ontologyHierarchy := CASE WHEN hierarchy OR parents + children_ > 0 THEN struct_pack(
            termId := identifier_s,
            ontologyPrefix := CASE
                WHEN slug IN ('cvterm', 'cv_term') OR slug = '' THEN CASE
                    WHEN contains(identifier_s, '-') THEN split_part(identifier_s, '-', 1)
                    WHEN contains(identifier_s, ':') THEN split_part(identifier_s, ':', 1)
                    WHEN starts_with(identifier_s, 'KW') THEN 'KW'
                    ELSE slug END
                ELSE slug END,
            label := label_s,
            definition := list_filter(annotations_s, lambda a: a.term IN
                ('description', 'biolink:description', 'IAO:0000115') AND a.value <> '')[1].value,
            ontologyId := NULL::VARCHAR,
            childCount := children_,
            parentCount := parents) END,
        identifiersTotal := len(kept),
        identifiers := list_transform(kept, lambda k: struct_pack(
            entityPk := pk_s, identifier := k.id, identifierType := k.type, type := k.type,
            id := k.id, curie := NULL::VARCHAR, source := k.source, isCanonical := k.canon)),
        annotations := annotations_json,
        displayName := CASE
            WHEN kind IN ('chemicalentity', 'chemical', 'smallmolecule', 'compound', 'metabolite', 'drug', 'lipid')
            THEN CASE
                WHEN label_s <> '' AND NOT regexp_full_match(label_s, '[0-9]+')
                    AND NOT starts_with(lower(label_s), 'inchi=') THEN label_s
                ELSE coalesce(
                    list_filter(chemical_names, lambda v: NOT regexp_full_match(v, '[0-9]+'))[1],
                    chemical_names[1], nullif(label_s, ''), nullif(identifier_s, ''), pk_s) END
            WHEN kind IN ('protein', 'gene') THEN coalesce(
                nullif(label_s, ''),
                list_filter(kept, lambda k: k.type IN ('genesymbol', 'gene_symbol', 'gene_name_primary',
                    'symbol', 'om:0200'))[1].id,
                list_filter(kept, lambda k: k.type IN ('uniprot', 'uniprotkb'))[1].id,
                list_filter(kept, lambda k: k.type IN ('name', 'om:0202', 'entry_name', 'entry name')
                    OR ends_with(k.type, ':name') OR contains(k.type, ' entry name'))[1].id,
                nullif(identifier_s, ''), pk_s)
            ELSE coalesce(
                nullif(label_s, ''),
                list_filter(kept, lambda k: k.type IN ('name', 'om:0202', 'entry_name', 'entry name')
                    OR ends_with(k.type, ':name') OR contains(k.type, ' entry name'))[1].id,
                nullif(identifier_s, ''), pk_s)
        END,
        taxonomyName := NULL::VARCHAR
    ) AS summary
    FROM shape_5
), relationship_items AS (
    SELECT s.position, struct_pack(
        relation := op_relation_json(s, coalesce(a.annotations, [])),
        subjectEntity := subj.summary,
        objectEntity := obj.summary,
        annotations := coalesce(a.annotations, [])) AS item
    FROM structural s
    LEFT JOIN structural_annotations a USING (_resource, relation_id)
    JOIN shaped subj ON subj.item = 'subject' AND subj.position = s.position
    JOIN shaped obj ON obj.item = 'object' AND obj.position = s.position
)
SELECT to_json(struct_pack(
    -- the main entity with its evidence count; merged as JSON, since DuckDB 2.0 drops
    -- the name of struct_insert's field inside a macro
    entity := json_merge_patch(to_json((SELECT summary FROM shaped WHERE item = 'main')),
        to_json(struct_pack(molecularEvidenceTotal := (SELECT evidence_total FROM merged)))),
    relationships := coalesce((SELECT list(item ORDER BY position) FROM relationship_items), []),
    relationshipsTotal := coalesce((SELECT any_value(total) FROM structural_page), 0),
    summary := struct_pack(interactionCount := (SELECT interactions FROM counts)),
    annotations := coalesce((SELECT list(struct_pack(relationPk := r.relation_key, predicate := r.predicate)
            ORDER BY r.relation_key, r._resource || '/' || r._resource_version)
        FROM association_page p JOIN page_relations r USING (_resource, relation_id)), [])
)) AS details
WHERE EXISTS (SELECT 1 FROM copies);
