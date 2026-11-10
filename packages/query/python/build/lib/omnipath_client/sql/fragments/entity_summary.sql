-- Entity rows as the API's entity summaries (shape/entity.py and display_name.py):
-- reads {{input}} with columns item, position, pk, raw_namespace, raw_identifier,
-- raw_label, raw_type, taxon, has_hierarchy, parent_count, child_count, relation_count,
-- gene_reference_keys, reference_entity_key, identifiers, annotations, source_keys;
-- gives item, position and summary.
WITH shape_1 AS (
    SELECT *,
        coalesce(pk, '') AS pk_s,
        trim(coalesce(nullif(raw_namespace, ''), 'unknown')) AS namespace_s,
        trim(coalesce(nullif(raw_identifier, ''), coalesce(pk, ''))) AS identifier_s,
        list_transform(coalesce(annotations, []), lambda a: struct_pack(
            term := trim(coalesce(a.term, '')), value := trim(coalesce(a.value, '')),
            quantity := a.quantity, source := trim(coalesce(a.source, '')),
            dataset := trim(coalesce(a.dataset, '')))) AS annotations_s
    FROM {{input}}
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
        taxonomyName := t.taxonomy_name
    ) AS summary
    FROM shape_5
    LEFT JOIN (SELECT taxon AS taxonomy_id, name AS taxonomy_name FROM op_taxonomy) t
        ON t.taxonomy_id = nullif(coalesce(shape_5.taxon, ''), '')
)
SELECT item, position, summary FROM shaped
