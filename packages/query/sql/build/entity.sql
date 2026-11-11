-- Each entity's copies as one entity: the richest copy's columns (the merged_entity
-- fragment, ranking copies by their identifier and annotation counts), its relations,
-- evidence, identifier and annotation counts, the resources that have it, and the
-- summary that search results, relation endpoints and group cards show (stored without
-- its identifiers, which op_summary adds). chemical: some copy is of a chemical type;
-- structure_group and gene_group: the structure block and Entrez gene reference all its
-- copies agree on, which group it in the grouped explorer. A relation end without an
-- entity row is an entity as its relations describe it (endpoint_stubs).
WITH rows_ AS (
    SELECT * FROM work.entity WHERE {{partition entity_key}}
    UNION ALL BY NAME
    SELECT entity_key, _resource, _resource_version, -1 AS entity_id,
        coalesce(nullif(stub.entity_type, ''), 'unknown') AS namespace,
        coalesce(nullif(stub.label, ''), entity_key) AS identifier, stub.label AS label,
        stub.entity_type AS entity_type, stub.taxon AS taxon, stub.reference_entity_key AS reference_entity_key
    FROM work.endpoint_stubs WHERE {{partition entity_key}}
), copies AS (
    SELECT * EXCLUDE (_resource_version), _resource || '/' || _resource_version AS resource,
        0 AS position,
        NULL::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
        NULL::STRUCT(term VARCHAR, "value" VARCHAR, quantity STRUCT(has_numeric_value DOUBLE,
            has_unit VARCHAR, has_unit_prefix VARCHAR, has_binary_relation VARCHAR,
            source_field VARCHAR, comparator VARCHAR), "source" VARCHAR)[] AS annotations
    FROM rows_
), merged AS (
    SELECT m.* REPLACE (coalesce(c.relation_count, 0) AS relation_count), n.entity_n
    FROM {{include:merged_entity copies}} m
    JOIN work.entity_numbers n ON n.entity_key = m.pk
    LEFT JOIN work.entity_relation_counts c USING (entity_n)
), shaped AS (
    SELECT * FROM {{include:entity_summary merged}}
), facts AS (
    SELECT entity_key,
        list(DISTINCT _resource ORDER BY _resource) AS resources,
        bool_or(coalesce(list_contains((SELECT list(entity_type) FROM op_chemical_types), entity_type), false)) AS chemical,
        CASE WHEN min(hash(group_connectivity)) FILTER (WHERE group_connectivity IS NOT NULL)
                = max(hash(group_connectivity)) FILTER (WHERE group_connectivity IS NOT NULL)
            THEN any_value(group_connectivity) FILTER (WHERE group_connectivity IS NOT NULL) END AS structure_group,
        CASE WHEN bool_and(coalesce(starts_with(reference_entity_key, 'entrez:'), false))
            THEN any_value(reference_entity_key) END AS gene_reference,
        min(hash(reference_entity_key)) = max(hash(reference_entity_key)) AS one_reference
    FROM rows_ GROUP BY entity_key
)
SELECT m.entity_n, m.pk AS entity_key, m.raw_namespace AS namespace, m.raw_identifier AS identifier,
    m.raw_label AS label, m.raw_type AS entity_type, m.taxon, m.has_hierarchy, m.parent_count,
    m.child_count, m.group_connectivity, m.reference_entity_key, m.gene_reference_keys,
    m.relation_count, m.evidence_total AS evidence_count,
    coalesce(l.identifier_count, 0) AS identifier_count, coalesce(l.annotation_count, 0) AS annotation_count,
    f.resources, f.chemical, f.structure_group,
    CASE WHEN f.one_reference THEN f.gene_reference END AS gene_group,
    struct_update(s.summary, identifiers := NULL) AS summary
FROM merged m
JOIN shaped s ON s.summary.entityPk = m.pk
JOIN facts f ON f.entity_key = m.pk
LEFT JOIN work.entity_list_counts l USING (entity_n)
