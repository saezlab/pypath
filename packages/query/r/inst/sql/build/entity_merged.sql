-- Every entity's copies merged into one row: the merged_entity fragment over all copies
-- (ranking them by their identifier and annotation counts), the summary that search
-- results, relation endpoints and group cards show, the merged identifier and annotation
-- lists an entity page shows, and the copies' resources, types, taxa and reference keys,
-- which filters test (an entity passes when some copy has the value). The lists come from
-- merged_identifiers and merged_annotations (staged in the build's scratch database).
WITH copies AS (
    SELECT * EXCLUDE (_resource_version), _resource || '/' || _resource_version AS resource,
        0 AS position,
        NULL::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
        NULL::STRUCT(term VARCHAR, "value" VARCHAR, quantity STRUCT(has_numeric_value DOUBLE,
            has_unit VARCHAR, has_unit_prefix VARCHAR, has_binary_relation VARCHAR,
            source_field VARCHAR, comparator VARCHAR), "source" VARCHAR, dataset VARCHAR)[] AS annotations
    FROM entity WHERE {{partition entity_key}}
), merged AS (
    SELECT * FROM {{include:merged_entity copies}}
), shaped AS (
    SELECT * FROM {{include:entity_summary merged}}
), facts AS (
    SELECT entity_key,
        list(DISTINCT _resource ORDER BY _resource) AS resources,
        list(DISTINCT entity_type ORDER BY entity_type) FILTER (WHERE entity_type IS NOT NULL) AS copy_types,
        list(DISTINCT taxon ORDER BY taxon) FILTER (WHERE taxon IS NOT NULL) AS copy_taxa,
        list(DISTINCT reference_entity_key ORDER BY reference_entity_key)
            FILTER (WHERE reference_entity_key IS NOT NULL) AS copy_reference_keys,
        count(*) AS copies
    FROM entity WHERE {{partition entity_key}} GROUP BY entity_key
)
SELECT m.pk AS entity_key, m.raw_namespace AS namespace, m.raw_identifier AS identifier,
    m.raw_label AS label, m.raw_type AS entity_type, m.taxon, m.has_hierarchy, m.parent_count,
    m.child_count, m.group_connectivity, m.relation_count, m.evidence_total,
    m.gene_reference_keys, m.reference_entity_key,
    f.resources, coalesce(f.copy_types, []) AS copy_types, coalesce(f.copy_taxa, []) AS copy_taxa,
    coalesce(f.copy_reference_keys, []) AS copy_reference_keys, f.copies,
    s.summary,
    coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations
FROM merged m
JOIN shaped s ON s.summary.entityPk = m.pk
JOIN facts f ON f.entity_key = m.pk
LEFT JOIN merged_identifiers i ON i.entity_key = m.pk
LEFT JOIN merged_annotations a ON a.entity_key = m.pk
