-- Merged entities ({{input}}: entity_merged rows with a position) with their identifier
-- and annotation lists, as the entity_summary fragment reads them, and their evidence
-- total.
SELECT 'main' AS item, position, entity_key AS pk, namespace AS raw_namespace,
    identifier AS raw_identifier, label AS raw_label, entity_type AS raw_type, taxon,
    has_hierarchy, parent_count, child_count, group_connectivity, relation_count,
    gene_reference_keys, reference_entity_key, identifiers, annotations,
    [entity_key] AS source_keys, evidence_total
FROM {{input}}
