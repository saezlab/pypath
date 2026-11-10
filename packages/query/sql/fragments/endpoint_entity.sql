-- A relation's endpoint as the entity_summary fragment reads it: the endpoint's entity
-- row of the relation's own resource, or the relation's own columns when that resource
-- lacks it. Reads {{input}}: item, position, key, side_type, side_label, side_reference,
-- relation_taxon and the entity row's columns as e_entity_key, e_namespace,
-- e_identifier, e_label, e_entity_type, e_taxon, e_has_hierarchy, e_parent_count,
-- e_child_count, e_relation_count, e_gene_reference_keys, e_reference_entity_key.
SELECT item, position,
    coalesce(e_entity_key, key) AS pk,
    CASE WHEN e_entity_key IS NULL THEN coalesce(nullif(side_type, ''), 'unknown') ELSE e_namespace END AS raw_namespace,
    CASE WHEN e_entity_key IS NULL THEN coalesce(nullif(side_label, ''), key) ELSE e_identifier END AS raw_identifier,
    CASE WHEN e_entity_key IS NULL THEN side_label ELSE e_label END AS raw_label,
    CASE WHEN e_entity_key IS NULL THEN side_type ELSE e_entity_type END AS raw_type,
    CASE WHEN e_entity_key IS NULL THEN relation_taxon ELSE e_taxon END AS taxon,
    e_has_hierarchy AS has_hierarchy, e_parent_count AS parent_count, e_child_count AS child_count,
    e_relation_count AS relation_count,
    coalesce(e_gene_reference_keys, []) AS gene_reference_keys,
    CASE WHEN e_entity_key IS NULL THEN side_reference ELSE e_reference_entity_key END AS reference_entity_key,
    []::STRUCT(ns VARCHAR, id VARCHAR, is_canonical BOOLEAN, source VARCHAR)[] AS identifiers,
    NULL AS annotations,
    []::VARCHAR[] AS source_keys
FROM {{input}}
