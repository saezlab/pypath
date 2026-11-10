-- Relations as the API lists them (_relation_items): the relation summary and both
-- endpoints' entity summaries. Reads {{input}}: one row per relation with its relation
-- columns, _resource, position and annotations (those op_relation_json reads); and
-- __endpoint_rows, the endpoints' merged entities (entity_key, summary). An endpoint
-- without an entity row is summarized from the relation's own columns. Gives position
-- and item.
WITH ends AS (
    SELECT side AS item, r.position,
        CASE side WHEN 'subject' THEN r.subject_entity_key ELSE r.object_entity_key END AS key,
        CASE side WHEN 'subject' THEN r.subject_type ELSE r.object_type END AS side_type,
        CASE side WHEN 'subject' THEN r.subject_label ELSE r.object_label END AS side_label,
        CASE side WHEN 'subject' THEN r.subject_reference_entity_key
            ELSE r.object_reference_entity_key END AS side_reference,
        r._resource, r.taxon AS relation_taxon
    FROM {{input}} r, (VALUES ('subject'), ('object')) AS sides(side)
), joined AS (
    SELECT p.*, e.entity_key AS e_entity_key, e.summary AS e_summary
    FROM ends p LEFT JOIN __endpoint_rows e ON e.entity_key = p.key
), unknown AS (
    SELECT item, position, key, side_type, side_label, side_reference, relation_taxon,
        NULL::VARCHAR AS e_entity_key, NULL::VARCHAR AS e_namespace, NULL::VARCHAR AS e_identifier,
        NULL::VARCHAR AS e_label, NULL::VARCHAR AS e_entity_type, NULL::VARCHAR AS e_taxon,
        NULL::BOOLEAN AS e_has_hierarchy, NULL::INTEGER AS e_parent_count, NULL::INTEGER AS e_child_count,
        NULL::BIGINT AS e_relation_count, NULL::VARCHAR[] AS e_gene_reference_keys,
        NULL::VARCHAR AS e_reference_entity_key
    FROM joined WHERE e_entity_key IS NULL
), unknown_rows AS (
    SELECT * FROM {{include:endpoint_entity unknown}}
), shaped AS (
    SELECT item, position, e_summary AS summary FROM joined WHERE e_entity_key IS NOT NULL
    UNION ALL
    SELECT item, position, summary FROM {{include:entity_summary unknown_rows}}
)
SELECT r.position, struct_pack(
    relation := op_relation_json(r, r.annotations),
    subjectEntity := s.summary,
    objectEntity := o.summary) AS item
FROM {{input}} r
JOIN shaped s ON s.item = 'subject' AND s.position = r.position
JOIN shaped o ON o.item = 'object' AND o.position = r.position
