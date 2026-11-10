-- Relations as the API lists them (_relation_items): the relation summary and both
-- endpoints' entity summaries. Reads {{input}}: one row per relation with its relation
-- columns, _resource, position and annotations (those op_relation_json reads); and
-- __endpoint_rows, the endpoints' entity rows (entity_by_key) read by key. Gives
-- position and item.
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
    SELECT p.*, e.entity_key AS e_entity_key, e.namespace AS e_namespace, e.identifier AS e_identifier,
        e.label AS e_label, e.entity_type AS e_entity_type, e.taxon AS e_taxon,
        e.has_hierarchy AS e_has_hierarchy, e.parent_count AS e_parent_count,
        e.child_count AS e_child_count, e.relation_count AS e_relation_count,
        e.gene_reference_keys AS e_gene_reference_keys, e.reference_entity_key AS e_reference_entity_key
    FROM ends p LEFT JOIN __endpoint_rows e ON e._resource = p._resource AND e.entity_key = p.key
), endpoint_rows AS (
    SELECT * FROM {{include:endpoint_entity joined}}
), shaped AS (
    SELECT * FROM {{include:entity_summary endpoint_rows}}
)
SELECT r.position, struct_pack(
    relation := op_relation_json(r, r.annotations),
    subjectEntity := s.summary,
    objectEntity := o.summary) AS item
FROM {{input}} r
JOIN shaped s ON s.item = 'subject' AND s.position = r.position
JOIN shaped o ON o.item = 'object' AND o.position = r.position
