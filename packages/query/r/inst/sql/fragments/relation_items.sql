-- Relations as the API lists them (_relation_items): the relation summary and both
-- endpoints' entity summaries. Reads {{input}}: one row per relation with its relation
-- columns, position and annotations (those op_relation_json reads); and __endpoint_rows,
-- the endpoints' entities (entity_n, entity_key, summary). Gives position and item.
WITH joined AS (
    SELECT r.*, s.entity_key AS subject_entity_key, o.entity_key AS object_entity_key,
        s.summary AS subject_summary, o.summary AS object_summary
    FROM {{input}} r
    JOIN __endpoint_rows s ON s.entity_n = r.subject_n
    JOIN __endpoint_rows o ON o.entity_n = r.object_n
)
SELECT position, struct_pack(
    relation := op_relation_json(joined, joined.annotations),
    subjectEntity := subject_summary,
    objectEntity := object_summary) AS item
FROM joined
