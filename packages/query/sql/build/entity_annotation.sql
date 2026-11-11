-- Each entity's annotations over its copies, each once (by term, value, quantity and
-- source), numbered (ordinal) in the order each is first seen.
SELECT entity_n, row_number() OVER (PARTITION BY entity_n ORDER BY first)::INTEGER AS ordinal,
    term, "value", quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
    quantity_has_binary_relation, quantity_source_field, quantity_comparator, "source"
FROM (
    SELECT c.entity_n, a.term, a."value", a.quantity_has_numeric_value, a.quantity_has_unit,
        a.quantity_has_unit_prefix, a.quantity_has_binary_relation, a.quantity_source_field,
        a.quantity_comparator, a."source",
        min(struct_pack(r := c.resource, v := c.resource_version, e := c.entity_id, o := a.ordinal)) AS first
    FROM work.entity_annotation a JOIN work.entity_copies c USING (_eid)
    WHERE {{partition c.entity_key}}
    GROUP BY ALL)
