-- Each relation's annotations over its copies, each once (by every column but the
-- position), numbered (ordinal) in the order each is first seen (copies by global id,
-- each copy's in its order). In slices of relation keys: the full release's annotations
-- do not sort in 5 GB at once.
SELECT relation_n, row_number() OVER (PARTITION BY relation_n ORDER BY first)::INTEGER AS ordinal,
    term, "value", quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
    quantity_has_binary_relation, quantity_source_field, quantity_comparator, "source", dataset, "scope"
FROM (
    SELECT c.relation_n, a.term, a."value", a.quantity_has_numeric_value, a.quantity_has_unit,
        a.quantity_has_unit_prefix, a.quantity_has_binary_relation, a.quantity_source_field,
        a.quantity_comparator, a."source", a.dataset, a."scope",
        min((a._rid::HUGEINT << 32) + a.ordinal) AS first
    FROM work.relation_annotation a JOIN work.relation_copies c USING (_rid)
    WHERE {{partition c.relation_key}}
    GROUP BY ALL)
