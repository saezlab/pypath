-- Each entity's annotations over its copies, in resource and entity order; with several
-- copies an annotation they share appears once (first seen), a single copy's are kept as
-- they are; read by entity_merged. A stage of the build, in its scratch database; built in
-- slices of entity keys ({{partition ...}}), since list aggregates hold their groups in
-- memory.
SELECT entity_key, list(struct_pack(term := term, "value" := "value", quantity := quantity,
    "source" := "source", dataset := dataset) ORDER BY first) AS annotations
FROM (
    SELECT e.entity_key, a.term, a."value",
        op_quantity(a.quantity_has_numeric_value, a.quantity_has_unit, a.quantity_has_unit_prefix,
            a.quantity_has_binary_relation, a.quantity_source_field, a.quantity_comparator) AS quantity,
        a."source", a.dataset,
        min(struct_pack(r := e._resource || '/' || e._resource_version, e := e.entity_id, o := a.ordinal)) AS first
    FROM entity_annotation a JOIN entity e USING (_eid)
    JOIN (SELECT entity_key, count(*) AS copies FROM entity WHERE {{partition entity_key}} GROUP BY entity_key) c
        USING (entity_key)
    GROUP BY e.entity_key, a.term, a."value", quantity, a."source", a.dataset,
        CASE WHEN c.copies = 1 THEN a.ordinal END)
GROUP BY entity_key
