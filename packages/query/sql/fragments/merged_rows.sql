-- Entities ({{input}}: entity rows with a position) with their identifier and annotation
-- lists (from the tables __identifier_rows and __annotation_rows of their identifiers and
-- annotations), as the entity_summary fragment reads them, and their evidence total.
WITH identifiers AS (
    SELECT entity_n, list(struct_pack(ns := ns, id := id, is_canonical := is_canonical, "source" := "source")
        ORDER BY ordinal) AS identifiers
    FROM __identifier_rows GROUP BY entity_n
), annotations AS (
    SELECT entity_n, list(struct_pack(term := term, "value" := "value",
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        "source" := "source") ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY entity_n
)
SELECT 'main' AS item, e.position, e.entity_n, e.entity_key AS pk, e.namespace AS raw_namespace,
    e.identifier AS raw_identifier, e.label AS raw_label, e.entity_type AS raw_type, e.taxon,
    e.has_hierarchy, e.parent_count, e.child_count, e.group_connectivity, e.relation_count,
    e.gene_reference_keys, e.reference_entity_key,
    coalesce(i.identifiers, []) AS identifiers, coalesce(a.annotations, []) AS annotations,
    [e.entity_key] AS source_keys, e.evidence_count AS evidence_total
FROM {{input}} e LEFT JOIN identifiers i USING (entity_n) LEFT JOIN annotations a USING (entity_n)
