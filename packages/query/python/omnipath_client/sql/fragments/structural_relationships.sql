-- Steps: a page of the structural relations (has_input, has_output, enabled_by,
-- has_member, has_part) with an end among the entity keys {{rel_keys}} (a list), by
-- predicate and page order: {{rel_limit}} from {{rel_offset}}, stated by a resource of
-- {{resources}}. Gives the tables __rel_page (with total) and __rel_items (position,
-- item: the relation item with its annotations).

-- step: rel_ns AS list
SELECT entity_n FROM entity WHERE entity_key IN ({{rel_keys*}}) ORDER BY entity_n

-- step: rel_endpoints AS table
SELECT relation_n, predicate, page_rank FROM relation_endpoint
WHERE entity_n IN ({{rel_ns+}})
    AND predicate IN ('has_input', 'has_output', 'enabled_by', 'has_member', 'has_part')
    AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: rel_entries AS table
-- a relation between two of the entities once
SELECT DISTINCT * FROM __rel_endpoints

-- step: rel_page AS table
SELECT relation_n, count(*) OVER () AS total,
    row_number() OVER (ORDER BY predicate, page_rank) AS position
FROM __rel_entries
QUALIFY position > {{rel_offset}} AND position <= {{rel_offset}} + {{rel_limit}}

-- step: rel_ids AS list
SELECT relation_n FROM __rel_page ORDER BY relation_n

-- step: endpoint_relations AS table
SELECT * FROM relation WHERE relation_n IN ({{rel_ids*}})

-- step: rel_annotation_rows AS table
SELECT * FROM relation_annotation WHERE relation_n IN ({{rel_ids+}})

-- steps: endpoint_rows

-- step: rel_items AS table
WITH annotations AS (
    SELECT relation_n, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __rel_annotation_rows GROUP BY ALL
), relations AS (
    SELECT r.*, p.position, coalesce(a.annotations, []) AS annotations
    FROM __endpoint_relations r JOIN __rel_page p USING (relation_n) LEFT JOIN annotations a USING (relation_n)
), items AS (
    SELECT * FROM {{include:relation_items relations}}
)
SELECT i.position, op_json_patch(to_json(i.item), to_json(struct_pack(annotations := r.annotations))) AS item
FROM items i JOIN relations r USING (position)
