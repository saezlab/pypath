-- Steps: a page of the structural relations (has_input, has_output, enabled_by,
-- has_member, has_part) with an endpoint among the entity keys {{rel_keys}} (a list),
-- by predicate, resource and relation id: {{rel_limit}} from {{rel_offset}}, in the
-- resources {{resources}}. Gives the tables __rel_page (with total) and __rel_items
-- (position, item: the relation item with its annotations).

-- step: rel_entries AS table
SELECT DISTINCT _rid, resource, relation_id, predicate FROM relation_endpoint
WHERE key IN ({{rel_keys*}}) AND key_kind = 'entity'
    AND predicate IN ('has_input', 'has_output', 'enabled_by', 'has_member', 'has_part')
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: rel_page AS table
SELECT _rid, count(*) OVER () AS total,
    row_number() OVER (ORDER BY predicate, resource, relation_id) AS position
FROM __rel_entries
QUALIFY position > {{rel_offset}} AND position <= {{rel_offset}} + {{rel_limit}}

-- step: rel_ids AS list
SELECT _rid FROM __rel_page ORDER BY _rid

-- step: endpoint_relations AS table
SELECT * FROM relation WHERE _rid IN ({{rel_ids*}})

-- step: rel_annotation_rows AS table
SELECT * FROM relation_annotation WHERE _rid IN ({{rel_ids*}})

-- steps: endpoint_rows

-- step: rel_items AS table
WITH annotations AS (
    SELECT _rid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __rel_annotation_rows GROUP BY ALL
), relations AS (
    SELECT r.*, p.position, coalesce(a.annotations, []) AS annotations
    FROM __endpoint_relations r JOIN __rel_page p USING (_rid) LEFT JOIN annotations a USING (_rid)
), items AS (
    SELECT * FROM {{include:relation_items relations}}
)
SELECT i.position, op_json_patch(to_json(i.item), to_json(struct_pack(annotations := r.annotations))) AS item
FROM items i JOIN relations r USING (position)
