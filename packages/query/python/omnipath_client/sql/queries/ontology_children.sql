-- query: ontology_children
-- A page of an ontology term's direct children within its namespace, by entity key.
-- param: term_id VARCHAR = ""
-- param: ontology_id VARCHAR = NULL
-- param: limit INTEGER = 100
-- param: offset INTEGER = 0
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"term_id": "GO:0008150", "limit": 5}

-- step: checked AS value
SELECT CASE WHEN {{limit}} BETWEEN 1 AND 1000 AND {{offset}} >= 0 THEN true
    ELSE error('invalid: Invalid ontology child page') END

-- steps: ontology_seed

-- steps: ontology_children

-- step: child_keys AS list
SELECT DISTINCT child_key FROM __child_rows ORDER BY 1

-- step: children AS table
SELECT entity_key, identifier, label, namespace, row_number() OVER (ORDER BY entity_key) AS position
FROM entity
WHERE entity_key IN ({{child_keys*}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))

-- step: response AS value
SELECT to_json(struct_pack(
    termId := {{term_id}},
    children := coalesce((SELECT list(struct_pack(id := identifier, termId := identifier,
            name := coalesce(nullif(label, ''), identifier), label := coalesce(nullif(label, ''), identifier),
            entityPk := entity_key, ontologyId := namespace) ORDER BY position)
        FROM __children WHERE position > {{offset}} AND position <= {{offset}} + {{limit}}), []),
    total := (SELECT count(*) FROM __children),
    nextCursor := CASE WHEN {{offset}} + {{limit}} < (SELECT count(*) FROM __children) THEN {{offset}} + {{limit}} END))
