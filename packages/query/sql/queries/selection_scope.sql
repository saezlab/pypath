-- query: selection_scope
-- An entity selection: the seed entities' keys (resolved from names, identifiers and
-- keys) and, unless include_associated is false, the entities at the other end of
-- their relations, in the order the relations list them.
-- param: entity_pks VARCHAR[] = NULL
-- param: include_associated BOOLEAN = true
-- param: annotation_term_ids VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_pks": ["TP53"]}

-- step: scope_entity_pks AS value
SELECT {{entity_pks}}

-- step: scope_include_associated AS value
SELECT {{include_associated}}

-- steps: selection_scope

-- step: response AS value
SELECT to_json(struct_pack(
    entityPks := {{scope_keys}},
    seedEntityPks := coalesce({{seeds}}, []::VARCHAR[]),
    termEntityPks := []::VARCHAR[],
    ontologyTermIds := coalesce({{annotation_term_ids}}, []::VARCHAR[]),
    criteriaCount := len(coalesce({{seeds}}, [])),
    expandedEntityCount := greatest(0, len({{scope_keys}}) - len(coalesce({{seeds}}, [])))))
