-- Steps: the relations with the keys {{relation_keys}} (a list) that a resource of
-- {{resources}} states, by key (the relation table is sorted by key). Gives the table
-- __found: relation_key, relation_n.

-- step: found AS table
SELECT relation_key, relation_n FROM relation
WHERE relation_key IN ({{relation_keys+}}) AND (op_unset({{resources}}) OR list_has_any(resources, {{resources}}))
