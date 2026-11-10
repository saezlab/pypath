-- query: stats_interaction_types
-- Relations per predicate and interaction class, most first.
-- param: resources VARCHAR[] = NULL
-- output: json response
-- output: table counts
-- example: {}

-- step: counts AS table
SELECT predicate, interaction_class, count(*) AS n FROM relation
WHERE op_unset({{resources}}) OR list_contains({{resources}}, _resource)
GROUP BY 1, 2 ORDER BY n DESC, predicate, interaction_class

-- step: response AS value
SELECT to_json(coalesce(list(struct_pack(interactionType := predicate, interactionClass := interaction_class,
    count := n) ORDER BY n DESC, predicate, interaction_class), []))
FROM __counts
