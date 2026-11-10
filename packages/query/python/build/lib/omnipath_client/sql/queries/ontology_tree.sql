-- query: ontology_tree
-- An ontology term's ancestry as a tree: every path from its roots down to the term
-- (at most 10 levels up, within its namespace), with the term's first 100 children.
-- param: term_id VARCHAR = ""
-- param: ontology_id VARCHAR = NULL
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"term_id": "GO:0006915"}

-- steps: ontology_seed

-- steps: ontology_children

-- step: has_more AS value
SELECT (SELECT count(*) FROM __child_rows) > 100

-- step: seed_children AS table
SELECT * FROM __child_rows ORDER BY child_key LIMIT 100

-- step: frontier AS value
SELECT CASE WHEN list_contains({{known}}, {{seed_key}}) THEN [{{seed_key}}] ELSE []::VARCHAR[] END

-- step: seen AS value
SELECT coalesce([{{seed_key}}], []::VARCHAR[])

-- step: tree_edges AS table
SELECT NULL::VARCHAR AS parent_key, NULL::VARCHAR AS parent_label, NULL::VARCHAR AS child_key,
    NULL::VARCHAR AS child_label WHERE false

-- repeat: level = 1 next {{level}} + 1 until walk_done max 10

-- step: found AS table
SELECT DISTINCT child_key, child_label, predicate, parent_key, parent_label FROM ontology_edge
WHERE child_key IN ({{frontier*}})
    AND (op_unset({{resources}}) OR list_contains({{resources}}, split_part(resource, '/', 1)))

-- step: parent_candidates AS list
SELECT DISTINCT parent_key FROM __found WHERE NOT list_contains({{known}}, parent_key) ORDER BY 1

-- step: known AS value
SELECT {{known}} || coalesce((SELECT list(DISTINCT entity_key) FROM entity_by_key
    WHERE entity_key IN ({{parent_candidates*}}) AND namespace = {{seed_namespace}}
        AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))), []::VARCHAR[])

-- step: walked AS table
SELECT * FROM __tree_edges
UNION
SELECT parent_key, parent_label, child_key, child_label FROM __found WHERE list_contains({{known}}, parent_key)

-- step: tree_edges AS table
SELECT * FROM __walked

-- step: frontier AS value
SELECT coalesce(list(DISTINCT parent_key ORDER BY parent_key), []::VARCHAR[]) FROM __found
WHERE list_contains({{known}}, parent_key) AND NOT list_contains({{seen}}, parent_key)

-- step: seen AS value
SELECT {{seen}} || {{frontier}}

-- step: walk_done AS value
SELECT len({{frontier}}) = 0

-- end

-- step: tree_keys AS list
SELECT DISTINCT k FROM (
    SELECT {{seed_key}} AS k UNION ALL SELECT parent_key FROM __tree_edges UNION ALL SELECT child_key FROM __tree_edges
    UNION ALL SELECT child_key FROM __seed_children)
WHERE k IS NOT NULL ORDER BY 1

-- step: meta AS table
-- each key's id and name: of its last copy, as a dict keeps them
SELECT entity_key, arg_max(struct_pack(
        id := coalesce(nullif(identifier, ''), nullif(label, ''), entity_key),
        name := coalesce(nullif(label, ''), coalesce(nullif(identifier, ''), nullif(label, ''), entity_key))),
    struct_pack(r := _resource, i := entity_id)) AS m
FROM entity_by_key
WHERE entity_key IN ({{tree_keys*}}) AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
GROUP BY entity_key

-- step: adjacency AS table
SELECT DISTINCT parent_key, child_key FROM __tree_edges
UNION SELECT {{seed_key}}, child_key FROM __seed_children

-- step: roots AS table
SELECT DISTINCT parent_key AS key FROM __adjacency
WHERE parent_key NOT IN (SELECT child_key FROM __tree_edges)

-- step: paths AS table
-- every path from a root (or the seed, without roots), a node met again on its path
-- is a leaf
WITH RECURSIVE walk(key, path, depth) AS (
    SELECT key, [key], 1 FROM (SELECT key FROM __roots UNION ALL
        SELECT {{seed_key}} WHERE (SELECT count(*) FROM __roots) = 0)
    UNION ALL
    SELECT a.child_key, w.path || [a.child_key], w.depth + 1
    FROM walk w JOIN __adjacency a ON a.parent_key = w.key
    WHERE NOT list_contains(w.path[1:-2], w.key) AND w.depth < 64
)
SELECT w.key, w.path, w.depth, coalesce(m.m.id, w.key) AS id, coalesce(m.m.name, w.key) AS name
FROM walk w LEFT JOIN __meta m ON m.entity_key = w.key

-- step: depth AS value
SELECT coalesce(max(depth), 0) FROM __paths

-- step: built AS table
SELECT NULL::VARCHAR[] AS path, NULL::VARCHAR AS name, NULL::VARCHAR AS key, NULL::JSON AS node WHERE false

-- repeat: d = {{depth}} next {{d}} - 1 until built_done max 64

-- step: built_next AS table
SELECT p.path, p.name, p.key, to_json(struct_pack(
        id := p.id, name := p.name,
        children := coalesce((SELECT list(b.node ORDER BY b.name, b.key) FROM __built b
            WHERE b.path[1:-2] = p.path AND len(b.path) = len(p.path) + 1), []::JSON[]),
        distance := CASE WHEN p.id = coalesce(nullif({{seed_identifier}}, ''), nullif({{seed_label}}, ''), {{term_id}})
            OR p.key = {{seed_key}} THEN 0 ELSE 1 END,
        hasMoreChildren := p.key = {{seed_key}} AND {{has_more}},
        nextCursor := CASE WHEN p.key = {{seed_key}} AND {{has_more}} THEN 100 END,
        ontologyId := {{seed_namespace}})) AS node
FROM __paths p WHERE p.depth = {{d}}
UNION ALL SELECT * FROM __built

-- step: built AS table
SELECT * FROM __built_next

-- step: built_done AS value
SELECT {{d}} <= 1

-- end

-- step: response AS value
SELECT CASE
    WHEN {{seed_key}} IS NULL THEN to_json(struct_pack(root := struct_pack(id := {{term_id}}, name := {{term_id}},
        children := []::VARCHAR[], distance := 0)))
    WHEN (SELECT count(*) FROM __tree_edges) = 0 AND (SELECT count(*) FROM __seed_children) = 0
        THEN to_json(struct_pack(root := NULL::JSON))
    WHEN (SELECT count(*) FROM __built WHERE len(path) = 1) = 1
        THEN to_json(struct_pack(root := (SELECT node FROM __built WHERE len(path) = 1 LIMIT 1)))
    ELSE to_json(struct_pack(root := struct_pack(id := 'root', name := 'Available ancestry (display root)',
        children := (SELECT list(node ORDER BY name, key) FROM __built WHERE len(path) = 1), distance := 1)))
END
