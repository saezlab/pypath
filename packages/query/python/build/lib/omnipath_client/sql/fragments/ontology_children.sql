-- Steps: the children of the seed (ontology_seed) within its namespace: the table
-- __child_rows (child_key, child_label; one row per edge, by child key). The value
-- known lists the keys found in the namespace (the seed's, the children's).

-- step: known AS value
SELECT coalesce(list(DISTINCT entity_key), []::VARCHAR[]) FROM entity_by_key
WHERE {{seed_key}} IS NOT NULL AND entity_key = {{seed_key}} AND namespace = {{seed_namespace}}
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))

-- step: child_edges AS table when len({{known}}) > 0
SELECT DISTINCT child_key, child_label, predicate, parent_key, parent_label FROM ontology_edge
WHERE parent_key = {{seed_key}}
    AND (op_unset({{resources}}) OR list_contains({{resources}}, split_part(resource, '/', 1)))

-- step: child_candidates AS list
SELECT DISTINCT child_key FROM __child_edges ORDER BY 1

-- step: known AS value
SELECT {{known}} || coalesce((SELECT list(DISTINCT entity_key) FROM entity_by_key
    WHERE entity_key IN ({{child_candidates*}}) AND namespace = {{seed_namespace}}
        AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))), []::VARCHAR[])

-- step: child_rows AS table
SELECT child_key, child_label FROM __child_edges WHERE list_contains({{known}}, child_key)
ORDER BY child_key
