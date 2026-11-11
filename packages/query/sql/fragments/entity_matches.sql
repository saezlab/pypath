-- Steps: what a search matches in its first phase with matches (all, key, prefix,
-- contains: the phases value), in the resources {{match_scope}}, passing the entity
-- filters (type_names, taxonomy_ids, reference_entity_keys and the entity key list
-- match_keys, NULL for none): the value match_phase and, for prefixes, the matching
-- entities (__match_prefix). The entity_match_rows fragment reads every matching entity
-- (all hits, not a page) without copying them.

-- step: match_hits AS table when list_contains({{phases}}, 'prefix')
-- each entity a term names, with the namespaces it must have (a CURIE's), or any
SELECT entity_n, bool_or(namespace IS NULL) AS any_namespace,
    list(DISTINCT namespace) FILTER (WHERE namespace IS NOT NULL) AS namespaces
FROM (
    SELECT entity_n, NULL::VARCHAR AS namespace FROM entity_term WHERE term >= {{ql}} AND term < {{q_upper}}
    UNION ALL
    SELECT entity_n, {{curie_ns}} FROM entity_term WHERE {{curie_ns}} IS NOT NULL AND term = lower({{curie_id}})
    UNION ALL
    SELECT entity_n, NULL FROM entity_group
    WHERE {{curie_ns}} IS NOT NULL AND kind = 'reference' AND group_key = {{curie_ns}} || ':' || {{curie_id}})
GROUP BY entity_n

-- step: match_prefix AS table when list_contains({{phases}}, 'prefix')
SELECT e.entity_n, e.entity_key, e.namespace, e.label, e.entity_type, e.taxon, e.resources,
    e.relation_count, e.chemical, e.structure_group, e.gene_group
FROM entity e JOIN __match_hits h USING (entity_n)
WHERE (h.any_namespace OR list_contains(h.namespaces, lower(coalesce(e.namespace, ''))))
    AND op_entity_ok(e, {{match_scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{match_keys}})

-- step: match_phase AS value
SELECT CASE
    WHEN list_contains({{phases}}, 'prefix') AND (SELECT count(*) FROM __match_prefix) > 0 THEN 'prefix'
    WHEN list_contains({{phases}}, 'contains') THEN 'contains'
    ELSE {{phases}}[1]
END
