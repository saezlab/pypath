-- Steps: what a search matches in its first phase with matches (all, key, prefix,
-- contains: the phases value), in the resources {{match_scope}}, passing the entity
-- filters (type_names, taxonomy_ids, reference_entity_keys and the entity key list
-- match_keys, NULL for none): the value match_phase and, for prefixes, the matching rows
-- (__match_prefix). The entity_match_rows fragment reads every matching row (all
-- hits, not a page) without copying them.

-- step: match_hits AS table when list_contains({{phases}}, 'prefix')
SELECT DISTINCT _eid, namespace FROM (
    SELECT _eid, NULL::VARCHAR AS namespace FROM entity_term
    WHERE term >= {{ql}} AND term < {{q_upper}} AND (op_unset({{match_scope}}) OR list_contains({{match_scope}}, _resource))
    UNION ALL
    SELECT _eid, {{curie_ns}} FROM entity_term
    WHERE {{curie_ns}} IS NOT NULL AND term = lower({{curie_id}})
        AND (op_unset({{match_scope}}) OR list_contains({{match_scope}}, _resource))
    UNION ALL
    SELECT _eid, NULL FROM entity_group
    WHERE {{curie_ns}} IS NOT NULL AND kind = 'reference' AND group_key = {{curie_ns}} || ':' || {{curie_id}}
        AND (op_unset({{match_scope}}) OR list_contains({{match_scope}}, _resource)))

-- step: match_prefix AS table when list_contains({{phases}}, 'prefix')
SELECT DISTINCT e.* FROM entity_by_key e JOIN __match_hits h USING (_eid)
WHERE (h.namespace IS NULL OR h.namespace = lower(coalesce(e.namespace, '')))
    AND op_entity_ok(e, {{match_scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{match_keys}})

-- step: match_phase AS value
SELECT CASE
    WHEN list_contains({{phases}}, 'prefix') AND (SELECT count(*) FROM __match_prefix) > 0 THEN 'prefix'
    WHEN list_contains({{phases}}, 'contains') THEN 'contains'
    ELSE {{phases}}[1]
END
