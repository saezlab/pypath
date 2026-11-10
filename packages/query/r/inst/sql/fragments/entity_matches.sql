-- Steps: the entity rows a search matches in its first phase with matches (all, key,
-- prefix, contains: the phases value), in the resources {{match_scope}}, passing the
-- entity filters (type_names, taxonomy_ids, reference_entity_keys and the entity key
-- list match_keys, NULL for none): the table __matches. Every hit, not a page.

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
SELECT DISTINCT e.* FROM entity e JOIN __match_hits h USING (_eid)
WHERE (h.namespace IS NULL OR h.namespace = lower(coalesce(e.namespace, '')))
    AND op_entity_ok(e, {{match_scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{match_keys}})

-- step: match_phase AS value
SELECT CASE
    WHEN list_contains({{phases}}, 'prefix') AND (SELECT count(*) FROM __match_prefix) > 0 THEN 'prefix'
    WHEN list_contains({{phases}}, 'contains') THEN 'contains'
    ELSE {{phases}}[1]
END

-- step: matches AS table
SELECT * FROM __match_prefix WHERE {{match_phase}} = 'prefix'
UNION ALL BY NAME
SELECT * FROM entity
WHERE {{match_phase}} <> 'prefix'
    AND op_entity_ok(entity, {{match_scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, {{match_keys}})
    AND CASE {{match_phase}}
        WHEN 'key' THEN entity_key = {{q}}
        WHEN 'contains' THEN contains(lower(label), {{ql}})
        ELSE true END
