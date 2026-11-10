-- Steps: the entity keys of {{key_tokens}} (a list) in the resources {{key_scope}}
-- (NULL: all): labels, identifiers, ns:id / ns|id CURIEs, reference keys, entity keys
-- and group keys (gene:<reference>, connectivity:<block>, entity:<key>), by exact
-- case-insensitive match. Gives the list entity_keys: the entity keys given, then the
-- groups' members, then the matches, sorted, without repeats.

-- step: key_classes AS table
SELECT token, position,
    CASE
        WHEN starts_with(token, 'gene:') OR starts_with(token, 'connectivity:') OR starts_with(token, 'entity:') THEN 'group'
        WHEN regexp_full_match(token, '[0-9a-fA-F]{64}') THEN 'key'
        ELSE 'other'
    END AS kind,
    nullif(lower(regexp_extract(token, '^([A-Za-z][\w.-]*)[:|](.+)$', 1)), '') AS curie_ns,
    nullif(regexp_extract(token, '^([A-Za-z][\w.-]*)[:|](.+)$', 2), '') AS curie_id
FROM (SELECT trim(t) AS token, i AS position
      FROM (SELECT unnest({{key_tokens}}) AS t, generate_subscripts({{key_tokens}}, 1) AS i)
      WHERE trim(coalesce(t, '')) <> '')

-- step: key_terms AS list
SELECT DISTINCT t FROM (
    SELECT lower(token) AS t FROM __key_classes WHERE kind = 'other'
    UNION SELECT lower(curie_id) FROM __key_classes WHERE kind = 'other' AND curie_ns IS NOT NULL)
ORDER BY t

-- step: key_term_rows AS table
SELECT term, _resource, _eid FROM entity_term
WHERE term IN ({{key_terms*}}) AND (op_unset({{key_scope}}) OR list_contains({{key_scope}}, _resource))

-- step: key_groups AS list
SELECT DISTINCT g FROM (
    SELECT curie_ns || ':' || curie_id AS g FROM __key_classes WHERE kind = 'other' AND curie_ns IS NOT NULL
    UNION SELECT substr(token, strpos(token, ':') + 1) FROM __key_classes
    WHERE kind = 'group' AND NOT starts_with(token, 'entity:'))
ORDER BY g

-- step: key_group_rows AS table
SELECT kind, group_key, _resource, _eid FROM entity_group
WHERE group_key IN ({{key_groups*}}) AND (op_unset({{key_scope}}) OR list_contains({{key_scope}}, _resource))

-- step: key_others AS list
SELECT DISTINCT token FROM __key_classes WHERE kind = 'other' ORDER BY token

-- step: key_rows AS table
SELECT entity_key, _resource, _eid FROM entity_by_key
WHERE entity_key IN ({{key_others*}}) AND (op_unset({{key_scope}}) OR list_contains({{key_scope}}, _resource))

-- step: key_hits AS table
-- (global id, the namespace the entity must have or NULL for any, from a group)
SELECT _eid, NULL::VARCHAR AS namespace, false AS member FROM __key_rows
UNION ALL
SELECT t._eid, NULL, false FROM __key_term_rows t
WHERE t.term IN (SELECT lower(token) FROM __key_classes WHERE kind = 'other')
UNION ALL
SELECT t._eid, c.curie_ns, false FROM __key_term_rows t
JOIN __key_classes c ON c.kind = 'other' AND c.curie_ns IS NOT NULL AND t.term = lower(c.curie_id)
UNION ALL
SELECT g._eid, NULL, false FROM __key_group_rows g
JOIN __key_classes c ON c.kind = 'other' AND c.curie_ns IS NOT NULL
    AND g.kind = 'reference' AND g.group_key = c.curie_ns || ':' || c.curie_id
UNION ALL
SELECT g._eid, NULL, true FROM __key_group_rows g
JOIN __key_classes c ON c.kind = 'group' AND NOT starts_with(c.token, 'entity:')
    AND g.kind = CASE WHEN starts_with(c.token, 'gene:') THEN 'reference' ELSE 'connectivity' END
    AND g.group_key = substr(c.token, strpos(c.token, ':') + 1)

-- step: key_hit_ids AS list
SELECT DISTINCT _eid FROM __key_hits ORDER BY _eid

-- step: key_hit_entities AS table
SELECT _eid, entity_key, namespace FROM entity WHERE _eid IN ({{key_hit_ids*}})

-- step: entity_keys AS list
WITH ordered AS (
    SELECT token AS k, 1 AS part, position AS n FROM __key_classes WHERE kind = 'key'
    UNION ALL
    SELECT substr(token, 8), 2, position FROM __key_classes WHERE kind = 'group' AND starts_with(token, 'entity:')
    UNION ALL
    SELECT DISTINCT e.entity_key, 3, 0 FROM __key_hits h JOIN __key_hit_entities e USING (_eid) WHERE h.member
    UNION ALL
    SELECT DISTINCT e.entity_key, 4, 0 FROM __key_hits h JOIN __key_hit_entities e USING (_eid)
    WHERE NOT h.member AND (h.namespace IS NULL OR h.namespace = lower(coalesce(e.namespace, '')))
)
SELECT k FROM ordered GROUP BY k ORDER BY min(struct_pack(p := part, n := n, k := k))
