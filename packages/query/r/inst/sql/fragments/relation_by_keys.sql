-- Steps: the selected resources' relations with the keys {{relation_keys}} (a list),
-- by the keys' first 60 bits (relation_key_lookup), then the full key on the rows
-- found; a key that is no sha256 hex is looked for in the relation table itself.
-- Gives the table __found: relation_key, _rid, _resource, resource, relation_id, in
-- resource order.

-- step: key_prefixes AS list
SELECT DISTINCT try_cast('0x' || substr(k, 1, 15) AS BIGINT) AS p
FROM unnest({{relation_keys}}) AS u(k) WHERE regexp_full_match(k, '[0-9a-f]{64}') ORDER BY p

-- step: key_candidates AS list
SELECT _rid FROM relation_key_lookup WHERE key_prefix IN ({{key_prefixes*}}) ORDER BY _rid

-- step: key_rows AS table
SELECT relation_key, _rid, _resource, _resource_version, relation_id FROM relation WHERE _rid IN ({{key_candidates*}})

-- step: other_keys AS list
SELECT k FROM unnest({{relation_keys}}) AS u(k) WHERE NOT regexp_full_match(k, '[0-9a-f]{64}') ORDER BY k

-- step: other_rows AS table when {{#other_keys}} > 0
SELECT relation_key, _rid, _resource, _resource_version, relation_id FROM relation WHERE relation_key IN ({{other_keys*}})

-- step: found AS table
SELECT DISTINCT relation_key, _rid, _resource, _resource || '/' || _resource_version AS resource, relation_id
FROM (SELECT * FROM __key_rows UNION ALL SELECT * FROM __other_rows)
WHERE list_contains({{relation_keys}}, relation_key)
    AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY resource, relation_id
