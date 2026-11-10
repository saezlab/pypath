-- Steps: the ontology term {{term_id}} (an identifier, label or entity key, or ns|id) in
-- the namespace {{ontology_id}} (NULL: any): the value seed (entity_key, label,
-- identifier, namespace; NULL when none), found among the sorted terms and by key. Two
-- entities with that name are an error, as is a namespace the ns|id contradicts.

-- step: seed_parts AS value
SELECT CASE
    WHEN contains({{term_id}}, '|') AND {{ontology_id}} IS NOT NULL
        AND lower({{ontology_id}}) <> lower(split_part({{term_id}}, '|', 1))
        THEN error('invalid: Term namespace contradicts ontologyId')
    WHEN contains({{term_id}}, '|') THEN [split_part({{term_id}}, '|', 1), substr({{term_id}}, strpos({{term_id}}, '|') + 1)]
    ELSE [{{ontology_id}}, {{term_id}}]
END

-- step: seed_hits AS list
SELECT DISTINCT _eid FROM entity_term
WHERE term = lower({{seed_parts}}[2]) AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource))
ORDER BY 1

-- step: seed_rows AS table
SELECT DISTINCT entity_key, label, identifier, namespace FROM (
    SELECT entity_key, label, identifier, namespace FROM entity WHERE _eid IN ({{seed_hits*}})
    UNION ALL
    SELECT entity_key, label, identifier, namespace FROM entity_by_key
    WHERE entity_key = {{seed_parts}}[2] AND (op_unset({{resources}}) OR list_contains({{resources}}, _resource)))
WHERE {{seed_parts}}[2] IN (entity_key, identifier, label)
    AND ({{seed_parts}}[1] IS NULL OR lower(coalesce(namespace, '')) = lower({{seed_parts}}[1]))

-- step: seed_checked AS value
SELECT CASE WHEN (SELECT count(DISTINCT entity_key) FROM __seed_rows) > 1
    THEN error('invalid: Ambiguous ontology term; use ontologyId and an exact identifier or entity key')
    ELSE true END

-- step: seed_first AS table
SELECT * FROM __seed_rows
ORDER BY coalesce(entity_key, ''), coalesce(label, ''), coalesce(identifier, ''), coalesce(namespace, '') LIMIT 1

-- step: seed_key AS value
SELECT entity_key FROM __seed_first

-- step: seed_label AS value
SELECT label FROM __seed_first

-- step: seed_identifier AS value
SELECT identifier FROM __seed_first

-- step: seed_namespace AS value
SELECT namespace FROM __seed_first
