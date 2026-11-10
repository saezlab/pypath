-- The entity rows the entity_matches steps found (their phase's), as a subquery: the
-- prefix matches, or the rows of the other phases read in place.
SELECT * FROM __match_prefix WHERE {{match_phase}} = 'prefix'
UNION ALL BY NAME
SELECT * FROM entity_by_key
WHERE {{match_phase}} <> 'prefix'
    AND ({{match_keys}} IS NULL OR entity_key IN ({{match_keys*}}))
    AND op_entity_ok(entity_by_key, {{match_scope}}, {{type_names}}, {{taxonomy_ids}}, {{reference_entity_keys}}, NULL)
    AND CASE {{match_phase}}
        WHEN 'key' THEN entity_key = {{q}}
        WHEN 'contains' THEN contains(lower(label), {{ql}})
        ELSE true END
