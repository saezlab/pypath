-- Steps: the parts of a relation selection by the relation filters (the parameters of
-- the relations query): Biolink names, sign values, forms (molecular form filters
-- given), keyed (endpoint filters given), scope (the resources read), the entity and
-- scope keys, the candidate relations of the endpoint filters (__candidate_rows) and,
-- under form filters, the relations with matching evidence (form_ids), unless the
-- including query's value read_candidates is false. The relation_conditions fragment
-- applies the other filters.

-- step: predicate_names AS value
SELECT (SELECT list(name) FROM op_biolink('predicate', {{predicates}}))

-- step: type_names AS value
SELECT (SELECT list(name) FROM op_biolink('class',
    CASE WHEN op_unset({{interaction_types}}) THEN {{entity_types}} ELSE {{interaction_types}} END))

-- step: partner_names AS value
SELECT (SELECT list(name) FROM op_biolink('class', {{partner_types}}))

-- step: sign_values AS value
SELECT CASE WHEN op_unset({{signs}}) AND {{sign}} IS NOT NULL THEN [{{sign}}] ELSE {{signs}} END

-- step: forms AS value
SELECT NOT (op_unset({{protein_entity_keys}}) AND op_unset({{transcript_entity_keys}})
    AND op_unset({{isoform_identifiers}}))

-- step: keyed AS value
SELECT NOT (op_unset({{entities}}) AND op_unset({{scope_entity_ids}}) AND op_unset({{reference_entity_keys}}))

-- step: scope AS value
-- the resources read: the selected ones and the source filter's; [''] selects none
SELECT CASE
    WHEN op_unset({{resources}}) AND op_unset({{sources}}) THEN NULL
    WHEN op_unset({{sources}}) THEN {{resources}}
    WHEN op_unset({{resources}}) THEN {{sources}}
    ELSE coalesce(nullif(list_intersect({{resources}}, {{sources}}), []), [''])
END

-- step: key_tokens AS list
SELECT unnest(coalesce({{entities}}, []))

-- step: key_scope AS value
SELECT {{scope}}

-- steps: entity_keys

-- step: entity_key_list AS value
SELECT {{entity_keys}}

-- step: key_tokens AS list
SELECT unnest(coalesce({{scope_entity_ids}}, []))

-- steps: entity_keys

-- step: scope_key_list AS value
SELECT {{entity_keys}}

-- step: selection_keys AS list
SELECT DISTINCT k FROM (
    SELECT unnest(coalesce({{entity_key_list}}, [])) AS k
    UNION ALL SELECT unnest(coalesce({{scope_key_list}}, []))
    UNION ALL SELECT unnest(coalesce({{reference_entity_keys}}, []))) ORDER BY k

-- step: selection_endpoints AS table
SELECT key, key_kind, side, _rid FROM relation_endpoint
WHERE key IN ({{selection_keys*}}) AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))

-- step: candidate_ids AS list
WITH by_entity AS (
    SELECT _rid FROM __selection_endpoints
    WHERE key_kind = 'entity' AND list_contains(coalesce({{entity_key_list}}, []), key)
        AND ({{gene_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{gene_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY _rid
    HAVING {{gene_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
), by_scope AS (
    SELECT _rid FROM __selection_endpoints
    WHERE key_kind = 'entity' AND list_contains(coalesce({{scope_key_list}}, []), key)
        AND ({{scope_endpoint_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{scope_endpoint_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY _rid
    HAVING {{scope_endpoint_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
), by_reference AS (
    SELECT DISTINCT _rid FROM __selection_endpoints
    WHERE key_kind = 'reference' AND list_contains(coalesce({{reference_entity_keys}}, []), key)
)
SELECT DISTINCT _rid FROM __selection_endpoints
WHERE (op_unset({{entities}}) OR _rid IN (SELECT _rid FROM by_entity))
    AND (op_unset({{scope_entity_ids}}) OR _rid IN (SELECT _rid FROM by_scope))
    AND (op_unset({{reference_entity_keys}}) OR _rid IN (SELECT _rid FROM by_reference))
ORDER BY _rid

-- step: candidate_rows AS table when {{read_candidates}}
SELECT * FROM relation WHERE _rid IN ({{candidate_ids*}})

-- step: candidate_evidence AS table when {{forms}} AND {{keyed}} AND {{read_candidates}}
SELECT * FROM relation_evidence WHERE _rid IN ({{candidate_ids*}})

-- step: form_ids AS list when {{forms}}
-- relations with an evidence row matching the molecular form filters
SELECT DISTINCT _rid FROM (
    SELECT * FROM __candidate_evidence
    UNION ALL BY NAME
    SELECT * FROM relation_evidence
    WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))) ev
WHERE op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
ORDER BY _rid
