-- Steps: the parts of a relation selection by the relation filters (the parameters of
-- the relations query): Biolink names, sign values, forms (molecular form filters
-- given), keyed (endpoint filters given), scope (the resources read), the entity and
-- scope keys, the candidate relations of the endpoint filters (__candidate_rows) and,
-- under form filters, the relations with matching evidence (form_ids), unless the
-- including query's value read_candidates is false. The including query's value
-- selection_sources (a source filter, or NULL) narrows the resources read. The relation_conditions fragment
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
    WHEN op_unset({{resources}}) AND op_unset({{selection_sources}}) THEN NULL
    WHEN op_unset({{selection_sources}}) THEN {{resources}}
    WHEN op_unset({{resources}}) THEN {{selection_sources}}
    ELSE coalesce(nullif(list_intersect({{resources}}, {{selection_sources}}), []), [''])
END

-- step: key_tokens AS list
SELECT unnest(coalesce({{entities}}, []))

-- step: key_scope AS value
SELECT {{scope}}

-- steps: entity_keys

-- step: entity_key_list AS value
SELECT {{entity_keys}}

-- step: entity_n_list AS list
SELECT entity_n FROM entity WHERE entity_key IN ({{entity_key_list*}}) ORDER BY entity_n

-- step: key_tokens AS list
SELECT unnest(coalesce({{scope_entity_ids}}, []))

-- steps: entity_keys

-- step: scope_key_list AS value
SELECT {{entity_keys}}

-- step: scope_n_list AS list
SELECT entity_n FROM entity WHERE entity_key IN ({{scope_key_list*}}) ORDER BY entity_n

-- step: selection_ns AS list
SELECT DISTINCT n FROM (SELECT unnest({{entity_n_list}}::INTEGER[]) AS n
    UNION ALL SELECT unnest({{scope_n_list}}::INTEGER[])) ORDER BY n

-- step: selection_endpoints AS table when {{read_candidates}}
SELECT entity_n, side, relation_n FROM relation_endpoint
WHERE entity_n IN ({{selection_ns+}}) AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))

-- step: selection_references AS table when {{read_candidates}} AND NOT op_unset({{reference_entity_keys}})
SELECT relation_n FROM relation_reference
WHERE reference_key IN ({{reference_entity_keys+}}) AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))

-- step: candidate_ids AS list when {{read_candidates}}
WITH by_entity AS (
    SELECT relation_n FROM __selection_endpoints
    WHERE list_contains({{entity_n_list}}, entity_n)
        AND ({{gene_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{gene_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY relation_n
    HAVING {{gene_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
), by_scope AS (
    SELECT relation_n FROM __selection_endpoints
    WHERE list_contains({{scope_n_list}}, entity_n)
        AND ({{scope_endpoint_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{scope_endpoint_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY relation_n
    HAVING {{scope_endpoint_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
)
SELECT DISTINCT relation_n FROM (
    SELECT relation_n FROM __selection_endpoints UNION ALL SELECT relation_n FROM __selection_references)
WHERE (op_unset({{entities}}) OR relation_n IN (SELECT relation_n FROM by_entity))
    AND (op_unset({{scope_entity_ids}}) OR relation_n IN (SELECT relation_n FROM by_scope))
    AND (op_unset({{reference_entity_keys}}) OR relation_n IN (SELECT relation_n FROM __selection_references))
ORDER BY relation_n

-- step: candidate_rows AS table when {{read_candidates}}
SELECT * FROM relation WHERE relation_n IN ({{candidate_ids*}})

-- step: candidate_evidence AS table when {{forms}} AND {{keyed}} AND {{read_candidates}}
SELECT * FROM relation_evidence WHERE relation_n IN ({{candidate_ids*}})

-- step: form_ids AS list when {{forms}}
-- relations with an evidence record matching the molecular form filters
SELECT DISTINCT relation_n FROM (
    SELECT * FROM __candidate_evidence
    UNION ALL BY NAME
    SELECT e.* FROM relation_evidence e
    WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, e.resource))) ev
WHERE op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
ORDER BY relation_n
