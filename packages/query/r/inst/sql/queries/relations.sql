-- query: relations
-- Relations as the explorer exports them: one row per relation with its annotations,
-- evidence and the product records its evidence names. Endpoint filters (entities,
-- scope entities, reference keys) read the endpoint table's ranges; the other filters
-- are conditions on the relations. At most max_relations relations, in global id order.
-- A relation lists only its own resource among its sources, so a source filter reads
-- only those resources.
-- param: entities VARCHAR[] = NULL
-- param: gene_mode VARCHAR = "any"
-- param: scope_entity_ids VARCHAR[] = NULL
-- param: scope_endpoint_mode VARCHAR = "any"
-- param: reference_entity_keys VARCHAR[] = NULL
-- param: sources VARCHAR[] = NULL
-- param: predicates VARCHAR[] = NULL
-- param: relation_categories VARCHAR[] = NULL
-- param: taxonomy_ids VARCHAR[] = NULL
-- param: interaction_types VARCHAR[] = NULL
-- param: entity_types VARCHAR[] = NULL
-- param: partner_types VARCHAR[] = NULL
-- param: is_directed BOOLEAN = NULL
-- param: sign INTEGER = NULL
-- param: signs INTEGER[] = NULL
-- param: min_evidence_count INTEGER = NULL
-- param: object_aspect_qualifier VARCHAR[] = NULL
-- param: object_direction_qualifier VARCHAR[] = NULL
-- param: causal_mechanism_qualifier VARCHAR[] = NULL
-- param: protein_entity_keys VARCHAR[] = NULL
-- param: transcript_entity_keys VARCHAR[] = NULL
-- param: isoform_identifiers VARCHAR[] = NULL
-- param: molecular_endpoint_mode VARCHAR = "any"
-- param: max_relations INTEGER = 100000
-- param: resources VARCHAR[] = NULL
-- output: table rows
-- example: {"entities": ["TP53"]}

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

-- step: endpoint_keys AS list
SELECT DISTINCT k FROM (
    SELECT unnest(coalesce({{entity_key_list}}, [])) AS k
    UNION ALL SELECT unnest(coalesce({{scope_key_list}}, []))
    UNION ALL SELECT unnest(coalesce({{reference_entity_keys}}, []))) ORDER BY k

-- step: endpoints AS table
SELECT key, key_kind, side, _rid FROM relation_endpoint
WHERE key IN ({{endpoint_keys*}}) AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))

-- step: candidate_ids AS list
WITH by_entity AS (
    SELECT _rid FROM __endpoints
    WHERE key_kind = 'entity' AND list_contains(coalesce({{entity_key_list}}, []), key)
        AND ({{gene_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{gene_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY _rid
    HAVING {{gene_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
), by_scope AS (
    SELECT _rid FROM __endpoints
    WHERE key_kind = 'entity' AND list_contains(coalesce({{scope_key_list}}, []), key)
        AND ({{scope_endpoint_mode}} NOT IN ('source', 'source_only') OR side = 'subject')
        AND ({{scope_endpoint_mode}} NOT IN ('target', 'target_only') OR side = 'object')
    GROUP BY _rid
    HAVING {{scope_endpoint_mode}} <> 'both' OR (bool_or(side = 'subject') AND bool_or(side = 'object'))
), by_reference AS (
    SELECT DISTINCT _rid FROM __endpoints
    WHERE key_kind = 'reference' AND list_contains(coalesce({{reference_entity_keys}}, []), key)
)
SELECT DISTINCT _rid FROM __endpoints
WHERE (op_unset({{entities}}) OR _rid IN (SELECT _rid FROM by_entity))
    AND (op_unset({{scope_entity_ids}}) OR _rid IN (SELECT _rid FROM by_scope))
    AND (op_unset({{reference_entity_keys}}) OR _rid IN (SELECT _rid FROM by_reference))
ORDER BY _rid

-- step: candidate_rows AS table
SELECT * FROM relation WHERE _rid IN ({{candidate_ids*}})

-- step: candidate_evidence AS table when {{forms}} AND {{keyed}}
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

-- step: selected AS table
SELECT r.* FROM (
    SELECT * FROM __candidate_rows
    UNION ALL BY NAME
    SELECT * FROM relation WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))
) r
WHERE (NOT {{forms}} OR r._rid IN ({{form_ids*}}))
    AND (op_unset({{relation_categories}}) OR list_contains({{relation_categories}}, r.category))
    AND (op_unset({{predicate_names}}) OR list_contains({{predicate_names}}, r.predicate))
    AND (op_unset({{object_aspect_qualifier}})
        OR list_has_any(r.object_aspect_qualifier, {{object_aspect_qualifier}}))
    AND (op_unset({{object_direction_qualifier}})
        OR list_has_any(r.object_direction_qualifier, {{object_direction_qualifier}}))
    AND (op_unset({{causal_mechanism_qualifier}})
        OR list_has_any(r.causal_mechanism_qualifier, {{causal_mechanism_qualifier}}))
    AND ({{is_directed}} IS NULL OR r.is_directed = {{is_directed}})
    AND (op_unset({{sign_values}}) OR list_contains({{sign_values}}, r.sign))
    AND (coalesce({{min_evidence_count}}, 0) <= 1 OR r.evidence_count >= {{min_evidence_count}})
    AND (op_unset({{type_names}})
        OR (op_unset({{partner_names}})
            AND (list_contains({{type_names}}, r.subject_type) OR list_contains({{type_names}}, r.object_type)))
        OR (NOT op_unset({{partner_names}})
            AND ((list_contains({{type_names}}, r.subject_type) AND list_contains({{partner_names}}, r.object_type))
              OR (list_contains({{type_names}}, r.object_type) AND list_contains({{partner_names}}, r.subject_type)))))
    AND (op_unset({{taxonomy_ids}}) OR list_contains({{taxonomy_ids}}, r.taxon))
    AND (op_unset({{sources}}) OR list_has_any(r.sources, {{sources}}))
ORDER BY r._rid
LIMIT {{max_relations}}

-- step: selected_ids AS list
SELECT _rid FROM __selected ORDER BY _rid

-- step: annotation_rows AS table
SELECT * FROM relation_annotation WHERE _rid IN ({{selected_ids*}})

-- step: evidence_rows AS table
SELECT * FROM relation_evidence WHERE _rid IN ({{selected_ids*}})

-- step: evidence AS table
SELECT _rid, list(struct_pack(
    source := source, dataset := dataset, row_id := row_id,
    upstream_id := upstream_id, annotations := annotations,
    subject_molecular_form := subject_molecular_form,
    object_molecular_form := object_molecular_form) ORDER BY ordinal) AS evidence
FROM __evidence_rows ev
WHERE NOT {{forms}} OR op_evidence_ok(ev, {{protein_entity_keys}}, {{transcript_entity_keys}},
    {{isoform_identifiers}}, {{molecular_endpoint_mode}})
GROUP BY ALL

-- step: product_keys AS list
SELECT DISTINCT k FROM (
    SELECT unnest([f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key]) AS k
    FROM __evidence, unnest(__evidence.evidence) AS occurrences(f))
WHERE k IS NOT NULL ORDER BY k

-- step: product_rows AS table
SELECT * EXCLUDE (key_rank, _eid, _resource, _resource_version), _resource || '/' || _resource_version AS resource
FROM entity_by_key
WHERE entity_key IN ({{product_keys*}}) AND (op_unset({{scope}}) OR list_contains({{scope}}, _resource))

-- step: rows AS table
WITH annotations AS (
    SELECT _rid, list(struct_pack(term := term, value := value,
        quantity := op_quantity(quantity_has_numeric_value, quantity_has_unit, quantity_has_unit_prefix,
            quantity_has_binary_relation, quantity_source_field, quantity_comparator),
        source := source, dataset := dataset, scope := scope) ORDER BY ordinal) AS annotations
    FROM __annotation_rows GROUP BY ALL
), nested AS (
    SELECT s.* EXCLUDE (evidence_count),
        coalesce(a.annotations, []) AS annotations,
        coalesce(e.evidence, []) AS evidence,
        CASE WHEN {{forms}} THEN len(coalesce(e.evidence, [])) ELSE s.evidence_count END AS evidence_count
    FROM __selected s LEFT JOIN annotations a USING (_rid) LEFT JOIN __evidence e USING (_rid)
), product_keys AS (
    SELECT DISTINCT n._rid, unnest([
        f.subject_molecular_form.protein_entity_key, f.subject_molecular_form.transcript_entity_key,
        f.object_molecular_form.protein_entity_key, f.object_molecular_form.transcript_entity_key
    ]) AS entity_key
    FROM nested n, unnest(n.evidence) AS occurrences(f)
), products AS (
    SELECT k._rid, list(pe ORDER BY pe.entity_key, pe.resource) AS referenced_product_records
    FROM product_keys k JOIN __product_rows pe USING (entity_key)
    GROUP BY ALL
)
SELECT nested.* EXCLUDE (_rid, _resource, _resource_version, relation_id),
    coalesce(products.referenced_product_records, []) AS referenced_product_records
FROM nested LEFT JOIN products USING (_rid)
ORDER BY nested._rid
