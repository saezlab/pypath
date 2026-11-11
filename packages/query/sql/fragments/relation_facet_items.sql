-- Steps: the relation filter options a selection leaves (__facet_items: facet, value,
-- category, n, label, position), most relations first. The scope (endpoint filters, molecular
-- forms, direction and a single sign) selects the relations, as combinations of their
-- facet columns; each filter group counts its values under the other groups' filters,
-- not its own, so its other values stay listed. Without an endpoint or form scope the
-- release's relation_facet table holds the combinations.

-- step: selection_sources AS value
SELECT NULL::VARCHAR[]

-- steps: relation_fast

-- step: read_candidates AS value
SELECT NOT {{fast}}

-- steps: relation_candidates

-- steps: relation_range

-- step: entity_type_names AS value
SELECT (SELECT list(name) FROM op_biolink('class', {{entity_types}}))

-- step: combinations AS table
WITH chosen AS (
    SELECT category, predicate, subject_type, object_type, taxon, sign, is_directed, sources,
        object_aspect_qualifier, object_direction_qualifier, causal_mechanism_qualifier,
        resources, evidence_count
    FROM __range_selected WHERE {{fast}}
    UNION ALL BY NAME
    SELECT category, predicate, subject_type, object_type, taxon, sign, is_directed, sources,
        object_aspect_qualifier, object_direction_qualifier, causal_mechanism_qualifier,
        resources, evidence_count
    FROM __candidate_rows WHERE NOT {{fast}} AND {{keyed}} AND (NOT {{forms}} OR relation_n IN ({{form_ids*}}))
    UNION ALL BY NAME
    SELECT category, predicate, subject_type, object_type, taxon, sign, is_directed, sources,
        object_aspect_qualifier, object_direction_qualifier, causal_mechanism_qualifier,
        resources, evidence_count
    FROM relation WHERE NOT {{keyed}} AND {{forms}} AND relation_n IN ({{form_ids*}})
        AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))
)
SELECT category, predicate, subject_type, object_type, taxon, sign, is_directed, sources,
    object_aspect_qualifier, object_direction_qualifier, causal_mechanism_qualifier, resources,
    max(evidence_count) AS evidence_count, count(*) AS relation_count
FROM chosen
WHERE ({{is_directed}} IS NULL OR is_directed = {{is_directed}}) AND ({{sign}} IS NULL OR sign = {{sign}})
GROUP BY ALL
UNION ALL BY NAME
SELECT category, predicate, subject_type, object_type, taxon, sign, is_directed, sources,
    object_aspect_qualifier, object_direction_qualifier, causal_mechanism_qualifier, resources,
    evidence_count, relation_count
FROM relation_facet
WHERE NOT ({{keyed}} OR {{forms}}) AND (op_unset({{scope}}) OR list_has_any(resources, {{scope}}))
    AND ({{is_directed}} IS NULL OR is_directed = {{is_directed}}) AND ({{sign}} IS NULL OR sign = {{sign}})

-- step: facet_rows AS table
-- each combination with whether each filter group's filters keep it
SELECT *,
        (op_unset({{relation_categories}}) OR list_contains({{relation_categories}}, category))
            AND (op_unset({{predicate_names}}) OR list_contains({{predicate_names}}, predicate)) AS g_type,
        op_unset({{type_names}})
            OR (op_unset({{partner_names}})
                AND (list_contains({{type_names}}, subject_type) OR list_contains({{type_names}}, object_type)))
            OR (NOT op_unset({{partner_names}})
                AND ((list_contains({{type_names}}, subject_type) AND list_contains({{partner_names}}, object_type))
                  OR (list_contains({{type_names}}, object_type) AND list_contains({{partner_names}}, subject_type)))) AS g_part,
        -- the participant facets leave out the type filters; an entity type filter is
        -- the fallback that stays
        op_unset({{entity_type_names}}) OR list_contains({{entity_type_names}}, subject_type)
            OR list_contains({{entity_type_names}}, object_type) AS g_part_fallback,
        op_unset({{sources}}) OR list_has_any(sources, {{sources}}) AS g_source,
        op_unset({{taxonomy_ids}}) OR list_contains({{taxonomy_ids}}, taxon) AS g_taxon,
        op_unset({{sign_values}}) OR list_contains({{sign_values}}, sign) AS g_sign,
        op_unset({{object_aspect_qualifier}}) OR list_has_any(object_aspect_qualifier, {{object_aspect_qualifier}}) AS g_q1,
        op_unset({{object_direction_qualifier}}) OR list_has_any(object_direction_qualifier, {{object_direction_qualifier}}) AS g_q2,
        op_unset({{causal_mechanism_qualifier}}) OR list_has_any(causal_mechanism_qualifier, {{causal_mechanism_qualifier}}) AS g_q3,
        ({{is_directed}} IS NULL OR is_directed = {{is_directed}})
            AND (coalesce({{min_evidence_count}}, 0) <= 1 OR evidence_count >= {{min_evidence_count}}) AS g_rest
    FROM __combinations

-- step: facet_items AS table
WITH c AS (SELECT * FROM __facet_rows), values_ AS (
    SELECT 'category' AS facet, category AS value, NULL::VARCHAR AS category_, relation_count
    FROM c WHERE g_rest AND g_part AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'predicate', predicate, category, relation_count
    FROM c WHERE g_rest AND g_part AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'participant_type', unnest(list_distinct([subject_type, object_type])), NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part_fallback AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'participant_pair', CASE WHEN subject_type <> '' AND object_type <> ''
        THEN array_to_string(list_sort([subject_type, object_type]), '|') END, NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part_fallback AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'source', unnest(list_distinct(sources)), NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_taxon AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'taxonomy_id', taxon, NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_source AND g_sign AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'sign', sign::VARCHAR, NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_source AND g_taxon AND g_q1 AND g_q2 AND g_q3
    UNION ALL
    SELECT 'object_aspect_qualifier', unnest(list_distinct(object_aspect_qualifier)), NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_source AND g_taxon AND g_sign AND g_q2 AND g_q3
    UNION ALL
    SELECT 'object_direction_qualifier', unnest(list_distinct(object_direction_qualifier)), NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q3
    UNION ALL
    SELECT 'causal_mechanism_qualifier', unnest(list_distinct(causal_mechanism_qualifier)), NULL, relation_count
    FROM c WHERE g_rest AND g_type AND g_part AND g_source AND g_taxon AND g_sign AND g_q1 AND g_q2
)
, counted AS (
    SELECT facet, value, category_ AS category, sum(relation_count) AS n
    FROM values_ WHERE value IS NOT NULL AND value <> ''
    GROUP BY facet, value, category_
)
-- taxa by name, joined on the counted values (an equality, not a loop over every row)
SELECT c.*, t.name AS label, row_number() OVER (ORDER BY c.n DESC, c.facet, c.value, coalesce(c.category, '')) AS position
FROM counted c LEFT JOIN op_taxonomy t ON t.taxon = CASE WHEN c.facet = 'taxonomy_id' THEN c.value END
ORDER BY position
