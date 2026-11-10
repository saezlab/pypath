-- query: relation_facets
-- The relation filter options a selection leaves, most relations first (each group's
-- values under the other groups' filters). With taxonomy_limit, only that many taxa
-- (those matching taxonomy_query, by id or name) and the selected ones.
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
-- param: resources VARCHAR[] = NULL
-- param: taxonomy_limit INTEGER = NULL
-- param: taxonomy_query VARCHAR = ""
-- output: json response
-- output: table items
-- example: {"entities": ["TP53"]}

-- steps: relation_facet_items

-- step: items AS table
WITH listed AS (
    SELECT facet, value, category, label, position FROM __facet_items
), missing AS (
    -- the selected taxa no relation has, after the others
    SELECT 'taxonomy_id' AS facet, t AS value, NULL::VARCHAR AS category, x.name AS label,
        1000000 + row_number() OVER (ORDER BY t) AS position
    FROM (SELECT DISTINCT unnest(coalesce({{taxonomy_ids}}, [])) AS t) s
    LEFT JOIN op_taxonomy x ON x.taxon = s.t
    WHERE {{taxonomy_limit}} IS NOT NULL AND t NOT IN (SELECT value FROM listed WHERE facet = 'taxonomy_id')
), everything AS (SELECT * FROM listed UNION ALL SELECT * FROM missing),
visible AS (
    SELECT value FROM (
        SELECT value, position FROM everything
        WHERE facet = 'taxonomy_id' AND (coalesce({{taxonomy_query}}, '') = ''
            OR contains(lower(value), lower(trim({{taxonomy_query}})))
            OR contains(lower(coalesce(label, '')), lower(trim({{taxonomy_query}}))))
        ORDER BY position LIMIT greatest(1, least(coalesce({{taxonomy_limit}}, 0), 10000)))
    UNION SELECT unnest(coalesce({{taxonomy_ids}}, []))
)
SELECT facet, value, category, label, position FROM everything
WHERE facet <> 'taxonomy_id' OR {{taxonomy_limit}} IS NULL OR value IN (SELECT value FROM visible)
ORDER BY position

-- step: response AS value
SELECT to_json(coalesce(list(CASE WHEN facet = 'taxonomy_id'
        THEN to_json(struct_pack(facetName := facet, facetValue := value, facetCategory := category, facetLabel := label))
        ELSE to_json(struct_pack(facetName := facet, facetValue := value, facetCategory := category)) END
    ORDER BY position), []))
FROM __items
