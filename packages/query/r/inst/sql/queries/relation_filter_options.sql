-- query: relation_filter_options
-- The values the relation filters offer within a selection: predicates by category,
-- sources and participant types.
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
-- output: json response
-- example: {}

-- steps: relation_facet_items

-- step: response AS value
SELECT to_json(struct_pack(
    predicatesByCategory := coalesce((
        SELECT map_from_entries(list(struct_pack(k := c, v := predicates) ORDER BY c))
        FROM (SELECT coalesce(category, 'uncategorized') AS c, list(DISTINCT value ORDER BY value) AS predicates
              FROM __facet_items WHERE facet = 'predicate' GROUP BY 1)), map([]::VARCHAR[], []::VARCHAR[][])),
    sources := coalesce((SELECT list(DISTINCT value ORDER BY value) FROM __facet_items WHERE facet = 'source'), []::VARCHAR[]),
    interactionTypes := coalesce((SELECT list(DISTINCT value ORDER BY value) FROM __facet_items
        WHERE facet = 'participant_type'), []::VARCHAR[])))
