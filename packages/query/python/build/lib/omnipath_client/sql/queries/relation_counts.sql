-- query: relation_counts
-- The matching relations counted by category, by predicate (the 25 most) and by sign.
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
-- example: {"entities": ["TP53"]}

-- step: selection_sources AS value
SELECT {{sources}}

-- steps: relation_fast

-- step: read_candidates AS value
SELECT NOT {{fast}}

-- steps: relation_candidates

-- steps: relation_range

-- step: counts AS table
SELECT category, predicate, sign, count(*) AS n FROM (
    SELECT r.category, r.predicate, r.sign FROM __range_selected r
    WHERE {{fast}} AND {{include:relation_conditions r}}
    UNION ALL
    SELECT r.category, r.predicate, r.sign FROM __candidate_rows r
    WHERE NOT {{fast}} AND {{keyed}} AND {{include:relation_conditions r}}
    UNION ALL
    SELECT r.category, r.predicate, r.sign FROM relation r
    WHERE NOT {{keyed}} AND (op_unset({{scope}}) OR list_contains({{scope}}, r._resource))
        AND {{include:relation_conditions r}})
GROUP BY ALL

-- step: response AS value
WITH categories AS (
    SELECT category AS k, sum(n) AS n FROM __counts WHERE coalesce(category, '') <> '' GROUP BY 1
), predicates AS (
    SELECT predicate AS k, sum(n) AS n FROM __counts WHERE coalesce(predicate, '') <> '' GROUP BY 1
    ORDER BY n DESC, k LIMIT 25
), signs AS (
    SELECT sign::VARCHAR AS k, sum(n) AS n FROM __counts WHERE sign IS NOT NULL GROUP BY 1
)
SELECT to_json(struct_pack(
    categories := coalesce((SELECT map_from_entries(list(struct_pack(k := k, v := n) ORDER BY n DESC, k)) FROM categories), map([]::VARCHAR[], []::BIGINT[])),
    predicates := coalesce((SELECT map_from_entries(list(struct_pack(k := k, v := n) ORDER BY n DESC, k)) FROM predicates), map([]::VARCHAR[], []::BIGINT[])),
    signs := coalesce((SELECT map_from_entries(list(struct_pack(k := k, v := n) ORDER BY k)) FROM signs), map([]::VARCHAR[], []::BIGINT[]))))
