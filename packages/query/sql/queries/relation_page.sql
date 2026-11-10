-- query: relation_page
-- A page of the relations matching the explorer's filters, in page order (category,
-- predicate, labels, resource, id), with both endpoints, the total and the next page's
-- offset. A plain entity or reference key selection is one range of the endpoint
-- table, in page order already; other selections read the candidates' relations.
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
-- param: limit INTEGER = 50
-- param: offset INTEGER = 0
-- output: json response
-- output: table items
-- example: {"entities": ["TP53"], "limit": 20}

-- step: details AS value
SELECT false

-- steps: relation_page_steps

-- step: response AS value
SELECT to_json(struct_pack(
    relations := coalesce((SELECT list(item.relation ORDER BY position) FROM __items), []),
    rows := coalesce((SELECT list(item ORDER BY position) FROM __items), []),
    total := {{total}},
    nextCursor := CASE WHEN {{offset}} + {{#page_ids}} < {{total}} THEN ({{offset}} + {{#page_ids}})::VARCHAR END))
