-- query: entity_relationships
-- A page of the structural relationships of entities (by key): each relation with both
-- endpoints and its annotations, the total and the next page's offset.
-- param: entity_keys VARCHAR[] = NULL
-- param: limit INTEGER = 50
-- param: offset INTEGER = 0
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_keys": ["missing"]}

-- step: rel_keys AS list
SELECT DISTINCT unnest(coalesce({{entity_keys}}, [])) AS k ORDER BY k

-- step: rel_limit AS value
SELECT {{limit}}

-- step: rel_offset AS value
SELECT {{offset}}

-- steps: structural_relationships

-- step: response AS value
WITH page AS (SELECT count(*) AS n, any_value(total) AS total FROM __rel_page)
SELECT to_json(struct_pack(
    relationships := coalesce((SELECT list(item ORDER BY position) FROM __rel_items), []),
    relationshipsTotal := CASE WHEN n = 0 THEN 0 ELSE total END,
    nextCursor := CASE WHEN n > 0 AND {{offset}} + n < total THEN ({{offset}} + n)::VARCHAR END))
FROM page
