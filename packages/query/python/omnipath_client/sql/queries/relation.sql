-- query: relation
-- One relation by its key, as the API's relation item: the relation summary and both
-- endpoints; NULL when no selected resource states it.
-- param: relation_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"relation_key": "missing"}

-- step: relation_keys AS list
SELECT {{relation_key}}

-- steps: relation_by_keys

-- step: first_id AS list
SELECT relation_n FROM __found LIMIT 1

-- step: endpoint_relations AS table
SELECT *, 1 AS position FROM relation WHERE relation_n IN ({{first_id*}})

-- steps: endpoint_rows

-- step: items AS table
WITH relations AS (SELECT *, op_qualifier_annotations(r) AS annotations FROM __endpoint_relations r)
SELECT * FROM {{include:relation_items relations}}

-- step: response AS value
SELECT to_json(item) FROM __items
