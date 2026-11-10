-- query: relation
-- One relation by its key, as the API's relation item: the relation summary and both
-- endpoints, from the first selected resource stating it; NULL when none does.
-- param: relation_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"relation_key": "missing"}

-- step: relation_keys AS list
SELECT {{relation_key}}

-- steps: relation_by_keys

-- step: first_id AS list
SELECT _rid FROM __found ORDER BY resource, relation_id LIMIT 1

-- step: rows AS table
SELECT *, 1 AS position FROM relation WHERE _rid IN ({{first_id*}})

-- step: endpoint_keys AS list
SELECT DISTINCT k FROM (SELECT subject_entity_key AS k FROM __rows UNION SELECT object_entity_key FROM __rows)
WHERE k IS NOT NULL ORDER BY k

-- step: endpoint_rows AS table
SELECT * FROM entity_by_key WHERE entity_key IN ({{endpoint_keys*}})

-- step: items AS table
WITH relations AS (SELECT *, op_qualifier_annotations(r) AS annotations FROM __rows r)
SELECT * FROM {{include:relation_items relations}}

-- step: response AS value
SELECT to_json(item) FROM __items
