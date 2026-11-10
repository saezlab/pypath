-- query: shouted
-- param: limit INTEGER = 3
-- output: json response
-- example: {"limit": 2}

-- step: rows AS table
SELECT * FROM item WHERE _eid IN (1, 2, 3, 4)

-- step: response AS value
SELECT to_json(list(shout ORDER BY key)) FROM {{include:shout __rows}}
