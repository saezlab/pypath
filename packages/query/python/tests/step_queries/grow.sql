-- query: grow
-- param: target INTEGER = 5
-- output: json response
-- example: {"target": 5}
-- example: {"target": 100}

-- repeat: n = 1 next {{n}} * 2 until enough max 4
-- step: found AS table
SELECT _eid FROM item ORDER BY _eid LIMIT {{n}}
-- step: enough AS value
SELECT count(*) >= {{target}} FROM __found
-- end

-- step: response AS value
SELECT to_json({'n': {{n}}, 'found': (SELECT count(*) FROM __found), 'enough': {{enough}}})
