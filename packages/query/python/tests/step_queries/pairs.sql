-- query: pairs
-- param: labels VARCHAR[] = NULL
-- param: mode VARCHAR = "any"
-- output: table matched
-- output: json response
-- example: {"labels": ["Alpha", "Gamma"]}
-- example: {"labels": ["Alpha"], "mode": "strict"}
-- example: {"labels": ["Alpha"], "mode": "other"}

-- step: checked AS value
SELECT CASE WHEN {{mode}} IN ('any', 'strict') THEN true ELSE error('invalid: mode must be any or strict') END

-- step: wanted AS rows
SELECT key, label FROM item WHERE list_contains({{labels}}, label) ORDER BY key

-- step: matched AS table
SELECT w.key, w.label, i._eid FROM {{table:wanted}} w JOIN item i USING (key)
WHERE {{mode}} = 'any' OR i._eid > 1 ORDER BY w.key

-- step: response AS value
SELECT to_json({'matched': (SELECT coalesce(list(key ORDER BY key), []) FROM __matched), 'mode': {{mode}}})
