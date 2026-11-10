-- query: lookup
-- param: keys VARCHAR[] = NULL
-- param: offset INTEGER = 0
-- param: limit INTEGER = 2
-- output: table page
-- output: json response
-- example: {"keys": ["b", "a", "zz"]}
-- example: {"keys": []}

-- step: ids AS list
SELECT _eid FROM item WHERE key IN ({{keys*}}) ORDER BY _eid

-- step: rows AS table when {{#ids}} > 0
SELECT * FROM item WHERE _eid IN ({{ids*}})

-- step: page AS table
SELECT key, label FROM __rows WHERE _eid IN ({{ids[offset:limit]}}) ORDER BY key

-- step: response AS value
SELECT to_json({'items': coalesce(list({'key': key, 'label': label} ORDER BY key), []), 'total': {{#ids}}}) FROM __page
