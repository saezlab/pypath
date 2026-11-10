-- query: resolve_keys
-- The entity keys of labels, identifiers, CURIEs, reference keys, entity keys and the
-- explorer's group keys, by exact case-insensitive match: the keys given, then the
-- groups' members, then the matches (sorted), without repeats.
-- param: tokens VARCHAR[] = NULL
-- param: resources VARCHAR[] = NULL
-- output: table keys
-- output: json response
-- example: {"tokens": ["TP53", "entrez:7157"]}

-- step: key_tokens AS list
SELECT unnest(coalesce({{tokens}}, []))

-- step: key_scope AS value
SELECT {{resources}}

-- steps: entity_keys

-- step: keys AS table
SELECT unnest({{entity_keys}}::VARCHAR[]) AS entity_key

-- step: response AS value
SELECT to_json(coalesce({{entity_keys}}::VARCHAR[], []::VARCHAR[]))
