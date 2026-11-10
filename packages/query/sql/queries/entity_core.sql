-- query: entity_core
-- An entity merged from its resource copies, as the API's {"entity": ...}; NULL when no
-- selected resource has it.
-- param: entity_key VARCHAR
-- param: resources VARCHAR[] = NULL
-- output: json response
-- example: {"entity_key": "missing"}

-- steps: entity_main

-- step: response AS value when {{entity_json}} IS NOT NULL
SELECT to_json(struct_pack(entity := {{entity_json}}::JSON))
