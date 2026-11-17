-- query: metsigdb_sets
-- The sets of the MetSigDB subset without their members: resource, id, label, type,
-- sub-type, organism and size, one row per set, in the order of the table's key.
-- param: resource VARCHAR[] = NULL
-- param: set_type VARCHAR[] = NULL
-- param: set_sub_type VARCHAR[] = NULL
-- param: organism BIGINT = NULL
-- param: sets VARCHAR[] = NULL
-- param: limit INTEGER = 1000
-- param: offset INTEGER = 0
-- param: total BOOLEAN = false
-- output: json response
-- example: {"set_type": ["pathway"], "limit": 5}

-- step: present AS value
SELECT CASE WHEN {{include:metsigdb_present x}} THEN true
    ELSE error('invalid: this release has no MetSigDB subset') END

-- step: page AS table
SELECT s.*, row_number() OVER (ORDER BY s.resource, s.set_source_id) AS n FROM (
    SELECT m.resource, m.set_source_id, any_value(m.set_entity_key) AS set_entity_key,
        any_value(m.set_label) AS set_label, any_value(m.set_type) AS set_type,
        any_value(m.set_sub_type) AS set_sub_type, any_value(m.organism) AS organism,
        any_value(m.set_size) AS set_size
    FROM metsigdb_membership m
    WHERE (op_unset({{resource}}) OR list_contains({{resource}}, m.resource))
        AND (op_unset({{set_type}}) OR list_contains({{set_type}}, m.set_type))
        AND (op_unset({{set_sub_type}}) OR list_contains({{set_sub_type}}, m.set_sub_type))
        AND ({{organism}} IS NULL OR m.organism = {{organism}})
        AND (op_unset({{sets}}) OR list_contains({{sets}}, m.set_source_id))
    GROUP BY m.resource, m.set_source_id
    ORDER BY m.resource, m.set_source_id
    LIMIT {{limit}} + 1 OFFSET {{offset}}
) s

-- step: matched AS value when {{total}}
SELECT count(*) FROM (SELECT 1 FROM metsigdb_membership m
    WHERE (op_unset({{resource}}) OR list_contains({{resource}}, m.resource))
        AND (op_unset({{set_type}}) OR list_contains({{set_type}}, m.set_type))
        AND (op_unset({{set_sub_type}}) OR list_contains({{set_sub_type}}, m.set_sub_type))
        AND ({{organism}} IS NULL OR m.organism = {{organism}})
        AND (op_unset({{sets}}) OR list_contains({{sets}}, m.set_source_id))
    GROUP BY m.resource, m.set_source_id)

-- step: response AS value
SELECT to_json(struct_pack(
    count := count(*) FILTER (WHERE n <= {{limit}}),
    has_more := count(*) > {{limit}},
    sets := coalesce(list(struct_pack(
        resource := resource, set := set_source_id, set_entity := set_entity_key,
        set_label := set_label, set_type := set_type, set_sub_type := set_sub_type,
        organism := organism, set_size := set_size) ORDER BY n) FILTER (WHERE n <= {{limit}}), []),
    total := {{matched}}))
FROM __page
