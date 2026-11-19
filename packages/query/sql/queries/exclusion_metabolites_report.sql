-- query: exclusion_metabolites_report
-- What the exclusion build measured, one row per finding: the size of each class and where
-- its members came from, how many seeds resolved, and the seeds and derived members the two
-- ways of listing disagree on. `detail` is JSON text.
-- param: kinds VARCHAR[] = NULL
-- param: classes VARCHAR[] = NULL
-- param: limit INTEGER = 1000
-- param: offset INTEGER = 0
-- param: total BOOLEAN = false
-- output: json response
-- example: {"kinds": ["class_size"]}

-- step: present AS value
SELECT CASE WHEN {{include:exclusion_present x}} THEN true
    ELSE error('invalid: this release has no exclusion metabolite subset') END

-- step: page AS table
SELECT r.*, row_number() OVER w AS n FROM (
    SELECT r.* FROM exclusion_metabolite_report r
    WHERE (op_unset({{kinds}}) OR list_contains({{kinds}}, r.kind))
        AND (op_unset({{classes}}) OR list_contains({{classes}}, r.exclusion_class))
    ORDER BY r.kind, coalesce(r.exclusion_class::VARCHAR, ''), coalesce(r.namespace::VARCHAR, ''), r.detail
    LIMIT {{limit}} + 1 OFFSET {{offset}}
) r
WINDOW w AS (ORDER BY r.kind, coalesce(r.exclusion_class::VARCHAR, ''), coalesce(r.namespace::VARCHAR, ''), r.detail)

-- step: matched AS value when {{total}}
SELECT count(*) FROM exclusion_metabolite_report r
WHERE (op_unset({{kinds}}) OR list_contains({{kinds}}, r.kind))
    AND (op_unset({{classes}}) OR list_contains({{classes}}, r.exclusion_class))

-- step: response AS value
SELECT to_json(struct_pack(
    count := count(*) FILTER (WHERE n <= {{limit}}),
    has_more := count(*) > {{limit}},
    rows := coalesce(list(struct_pack(
        kind := kind, exclusion_class := exclusion_class, namespace := namespace,
        resolved := resolved, fallback_only := fallback_only, unresolved := unresolved,
        detail := detail::JSON, build_id := build_id) ORDER BY n) FILTER (WHERE n <= {{limit}}), []),
    total := {{matched}}))
FROM __page
