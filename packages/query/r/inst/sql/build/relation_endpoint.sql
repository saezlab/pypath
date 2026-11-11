-- Each relation under both its ends' entity numbers, sorted by entity and then page
-- order, with its filter and facet columns: an entity's relations are one range, in
-- page order.
SELECT entity_n, side, relation_n, page_rank, {{relation_columns}}
FROM (
    SELECT subject_n AS entity_n, 'subject' AS side, * FROM omnipath.relation
    UNION ALL
    SELECT object_n, 'object', * FROM omnipath.relation)
