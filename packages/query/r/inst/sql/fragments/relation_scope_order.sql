-- Steps: the relations of the resources read (the including query's value scope) that
-- the relation filters keep, without endpoint filters: the table __scope_order
-- (relation_n, page_rank, category, predicate, sign). The relations of a resource are
-- one range of the table relation_resource (sorted by resource), where a list of
-- resources on the relation table is read in every row. Empty unless a scope is read.
-- Needs the relation_candidates steps.

-- step: scope_order AS table when NOT {{keyed}} AND NOT op_unset({{scope}})
SELECT r.relation_n, r.page_rank, r.category, r.predicate, r.sign FROM relation r
WHERE r.relation_n IN (SELECT relation_n FROM relation_resource WHERE resource IN ({{scope*}}))
    AND {{include:relation_conditions r}}
