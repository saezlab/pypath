-- Each relation under the reference keys its copies give its ends (a gene's Entrez key,
-- for instance), sorted by key and then page order, with its filter and facet columns.
WITH refs AS (
    SELECT DISTINCT c.relation_n, x.side, x.k
    FROM work.relation r JOIN work.relation_copies c USING (_rid),
        LATERAL (VALUES ('subject', r.subject_reference_entity_key), ('object', r.object_reference_entity_key)) AS x(side, k)
    WHERE x.k IS NOT NULL
)
SELECT refs.k AS reference_key, refs.side, rel.relation_n, rel.page_rank, {{relation_columns}}
FROM refs JOIN omnipath.relation rel USING (relation_n)
