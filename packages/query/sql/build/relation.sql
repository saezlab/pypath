-- The merged relations with their place in the page order (page_rank) and their number
-- of merged annotations.
SELECT m.*, k.page_rank, coalesce(a.n, 0)::BIGINT AS annotation_count
FROM work.relation_merged m
JOIN work.relation_rank k USING (relation_n)
LEFT JOIN (SELECT relation_n, count(*) AS n FROM omnipath.relation_annotation GROUP BY relation_n) a USING (relation_n)
