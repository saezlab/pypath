-- Each relation under the reference keys its copies give its ends (a gene's Entrez key,
-- for instance), sorted by key and then page order, with its filter and facet columns.
SELECT refs.reference_key, refs.side, rel.relation_n, rel.page_rank, {{relation_columns}}
FROM work.relation_refs refs JOIN omnipath.relation rel USING (relation_n)
