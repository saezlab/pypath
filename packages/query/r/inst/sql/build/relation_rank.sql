-- Each relation's place in the page order of every relation: category, predicate, the
-- ends' labels, key.
SELECT r.relation_n, row_number() OVER (ORDER BY r.category, r.predicate, s.label, o.label, r.relation_key)::INTEGER AS page_rank
FROM work.relation_merged r
JOIN omnipath.entity s ON s.entity_n = r.subject_n
JOIN omnipath.entity o ON o.entity_n = r.object_n
