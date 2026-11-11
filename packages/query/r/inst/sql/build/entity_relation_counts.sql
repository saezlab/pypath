-- Each entity's number of relations (a relation from an entity to itself counts once).
SELECT entity_n, count(*)::BIGINT AS relation_count FROM (
    SELECT subject_n AS entity_n FROM work.relation_merged
    UNION ALL SELECT object_n FROM work.relation_merged WHERE object_n <> subject_n)
GROUP BY entity_n
