-- Each entity's numbers of merged identifiers and annotations.
SELECT entity_n, sum(i)::BIGINT AS identifier_count, sum(a)::BIGINT AS annotation_count FROM (
    SELECT entity_n, 1 AS i, 0 AS a FROM omnipath.entity_identifier
    UNION ALL SELECT entity_n, 0, 1 FROM omnipath.entity_annotation)
GROUP BY entity_n
