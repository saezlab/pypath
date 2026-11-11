-- The reference keys each relation's copies give its ends, once per relation, end and key.
SELECT DISTINCT c.relation_n, x.side, x.k AS reference_key
FROM work.relation r JOIN work.relation_copies c USING (_rid),
    LATERAL (VALUES ('subject', r.subject_reference_entity_key), ('object', r.object_reference_entity_key)) AS x(side, k)
WHERE x.k IS NOT NULL
