-- The lowercase labels and identifiers of every copy, sorted by term: a prefix is one
-- range of entity numbers.
SELECT DISTINCT t.term, t.kind, c.entity_n
FROM work.entity_term t JOIN work.entity_copies c USING (_eid)
