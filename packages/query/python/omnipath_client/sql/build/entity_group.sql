-- The entities under each group key their copies give (reference, gene, connectivity).
SELECT DISTINCT g.group_key, g.kind, c.entity_n
FROM work.entity_group g JOIN work.entity_copies c USING (_eid)
