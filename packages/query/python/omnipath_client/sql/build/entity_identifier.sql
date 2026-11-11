-- Each entity's identifiers over its copies, united by (ns, id), numbered (ordinal) in
-- the order each is first seen: copies by resource and id (their global id), each
-- copy's in its order. The first sighting is one integer, so that the grouping spills
-- to disk instead of needing memory for all of it.
SELECT entity_n, row_number() OVER (PARTITION BY entity_n ORDER BY first)::INTEGER AS ordinal,
    ns, id, is_canonical, "source"
FROM (
    SELECT c.entity_n, i.ns, i.id, min((i._eid::HUGEINT << 32) + i.ordinal) AS first,
        arg_min(i.is_canonical, (i._eid::HUGEINT << 32) + i.ordinal) AS is_canonical,
        arg_min(i."source", (i._eid::HUGEINT << 32) + i.ordinal) AS "source"
    FROM work.entity_identifier i JOIN work.entity_copies c USING (_eid)
    GROUP BY c.entity_n, i.ns, i.id)
