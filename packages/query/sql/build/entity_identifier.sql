-- Each entity's identifiers over its copies, united by (ns, id), numbered (ordinal) in
-- the order each is first seen: copies by resource and id, each copy's in its order.
SELECT entity_n, row_number() OVER (PARTITION BY entity_n ORDER BY first)::INTEGER AS ordinal,
    ns, id, is_canonical, "source"
FROM (
    SELECT c.entity_n, i.ns, i.id,
        min(struct_pack(r := c.resource, v := c.resource_version, e := c.entity_id, o := i.ordinal)) AS first,
        arg_min(i.is_canonical, struct_pack(r := c.resource, v := c.resource_version, e := c.entity_id, o := i.ordinal)) AS is_canonical,
        arg_min(i."source", struct_pack(r := c.resource, v := c.resource_version, e := c.entity_id, o := i.ordinal)) AS "source"
    FROM work.entity_identifier i JOIN work.entity_copies c USING (_eid)
    GROUP BY c.entity_n, i.ns, i.id)
