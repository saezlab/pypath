-- Each entity's identifiers over its copies, united by (ns, id) and ordered by where each
-- is first seen (copies in resource and entity order, each copy's in its order); read by
-- entity_merged. A stage of the build, in its scratch database; built in slices of
-- entity keys ({{partition ...}}), since list aggregates hold their groups in memory.
SELECT entity_key, list(struct_pack(ns := ns, id := id, is_canonical := is_canonical, source := source)
    ORDER BY first) AS identifiers
FROM (
    SELECT e.entity_key, i.ns, i.id,
        min(struct_pack(r := e._resource || '/' || e._resource_version, e := e.entity_id, o := i.ordinal)) AS first,
        arg_min(i.is_canonical, struct_pack(r := e._resource || '/' || e._resource_version, e := e.entity_id, o := i.ordinal)) AS is_canonical,
        arg_min(i.source, struct_pack(r := e._resource || '/' || e._resource_version, e := e.entity_id, o := i.ordinal)) AS source
    FROM entity_identifier i JOIN entity e USING (_eid)
    WHERE {{partition e.entity_key}}
    GROUP BY e.entity_key, i.ns, i.id)
GROUP BY entity_key
