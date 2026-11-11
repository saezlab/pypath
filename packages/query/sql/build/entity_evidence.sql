-- Every evidence record of each entity, with the resource (and the resource's own
-- entity id) that states it.
SELECT c.entity_n, c.resource, c.resource_version, c.entity_id, e.ordinal, e."source", e.dataset,
    e.row_id, e.upstream_id, e.annotations, e.molecular_form
FROM work.entity_evidence e JOIN work.entity_copies c USING (_eid)
