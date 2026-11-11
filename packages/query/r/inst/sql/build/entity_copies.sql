-- Each resource's entity row (copy) with its entity's number, by global id: how the
-- resources' identifier, annotation, evidence, term and group rows find their entity.
SELECT e._eid, n.entity_n, e.entity_key, e._resource AS resource, e._resource_version AS resource_version,
    e.entity_id
FROM work.entity e JOIN work.entity_numbers n USING (entity_key)
