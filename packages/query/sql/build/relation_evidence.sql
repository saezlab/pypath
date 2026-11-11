-- Every evidence record of each relation, with the resource (and the resource's own
-- relation id) that states it.
SELECT c.relation_n, c.resource, c.resource_version, c.relation_id, e.ordinal, e."source", e.dataset,
    e.row_id, e.upstream_id, e.annotations, e.subject_molecular_form, e.object_molecular_form
FROM work.relation_evidence e JOIN work.relation_copies c USING (_rid)
WHERE {{partition c.relation_key}}
