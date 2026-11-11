-- Each resource's relation row (copy) with its relation's number, by global id.
SELECT r._rid, n.relation_n, r.relation_key, r._resource AS resource, r._resource_version AS resource_version,
    r.relation_id
FROM work.relation r JOIN work.relation_numbers n USING (relation_key)
