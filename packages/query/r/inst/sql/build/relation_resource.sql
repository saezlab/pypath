-- The resources that state each relation, one row per relation and resource, sorted by
-- resource: the relations of a resource are one range of rows (the relation table's
-- list of resources is read in every row).
SELECT unnest(resources) AS resource, relation_n FROM omnipath.relation
