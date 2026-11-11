-- Steps: the endpoint entities of the relations in __endpoint_relations: the table
-- __endpoint_rows (entity_n, entity_key, label, summary), which the relation_items fragment
-- reads.

-- step: endpoint_ids AS list
SELECT DISTINCT n FROM (
    SELECT subject_n AS n FROM __endpoint_relations UNION SELECT object_n FROM __endpoint_relations)
WHERE n IS NOT NULL ORDER BY n

-- step: endpoint_rows AS table
SELECT entity_n, entity_key, label, op_summary(summary) AS summary FROM entity WHERE entity_n IN ({{endpoint_ids*}})
