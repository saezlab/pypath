-- Steps: the endpoint entities of the relations in __endpoint_relations, as merged at
-- build: the table __endpoint_rows (entity_key, summary), which the relation_items
-- fragment reads. One indexed lookup per key.

-- step: endpoint_keys AS list
SELECT DISTINCT k FROM (
    SELECT subject_entity_key AS k FROM __endpoint_relations UNION SELECT object_entity_key FROM __endpoint_relations)
WHERE k IS NOT NULL ORDER BY k

-- step: endpoint_rows AS table
SELECT entity_key, op_summary(summary) AS summary FROM entity_merged WHERE entity_key IN ({{endpoint_keys*}})
