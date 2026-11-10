-- Steps: the endpoint entities of the relations in __endpoint_relations, each from its
-- relation's own resource, with the columns an endpoint's summary shows: the table
-- __endpoint_rows (which the relation_items fragment reads). The keys are found in the
-- key-sorted copy reading three columns, then the rows by global id (an index).

-- step: endpoint_keys AS list
SELECT DISTINCT k FROM (
    SELECT subject_entity_key AS k FROM __endpoint_relations UNION SELECT object_entity_key FROM __endpoint_relations)
WHERE k IS NOT NULL ORDER BY k

-- step: endpoint_candidates AS table
SELECT entity_key, _resource, _eid FROM entity_by_key WHERE entity_key IN ({{endpoint_keys*}})

-- step: endpoint_ids AS list
SELECT DISTINCT c._eid FROM __endpoint_candidates c
JOIN (SELECT _resource, subject_entity_key AS k FROM __endpoint_relations
      UNION SELECT _resource, object_entity_key FROM __endpoint_relations) e
    ON e._resource = c._resource AND e.k = c.entity_key
ORDER BY 1

-- step: endpoint_rows AS table
SELECT entity_key, _resource, namespace, identifier, label, entity_type, taxon, has_hierarchy,
    parent_count, child_count, relation_count, gene_reference_keys, reference_entity_key
FROM entity WHERE _eid IN ({{endpoint_ids*}})
