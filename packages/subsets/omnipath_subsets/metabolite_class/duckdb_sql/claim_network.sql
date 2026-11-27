-- Signal 3: the chemicals that a reaction of the named resources takes as input or gives
-- as output. The reaction is the subject of the relation and the chemical its object.
--
-- The build writes the resources in as literals, not as a list parameter: a literal
-- condition on the resource lets DuckDB read only those resources' files. The type of the
-- chemical is the one the reaction's own resource gives it.
INSERT INTO mc_claim
SELECT DISTINCT r.object_entity_key, $metabolite_class, $priority, r.resource
FROM release_relation r
JOIN resource_entity c ON c.resource = r.resource AND c.entity_key = r.object_entity_key
WHERE r.resource IN ({resources})
  AND c.resource IN ({resources})
  AND list_contains($predicates, r.predicate)
  AND list_contains($entity_types, c.entity_type);
