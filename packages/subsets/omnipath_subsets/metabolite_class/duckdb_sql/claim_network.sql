-- Signal 3: the chemicals that a reaction of the named resources takes as input or gives
-- as output. The reaction is the subject of the relation and the chemical its object.
INSERT INTO mc_claim
SELECT DISTINCT r.object_entity_key, $metabolite_class, $priority, r.resource
FROM release_relation r
JOIN release_entity c ON c.entity_key = r.object_entity_key
WHERE list_contains($resources, r.resource)
  AND list_contains($predicates, r.predicate)
  AND list_contains($entity_types, c.entity_type);
