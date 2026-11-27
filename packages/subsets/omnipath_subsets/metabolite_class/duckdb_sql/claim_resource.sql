-- Signal 1: the chemicals of one release resource. A cross-reference from another resource
-- does not count; only the resource's own records do.
INSERT INTO mc_claim
SELECT DISTINCT entity_key, $metabolite_class, $priority, resource
FROM resource_entity
WHERE resource = $resource
  AND list_contains($entity_types, entity_type);
