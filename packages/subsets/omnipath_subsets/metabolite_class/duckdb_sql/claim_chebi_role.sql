-- Signal 2: the chemicals that ChEBI gives a role, or a role below it.
--
-- The walk starts from ChEBI's own copy of the role term and descends ChEBI's
-- `subclass_of` edges, as the cofactor walk of the exclusion product does. A carrier is
-- the subject of a ChEBI `has_chemical_role` edge to one of those roles. The role terms
-- themselves are not carriers.
INSERT INTO mc_claim
WITH RECURSIVE edge AS (
  SELECT DISTINCT subject_entity_key AS child, object_entity_key AS parent
  FROM release_relation
  WHERE resource = 'chebi'
    AND statement_kind IN ('ontology', 'ontology_axiom')
    AND predicate = 'subclass_of'
),
root AS (
  SELECT entity_key FROM resource_entity WHERE resource = 'chebi' AND identifier = $root
  UNION
  SELECT e.entity_key
  FROM resource_identifier i
  JOIN resource_entity e ON e.resource = i.resource AND e.entity_id = i.entity_id
  WHERE i.resource = 'chebi' AND i.ns = 'chebi' AND i.id = $root
),
walk (entity_key, depth) AS (
  SELECT entity_key, 0 FROM root
  UNION
  SELECT e.child, w.depth + 1
  FROM walk w
  JOIN edge e ON e.parent = w.entity_key
  WHERE w.depth < 20
)
SELECT DISTINCT r.subject_entity_key, $metabolite_class, $priority, 'chebi'
FROM release_relation r
JOIN (SELECT DISTINCT entity_key FROM walk) role ON role.entity_key = r.object_entity_key
JOIN release_entity c ON c.entity_key = r.subject_entity_key
WHERE r.resource = 'chebi'
  AND r.statement_kind IN ('ontology', 'ontology_axiom')
  AND r.predicate = 'has_chemical_role'
  AND list_contains($entity_types, c.entity_type);
