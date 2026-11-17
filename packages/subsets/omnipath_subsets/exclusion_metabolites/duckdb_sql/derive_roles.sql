-- Cofactors from the chemical ontology: every chemical that carries the ChEBI
-- cofactor role, or a role below it, by a `has_chemical_role` edge.
--
-- A role is an entity of its own in this release, joined to a chemical by an
-- edge. The role term and the roles under it are found by the same bounded walk
-- as the ion closure. The role term itself is not a cofactor and is not listed.
--
-- The role misses partners that are cofactors in practice (ADP, CoA, GTP on the
-- current build), which is what the reviewed core list is for.

INSERT INTO exclusion_metabolite
WITH RECURSIVE walk (exclusion_class, entity_key, curie, depth) AS (
  SELECT exclusion_class, entity_key, curie, 0 FROM exclusion_role_root
  UNION
  SELECT w.exclusion_class, e.child, w.curie, w.depth + 1
  FROM walk w
  JOIN exclusion_subclass_edge e ON e.parent = w.entity_key
  WHERE w.depth < 20
),
role AS (
  SELECT DISTINCT ON (exclusion_class, entity_key) exclusion_class, entity_key, curie, depth
  FROM walk
  ORDER BY exclusion_class, entity_key, depth, curie
),
carrier AS (
  SELECT DISTINCT ON (r.subject_entity_key, role.exclusion_class)
         r.subject_entity_key AS entity_key, role.exclusion_class, role.curie, role.depth
  FROM release_relation r
  JOIN role ON role.entity_key = r.object_entity_key
  WHERE r.resource = 'chebi'
    AND r.statement_kind IN ('ontology', 'ontology_axiom')
    AND r.predicate = 'has_chemical_role'
  ORDER BY r.subject_entity_key, role.exclusion_class, role.depth, role.curie
)
SELECT c.entity_key,
       c.exclusion_class,
       'chebi_role',
       'derived',
       'entity',
       NULL,
       NULL,
       CASE WHEN c.depth = 0 THEN 'has_chemical_role ' || c.curie
            ELSE 'has_chemical_role a role below ' || c.curie || ' at depth ' || c.depth END,
       NULL,
       $list_version,
       $build_id
FROM carrier c
JOIN release_entity e ON e.entity_key = c.entity_key
WHERE list_contains($entity_types, e.entity_type);
