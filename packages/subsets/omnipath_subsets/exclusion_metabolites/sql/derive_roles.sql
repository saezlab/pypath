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
  (entity_id, exclusion_class, assignment_source, origin, match_kind, seed_namespace,
   seed_identifier, evidence, via_entity_id, list_version, build_id)
WITH RECURSIVE root (exclusion_class, entity_id, curie) AS (
  SELECT * FROM unnest(%(role_classes)s::text[], %(role_entities)s::uuid[], %(role_curies)s::text[])
),
walk (exclusion_class, entity_id, curie, depth) AS (
  SELECT exclusion_class, entity_id, curie, 0 FROM root
  UNION
  SELECT w.exclusion_class, e.child, w.curie, w.depth + 1
  FROM walk w
  JOIN exclusion_subclass_edge e ON e.parent = w.entity_id
  WHERE w.depth < 20
),
role AS (
  SELECT DISTINCT ON (exclusion_class, entity_id) exclusion_class, entity_id, curie, depth
  FROM walk
  ORDER BY exclusion_class, entity_id, depth, curie
),
carrier AS (
  SELECT DISTINCT ON (r.subject_entity_id, role.exclusion_class)
         r.subject_entity_id AS entity_id, role.exclusion_class, role.curie, role.depth
  FROM entity_ontology_relation r
  JOIN data_source s ON s.source_id = r.source_id AND s.name = 'chebi'
  JOIN vocab_relation_predicate p
    ON p.relation_predicate_id = r.predicate_id AND p.name = 'has_chemical_role'
  JOIN role ON role.entity_id = r.object_entity_id
  ORDER BY r.subject_entity_id, role.exclusion_class, role.depth, role.curie
)
SELECT c.entity_id,
       c.exclusion_class,
       'chebi_role',
       'derived',
       'entity',
       NULL,
       NULL,
       CASE WHEN c.depth = 0 THEN 'has_chemical_role ' || c.curie
            ELSE 'has_chemical_role a role below ' || c.curie || ' at depth ' || c.depth END,
       NULL,
       %(list_version)s,
       %(build_id)s
FROM carrier c
JOIN entity e ON e.entity_id = c.entity_id
JOIN vocab_entity_type vt
  ON vt.entity_type_id = e.entity_type_id AND vt.name = ANY(%(entity_types)s::text[]);
