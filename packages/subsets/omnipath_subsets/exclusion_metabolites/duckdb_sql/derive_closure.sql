-- Ions and atoms from the chemical ontology: every `subclass_of` descendant of
-- the class's root, the root included.
--
-- The seeds of MetaProViz are a floor. The ontology knows more, and it moves
-- with each release, so the descendants are read here and not typed in. A member
-- the seeds also name gets its own row for this source.
--
-- The depth guard stops a cycle in the source from running away. The walk is
-- a UNION, not a UNION ALL: a row already seen at the same depth is not walked
-- again, so a densely connected cycle costs entities times depth and not the
-- number of paths through it. The shortest path to each entity wins.

INSERT INTO exclusion_metabolite
WITH RECURSIVE walk (exclusion_class, entity_key, curie, depth) AS (
  SELECT exclusion_class, entity_key, curie, 0 FROM exclusion_root
  UNION
  SELECT w.exclusion_class, e.child, w.curie, w.depth + 1
  FROM walk w
  JOIN exclusion_subclass_edge e ON e.parent = w.entity_key
  WHERE w.depth < 20
),
nearest AS (
  SELECT DISTINCT ON (exclusion_class, entity_key)
         exclusion_class, entity_key, curie, depth
  FROM walk
  ORDER BY exclusion_class, entity_key, depth, curie
)
SELECT n.entity_key,
       n.exclusion_class,
       'chebi_hierarchy',
       'derived',
       'entity',
       NULL,
       NULL,
       CASE WHEN n.depth = 0 THEN 'the root ' || n.curie
            ELSE 'subclass_of descendant of ' || n.curie || ' at depth ' || n.depth END,
       NULL,
       $list_version,
       $build_id
FROM nearest n
JOIN release_entity e ON e.entity_key = n.entity_key
WHERE list_contains($entity_types, e.entity_type);
