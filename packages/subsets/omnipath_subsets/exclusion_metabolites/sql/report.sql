-- The coverage rows of the report: for each class and namespace, how many seeds
-- found an entity, how many found only a name-identified one, and which found
-- nothing.
--
-- The seed counts come from the resolution, which runs in Python so that one
-- spelling rule serves both the build and its tests. This file writes them.

INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, namespace, resolved, fallback_only, unresolved, detail, build_id)
SELECT 'seed_coverage',
       c.exclusion_class,
       c.namespace,
       c.resolved,
       c.fallback_only,
       c.unresolved,
       c.detail,
       %(build_id)s
FROM jsonb_to_recordset(%(coverage)s::jsonb)
  AS c(exclusion_class text, namespace text, resolved integer, fallback_only integer,
       unresolved integer, detail jsonb)
ORDER BY c.exclusion_class, c.namespace;

-- How big each class is, and how its members were found.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
SELECT 'class_size',
       c.cls,
       count(DISTINCT x.entity_id),
       jsonb_build_object(
         'entity_matches', count(DISTINCT x.entity_id) FILTER (WHERE x.match_kind = 'entity'),
         'structure_key_matches',
           count(DISTINCT x.entity_id) FILTER (WHERE x.match_kind = 'structure_key'),
         'by_source', coalesce((
           SELECT jsonb_object_agg(q.assignment_source, q.n)
           FROM (SELECT assignment_source, count(*) AS n FROM exclusion_metabolite
                 WHERE exclusion_class = c.cls GROUP BY 1) q
         ), '{}'::jsonb)
       ),
       %(build_id)s
FROM unnest(ARRAY['ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors', 'cofactors_extended']) AS c(cls)
LEFT JOIN exclusion_metabolite x ON x.exclusion_class = c.cls
GROUP BY c.cls;

-- Seeds the closure does not reach, and closure members the seeds lack. Only
-- for a class whose root resolved: without a root the whole class would read as
-- missed.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH missed AS (
  SELECT s.exclusion_class, s.entity_id, min(s.seed_identifier) AS identifier
  FROM exclusion_metabolite s
  WHERE s.assignment_source = 'metaproviz_seed' AND s.match_kind = 'entity'
    AND s.exclusion_class = ANY(%(derived_classes)s::text[])
    AND NOT EXISTS (
      SELECT 1 FROM exclusion_metabolite d
      WHERE d.entity_id = s.entity_id AND d.exclusion_class = s.exclusion_class
        AND d.assignment_source = 'chebi_hierarchy')
  GROUP BY s.exclusion_class, s.entity_id
)
SELECT 'seed_not_derived',
       c.cls,
       count(m.entity_id),
       jsonb_build_object(
         'entities', coalesce(jsonb_agg(m.entity_id ORDER BY m.entity_id)
                              FILTER (WHERE m.entity_id IS NOT NULL), '[]'::jsonb),
         'identifiers', coalesce(jsonb_agg(m.identifier ORDER BY m.entity_id)
                                 FILTER (WHERE m.entity_id IS NOT NULL), '[]'::jsonb)
       ),
       %(build_id)s
FROM unnest(%(derived_classes)s::text[]) AS c(cls)
LEFT JOIN missed m ON m.exclusion_class = c.cls
GROUP BY c.cls;

INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH extra AS (
  SELECT d.exclusion_class, d.entity_id,
         coalesce(e.label, e.canonical_identifier) AS label,
         row_number() OVER (PARTITION BY d.exclusion_class
                            ORDER BY coalesce(e.label, e.canonical_identifier), d.entity_id) AS n
  FROM exclusion_metabolite d
  JOIN entity e ON e.entity_id = d.entity_id
  WHERE d.assignment_source = 'chebi_hierarchy' AND d.match_kind = 'entity'
    AND NOT EXISTS (
      SELECT 1 FROM exclusion_metabolite s
      WHERE s.entity_id = d.entity_id AND s.exclusion_class = d.exclusion_class
        AND s.origin = 'curated' AND s.match_kind = 'entity')
)
SELECT 'derived_not_seeded',
       c.cls,
       count(x.entity_id),
       jsonb_build_object(
         'entities', coalesce(jsonb_agg(x.entity_id ORDER BY x.n)
                              FILTER (WHERE x.n <= 1000), '[]'::jsonb),
         'labels', coalesce(jsonb_agg(x.label ORDER BY x.n)
                            FILTER (WHERE x.n <= 1000), '[]'::jsonb),
         'truncated', count(x.entity_id) > 1000
       ),
       %(build_id)s
FROM unnest(%(derived_classes)s::text[]) AS c(cls)
LEFT JOIN extra x ON x.exclusion_class = c.cls
GROUP BY c.cls;

-- Keys that pulled in variants, so a person can see what a key match added and
-- remove a wanted metabolite with a veto.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH keyed AS (
  SELECT x.exclusion_class, g.group_key,
         min(x.via_entity_id::text) AS via,
         array_agg(DISTINCT coalesce(e.label, e.canonical_identifier)
                   ORDER BY coalesce(e.label, e.canonical_identifier)) AS labels
  FROM exclusion_metabolite x
  JOIN chemical_resolution_group_member g ON g.entity_id = x.entity_id AND g.level_id = 1
  JOIN entity e ON e.entity_id = x.entity_id
  WHERE x.match_kind = 'structure_key'
  GROUP BY x.exclusion_class, g.group_key
),
numbered AS (
  SELECT k.*, row_number() OVER (PARTITION BY exclusion_class ORDER BY group_key) AS n
  FROM keyed k
)
SELECT 'shared_key',
       exclusion_class,
       count(*),
       jsonb_build_object(
         'keys', jsonb_agg(jsonb_build_object('key', group_key, 'via', via, 'labels', to_jsonb(labels))
                           ORDER BY n) FILTER (WHERE n <= 500),
         'truncated', count(*) > 500
       ),
       %(build_id)s
FROM numbered
GROUP BY exclusion_class;

-- Listed members that carry no structure key: nothing can be pulled in for them.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH unkeyed AS (
  SELECT x.exclusion_class, x.entity_id,
         min(coalesce(e.label, e.canonical_identifier)) AS label
  FROM exclusion_metabolite x
  JOIN entity e ON e.entity_id = x.entity_id
  WHERE x.match_kind = 'entity'
    AND NOT EXISTS (
      SELECT 1 FROM chemical_resolution_group_member g
      WHERE g.entity_id = x.entity_id AND g.level_id = 1)
  GROUP BY x.exclusion_class, x.entity_id
),
numbered AS (
  SELECT u.*, row_number() OVER (PARTITION BY exclusion_class ORDER BY label, entity_id) AS n
  FROM unkeyed u
)
SELECT 'no_key',
       c.cls,
       count(u.entity_id),
       jsonb_build_object(
         'labels', coalesce(jsonb_agg(u.label ORDER BY u.n) FILTER (WHERE u.n <= 1000),
                            '[]'::jsonb),
         'truncated', count(u.entity_id) > 1000
       ),
       %(build_id)s
FROM unnest(ARRAY['ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors', 'cofactors_extended']) AS c(cls)
LEFT JOIN numbered u ON u.exclusion_class = c.cls
GROUP BY c.cls;

-- Every veto, with what it removed. A veto that removed a derived or key-matched
-- row is the case a reader needs to see; one that removed nothing is here too,
-- because it may be stale.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, namespace, resolved, detail, build_id)
SELECT 'veto_applied',
       v.exclusion_class,
       v.namespace,
       count(r.entity_id),
       jsonb_build_object(
         'identifier', v.identifier,
         'note', v.note,
         'resolved', v.resolved,
         'entities', v.entities,
         'removed', coalesce(
           jsonb_agg(jsonb_build_object('entity', r.entity_id,
                                        'assignment_source', r.assignment_source,
                                        'match_kind', r.match_kind)
                     ORDER BY r.entity_id, r.assignment_source, r.match_kind)
             FILTER (WHERE r.entity_id IS NOT NULL), '[]'::jsonb)
       ),
       %(build_id)s
FROM jsonb_to_recordset(%(vetoes)s::jsonb)
  AS v(exclusion_class text, namespace text, identifier text, note text, resolved boolean,
       entities jsonb)
LEFT JOIN exclusion_veto_removed r
  ON r.veto_class = v.exclusion_class AND r.namespace = v.namespace
 AND r.identifier = v.identifier
GROUP BY v.exclusion_class, v.namespace, v.identifier, v.note, v.resolved, v.entities;
