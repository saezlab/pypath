-- How big each class is, and how its members were found.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH by_source AS (
  SELECT exclusion_class, json_group_object(assignment_source, n) AS sources
  FROM (SELECT exclusion_class, assignment_source, count(*) AS n
        FROM exclusion_metabolite GROUP BY 1, 2)
  GROUP BY 1
)
SELECT 'class_size',
       c.cls,
       count(DISTINCT x.entity_key)::INTEGER,
       json_object(
         'entity_matches', count(DISTINCT x.entity_key) FILTER (WHERE x.match_kind = 'entity'),
         'structure_key_matches',
           count(DISTINCT x.entity_key) FILTER (WHERE x.match_kind = 'structure_key'),
         'by_source', coalesce(any_value(s.sources), '{}'::JSON)
       )::VARCHAR,
       $build_id
FROM unnest(['ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors', 'cofactors_extended'])
  AS c(cls)
LEFT JOIN exclusion_metabolite x ON x.exclusion_class = c.cls
LEFT JOIN by_source s ON s.exclusion_class = c.cls
GROUP BY c.cls;
