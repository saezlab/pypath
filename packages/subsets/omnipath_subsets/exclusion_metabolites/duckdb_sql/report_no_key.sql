-- Listed members that carry no structure key: nothing can be pulled in for them.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH unkeyed AS (
  SELECT x.exclusion_class, x.entity_key,
         min(coalesce(e.label, e.identifier)) AS label
  FROM exclusion_metabolite x
  JOIN release_entity e ON e.entity_key = x.entity_key
  WHERE x.match_kind = 'entity'
    AND NOT EXISTS (SELECT 1 FROM structure_group g WHERE g.entity_key = x.entity_key)
  GROUP BY x.exclusion_class, x.entity_key
),
numbered AS (
  SELECT u.*, row_number() OVER (PARTITION BY exclusion_class ORDER BY label, entity_key) AS n
  FROM unkeyed u
)
SELECT 'no_key',
       c.cls,
       count(u.entity_key)::INTEGER,
       json_object(
         'labels', coalesce(to_json(list(u.label ORDER BY u.n) FILTER (WHERE u.n <= 1000)),
                            '[]'::JSON),
         'truncated', count(u.entity_key) > 1000
       )::VARCHAR,
       $build_id
FROM unnest(['ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors', 'cofactors_extended'])
  AS c(cls)
LEFT JOIN numbered u ON u.exclusion_class = c.cls
GROUP BY c.cls;
