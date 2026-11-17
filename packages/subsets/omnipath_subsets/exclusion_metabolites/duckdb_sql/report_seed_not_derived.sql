-- Seeds the closure does not reach. Only for a class whose root resolved:
-- without a root the whole class would read as missed.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH missed AS (
  SELECT s.exclusion_class, s.entity_key, min(s.seed_identifier) AS identifier
  FROM exclusion_metabolite s
  WHERE s.assignment_source = 'metaproviz_seed' AND s.match_kind = 'entity'
    AND list_contains($derived_classes, s.exclusion_class)
    AND NOT EXISTS (
      SELECT 1 FROM exclusion_metabolite d
      WHERE d.entity_key = s.entity_key AND d.exclusion_class = s.exclusion_class
        AND d.assignment_source = 'chebi_hierarchy')
  GROUP BY s.exclusion_class, s.entity_key
)
SELECT 'seed_not_derived',
       c.cls,
       count(m.entity_key)::INTEGER,
       json_object(
         'entities', coalesce(to_json(list(m.entity_key ORDER BY m.entity_key)
                                      FILTER (WHERE m.entity_key IS NOT NULL)), '[]'::JSON),
         'identifiers', coalesce(to_json(list(m.identifier ORDER BY m.entity_key)
                                         FILTER (WHERE m.entity_key IS NOT NULL)), '[]'::JSON)
       )::VARCHAR,
       $build_id
FROM unnest($derived_classes) AS c(cls)
LEFT JOIN missed m ON m.exclusion_class = c.cls
GROUP BY c.cls;
