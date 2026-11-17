-- Closure members the curated seeds lack.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH extra AS (
  SELECT d.exclusion_class, d.entity_key,
         coalesce(e.label, e.identifier) AS label,
         row_number() OVER (PARTITION BY d.exclusion_class
                            ORDER BY coalesce(e.label, e.identifier), d.entity_key) AS n
  FROM exclusion_metabolite d
  JOIN release_entity e ON e.entity_key = d.entity_key
  WHERE d.assignment_source = 'chebi_hierarchy' AND d.match_kind = 'entity'
    AND NOT EXISTS (
      SELECT 1 FROM exclusion_metabolite s
      WHERE s.entity_key = d.entity_key AND s.exclusion_class = d.exclusion_class
        AND s.origin = 'curated' AND s.match_kind = 'entity')
)
SELECT 'derived_not_seeded',
       c.cls,
       count(x.entity_key)::INTEGER,
       json_object(
         'entities', coalesce(to_json(list(x.entity_key ORDER BY x.n)
                                      FILTER (WHERE x.n <= 1000)), '[]'::JSON),
         'labels', coalesce(to_json(list(x.label ORDER BY x.n)
                                    FILTER (WHERE x.n <= 1000)), '[]'::JSON),
         'truncated', count(x.entity_key) > 1000
       )::VARCHAR,
       $build_id
FROM unnest($derived_classes) AS c(cls)
LEFT JOIN extra x ON x.exclusion_class = c.cls
GROUP BY c.cls;
