-- Keys that pulled in variants, so a person can see what a key match added and
-- remove a wanted metabolite with a veto.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH keyed AS (
  SELECT x.exclusion_class, g.group_key,
         min(x.via_entity_key) AS via,
         list(DISTINCT coalesce(e.label, e.identifier) ORDER BY coalesce(e.label, e.identifier))
           AS labels
  FROM exclusion_metabolite x
  JOIN structure_group g ON g.entity_key = x.entity_key
  JOIN release_entity e ON e.entity_key = x.entity_key
  WHERE x.match_kind = 'structure_key'
  GROUP BY x.exclusion_class, g.group_key
),
numbered AS (
  SELECT k.*, row_number() OVER (PARTITION BY exclusion_class ORDER BY group_key) AS n
  FROM keyed k
)
SELECT 'shared_key',
       exclusion_class,
       count(*)::INTEGER,
       json_object(
         'keys', to_json(list(json_object('key', group_key, 'via', via, 'labels', labels)
                              ORDER BY n) FILTER (WHERE n <= 500)),
         'truncated', count(*) > 500
       )::VARCHAR,
       $build_id
FROM numbered
GROUP BY exclusion_class;
