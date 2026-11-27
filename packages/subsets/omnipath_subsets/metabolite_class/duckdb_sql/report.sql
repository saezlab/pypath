-- The figures of the build, one row for each kind and class.

-- The molecules of each class, zero included.
INSERT INTO metabolite_class_report
SELECT 'class_size', c.name, count(m.entity_key)::INTEGER, NULL, $build_id
FROM unnest($classes) AS c(name)
LEFT JOIN metabolite_class m ON m.metabolite_class = c.name
GROUP BY c.name;

-- The molecules each combination of signals claims, before the collapse to one class.
INSERT INTO metabolite_class_report
WITH combination AS (
  SELECT entity_key, list(DISTINCT metabolite_class ORDER BY metabolite_class) AS classes
  FROM mc_claim
  GROUP BY entity_key
)
SELECT 'overlap', NULL, count(*)::INTEGER, json_object('classes', classes)::VARCHAR, $build_id
FROM combination
GROUP BY classes;

-- The connectivity groups that pulled in variants, and the variants, as the exclusion
-- product reports them.
INSERT INTO metabolite_class_report
WITH keyed AS (
  SELECT m.metabolite_class, a.group_key,
         min(m.via_entity_key) AS via,
         list(DISTINCT m.label ORDER BY m.label) FILTER (WHERE m.label IS NOT NULL) AS labels
  FROM metabolite_class m
  JOIN mc_assigned a USING (entity_key)
  WHERE m.match_kind = 'structure_key'
  GROUP BY m.metabolite_class, a.group_key
),
numbered AS (
  SELECT k.*, row_number() OVER (PARTITION BY metabolite_class ORDER BY group_key) AS n
  FROM keyed k
)
SELECT 'shared_key', c.name, count(u.group_key)::INTEGER,
       json_object(
         'keys', coalesce(to_json(list(json_object('key', u.group_key, 'via', u.via,
                                                   'labels', u.labels) ORDER BY u.n)
                                  FILTER (WHERE u.n <= 500)), '[]'::JSON),
         'truncated', count(u.group_key) > 500
       )::VARCHAR,
       $build_id
FROM unnest($classes) AS c(name)
LEFT JOIN numbered u ON u.metabolite_class = c.name
GROUP BY c.name;

-- The molecules a source claims that have no structure key, so no variant can join them:
-- lipid names, native fallbacks and records without a structure.
INSERT INTO metabolite_class_report
WITH unkeyed AS (
  SELECT m.metabolite_class, m.entity_key, coalesce(m.label, m.entity_key) AS label
  FROM metabolite_class m
  JOIN mc_assigned a USING (entity_key)
  WHERE m.match_kind = 'entity' AND a.group_key IS NULL
),
numbered AS (
  SELECT u.*, row_number() OVER (PARTITION BY metabolite_class ORDER BY label, entity_key) AS n
  FROM unkeyed u
)
SELECT 'no_key', c.name, count(u.entity_key)::INTEGER,
       json_object(
         'labels', coalesce(to_json(list(u.label ORDER BY u.n) FILTER (WHERE u.n <= 1000)),
                            '[]'::JSON),
         'truncated', count(u.entity_key) > 1000
       )::VARCHAR,
       $build_id
FROM unnest($classes) AS c(name)
LEFT JOIN numbered u ON u.metabolite_class = c.name
GROUP BY c.name;
