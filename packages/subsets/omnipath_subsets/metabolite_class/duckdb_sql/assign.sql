-- The class of each molecule, after the connectivity groups.
--
-- A molecule's own class is its first claim. A connectivity group takes the best own
-- class of its members, and every chemical of the group gets it: the members that a
-- source claims and the ones it does not. A molecule without a structure key (a lipid
-- name, a native fallback) keeps its own class. No entity is merged with another.
CREATE TEMP TABLE mc_assigned AS
WITH own AS (
  SELECT c.entity_key, min(k.precedence) AS precedence
  FROM mc_claim c
  JOIN mc_class k ON k.name = c.metabolite_class
  GROUP BY c.entity_key
),
grouped AS (
  SELECT entity_key, group_connectivity AS group_key
  FROM release_entity
  WHERE group_connectivity IS NOT NULL
    AND list_contains($entity_types, entity_type)
),
best AS (
  SELECT g.group_key, min(o.precedence) AS precedence
  FROM own o
  JOIN grouped g USING (entity_key)
  GROUP BY g.group_key
),
via AS (
  SELECT g.group_key, min(o.entity_key) AS via_entity_key
  FROM own o
  JOIN grouped g USING (entity_key)
  JOIN best b ON b.group_key = g.group_key AND b.precedence = o.precedence
  GROUP BY g.group_key
),
member AS (
  SELECT entity_key FROM own
  UNION
  SELECT g.entity_key FROM grouped g JOIN best b USING (group_key)
)
SELECT m.entity_key,
       g.group_key,
       k.name AS metabolite_class,
       CASE WHEN o.precedence = coalesce(b.precedence, o.precedence) THEN 'entity'
            ELSE 'structure_key' END AS match_kind,
       CASE WHEN o.precedence = coalesce(b.precedence, o.precedence) THEN NULL
            ELSE v.via_entity_key END AS via_entity_key
FROM member m
LEFT JOIN own o USING (entity_key)
LEFT JOIN grouped g USING (entity_key)
LEFT JOIN best b ON b.group_key = g.group_key
LEFT JOIN via v ON v.group_key = g.group_key
JOIN mc_class k ON k.precedence = coalesce(b.precedence, o.precedence);
