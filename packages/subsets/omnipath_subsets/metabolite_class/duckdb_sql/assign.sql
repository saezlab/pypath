-- The class of each molecule, after the connectivity groups.
--
-- A molecule's own class is its first claim. A connectivity group takes the best own
-- class of its members, and every chemical of the group gets it: the members that a
-- source claims and the ones it does not. A molecule without a structure key (a lipid
-- name, a native fallback) keeps its own class. No entity is merged with another.
--
-- Each step is a table of its own. As one statement of common table expressions, DuckDB
-- planned it into a join that ran past the memory limit on a full release; as steps, each
-- takes a second.

-- The best claim of each claimed molecule.
CREATE TEMP TABLE mc_own AS
SELECT c.entity_key, min(k.precedence) AS precedence
FROM mc_claim c
JOIN mc_class k ON k.name = c.metabolite_class
GROUP BY c.entity_key;

-- The structure key of each chemical, as release_entity gives it (the first resource in
-- name order that has one), read from the chemicals' rows only: the merged view would
-- aggregate every entity of the release.
CREATE TEMP TABLE mc_grouped AS
SELECT entity_key, arg_min(group_connectivity, resource) AS group_key
FROM resource_entity
WHERE group_connectivity IS NOT NULL
  AND list_contains($entity_types, entity_type)
GROUP BY entity_key;

-- The best class of each group that has a claimed member, and the lowest key of the
-- members that give it.
CREATE TEMP TABLE mc_best AS
SELECT g.group_key,
       min(o.precedence) AS precedence,
       arg_min(o.entity_key, (o.precedence, o.entity_key)) AS via_entity_key
FROM mc_own o
JOIN mc_grouped g USING (entity_key)
GROUP BY g.group_key;

-- Every claimed molecule, and every chemical of a group with a claimed member.
CREATE TEMP TABLE mc_assigned AS
WITH member AS (
  SELECT entity_key FROM mc_own
  UNION
  SELECT g.entity_key FROM mc_grouped g JOIN mc_best b USING (group_key)
),
placed AS (
  SELECT m.entity_key,
         g.group_key,
         o.precedence AS own_precedence,
         coalesce(b.precedence, o.precedence) AS precedence,
         b.via_entity_key
  FROM member m
  LEFT JOIN mc_own o USING (entity_key)
  LEFT JOIN mc_grouped g USING (entity_key)
  LEFT JOIN mc_best b ON b.group_key = g.group_key
)
SELECT p.entity_key,
       p.group_key,
       k.name AS metabolite_class,
       CASE WHEN p.own_precedence = p.precedence THEN 'entity' ELSE 'structure_key' END
         AS match_kind,
       CASE WHEN p.own_precedence = p.precedence THEN NULL ELSE p.via_entity_key END
         AS via_entity_key
FROM placed p
JOIN mc_class k ON k.precedence = p.precedence;
