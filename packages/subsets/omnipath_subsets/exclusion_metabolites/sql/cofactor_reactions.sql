-- How many reactions each listed cofactor takes part in. A review aid: a person
-- reading the class sees which members are everywhere (ATP) and which barely
-- appear. It never adds a member; the count is read from the members already
-- listed.

INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH counted AS (
  SELECT m.entity_id,
         coalesce(e.label, e.canonical_identifier) AS label,
         count(DISTINCT p.interaction_id) AS reactions
  FROM (SELECT DISTINCT entity_id FROM exclusion_metabolite
        WHERE exclusion_class = 'cofactors') m
  JOIN entity e ON e.entity_id = m.entity_id
  JOIN interaction_party p ON p.entity_id = m.entity_id
  GROUP BY m.entity_id, coalesce(e.label, e.canonical_identifier)
),
numbered AS (
  SELECT c.*, row_number() OVER (ORDER BY reactions DESC, label, entity_id) AS n FROM counted c
)
SELECT 'cofactor_reactions',
       'cofactors',
       count(*),
       jsonb_build_object(
         'cofactors', coalesce(
           jsonb_agg(jsonb_build_object('entity', entity_id, 'label', label, 'reactions', reactions)
                     ORDER BY n) FILTER (WHERE n <= 1000), '[]'::jsonb),
         'truncated', count(*) > 1000
       ),
       %(build_id)s
FROM numbered
HAVING count(*) > 0;
