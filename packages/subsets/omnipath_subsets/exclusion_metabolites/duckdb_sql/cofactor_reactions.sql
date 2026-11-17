-- How many reactions each listed cofactor takes part in. A review aid: a person
-- reading the class sees which members are everywhere (ATP) and which barely
-- appear. It never adds a member; the count is read from the members already
-- listed.
--
-- A reaction is an entity of its own (a molecular activity) with one
-- participant edge to each member.

INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, resolved, detail, build_id)
WITH counted AS (
  SELECT m.entity_key,
         coalesce(e.label, e.identifier) AS label,
         count(DISTINCT p.subject_entity_key) AS reactions
  FROM (SELECT DISTINCT entity_key FROM exclusion_metabolite
        WHERE exclusion_class = 'cofactors') m
  JOIN release_entity e ON e.entity_key = m.entity_key
  JOIN release_relation p ON p.object_entity_key = m.entity_key
   AND p.predicate IN ('has_participant', 'has_input', 'has_output')
  JOIN release_entity reaction
    ON reaction.entity_key = p.subject_entity_key AND reaction.entity_type = 'molecular_activity'
  GROUP BY m.entity_key, coalesce(e.label, e.identifier)
),
numbered AS (
  SELECT c.*, row_number() OVER (ORDER BY reactions DESC, label, entity_key) AS n FROM counted c
)
SELECT 'cofactor_reactions',
       'cofactors',
       count(*)::INTEGER,
       json_object(
         'cofactors', coalesce(
           to_json(list(json_object('entity', entity_key, 'label', label, 'reactions', reactions)
                        ORDER BY n) FILTER (WHERE n <= 1000)), '[]'::JSON),
         'truncated', count(*) > 1000
       )::VARCHAR,
       $build_id
FROM numbered
HAVING count(*) > 0;
