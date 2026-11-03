-- A veto removes the entity from its class, whatever added it: a seed, the
-- ontology, an addition or a structure-key match. What it removes is recorded
-- so the report can name it.
--
-- It runs twice. Before the structure-key expansion, so a vetoed member does not
-- pull its variants in. After it, so a vetoed entity that a different member
-- would pull in as a variant stays out.

WITH gone AS (
  DELETE FROM exclusion_metabolite x
  USING exclusion_veto v
  WHERE x.entity_id = v.entity_id AND x.exclusion_class = v.exclusion_class
  RETURNING v.exclusion_class AS veto_class, v.namespace, v.identifier, x.entity_id,
            x.assignment_source, x.match_kind
)
INSERT INTO exclusion_veto_removed
SELECT veto_class, namespace, identifier, entity_id, assignment_source, match_kind FROM gone;
