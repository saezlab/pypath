-- A veto removes the entity from its class, whatever added it: a seed, the
-- ontology, an addition or a structure-key match. What it removes is recorded
-- so the report can name it.
--
-- It runs twice. Before the structure-key expansion, so a vetoed member does not
-- pull its variants in. After it, so a vetoed entity that a different member
-- would pull in as a variant stays out.
INSERT INTO exclusion_veto_removed
SELECT v.exclusion_class, v.namespace, v.identifier, x.entity_key, x.assignment_source,
       x.match_kind
FROM exclusion_metabolite x
JOIN exclusion_veto v ON x.entity_key = v.entity_key AND x.exclusion_class = v.exclusion_class;
