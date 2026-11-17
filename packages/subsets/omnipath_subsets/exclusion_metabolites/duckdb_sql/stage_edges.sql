-- ChEBI's own subclass_of hierarchy, staged once for the two derivations that
-- walk it (the ion and atom closure and the cofactor role).
--
-- Only ChEBI's own hierarchy: the roots are ChEBI terms, and a release holds
-- far more subclass_of rows from other ontologies.
CREATE OR REPLACE TEMP TABLE exclusion_subclass_edge AS
SELECT DISTINCT subject_entity_key AS child, object_entity_key AS parent
FROM release_relation
WHERE resource = 'chebi'
  AND statement_kind IN ('ontology', 'ontology_axiom')
  AND predicate = 'subclass_of';
