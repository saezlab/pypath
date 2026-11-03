-- ChEBI's own subclass_of hierarchy, staged once for the two derivations that
-- walk it (the ion and atom closure and the cofactor role).
--
DROP TABLE IF EXISTS exclusion_subclass_edge;

-- Only ChEBI's own hierarchy: the roots are ChEBI terms, and a release holds
-- about six times as many subclass_of rows from other ontologies.
CREATE TEMP TABLE exclusion_subclass_edge AS
SELECT DISTINCT r.subject_entity_id AS child, r.object_entity_id AS parent
FROM entity_ontology_relation r
JOIN data_source s ON s.source_id = r.source_id AND s.name = 'chebi'
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = r.predicate_id AND p.name = 'subclass_of';

CREATE INDEX ON exclusion_subclass_edge (parent);
ANALYZE exclusion_subclass_edge;
