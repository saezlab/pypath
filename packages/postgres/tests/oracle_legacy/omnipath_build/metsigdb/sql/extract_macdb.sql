-- MACdb membership extraction: one hop, with both ends on the evidence side.
--
-- MACdb is the only resource whose set does not sit on the entity side of the
-- relation. The trait is an evidence, so it needs the same resolution join as
-- the metabolite. Both ends must resolve, because both are published: the
-- trait becomes the set identity, and the metabolite becomes the anchor.
--
-- MACdb is also the only resource contributing the `disease` semantic. It
-- publishes metabolite-to-cancer-trait associations, and its traits resolve to
-- controlled-vocabulary terms.

DROP TABLE IF EXISTS metsigdb_stage;

CREATE TEMP TABLE metsigdb_stage AS
SELECT DISTINCT ON (st.canonical_identifier, cpd.entity_id)
       st.canonical_identifier AS set_source_id,
       st.entity_id            AS set_entity_id,
       cpd.entity_id           AS metabolite_entity_id,
       jsonb_build_object(
         'source_id',  re.source_id,
         'dataset_id', re.dataset_id,
         'row_id',     re.row_id
       ) AS provenance_record,
       -- MACdb records a type for every trait, and pypath models it as a
       -- `part_of` relation. Only 153 of the 269 traits are diseases: the rest
       -- are interventions, phenotypes, genotypes and gene abnormalities, all
       -- carrying `set_type = disease` because a resource contributes one
       -- semantic. The sub-type is where the difference lives, and it comes
       -- from the source, not from reading the label.
       tt.trait_type AS set_sub_type,
       NULL::jsonb AS set_context
FROM relation_evidence re
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = re.predicate_id
 AND p.name = 'associated_with'
JOIN entity_evidence_resolution trait
  ON trait.source_id = re.source_id
 AND trait.entity_evidence_id = re.subject_entity_evidence_id
 AND trait.status_id = 1
JOIN entity st
  ON st.entity_id = trait.entity_id
 AND st.entity_type_id = %(set_entity_type_id)s
JOIN entity_evidence_resolution cpd
  ON cpd.source_id = re.source_id
 AND cpd.entity_evidence_id = re.object_entity_evidence_id
 AND cpd.status_id = 1
JOIN entity met
  ON met.entity_id = cpd.entity_id
 AND met.entity_type_id = %(chemical_entity_type_id)s
LEFT JOIN (
  SELECT tr.subject_entity_id AS entity_id,
         min(o.canonical_identifier) AS trait_type
  FROM entity_ontology_relation tr
  JOIN vocab_relation_predicate tp
    ON tp.relation_predicate_id = tr.predicate_id AND tp.name = 'part_of'
  JOIN entity o ON o.entity_id = tr.object_entity_id
  WHERE tr.source_id = %(source_id)s
  GROUP BY 1
) tt ON tt.entity_id = trait.entity_id
WHERE re.source_id = %(source_id)s
ORDER BY st.canonical_identifier, cpd.entity_id, re.row_id
LIMIT %(max_records)s;
