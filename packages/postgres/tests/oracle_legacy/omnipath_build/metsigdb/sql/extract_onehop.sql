-- One-hop membership extraction: a set entity is `associated_with` a chemical.
--
-- Reactome, WikiPathways and ClassyFire's direct assignments all have this
-- shape. The set sits on the entity side of the relation, the metabolite on the
-- evidence side, and one join through `entity_evidence_resolution` turns the
-- evidence into a canonical entity.
--
-- Only `status_id = 1` publishes. Status 2 carries a placeholder entity minted
-- for an unresolved key, not a canonical one, and the row contract excludes it.
--
-- The source names the same metabolite in a set several times, so DISTINCT ON
-- folds the evidence rows to one row per pair and keeps the lowest-numbered
-- evidence as the provenance record.

DROP TABLE IF EXISTS metsigdb_stage;

CREATE TEMP TABLE metsigdb_stage AS
SELECT DISTINCT ON (st.canonical_identifier, res.entity_id)
       st.canonical_identifier AS set_source_id,
       st.entity_id            AS set_entity_id,
       res.entity_id           AS metabolite_entity_id,
       jsonb_build_object(
         'source_id',  re.source_id,
         'dataset_id', re.dataset_id,
         'row_id',     re.row_id
       ) AS provenance_record,
       -- Every extraction produces the same stage shape, so one publication
       -- path serves all five resources. Only ClassyFire fills this column.
       NULL::text  AS set_sub_type,
       NULL::jsonb AS set_context
FROM relation_evidence re
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = re.predicate_id
 AND p.name = 'associated_with'
JOIN entity st
  ON st.entity_id = re.subject_entity_id
 AND st.entity_type_id = %(set_entity_type_id)s
JOIN entity_evidence_resolution res
  ON res.source_id = re.source_id
 AND res.entity_evidence_id = re.object_entity_evidence_id
 AND res.status_id = 1
JOIN entity met
  ON met.entity_id = res.entity_id
 AND met.entity_type_id = %(chemical_entity_type_id)s
WHERE re.source_id = %(source_id)s
ORDER BY st.canonical_identifier, res.entity_id, re.row_id
LIMIT %(max_records)s;
