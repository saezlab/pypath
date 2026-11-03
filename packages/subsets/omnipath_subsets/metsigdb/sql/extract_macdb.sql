-- Main MACdb set/chemical publication, with explicit source trait attribute.
DROP TABLE IF EXISTS metsigdb_stage;
CREATE TEMP TABLE metsigdb_stage AS
WITH eligible AS (
  SELECT re.*, se.entity_type_id AS subject_type, oe.entity_type_id AS object_type,
         p.name AS predicate
  FROM relation_evidence re
  JOIN vocab_relation_predicate p ON p.relation_predicate_id = re.predicate_id
  JOIN entity se ON se.entity_id = re.subject_entity_id
  JOIN entity oe ON oe.entity_id = re.object_entity_id
  WHERE re.source_id = %(source_id)s
), pair AS (
  SELECT e.*,
    CASE WHEN subject_type = %(set_entity_type_id)s
         THEN subject_entity_id ELSE object_entity_id END AS set_entity_id,
    CASE WHEN subject_type = %(set_entity_type_id)s
         THEN object_entity_id ELSE subject_entity_id END AS metabolite_entity_id
  FROM eligible e
  WHERE (
    subject_type = %(set_entity_type_id)s
    AND object_type = ANY(%(chemical_entity_type_ids)s)
    AND (predicate = 'associated_with'
      OR (%(is_pathway)s AND predicate IN ('has_part','has_member','has_participant')))
  ) OR (
    object_type = %(set_entity_type_id)s
    AND subject_type = ANY(%(chemical_entity_type_ids)s)
    AND (predicate = 'associated_with'
      OR (%(is_pathway)s AND predicate IN ('part_of','member_of','participates_in')))
  )
)
SELECT DISTINCT ON (st.canonical_identifier, pair.metabolite_entity_id)
       st.canonical_identifier AS set_source_id, st.entity_id AS set_entity_id,
       pair.metabolite_entity_id,
       jsonb_build_object('source_id',pair.source_id,
                         'dataset_id',pair.dataset_id,'row_id',pair.row_id)
         AS provenance_record,
       tt.trait_type AS set_sub_type, NULL::jsonb AS set_context
FROM pair JOIN entity st ON st.entity_id = pair.set_entity_id
LEFT JOIN (
  SELECT er.entity_id, min(a.value) AS trait_type
  FROM entity_evidence_resolution er
  JOIN entity_evidence_annotation ea
    ON ea.source_id=er.source_id AND ea.entity_evidence_id=er.entity_evidence_id
  JOIN annotation a ON a.annotation_key=ea.annotation_key
  WHERE er.source_id=%(source_id)s AND a.term='macdb:trait_type'
  GROUP BY er.entity_id
) tt ON tt.entity_id=st.entity_id
ORDER BY st.canonical_identifier, pair.metabolite_entity_id, pair.row_id
LIMIT %(max_records)s;
