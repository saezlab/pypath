-- One-hop extraction with typed Biolink orientation and the published-chemical
-- eligibility; no resolver outcome is invented. Reactome and WikiPathways read it, and
-- MACdb with the trait type of its set.
--
-- The set is the endpoint of the set entity type, the metabolite the chemical one. A
-- pathway lists its members under `has_part`, `has_member` or `has_participant` and a
-- member lists its pathway under the inverse; every resource uses `associated_with`.

CREATE OR REPLACE TEMP TABLE metsigdb_stage AS
WITH eligible AS (
  SELECT r.*, se.entity_type AS subject_type, oe.entity_type AS object_type
  FROM release_relation_evidence r
  JOIN release_entity se ON se.entity_key = r.subject_entity_key
  JOIN release_entity oe ON oe.entity_key = r.object_entity_key
  WHERE r.resource = $source
    AND r.statement_kind NOT IN ('ontology', 'ontology_axiom')
), pair AS (
  SELECT e.*,
    CASE WHEN subject_type = $set_entity_type
         THEN subject_entity_key ELSE object_entity_key END AS set_entity_key,
    CASE WHEN subject_type = $set_entity_type
         THEN object_entity_key ELSE subject_entity_key END AS metabolite_entity_key
  FROM eligible e
  WHERE (
    subject_type = $set_entity_type
    AND list_contains($chemical_entity_types, object_type)
    AND (predicate = 'associated_with'
      OR ($is_pathway AND predicate IN ('has_part', 'has_member', 'has_participant')))
  ) OR (
    object_type = $set_entity_type
    AND list_contains($chemical_entity_types, subject_type)
    AND (predicate = 'associated_with'
      OR ($is_pathway AND predicate IN ('part_of', 'member_of', 'participates_in')))
  )
)
SELECT DISTINCT ON (st.identifier, pair.metabolite_entity_key)
       st.identifier AS set_source_id,
       st.entity_key AS set_entity_key,
       pair.metabolite_entity_key,
       json_object('source', pair.resource, 'dataset', pair.dataset,
                   'row_id', pair.row_id)::VARCHAR AS provenance_record,
       CASE WHEN $trait_types THEN tt.trait_type END AS set_sub_type,
       NULL::VARCHAR AS set_context
FROM pair
JOIN release_entity st ON st.entity_key = pair.set_entity_key
LEFT JOIN (
  SELECT entity_key, min(value) AS trait_type
  FROM release_entity_annotation
  WHERE resource = $source AND term = 'macdb:trait_type'
  GROUP BY entity_key
) tt ON tt.entity_key = st.entity_key
ORDER BY st.identifier, pair.metabolite_entity_key, pair.row_id;
