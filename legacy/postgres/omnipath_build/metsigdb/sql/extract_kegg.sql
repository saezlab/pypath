-- KEGG membership extraction: two hops, through the reaction.
--
-- KEGG publishes no compound-to-pathway edge. A pathway is `associated_with` a
-- reaction, and the reaction `has_participant` a compound, so a compound
-- reaches a pathway only through a reaction assigned to it.
--
-- The reaction is an internal join key and is never published, so this
-- extraction does not require it to resolve canonically. It requires the
-- compound end to resolve, because that is the end the substrate publishes.
-- The difference is large: requiring both ends gives 1,272 memberships, and
-- requiring the compound alone gives 4,969.
--
-- Seven in ten KEGG compound participations stay unresolved in the core build,
-- so KEGG publishes far fewer memberships than its raw edge count suggests.
-- That is an upstream gap, recorded in contracts/mapping-rules.md.
--
-- This file is interpolated with bind parameters, so it carries no bare percent
-- sign, not even inside a comment.

DROP TABLE IF EXISTS metsigdb_kegg_pathway_reaction;
DROP TABLE IF EXISTS metsigdb_kegg_reaction_compound;
DROP TABLE IF EXISTS metsigdb_stage;

-- Both hops read `entity_evidence_resolution` directly. The table is
-- partitioned by source, so one partition is scanned, and staging it first
-- would cost a 964,161-row temporary copy for no gain.

-- Hop one: pathway to reaction. The reaction end takes any resolution status,
-- because the reaction is an internal join key that is never published.
CREATE TEMP TABLE metsigdb_kegg_pathway_reaction AS
SELECT pw.canonical_identifier AS set_source_id,
       pw.entity_id            AS set_entity_id,
       rx.entity_id            AS reaction_entity_id
FROM relation_evidence re
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = re.predicate_id
 AND p.name = 'associated_with'
JOIN entity pw
  ON pw.entity_id = re.subject_entity_id
 AND pw.entity_type_id = %(set_entity_type_id)s
JOIN entity_evidence_resolution rx
  ON rx.source_id = re.source_id
 AND rx.entity_evidence_id = re.object_entity_evidence_id
 AND rx.entity_id IS NOT NULL
WHERE re.source_id = %(source_id)s;

CREATE INDEX ON metsigdb_kegg_pathway_reaction (reaction_entity_id);
ANALYZE metsigdb_kegg_pathway_reaction;

-- Hop two: reaction to compound. The compound end must be resolved.
CREATE TEMP TABLE metsigdb_kegg_reaction_compound AS
SELECT rx.entity_id  AS reaction_entity_id,
       cpd.entity_id AS metabolite_entity_id,
       re.dataset_id,
       re.row_id,
       re.source_id
FROM relation_evidence re
JOIN vocab_relation_predicate p
  ON p.relation_predicate_id = re.predicate_id
 AND p.name = 'has_participant'
JOIN entity_evidence_resolution rx
  ON rx.source_id = re.source_id
 AND rx.entity_evidence_id = re.subject_entity_evidence_id
 AND rx.entity_id IS NOT NULL
JOIN entity_evidence_resolution cpd
  ON cpd.source_id = re.source_id
 AND cpd.entity_evidence_id = re.object_entity_evidence_id
 AND cpd.status_id = 1
JOIN entity met
  ON met.entity_id = cpd.entity_id
 AND met.entity_type_id = %(chemical_entity_type_id)s
WHERE re.source_id = %(source_id)s;

CREATE INDEX ON metsigdb_kegg_reaction_compound (reaction_entity_id);
ANALYZE metsigdb_kegg_reaction_compound;

-- One row per pathway and compound. The provenance names the reaction the
-- compound travelled through, so the second hop stays visible in the data.
CREATE TEMP TABLE metsigdb_stage AS
SELECT DISTINCT ON (pr.set_source_id, rc.metabolite_entity_id)
       pr.set_source_id,
       pr.set_entity_id,
       rc.metabolite_entity_id,
       jsonb_build_object(
         'source_id',    rc.source_id,
         'dataset_id',   rc.dataset_id,
         'row_id',       rc.row_id,
         'via_reaction', rc.reaction_entity_id
       ) AS provenance_record,
       -- KEGG's global and overview maps are whole-metabolism diagrams, not
       -- pathways: eleven of them hold half of KEGG's memberships. The list is
       -- KEGG's own, declared in mapping.py, never guessed from a label.
       CASE WHEN pr.set_source_id = ANY(%(overview_maps)s)
            THEN 'overview_map' ELSE 'metabolic_map' END AS set_sub_type,
       NULL::jsonb AS set_context
FROM metsigdb_kegg_pathway_reaction pr
JOIN metsigdb_kegg_reaction_compound rc
  ON rc.reaction_entity_id = pr.reaction_entity_id
ORDER BY pr.set_source_id, rc.metabolite_entity_id, rc.row_id
LIMIT %(max_records)s;

DROP TABLE metsigdb_kegg_pathway_reaction;
DROP TABLE metsigdb_kegg_reaction_compound;
