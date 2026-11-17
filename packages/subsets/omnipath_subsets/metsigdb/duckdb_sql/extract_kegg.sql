-- KEGG membership extraction: two hops, through the reaction.
--
-- KEGG publishes no compound-to-pathway edge. A pathway is `associated_with` a
-- reaction, and the reaction has an input or an output compound, so a compound
-- reaches a pathway only through a reaction assigned to it.
--
-- The reaction is an internal join key and is never published. It is read by type, and
-- only the compound end must be a chemical.

-- Hop one: pathway to reaction.
CREATE OR REPLACE TEMP TABLE metsigdb_kegg_pathway_reaction AS
SELECT pw.identifier  AS set_source_id,
       pw.entity_key  AS set_entity_key,
       rx.entity_key  AS reaction_entity_key
FROM release_relation_evidence re
JOIN release_entity pw ON pw.entity_key IN (re.subject_entity_key, re.object_entity_key)
 AND pw.entity_type = $set_entity_type
JOIN release_entity rx ON rx.entity_key = CASE WHEN pw.entity_key = re.subject_entity_key
                                              THEN re.object_entity_key
                                              ELSE re.subject_entity_key END
 AND rx.entity_type = 'molecular_activity'
WHERE re.resource = $source
  AND re.statement_kind NOT IN ('ontology', 'ontology_axiom')
  AND (
    (pw.entity_key = re.subject_entity_key AND re.predicate IN
      ('associated_with', 'has_part', 'has_member', 'has_participant')) OR
    (pw.entity_key = re.object_entity_key AND re.predicate IN
      ('associated_with', 'part_of', 'member_of', 'participates_in')));

-- Hop two: reaction to compound. The compound end must be a chemical.
CREATE OR REPLACE TEMP TABLE metsigdb_kegg_reaction_compound AS
SELECT rx.entity_key  AS reaction_entity_key,
       cpd.entity_key AS metabolite_entity_key,
       re.resource, re.dataset, re.row_id
FROM release_relation_evidence re
JOIN release_entity rx ON rx.entity_key = re.subject_entity_key
 AND rx.entity_type = 'molecular_activity'
JOIN release_entity cpd ON cpd.entity_key = re.object_entity_key
 AND list_contains($chemical_entity_types, cpd.entity_type)
WHERE re.resource = $source
  AND re.statement_kind NOT IN ('ontology', 'ontology_axiom')
  AND re.predicate IN ('has_input', 'has_output');

-- One row per pathway and compound. The provenance names the reaction the compound
-- travelled through, so the second hop stays visible in the data.
CREATE OR REPLACE TEMP TABLE metsigdb_stage AS
SELECT DISTINCT ON (pr.set_source_id, rc.metabolite_entity_key)
       pr.set_source_id,
       pr.set_entity_key,
       rc.metabolite_entity_key,
       json_object('source', rc.resource, 'dataset', rc.dataset, 'row_id', rc.row_id,
                   'via_reaction', rc.reaction_entity_key)::VARCHAR AS provenance_record,
       -- KEGG's global and overview maps are whole-metabolism diagrams, not pathways. The
       -- list is KEGG's own, declared in mapping.py, never guessed from a label. A map is
       -- named by its number: the reference map rn01100 and an organism's hsa01100 are
       -- the same diagram.
       CASE WHEN list_contains($overview_map_numbers,
                               regexp_extract(pr.set_source_id, '^[a-z]+([0-9]{5})$', 1))
            THEN 'overview_map' ELSE 'metabolic_map' END AS set_sub_type,
       NULL::VARCHAR AS set_context
FROM metsigdb_kegg_pathway_reaction pr
JOIN metsigdb_kegg_reaction_compound rc
  ON rc.reaction_entity_key = pr.reaction_entity_key
ORDER BY pr.set_source_id, rc.metabolite_entity_key, rc.row_id;

DROP TABLE metsigdb_kegg_pathway_reaction;
DROP TABLE metsigdb_kegg_reaction_compound;
