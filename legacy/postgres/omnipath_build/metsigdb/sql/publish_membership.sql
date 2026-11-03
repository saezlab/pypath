-- Project the staged pairs onto the canonical layer, on both sides.
--
-- Every resource reaches the table through this file. The extraction differs
-- per resource; the publication does not.
--
-- Two temporary tables come out of it. `metsigdb_projection` carries the
-- metabolite side, `metsigdb_set_attrs` the set side.

DROP TABLE IF EXISTS metsigdb_projection;
DROP TABLE IF EXISTS metsigdb_set_attrs;

CREATE TEMP TABLE metsigdb_projection AS
SELECT e.entity_id,
       -- The row contract requires a label. Every metabolite has one today,
       -- and the fallback keeps the column honest if that ever changes.
       coalesce(e.label, e.canonical_identifier) AS metabolite_label,
       vt.name                                   AS metabolite_entity_type,
       -- Connectivity-level group key: one value per molecular skeleton, so
       -- charge, stereo and tautomer variants of one molecule agree on it.
       max(g.group_key)                          AS metabolite_structure_key,
       max(v.value) FILTER (WHERE it.name = 'Standard Inchi Key:MI:1101') AS inchikey,
       max(v.value) FILTER (WHERE it.name = 'Smiles:MI:0239')             AS smiles,
       -- HMDB ships two accession widths, HMDB00008 and HMDB0000001. Both
       -- reach the build, and a consumer joining on the column silently
       -- misses a tenth of it. The seven-digit form is the current one.
       max(
         CASE WHEN v.value ~ '^HMDB[0-9]+$'
              THEN 'HMDB' || lpad(substring(v.value FROM 5), 7, '0')
              ELSE v.value
         END
       ) FILTER (WHERE it.name = 'Hmdb:OM:0004')                          AS hmdb,
       max(v.value) FILTER (WHERE it.name = 'Pubchem Compound:OM:0002')   AS pubchem,
       max(v.value) FILTER (WHERE it.name = 'Chebi:MI:0474')              AS chebi,
       max(v.value) FILTER (WHERE it.name = 'Kegg Compound:MI:2012')      AS kegg
FROM (SELECT DISTINCT metabolite_entity_id FROM metsigdb_stage) s
JOIN entity e ON e.entity_id = s.metabolite_entity_id
JOIN vocab_entity_type vt ON vt.entity_type_id = e.entity_type_id
LEFT JOIN chemical_resolution_group_member g
  ON g.entity_id = e.entity_id AND g.level_id = 1
LEFT JOIN entity_identifier_lookup l ON l.entity_id = e.entity_id
LEFT JOIN identifier_evidence v ON v.identifier_id = l.identifier_id
LEFT JOIN vocab_identifier_type it ON it.identifier_type_id = v.identifier_type_id
GROUP BY e.entity_id, coalesce(e.label, e.canonical_identifier), vt.name;

CREATE INDEX ON metsigdb_projection (entity_id);

ANALYZE metsigdb_projection;

-- The set side: a readable name, and the species the source recorded.
--
-- The name joins on the set entity, never on the source string. MACdb trait
-- ids are bare integers and collide with ChEBI ids in
-- entity.canonical_identifier, so a string join returns other resources' sets.
--
-- The species is read from the taxonomy the source put on its set evidence.
-- No set entity in this build carries a taxonomy_id of its own, and the
-- identifier is not a reliable substitute: WP ids name no species, and
-- deriving one from an R-HSA prefix guesses where the data already answers.
-- Measured across all 4,569 set entities, no entity's evidences disagree.
--
-- The name has two sources and the ontology term wins. A curated term carries a
-- label; a set that is not an ontology term at all carries its name as a `Name`
-- identifier instead, which is how the source published it. WikiPathways
-- pathways are the second case: 745 of 818 carry `Name`, none is an ontology
-- term, and reading only the term published null for every one of them.
--
-- `starts_with` rather than `LIKE`, because psycopg2 interpolates this file and
-- a per-cent sign in it — escaped or not — is the bug that failed the cycle 010
-- load four times.
CREATE TEMP TABLE metsigdb_set_attrs AS
SELECT s.set_entity_id AS entity_id,
       coalesce(max(t.label), max(nm.name)) AS set_label,
       min(tax.taxonomy_id) AS organism
FROM (SELECT DISTINCT set_entity_id FROM metsigdb_stage) s
LEFT JOIN entity_ontology_term t ON t.term_entity_id = s.set_entity_id
LEFT JOIN (
  SELECT l.entity_id, min(ie.value) AS name
  FROM (SELECT DISTINCT set_entity_id FROM metsigdb_stage) staged
  JOIN entity_identifier_lookup l ON l.entity_id = staged.set_entity_id
  JOIN identifier_evidence ie ON ie.identifier_id = l.identifier_id
  JOIN vocab_identifier_type it ON it.identifier_type_id = ie.identifier_type_id
  WHERE starts_with(it.name, 'Name:')
  GROUP BY l.entity_id
) nm ON nm.entity_id = s.set_entity_id
LEFT JOIN (
  SELECT r.entity_id, ee.taxonomy_id
  FROM entity_evidence ee
  JOIN entity_evidence_resolution r
    ON r.source_id = ee.source_id
   AND r.entity_evidence_id = ee.entity_evidence_id
   AND r.status_id = 1
  WHERE ee.source_id = %(source_id)s AND ee.taxonomy_id IS NOT NULL
) tax ON tax.entity_id = s.set_entity_id
GROUP BY s.set_entity_id;

CREATE INDEX ON metsigdb_set_attrs (entity_id);

ANALYZE metsigdb_set_attrs;
