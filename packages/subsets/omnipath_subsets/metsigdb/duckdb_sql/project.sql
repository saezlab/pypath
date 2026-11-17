-- Project the staged pairs onto the release's entities, on both sides.
--
-- Every resource reaches the table through this file. The extraction differs per
-- resource; the publication does not. Two temporary tables come out of it:
-- `metsigdb_projection` carries the metabolite side, `metsigdb_set_attrs` the set side.

CREATE OR REPLACE TEMP TABLE metsigdb_projection AS
SELECT e.entity_key,
       coalesce(e.label, e.identifier)  AS metabolite_label,
       e.entity_type                    AS metabolite_entity_type,
       -- Connectivity-level group key: one value per molecular skeleton, so charge, stereo
       -- and tautomer variants of one molecule agree on it.
       any_value(e.group_connectivity)  AS metabolite_structure_key,
       max(i.id) FILTER (WHERE i.ns = 'inchikey') AS inchikey,
       max(i.id) FILTER (WHERE i.ns = 'smiles')   AS smiles,
       -- HMDB ships two accession widths, HMDB00008 and HMDB0000001, and sources spell the
       -- prefix in more than one way. A consumer joining on the column silently misses what
       -- differs. The seven-digit form is the current one.
       max(CASE WHEN regexp_matches(i.id, '^\s*(?i:HMDB)[\s:_-]*[0-9]+\s*$')
                THEN 'HMDB' || lpad(ltrim(regexp_extract(i.id, '([0-9]+)\s*$', 1), '0'), 7, '0')
                ELSE i.id END) FILTER (WHERE i.ns = 'hmdb') AS hmdb,
       max(i.id) FILTER (WHERE i.ns = 'pubchem')  AS pubchem,
       -- The namespace-local ChEBI accession; the release keeps the CURIE.
       max(CASE WHEN regexp_matches(i.id, '^(?i:CHEBI):[0-9]+$')
                THEN substr(i.id, 7) ELSE i.id END) FILTER (WHERE i.ns = 'chebi') AS chebi,
       max(i.id) FILTER (WHERE i.ns = 'kegg')     AS kegg
FROM (SELECT DISTINCT metabolite_entity_key FROM metsigdb_stage) s
JOIN release_entity e ON e.entity_key = s.metabolite_entity_key
LEFT JOIN release_identifier i ON i.entity_key = e.entity_key
GROUP BY e.entity_key, coalesce(e.label, e.identifier), e.entity_type;

-- The set side: a readable name, and the species the source recorded.
--
-- The name is the label the release holds for the set entity; a set that has none carries
-- its name as a `name` identifier, which is how the source published it.
--
-- The species is read from the taxon the source recorded on its own row of the set entity.
-- A set the resource declares organism-free (a reference map) has none, whatever a source
-- stamped on it.
CREATE OR REPLACE TEMP TABLE metsigdb_set_attrs AS
SELECT s.set_entity_key AS entity_key,
       coalesce(any_value(e.label), min(i.id) FILTER (WHERE i.ns = 'name')) AS set_label,
       CASE WHEN $organism_free_pattern IS NOT NULL
                 AND regexp_matches(any_value(e.identifier), $organism_free_pattern)
            THEN NULL
            ELSE any_value(try_cast(nullif(src.taxon, '') AS BIGINT)) END AS organism
FROM (SELECT DISTINCT set_entity_key FROM metsigdb_stage) s
JOIN release_entity e ON e.entity_key = s.set_entity_key
LEFT JOIN release_identifier i ON i.entity_key = s.set_entity_key
LEFT JOIN (
  SELECT entity_key, min(taxon) AS taxon
  FROM resource_entity
  WHERE resource = $source AND taxon IS NOT NULL AND taxon <> ''
  GROUP BY entity_key
) src ON src.entity_key = s.set_entity_key
GROUP BY s.set_entity_key;
