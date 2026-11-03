-- chemical_resolution_report.sql
--
-- Chemical entity resolution — read-only report over a built OmniPath
-- database. Promoted from spec 011's research/chemical-resolution-report.sql
-- (T129) into the standing analysis directory, alongside
-- gene_protein_resolution_report.sql -- run it before and after any change
-- to the chemical resolution path to see what moved.
--
-- Usage (build PG has 1 GB /dev/shm; a web PG may have less — this script
-- disables parallel gather so it runs anywhere):
--   PGPASSWORD=… psql -h localhost -p <port> -U omnipath -d omnipath \
--     -f analysis/chemical_resolution_report.sql
--
-- Every statement is read-only. Temp tables live for the session only.
-- Vocab ids are assigned SERIALLY per build — never hard-coded, resolved by
-- name below (same convention as gene_protein_resolution_report.sql).
--
-- Some sections need the identifier-translation (omnipath-utils) database as
-- well. Those are marked "utils" and carry the query to run there against a
-- copied list.

\set ON_ERROR_STOP on
\pset pager off
SET max_parallel_workers_per_gather = 0;   -- the build containers have a small /dev/shm
SET work_mem = '512MB';

-- Resolve the vocab ids we need by NAME (serial per build).
SELECT entity_type_id AS chem_type_id FROM vocab_entity_type WHERE name = 'Chemical:OM:0037' \gset
SELECT identifier_type_id AS chebi_tid FROM vocab_identifier_type WHERE name = 'Chebi:MI:0474' \gset
SELECT identifier_type_id AS pubchem_tid FROM vocab_identifier_type WHERE name = 'Pubchem Compound:OM:0002' \gset
SELECT identifier_type_id AS inchikey_tid FROM vocab_identifier_type WHERE name = 'Standard Inchi Key:MI:1101' \gset
SELECT identifier_type_id AS smiles_tid FROM vocab_identifier_type WHERE name = 'Smiles:MI:0239' \gset
SELECT identifier_type_id AS unresolved_tid FROM vocab_identifier_type WHERE name = 'omnipath:unresolved_entity_key' \gset

\echo '== 1. build identity =='
SELECT build_id, built_at, partial_build,
       (SELECT count(*) FROM entity)   AS entities,
       (SELECT count(*) FROM relation) AS relations
  FROM build_manifest;

\echo '== 2. chemical entities by canonical identifier type and mechanism =='
SELECT coalesce(vit.name, '(none)') AS canonical_type,
       e.resolution_mechanism,
       vrs.name AS status,
       count(*) AS entities
  FROM entity e
  JOIN vocab_resolution_status vrs ON vrs.resolution_status_id = e.resolution_status_id
  LEFT JOIN vocab_identifier_type vit ON vit.identifier_type_id = e.canonical_identifier_type_id
 WHERE e.entity_type_id = :chem_type_id
 GROUP BY 1, 2, 3
 ORDER BY 4 DESC;

\echo '== 3. chemical entities per source, by canonical identifier type =='
DROP TABLE IF EXISTS source_chemical;
CREATE TEMP TABLE source_chemical AS
SELECT DISTINCT ds.name AS source, r.entity_id
  FROM entity_evidence_resolution r
  JOIN entity e ON e.entity_id = r.entity_id AND e.entity_type_id = :chem_type_id
  JOIN data_source ds ON ds.source_id = r.source_id;
CREATE INDEX ON source_chemical (entity_id);
ANALYZE source_chemical;

SELECT sc.source,
       coalesce(vit.name, '(unresolved)') AS canonical_type,
       count(*) AS entities
  FROM source_chemical sc
  JOIN entity e USING (entity_id)
  LEFT JOIN vocab_identifier_type vit ON vit.identifier_type_id = e.canonical_identifier_type_id
 GROUP BY 1, 2
 ORDER BY 1, 3 DESC;

\echo '== 4. shared metabolites between resources =='
SELECT a.source, b.source, count(*) AS shared_entities
  FROM source_chemical a
  JOIN source_chemical b ON b.entity_id = a.entity_id AND a.source < b.source
 GROUP BY 1, 2
HAVING count(*) > 0
 ORDER BY 3 DESC
 LIMIT 60;

\echo '== 5. which identifiers each source supplies for chemicals =='
SELECT ds.name AS source, vit.name AS identifier_type, count(*) AS mention_rows
  FROM entity_evidence ee
  JOIN data_source ds ON ds.source_id = ee.source_id
  JOIN entity_evidence_identifier eei
    ON eei.source_id = ee.source_id AND eei.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
  JOIN vocab_identifier_type vit ON vit.identifier_type_id = ie.identifier_type_id
 WHERE ee.entity_type_id = :chem_type_id
 GROUP BY 1, 2
 ORDER BY 1, 3 DESC;

\echo '== 6. structure authorities: sources that supply SMILES or InChIKey =='
SELECT ds.name AS source,
       count(*) FILTER (WHERE ie.identifier_type_id = :smiles_tid) AS smiles,
       count(*) FILTER (WHERE ie.identifier_type_id = :inchikey_tid) AS inchikey
  FROM entity_evidence ee
  JOIN data_source ds ON ds.source_id = ee.source_id
  JOIN entity_evidence_identifier eei
    ON eei.source_id = ee.source_id AND eei.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
 WHERE ee.entity_type_id = :chem_type_id
   AND ie.identifier_type_id IN (:inchikey_tid, :smiles_tid)
 GROUP BY 1
HAVING count(*) FILTER (WHERE ie.identifier_type_id IN (:inchikey_tid, :smiles_tid)) > 0
 ORDER BY 2 DESC NULLS LAST, 3 DESC;

\echo '== 7. cross-reference consistency: every source against ChEBI =='
DROP TABLE IF EXISTS chebi_claim;
CREATE TEMP TABLE chebi_claim AS
SELECT ds.name AS source, ch.value AS chebi, upper(ik.value) AS inchikey
  FROM entity_evidence ee
  JOIN data_source ds ON ds.source_id = ee.source_id
  JOIN entity_evidence_identifier a
    ON a.source_id = ee.source_id AND a.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence ch ON ch.identifier_id = a.identifier_id AND ch.identifier_type_id = :chebi_tid
  JOIN entity_evidence_identifier b
    ON b.source_id = ee.source_id AND b.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence ik ON ik.identifier_id = b.identifier_id AND ik.identifier_type_id = :inchikey_tid
 WHERE ee.entity_type_id = :chem_type_id;

DROP TABLE IF EXISTS chebi_authority;
CREATE TEMP TABLE chebi_authority AS
SELECT chebi, min(inchikey) AS inchikey, count(DISTINCT inchikey) AS own_keys
  FROM chebi_claim WHERE source = 'chebi' GROUP BY 1;

SELECT 'ChEBI identifiers with more than one own InChIKey' AS check, count(*) AS n
  FROM chebi_authority WHERE own_keys > 1;

SELECT c.source,
       count(*) AS compared,
       count(*) FILTER (WHERE c.inchikey = a.inchikey) AS agree,
       count(*) FILTER (WHERE c.inchikey <> a.inchikey
                          AND left(c.inchikey, 14) = left(a.inchikey, 14)) AS layer_difference,
       count(*) FILTER (WHERE left(c.inchikey, 14) <> left(a.inchikey, 14)) AS different_structure
  FROM chebi_claim c
  JOIN chebi_authority a USING (chebi)
 WHERE c.source <> 'chebi'
 GROUP BY 1
 ORDER BY 2 DESC;

\echo '== 8. cross-reference pairs: ChEBI against PubChem, per resource =='
\echo '   Export the pairs, then resolve both sides in the utils database.'
DROP TABLE IF EXISTS xref_pair;
CREATE TEMP TABLE xref_pair AS
SELECT DISTINCT ds.name AS source, ch.value AS chebi, pc.value AS pubchem
  FROM entity_evidence ee
  JOIN data_source ds ON ds.source_id = ee.source_id
  JOIN entity_evidence_identifier a
    ON a.source_id = ee.source_id AND a.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence ch ON ch.identifier_id = a.identifier_id AND ch.identifier_type_id = :chebi_tid
  JOIN entity_evidence_identifier b
    ON b.source_id = ee.source_id AND b.entity_evidence_id = ee.entity_evidence_id
  JOIN identifier_evidence pc ON pc.identifier_id = b.identifier_id AND pc.identifier_type_id = :pubchem_tid
 WHERE ee.entity_type_id = :chem_type_id;
SELECT source, count(*) AS pairs FROM xref_pair GROUP BY 1 ORDER BY 2 DESC;
-- \copy (SELECT * FROM xref_pair) TO '/tmp/xref_pairs.csv' WITH (FORMAT csv)

\echo '== 9. the unresolved chemical tail: which identifiers it carries =='
SELECT vit.name AS identifier_type,
       count(*) AS mention_rows,
       count(DISTINCT r.entity_id) AS entities
  FROM entity e
  JOIN entity_evidence_resolution r ON r.entity_id = e.entity_id
  JOIN entity_evidence_identifier eei
    ON eei.source_id = r.source_id AND eei.entity_evidence_id = r.entity_evidence_id
  JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
  JOIN vocab_identifier_type vit ON vit.identifier_type_id = ie.identifier_type_id
 WHERE e.entity_type_id = :chem_type_id AND e.canonical_identifier_type_id = :unresolved_tid
 GROUP BY 1
 ORDER BY 2 DESC;

-- \copy (SELECT DISTINCT e.entity_id::text, vit.name, ie.value FROM entity e
--   JOIN entity_evidence_resolution r ON r.entity_id = e.entity_id
--   JOIN entity_evidence_identifier eei ON eei.source_id = r.source_id
--        AND eei.entity_evidence_id = r.entity_evidence_id
--   JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
--   JOIN vocab_identifier_type vit ON vit.identifier_type_id = ie.identifier_type_id
--  WHERE e.entity_type_id = :chem_type_id AND e.canonical_identifier_type_id = :unresolved_tid)
--   TO '/tmp/unresolved_chemical_ids.csv' WITH (FORMAT csv)

\echo '== 10. malformed canonical structure keys =='
SELECT count(*) FILTER (WHERE canonical_identifier LIKE 'InChIKey=%')                AS prefixed,
       count(*) FILTER (WHERE canonical_identifier ~ '^[A-Z]{14}-[A-Z]{10}-[A-Z]$')  AS well_formed,
       count(*)                                                                      AS total
  FROM entity
 WHERE entity_type_id = :chem_type_id AND canonical_identifier_type_id = :inchikey_tid;

\echo '== 11. dead resolver tables (expect 0 until the build writes them) =='
SELECT 'resolver_chemical_identifier_lookup' AS relation, count(*) FROM resolver_chemical_identifier_lookup
UNION ALL
SELECT 'resolver_chemical_identifier_lookup_ambiguous', count(*) FROM resolver_chemical_identifier_lookup_ambiguous;

\echo '== 12. per-resource coverage, from the standing table (spec 011 T126-T128) =='
SELECT ds.name AS source, vit.name AS namespace, crc.role,
       crc.mentions, crc.entities, crc.reached_structure, crc.reached_name,
       crc.unresolved, crc.conflicted
  FROM chemical_resolution_coverage crc
  JOIN data_source ds ON ds.source_id = crc.source_id
  JOIN vocab_identifier_type vit ON vit.identifier_type_id = crc.identifier_type_id
 ORDER BY ds.name, vit.name, crc.role;

-- ---------------------------------------------------------------------------
-- utils: run these against omnipath-utils, after copying the exported lists in.
-- ---------------------------------------------------------------------------
-- SET search_path = omnipath_utils, public;
--
-- -- source coverage of the identifier-to-structure projection
-- SELECT source_type,
--        count(*) AS rows,
--        count(DISTINCT source_id) AS ids,
--        count(*) FILTER (WHERE inchikey LIKE 'InChIKey=%') AS malformed
--   FROM resolver_chemical GROUP BY 1 ORDER BY 2 DESC;
--
-- -- how much of the ChEBI-canonical tail would resolve with a normalized key
-- CREATE TEMP TABLE cb (id text);
-- \copy cb FROM '/tmp/chebi_ids.txt'
-- SELECT count(*) AS resolvable, count(*) FILTER (WHERE n = 1) AS unambiguous
--   FROM (SELECT cb.id, count(DISTINCT r.inchikey) AS n
--           FROM cb JOIN resolver_chemical r
--             ON r.source_type = 'chebi' AND r.source_id = 'CHEBI:' || cb.id
--          GROUP BY 1) t;
