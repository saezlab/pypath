-- gene_protein_resolution_report.sql
--
-- Evaluate gene/protein entity resolution in a built OmniPath DB: how many
-- unique genes, how many protein mentions ground to a real gene vs a virtual
-- ("unknown") gene vs stay unresolved, WHY the unresolved ones fail (by the id
-- they carry), and the same by resource + the per-resource source id type /
-- translation path. Read-only; safe to run on a built instance.
--
-- Usage (build PG has 1 GB /dev/shm; a web PG may have less — this script
-- disables parallel gather so it runs anywhere):
--   PGPASSWORD=… psql -h localhost -p <port> -U omnipath -d omnipath \
--     -f analysis/gene_protein_resolution_report.sql
--
-- Key modelling facts (verify against the DB — vocab ids are assigned SERIALLY
-- per build, so NEVER hard-code them; we look them up by name below):
--   * entity_type Gene:MI:0250 / Protein:MI:0326 — the CANONICAL entity type.
--   * entity_evidence keeps the ORIGINAL mention type; the resolved/canonical
--     type is on `entity` (reached via entity_evidence_resolution.entity_id).
--   * molecular_type 'protein' on entity_evidence_resolution = the mention was
--     a protein; its outcome is the type of the entity it resolved to:
--       - a Gene entity, resolution_mechanism resolver|gene_anchor  -> GROUNDED
--       - a Gene entity, resolution_mechanism 'unknown_gene'        -> VIRTUAL gene
--       - still a Protein entity (status unresolved)                -> UNRESOLVED
--   * entity_evidence_resolution.reason_id is NOT populated — we derive the
--     unresolved reason from the identifier types the mention carries.

\set ON_ERROR_STOP on
SET max_parallel_workers_per_gather = 0;   -- avoid /dev/shm pressure on web PGs
SET work_mem = '512MB';

-- Resolve the vocab ids we need by NAME (serial per build).
SELECT entity_type_id AS gene_type_id   FROM vocab_entity_type WHERE name = 'Gene:MI:0250' \gset
SELECT entity_type_id AS protein_type_id FROM vocab_entity_type WHERE name = 'Protein:MI:0326' \gset
SELECT molecular_type_id AS prot_mol_id FROM vocab_molecular_type WHERE name = 'protein' \gset

-- One row per protein MENTION (evidence), tagged with resource, taxon presence
-- and outcome. This is the spine of the whole report.
DROP TABLE IF EXISTS _prot;
CREATE TEMP TABLE _prot AS
SELECT eer.source_id,
       eer.entity_evidence_id,
       ds.name AS resource,
       (ee.taxonomy_id IS NULL) AS no_taxon,
       CASE WHEN vrs.name = 'unresolved'                     THEN 'unresolved'
            WHEN ge.resolution_mechanism = 'unknown_gene'    THEN 'unknown_gene'
            ELSE 'grounded' END AS outcome
FROM entity_evidence_resolution eer
JOIN entity_evidence ee USING (source_id, entity_evidence_id)
JOIN data_source ds ON ds.source_id = eer.source_id
JOIN vocab_resolution_status vrs ON vrs.resolution_status_id = eer.status_id
LEFT JOIN entity ge ON ge.entity_id = eer.entity_id
WHERE eer.molecular_type_id = :prot_mol_id;

-- One row per mention with its "primary" gene-relevant id class + boolean flags
-- for the reason breakdown. A mention can carry several ids; primary_id picks
-- the highest-priority groundable one.
DROP TABLE IF EXISTS _prot_ids;
CREATE TEMP TABLE _prot_ids AS
SELECT p.source_id, p.entity_evidence_id, p.resource, p.outcome, p.no_taxon,
  bool_or(it.name = 'Uniprot:MI:1097')          AS uniprot,
  bool_or(it.name = 'Uniprot Trembl:MI:1099')   AS trembl,
  bool_or(it.name = 'Ensembl:MI:0476')          AS ensembl,
  bool_or(it.name = 'Gene Name Primary:OM:0200') AS genesymbol,
  bool_or(it.name = 'Entrez:MI:0477')           AS entrez,
  bool_or(it.name IN ('Smiles:MI:0239','Standard Inchi Key:MI:1101','Chembl Compound:MI:0967')) AS chem_id,
  CASE WHEN bool_or(it.name = 'Uniprot:MI:1097')          THEN 'uniprot'
       WHEN bool_or(it.name = 'Ensembl:MI:0476')          THEN 'ensembl'
       WHEN bool_or(it.name = 'Entrez:MI:0477')           THEN 'entrez'
       WHEN bool_or(it.name = 'Gene Name Primary:OM:0200') THEN 'genesymbol'
       ELSE 'other' END AS primary_id
FROM _prot p
JOIN entity_evidence_identifier eei USING (source_id, entity_evidence_id)
JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
JOIN vocab_identifier_type it ON it.identifier_type_id = ie.identifier_type_id
GROUP BY 1,2,3,4,5;

\echo
\echo ==================================================================
\echo == 1. UNIQUE GENE ENTITIES (by resolution mechanism)
\echo ==================================================================
SELECT resolution_mechanism, count(*) AS unique_genes
FROM entity WHERE entity_type_id = :gene_type_id
GROUP BY 1 ORDER BY 2 DESC;
SELECT count(*) AS unique_genes_total FROM entity WHERE entity_type_id = :gene_type_id;
SELECT count(*) AS unique_unresolved_protein_entities FROM entity WHERE entity_type_id = :protein_type_id;

\echo
\echo ==================================================================
\echo == 2. PROTEIN MENTIONS by outcome (grounding coverage)
\echo ==================================================================
SELECT outcome, count(*) AS mentions,
       round(100.0*count(*)/sum(count(*)) OVER (), 2) AS pct
FROM _prot GROUP BY 1 ORDER BY 2 DESC;

\echo
\echo ==================================================================
\echo == 3. BY RESOURCE: protein mentions by outcome
\echo ==================================================================
SELECT resource,
  count(*) AS total,
  count(*) FILTER (WHERE outcome='grounded')     AS grounded,
  count(*) FILTER (WHERE outcome='unknown_gene')  AS unknown_gene,
  count(*) FILTER (WHERE outcome='unresolved')    AS unresolved,
  round(100.0*count(*) FILTER (WHERE outcome='grounded')/count(*),1) AS pct_grounded
FROM _prot GROUP BY resource ORDER BY total DESC;

\echo
\echo ==================================================================
\echo == 4. WHY UNRESOLVED: best available id class x taxon presence
\echo ==================================================================
SELECT
  CASE WHEN uniprot OR trembl THEN '1 UniProt present, not in resolver (obsolete / TrEMBL / unmaterialized organism)'
       WHEN ensembl    THEN '2 Ensembl present, not bridged to Entrez'
       WHEN genesymbol THEN '3 gene symbol only'
       WHEN entrez     THEN '4 Entrez present, non-resolving (obsolete / outside gene space)'
       WHEN chem_id    THEN '5 only chemical ids (mention mis-typed as protein)'
       ELSE '6 no groundable id (internal accession / name only)' END AS reason,
  count(*) FILTER (WHERE NOT no_taxon) AS has_taxon,
  count(*) FILTER (WHERE no_taxon)     AS no_taxon,
  count(*) AS total
FROM _prot_ids WHERE outcome='unresolved' GROUP BY 1 ORDER BY 1;

\echo
\echo ==================================================================
\echo == 5. PER-RESOURCE primary source id type + grounding rate
\echo ==    (the original id type / translation path each resource uses)
\echo ==================================================================
SELECT resource, primary_id,
  count(*) AS mentions,
  round(100.0*count(*) FILTER (WHERE outcome='grounded')/count(*),1) AS pct_grounded
FROM _prot_ids GROUP BY resource, primary_id
HAVING count(*) >= 1000
ORDER BY resource, mentions DESC;

\echo
\echo ==================================================================
\echo == 6. Translation-path reference (source id type -> Entrez anchor)
\echo ==    A gene identifier reaches its Entrez anchor through gene space
\echo ==    (NCBI gene2ensembl + Ensembl BioMart), never by routing through
\echo ==    UniProt; the build collapses ENSG/ENSP into the 'Ensembl' id type.
\echo ==================================================================
-- uniprot     -> Entrez : UniProt idmapping (uniprot->entrez), direct
-- ensembl     -> Entrez : NCBI gene2ensembl (ensp/ensg->entrez, all transcripts,
--                          772 taxa) + Ensembl BioMart (ensp->ensg, ensg->symbol)
-- genesymbol  -> Entrez : Ensembl BioMart (symbol->ensg) -> gene2ensembl, or via UniProt
-- entrez      -> Entrez : identity
-- (see omnipath-utils resolver_protein.sql / resolver_gene)
