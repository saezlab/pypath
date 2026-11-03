-- resolver_gene FTP core (spec 007 R10 / Phase 3P, T064/T066).
--
-- The id_mapping_ftp-derived half of resolver_gene. This is the EXPENSIVE part
-- (the full-UniProt idmapping inversion, ~46 GB heap) and it changes ONLY on an
-- FTP reload, so create_resolver_views() rebuilds it only when id_mapping_ftp
-- changed (a full FTP load DROPs id_mapping_ftp CASCADE, dropping this table; its
-- absence is the rebuild trigger). Every additive/curated load reuses it as-is and
-- only rebuilds the cheap resolver_gene_curated delta.
--
-- resolver_gene (the name consumers read: omnipath-build DuckDB ATTACH, web/API)
-- is a VIEW = resolver_gene_ftp UNION ALL resolver_gene_curated (see
-- resolver_protein.sql). Columns are unchanged: (ncbi_tax_id, source_type,
-- source_id, entrez).
--
-- Perf (T063 diagnosis -> T066 fix): the old monolithic resolver_gene was fully
-- serial (~3 h) because the up_* CTEs over id_mapping_ftp were multiply-referenced
-- WITH fences that Postgres materialised once and then sorted-with-spill at
-- work_mem=64MB. Here we (1) SET LOCAL work_mem high + enable parallelism, (2)
-- pre-materialise each up_* extraction ONCE (a single parallel seq-scan of
-- id_mapping_ftp each, not one per reference) into an indexed UNLOGGED table, then
-- (3) do the bridge joins as cheap indexed joins and finish with one UNION dedup.
--
-- The up_* intermediates are UNLOGGED regular tables, NOT TEMP: parallel workers
-- cannot read a leader's TEMP tables, so a TEMP intermediate would force the final
-- union to run serially. They are created and dropped inside create_resolver_views'
-- single transaction, so a mid-script failure rolls them back automatically; the
-- explicit DROPs at the end clean up on success.

SET LOCAL work_mem = '4GB';
SET LOCAL max_parallel_workers_per_gather = 8;
SET LOCAL max_parallel_workers = 8;
SET LOCAL maintenance_work_mem = '4GB';
SET LOCAL max_parallel_maintenance_workers = 8;

CREATE SCHEMA IF NOT EXISTS omnipath_utils;

-- Drop whatever kind currently exists (view/matview from an earlier era, or the
-- table this file now builds) — IF EXISTS does not span a relkind mismatch. The
-- resolver_gene VIEW and resolver_gene_curated depend on this, so CASCADE.
DO $$
DECLARE k "char";
BEGIN
  SELECT c.relkind INTO k FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
   WHERE n.nspname = 'omnipath_utils' AND c.relname = 'resolver_gene_ftp';
  IF k = 'v' THEN EXECUTE 'DROP VIEW omnipath_utils.resolver_gene_ftp CASCADE';
  ELSIF k = 'm' THEN EXECUTE 'DROP MATERIALIZED VIEW omnipath_utils.resolver_gene_ftp CASCADE';
  ELSIF k = 'r' THEN EXECUTE 'DROP TABLE omnipath_utils.resolver_gene_ftp CASCADE';
  END IF;
END $$;

DROP TABLE IF EXISTS omnipath_utils._rgf_up_entrez;
DROP TABLE IF EXISTS omnipath_utils._rgf_up_ensg;
DROP TABLE IF EXISTS omnipath_utils._rgf_up_ensp;
DROP TABLE IF EXISTS omnipath_utils._rgf_up_symbol;
DROP TABLE IF EXISTS omnipath_utils._rgf_ensg_entrez;

-- ===== base extractions from id_mapping_ftp (materialised ONCE, indexed) =====
CREATE UNLOGGED TABLE omnipath_utils._rgf_up_entrez AS
SELECT DISTINCT m.ncbi_tax_id, m.source_id AS uniprot, m.target_id AS entrez
FROM omnipath_utils.id_mapping_ftp m
JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez';
CREATE INDEX ON omnipath_utils._rgf_up_entrez (ncbi_tax_id, uniprot);

CREATE UNLOGGED TABLE omnipath_utils._rgf_up_ensg AS
SELECT DISTINCT m.ncbi_tax_id, m.source_id AS uniprot, m.target_id AS ensg
FROM omnipath_utils.id_mapping_ftp m
JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensg';
CREATE INDEX ON omnipath_utils._rgf_up_ensg (ncbi_tax_id, uniprot);

CREATE UNLOGGED TABLE omnipath_utils._rgf_up_ensp AS
SELECT DISTINCT m.ncbi_tax_id, m.source_id AS uniprot, m.target_id AS ensp
FROM omnipath_utils.id_mapping_ftp m
JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensp';
CREATE INDEX ON omnipath_utils._rgf_up_ensp (ncbi_tax_id, uniprot);

CREATE UNLOGGED TABLE omnipath_utils._rgf_up_symbol AS
SELECT DISTINCT m.ncbi_tax_id, m.source_id AS uniprot, m.target_id AS genesymbol
FROM omnipath_utils.id_mapping_ftp m
JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'genesymbol'
WHERE m.target_id IS NOT NULL;
CREATE INDEX ON omnipath_utils._rgf_up_symbol (ncbi_tax_id, uniprot);

-- ensg -> entrez over a shared UniProt (the id_mapping_ftp bridge).
CREATE UNLOGGED TABLE omnipath_utils._rgf_ensg_entrez AS
SELECT DISTINCT g.ncbi_tax_id, g.ensg, e.entrez
FROM omnipath_utils._rgf_up_ensg g
JOIN omnipath_utils._rgf_up_entrez e
  ON e.ncbi_tax_id = g.ncbi_tax_id AND e.uniprot = g.uniprot;
CREATE INDEX ON omnipath_utils._rgf_ensg_entrez (ncbi_tax_id, ensg);

-- ===== resolver_gene_ftp: the id_mapping_ftp-derived UNION branches =====
-- Byte-for-byte the same branches as the old monolithic resolver_gene EXCEPT the
-- two secondary-UniProt-AC branches, which move to resolver_gene_curated (they
-- depend on the curated uniprot-sec map and are reformulated there to join this
-- table, so no id_mapping_ftp re-scan). Column names come from the first branch.
CREATE TABLE omnipath_utils.resolver_gene_ftp AS
-- entrez -> entrez (identity)
SELECT DISTINCT ncbi_tax_id, 'entrez'::varchar(64) AS source_type,
       entrez AS source_id, entrez
FROM omnipath_utils._rgf_up_entrez
UNION
-- ensg -> entrez (gene space)
SELECT DISTINCT ncbi_tax_id, 'ensg', ensg, entrez FROM omnipath_utils._rgf_ensg_entrez
UNION
-- ensp -> ensg -> entrez (protein -> gene -> gene; NOT via uniprot)
SELECT DISTINCT p.ncbi_tax_id, 'ensp', p.ensp, ge.entrez
FROM omnipath_utils._rgf_up_ensp p
JOIN omnipath_utils._rgf_up_ensg g
  ON g.ncbi_tax_id = p.ncbi_tax_id AND g.uniprot = p.uniprot
JOIN omnipath_utils._rgf_ensg_entrez ge
  ON ge.ncbi_tax_id = p.ncbi_tax_id AND ge.ensg = g.ensg
UNION
-- genesymbol -> ensg -> entrez
SELECT DISTINCT s.ncbi_tax_id, 'genesymbol', s.genesymbol, ge.entrez
FROM omnipath_utils._rgf_up_symbol s
JOIN omnipath_utils._rgf_up_ensg g
  ON g.ncbi_tax_id = s.ncbi_tax_id AND g.uniprot = s.uniprot
JOIN omnipath_utils._rgf_ensg_entrez ge
  ON ge.ncbi_tax_id = s.ncbi_tax_id AND ge.ensg = g.ensg
UNION
-- uniprot -> entrez (the accession itself; uniprot IS a protein, direct is fine)
SELECT DISTINCT ncbi_tax_id, 'uniprot', uniprot, entrez FROM omnipath_utils._rgf_up_entrez
UNION
-- uniprot -> ensg -> entrez (recover entries with an ensg but no GeneID)
SELECT DISTINCT g.ncbi_tax_id, 'uniprot', g.uniprot, ge.entrez
FROM omnipath_utils._rgf_up_ensg g
JOIN omnipath_utils._rgf_ensg_entrez ge
  ON ge.ncbi_tax_id = g.ncbi_tax_id AND ge.ensg = g.ensg
UNION
-- genesymbol -> entrez DIRECT via a single shared uniprot (supplement to the ensg
-- path): recovers symbols whose UniProt carries a GeneID but no ensg.
SELECT DISTINCT s.ncbi_tax_id, 'genesymbol', s.genesymbol, e.entrez
FROM omnipath_utils._rgf_up_symbol s
JOIN omnipath_utils._rgf_up_entrez e
  ON e.ncbi_tax_id = s.ncbi_tax_id AND e.uniprot = s.uniprot
UNION
-- uniprot -> its idmapping gene symbol -> entrez (NCBI gene_info). Recovers a
-- UniProt that the idmapping gives ONLY a gene symbol -- no GeneID and no Ensembl
-- gene -- which is the norm for the TrEMBL long tail of less-studied organisms, so
-- none of the arms above (direct entrez, ensg bridge) can reach its gene. Anchors
-- it when its symbol is a real NCBI gene symbol OR a synonym for the same organism
-- ('genesymbol' + 'genesymbol-syn' in gene_info); the synonym side closes part of
-- the UniProt-vs-NCBI symbol divergence. gene_info (curated id_mapping) is small
-- and, like id_mapping_ftp, is refreshed with this FTP core at each gene rebuild.
SELECT DISTINCT s.ncbi_tax_id, 'uniprot', s.uniprot, gi.entrez
FROM omnipath_utils._rgf_up_symbol s
JOIN (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS genesymbol, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id
     AND st.name IN ('genesymbol', 'genesymbol-syn')
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
) gi ON gi.ncbi_tax_id = s.ncbi_tax_id AND gi.genesymbol = s.genesymbol;

DROP TABLE omnipath_utils._rgf_up_entrez;
DROP TABLE omnipath_utils._rgf_up_ensg;
DROP TABLE omnipath_utils._rgf_up_ensp;
DROP TABLE omnipath_utils._rgf_up_symbol;
DROP TABLE omnipath_utils._rgf_ensg_entrez;

-- Keyed-lookup indexes (same shape as the old monolithic resolver_gene, so the
-- resolver_gene UNION-ALL view push-probes each child by (source_type, source_id)).
CREATE INDEX resolver_gene_ftp_st_si_idx
    ON omnipath_utils.resolver_gene_ftp (source_type, source_id)
    INCLUDE (ncbi_tax_id, entrez);
CREATE INDEX resolver_gene_ftp_key_idx
    ON omnipath_utils.resolver_gene_ftp (ncbi_tax_id, source_type, source_id);

ANALYZE omnipath_utils.resolver_gene_ftp;
