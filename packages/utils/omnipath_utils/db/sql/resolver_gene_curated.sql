-- resolver_gene curated delta + the resolver_gene UNION-ALL view
-- (spec 007 R10 / Phase 3P, T064/T067). Split out of resolver_protein.sql so
-- create_resolver_views can rebuild it independently (cheap: ~166s on the 5103
-- cumulative base). Depends on resolver_gene_ftp (the FTP core) existing.

-- ===== resolver_gene: curated delta + FTP core, unioned =====
-- Gene-anchor projection (spec 002, M-Genes / US7 / FR-026): maps any in-scope
-- source identifier to its NCBI Gene (Entrez) anchor, per taxon. The gene is the
-- canonical collapsing identity; entrez is the anchor (R17 — Entrez covers 50,018
-- taxa vs Ensembl's 240). entrez maps to itself.
--
-- Columns: (ncbi_tax_id, source_type, source_id, entrez). omnipath-build reads
-- resolver_gene via DuckDB ATTACH; web/API read it too. The name and columns are
-- UNCHANGED — resolver_gene is now a VIEW = resolver_gene_ftp UNION ALL
-- resolver_gene_curated (spec 007 R10 / Phase 3P, T064):
--   * resolver_gene_ftp     — the id_mapping_ftp-derived branches (expensive;
--     rebuilt only on an FTP reload — see resolver_gene_ftp.sql).
--   * resolver_gene_curated — the curated-id_mapping branches incl. all the 007
--     US1 anchors (gene_info / gene2accession / KEGG / Ensembl Genomes / NCBI
--     gene2ensembl / Ensembl BioMart). Cheap; rebuilt on EVERY additive load so a
--     newly loaded curated mapping is picked up without re-deriving the FTP core.
--
-- UNION ALL (not UNION): so the resolver_gene VIEW flattens to an appendrel and the
-- build's keyed join push-probes each child's (source_type, source_id) index,
-- instead of materialising + deduping the ~83M-row union on every read (a UNION
-- set-op is a join-pushdown barrier). Rows that both children emit therefore appear
-- twice in the view; that is inert — every consumer collapses with SELECT DISTINCT
-- / count(DISTINCT) / row_number (omnipath-build resolver_lookup gates, and the
-- combined resolver below), so a duplicate never changes a resolution or an
-- ambiguity count.

-- DROP + parallel CTAS (not INSERT ... SELECT, whose SELECT will not parallelise —
-- the CTAS-vs-INSERT gotcha). The DROP ... CASCADE takes out the resolver_gene view
-- and, transitively, resolver_gene_protein_combined (which reads the view). The view
-- is recreated below; combined is NOT rebuilt here — resolver_combined.sql rebuilds
-- it on a full run, and a gene-only rebuild (resolvers={'gene'}) intentionally leaves
-- it absent until the next full/promotion build (nothing in the live keyed path reads
-- it; only omnipath-build's pre-built parquet/duckdb resolver mode does).
-- The resolver_gene view (recreated at the end of this file) canonicalizes each
-- Entrez anchor's taxon to species level via taxon_species (strain->species, from
-- the NCBI taxdump; see load_taxonomy_species). Guarantee the table exists so the
-- view's dependency holds even when this file runs before the taxonomy loader (an
-- empty table just leaves every taxon unchanged).
CREATE TABLE IF NOT EXISTS omnipath_utils.taxon_species (
    tax_id integer PRIMARY KEY,
    species_tax_id integer NOT NULL
);

DROP TABLE IF EXISTS omnipath_utils.resolver_gene_curated CASCADE;
CREATE TABLE omnipath_utils.resolver_gene_curated AS
-- secondary -> primary UniProt AC (organism-agnostic, tax 0; ADR 0006). Lets a
-- resource that supplies a secondary accession still anchor to the gene.
WITH sec_pri AS (
    SELECT m.source_id AS sec, m.target_id AS pri
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot-sec'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'uniprot-pri'
),
-- AUTHORITATIVE gene space (curated id_mapping): NCBI gene2ensembl gives ensp->entrez
-- and ensg->entrez DIRECT (all transcripts, 772 taxa); Ensembl BioMart gives
-- ensp->ensg and ensg->genesymbol. These are the primary paths.
g2e_ensp AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensp, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensp'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
),
g2e_ensg AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensg, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensg'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
),
bm_ensp_ensg AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensp, m.target_id AS ensg
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensp'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensg'
),
bm_symbol_ensg AS (
    SELECT DISTINCT m.ncbi_tax_id, m.target_id AS genesymbol, m.source_id AS ensg
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensg'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'genesymbol'
),
-- ===== 007 US1: authoritative NCBI gene-space + KEGG anchors =====
-- gene_info: genesymbol / synonym -> entrez DIRECT, all organisms, independent of
-- UniProt (covers genes with no protein product). Both primary symbols and synonyms
-- are emitted under source_type 'genesymbol' so a symbol lookup that hits a synonym
-- still resolves; ambiguous synonyms are dropped by the combined resolver's gate.
gi_symbol AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS genesymbol, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id
     AND st.name IN ('genesymbol', 'genesymbol-syn')
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
),
-- gene2accession: RefSeq RNA/protein -> entrez DIRECT.
g2a_refseq AS (
    SELECT DISTINCT m.ncbi_tax_id, st.name AS source_type,
           m.source_id AS refseq, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id
     AND st.name IN ('refseqn', 'refseqp')
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
),
-- KEGG gene id -> entrez.
kg_entrez AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS kegg_gene, m.target_id AS entrez
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'kegg_gene'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
),
-- Ensembl Genomes divisions (BioMart genomes): gene -> genesymbol, and
-- protein/transcript -> gene. Anchored to entrez via the gi_symbol bridge below.
egg_symbol AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensgg, m.target_id AS genesymbol
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensgg'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'genesymbol'
),
egg_prot AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensgp, m.target_id AS ensgg
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensgp'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensgg'
),
egg_tx AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS ensgt, m.target_id AS ensgg
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'ensgt'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensgg'
),
-- 007 US3: RNAcentral bridges a miRBase accession (MI####/MIMAT####) to its
-- Ensembl gene and gene symbol via the shared URS sequence id (load_rnacentral_
-- mirbase_gene). Anchored to entrez below through the gene2ensembl ensg path and
-- the gene_info symbol path, so a miRBase mention reaches its host gene.
rc_mirbase_ensg AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS mirbase, m.target_id AS ensg
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'mirbase'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'ensg'
),
rc_mirbase_symbol AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS mirbase, m.target_id AS genesymbol
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'mirbase'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'genesymbol'
)
-- ===== AUTHORITATIVE gene-space paths (NCBI gene2ensembl + Ensembl BioMart) =====
-- Column names come from the first branch: (ncbi_tax_id, source_type, source_id, entrez).
-- ensp -> entrez DIRECT (gene2ensembl, every transcript, 772 taxa).
SELECT DISTINCT ncbi_tax_id, 'ensp'::varchar(64) AS source_type,
       ensp AS source_id, entrez
FROM g2e_ensp
UNION
-- ensp -> ensg (BioMart) -> entrez (gene2ensembl): catches ENSPs not directly in g2e.
SELECT DISTINCT b.ncbi_tax_id, 'ensp', b.ensp, ge.entrez
FROM bm_ensp_ensg b
JOIN g2e_ensg ge ON ge.ncbi_tax_id = b.ncbi_tax_id AND ge.ensg = b.ensg
UNION
-- ensg -> entrez DIRECT (gene2ensembl).
SELECT DISTINCT ncbi_tax_id, 'ensg', ensg, entrez FROM g2e_ensg
UNION
-- genesymbol -> ensg (BioMart) -> entrez (gene2ensembl): authoritative symbol path.
SELECT DISTINCT s.ncbi_tax_id, 'genesymbol', s.genesymbol, ge.entrez
FROM bm_symbol_ensg s
JOIN g2e_ensg ge ON ge.ncbi_tax_id = s.ncbi_tax_id AND ge.ensg = s.ensg
UNION
-- ===== 007 US1: DIRECT NCBI gene-space + KEGG anchors =====
-- genesymbol (and synonyms) -> entrez DIRECT (gene_info).
SELECT DISTINCT ncbi_tax_id, 'genesymbol', genesymbol, entrez FROM gi_symbol
UNION
-- refseqn / refseqp -> entrez DIRECT (gene2accession).
SELECT DISTINCT ncbi_tax_id, source_type, refseq, entrez FROM g2a_refseq
UNION
-- kegg_gene -> entrez DIRECT (kegg_gene).
SELECT DISTINCT ncbi_tax_id, 'kegg_gene', kegg_gene, entrez FROM kg_entrez
UNION
-- Ensembl Genomes gene -> genesymbol -> entrez (division bridge via gene_info).
SELECT DISTINCT e.ncbi_tax_id, 'ensgg', e.ensgg, g.entrez
FROM egg_symbol e
JOIN gi_symbol g ON g.ncbi_tax_id = e.ncbi_tax_id AND g.genesymbol = e.genesymbol
UNION
-- Ensembl Genomes protein -> gene -> genesymbol -> entrez.
SELECT DISTINCT p.ncbi_tax_id, 'ensgp', p.ensgp, g.entrez
FROM egg_prot p
JOIN egg_symbol e ON e.ncbi_tax_id = p.ncbi_tax_id AND e.ensgg = p.ensgg
JOIN gi_symbol g ON g.ncbi_tax_id = p.ncbi_tax_id AND g.genesymbol = e.genesymbol
UNION
-- Ensembl Genomes transcript -> gene -> genesymbol -> entrez.
SELECT DISTINCT t.ncbi_tax_id, 'ensgt', t.ensgt, g.entrez
FROM egg_tx t
JOIN egg_symbol e ON e.ncbi_tax_id = t.ncbi_tax_id AND e.ensgg = t.ensgg
JOIN gi_symbol g ON g.ncbi_tax_id = t.ncbi_tax_id AND g.genesymbol = e.genesymbol
UNION
-- miRBase -> ensg (RNAcentral) -> entrez (gene2ensembl).
SELECT DISTINCT r.ncbi_tax_id, 'mirbase', r.mirbase, ge.entrez
FROM rc_mirbase_ensg r
JOIN g2e_ensg ge ON ge.ncbi_tax_id = r.ncbi_tax_id AND ge.ensg = r.ensg
UNION
-- miRBase -> genesymbol (RNAcentral, HGNC) -> entrez (gene_info).
SELECT DISTINCT r.ncbi_tax_id, 'mirbase', r.mirbase, g.entrez
FROM rc_mirbase_symbol r
JOIN gi_symbol g ON g.ncbi_tax_id = r.ncbi_tax_id AND g.genesymbol = r.genesymbol
UNION
-- ===== secondary UniProt AC -> primary -> entrez (curated sec_pri x FTP core) =====
-- Reformulation of the two old up_entrez/up_ensg secondary-AC branches: the FTP
-- core already emits every ('uniprot', primary_ac -> entrez) row (direct AND via
-- the ensg bridge), so a secondary AC anchors by joining resolver_gene_ftp on its
-- primary, WITHOUT re-scanning id_mapping_ftp. Provably equal to the old
-- (up_entrez JOIN sec_pri) UNION (up_ensg JOIN sec_pri JOIN ensg_entrez).
SELECT DISTINCT rgf.ncbi_tax_id, 'uniprot', sp.sec, rgf.entrez
FROM sec_pri sp
JOIN omnipath_utils.resolver_gene_ftp rgf
  ON rgf.source_type = 'uniprot' AND rgf.source_id = sp.pri
WHERE sp.sec IS NOT NULL;

-- Keyed-lookup indexes (mirror resolver_gene_ftp so the UNION-ALL view probes both
-- children the same way).
CREATE INDEX resolver_gene_curated_st_si_idx
    ON omnipath_utils.resolver_gene_curated (source_type, source_id)
    INCLUDE (ncbi_tax_id, entrez);
CREATE INDEX resolver_gene_curated_key_idx
    ON omnipath_utils.resolver_gene_curated (ncbi_tax_id, source_type, source_id);
ANALYZE omnipath_utils.resolver_gene_curated;

-- resolver_gene: the unchanged name + columns consumers read — a UNION-ALL view over
-- the FTP core + the curated delta. The DROP ... CASCADE above already removed any
-- prior view; a one-time transition drop covers the pre-split case where resolver_gene
-- is still the monolithic materialised view / table.
DO $$
DECLARE k "char";
BEGIN
  SELECT c.relkind INTO k FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
   WHERE n.nspname = 'omnipath_utils' AND c.relname = 'resolver_gene';
  IF k = 'v' THEN EXECUTE 'DROP VIEW omnipath_utils.resolver_gene CASCADE';
  ELSIF k = 'm' THEN EXECUTE 'DROP MATERIALIZED VIEW omnipath_utils.resolver_gene CASCADE';
  ELSIF k = 'r' THEN EXECUTE 'DROP TABLE omnipath_utils.resolver_gene CASCADE';
  END IF;
END $$;
-- The taxon each Entrez anchor is emitted under is canonicalized to species level
-- via taxon_species (an Entrez GeneID belongs to one species; the id_mapping_ftp
-- taxon can be a strain, so the same GeneID otherwise surfaces under several NCBI
-- strain taxa -> a gene duplicated across taxa in the built graph). coalesce keeps
-- the original taxon when taxon_species has no row (species already, or genus and
-- above). The keyed lookup is unchanged: consumers probe by (source_type,
-- source_id) and READ the taxon, so canonicalizing it does not move the key.
CREATE VIEW omnipath_utils.resolver_gene AS
SELECT
  coalesce(ts.species_tax_id, g.ncbi_tax_id) AS ncbi_tax_id,
  g.source_type, g.source_id, g.entrez
FROM (
  SELECT ncbi_tax_id, source_type, source_id, entrez
  FROM omnipath_utils.resolver_gene_ftp
  UNION ALL
  SELECT ncbi_tax_id, source_type, source_id, entrez
  FROM omnipath_utils.resolver_gene_curated
) g
LEFT JOIN omnipath_utils.taxon_species ts ON ts.tax_id = g.ncbi_tax_id;
