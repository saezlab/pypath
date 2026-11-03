-- resolver_gene_protein_global (spec 002 T069/R25/US7; split out 007 R10 Phase 3P
-- T067). Pure id_mapping_ftp (uniprot/entrez -> entrez, taxon-agnostic), so it is
-- an FTP core: create_resolver_views rebuilds it only on an FTP reload.

-- Global (taxon-agnostic) UniProt/Entrez gene anchor (spec 002, T069 / R25 / US7).
--
-- A UniProt AC (and an Entrez id) uniquely determines its organism + gene, so this
-- slice carries NO taxon filter: it lets omnipath-build gene-anchor proteins that
-- arrive with no (or a mismatched) taxonomy — e.g. Rhea / Brenda / TCDB / ChEMBL
-- enzyme & target participants, which reference UniProt without an organism. The
-- resolved entity then inherits the organism *derived from the AC* (the row's
-- ncbi_tax_id) via the build's coalesce(rl.taxonomy_id, ee.taxonomy_id). Only
-- uniprot/entrez source types belong here; genesymbol/ensg/ensp stay per-taxon in
-- resolver_gene because a symbol repeats across organisms. The build's
-- taxonomy-optional-unambiguous-key gate keeps an AC that maps to >1 gene from
-- force-merging (it falls to the ambiguous lookup / multi-gene split instead).
--
-- Scale: ~24M rows over ~50k taxa (the full id_mapping_ftp uniprot->entrez map);
-- omnipath-build's needed_resolver_lookup filters this to the ACs actually
-- referenced, so resolution stays cheap. Read by omnipath-build via DuckDB ATTACH.
-- MATERIALIZED + indexed (2026-07-04): the build's keyed lookup probes this for
-- proteins that arrive with NO taxon (Rhea/Brenda/TCDB/ChEMBL reference UniProt
-- without an organism). As a plain view a keyed join seq-scans the 37M-row
-- derivation; materialise + a (source_type, source_id) index so it is an index
-- probe. The build emits its matches taxon-agnostically (taxonomy NULL) so the
-- no-taxon evidence matches — the taxon-bearing resolver_gene rows cannot.
DO $$
DECLARE k "char";
BEGIN
  SELECT c.relkind INTO k FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
   WHERE n.nspname = 'omnipath_utils' AND c.relname = 'resolver_gene_protein_global';
  IF k = 'v' THEN EXECUTE 'DROP VIEW omnipath_utils.resolver_gene_protein_global CASCADE';
  ELSIF k = 'm' THEN EXECUTE 'DROP MATERIALIZED VIEW omnipath_utils.resolver_gene_protein_global CASCADE';
  END IF;
END $$;
CREATE MATERIALIZED VIEW omnipath_utils.resolver_gene_protein_global AS
WITH up_entrez AS (
    SELECT DISTINCT m.ncbi_tax_id, m.source_id AS uniprot, m.target_id AS entrez
    FROM omnipath_utils.id_mapping_ftp m
    JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
    JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id AND tt.name = 'entrez'
)
-- uniprot -> entrez (the accession itself as a source id)
SELECT DISTINCT ncbi_tax_id, 'uniprot' AS source_type, uniprot AS source_id, entrez
FROM up_entrez
WHERE uniprot IS NOT NULL AND entrez IS NOT NULL
UNION
-- entrez -> entrez (identity)
SELECT DISTINCT ncbi_tax_id, 'entrez' AS source_type, entrez AS source_id, entrez
FROM up_entrez
WHERE entrez IS NOT NULL;

CREATE INDEX IF NOT EXISTS resolver_gene_protein_global_st_si_idx
    ON omnipath_utils.resolver_gene_protein_global (source_type, source_id)
    INCLUDE (entrez);
