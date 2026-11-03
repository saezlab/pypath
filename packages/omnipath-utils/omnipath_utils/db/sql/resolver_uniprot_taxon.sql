-- UniProt accession -> organism (NCBI taxon) lookup.
--
-- One row per UniProt accession giving the organism it belongs to. Consumers
-- (omnipath-build via DuckDB ATTACH) use it to recover the organism of a protein
-- mention that cited an accession but no organism: many reaction and binding-
-- affinity resources reference a UniProt accession without stating the species,
-- yet the accession itself names the organism. Recovering it lets such a mention
-- reach its gene through the per-organism gene map and, failing that, be typed as
-- a gene keyed by the accession (which requires a known organism).
--
-- Derived from the full-UniProt FTP idmapping (source rows keyed by accession,
-- each carrying its organism). Like the other id_mapping_ftp-derived cores this is
-- the expensive part (a single parallel scan of the ~46 GB table) and changes ONLY
-- on an FTP reload, so create_resolver_views() rebuilds it only when missing or
-- force_ftp_core (which populate_from_ftp passes after every FTP swap). Every
-- additive/curated load reuses it as-is.

SET LOCAL work_mem = '4GB';
SET LOCAL max_parallel_workers_per_gather = 8;
SET LOCAL max_parallel_workers = 8;
SET LOCAL maintenance_work_mem = '4GB';
SET LOCAL max_parallel_maintenance_workers = 8;

CREATE SCHEMA IF NOT EXISTS omnipath_utils;

-- Drop whatever kind currently exists (view/matview from an earlier era, or the
-- table this file builds) — IF EXISTS does not span a relkind mismatch.
DO $$
DECLARE k "char";
BEGIN
  SELECT c.relkind INTO k FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
   WHERE n.nspname = 'omnipath_utils' AND c.relname = 'resolver_uniprot_taxon';
  IF k = 'v' THEN EXECUTE 'DROP VIEW omnipath_utils.resolver_uniprot_taxon CASCADE';
  ELSIF k = 'm' THEN EXECUTE 'DROP MATERIALIZED VIEW omnipath_utils.resolver_uniprot_taxon CASCADE';
  ELSIF k = 'r' THEN EXECUTE 'DROP TABLE omnipath_utils.resolver_uniprot_taxon CASCADE';
  END IF;
END $$;

-- One organism per accession. An accession belongs to exactly one organism, so
-- the aggregate collapses to a single row; min() is a deterministic tie-break for
-- the rare accession seen against more than one taxon. The accession-format filter
-- keeps only well-formed UniProtKB accessions (both the 6- and 10-character forms).
CREATE TABLE omnipath_utils.resolver_uniprot_taxon AS
SELECT m.source_id AS uniprot, min(m.ncbi_tax_id) AS ncbi_tax_id
FROM omnipath_utils.id_mapping_ftp m
JOIN omnipath_utils.id_type st ON m.source_type_id = st.id AND st.name = 'uniprot'
WHERE m.ncbi_tax_id IS NOT NULL
  AND m.ncbi_tax_id <> 0
  AND m.source_id ~ '^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
GROUP BY m.source_id;

-- Keyed-lookup index: the build ships this shard's accessions to a temp table and
-- probes by accession (organism is the payload).
CREATE INDEX resolver_uniprot_taxon_key_idx
    ON omnipath_utils.resolver_uniprot_taxon (uniprot)
    INCLUDE (ncbi_tax_id);

ANALYZE omnipath_utils.resolver_uniprot_taxon;
