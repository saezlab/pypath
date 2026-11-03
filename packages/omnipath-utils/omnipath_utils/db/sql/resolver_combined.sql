-- resolver_gene_protein_combined (spec 002; split out for independent rebuild,
-- 007 R10 / Phase 3P T067). Depends on the resolver_gene view (resolver_gene_ftp
-- + resolver_gene_curated). Read only by omnipath-build's pre-built parquet/duckdb
-- resolver mode, NOT the live-PG keyed path — so iteration against a live utils DB
-- can skip it (create_resolver_views resolvers=...).

-- Combined gene/protein resolver. It emits one canonical target per
-- (taxon, source_type, source_id): Entrez if uniquely derivable, otherwise the
-- primary UniProt if uniquely derivable. Ambiguous keys are omitted and remain
-- unresolved in consumers.
DO $$
DECLARE k "char";
BEGIN
  SELECT c.relkind INTO k FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
   WHERE n.nspname = 'omnipath_utils'
     AND c.relname = 'resolver_gene_protein_combined';
  IF k = 'v' THEN
    EXECUTE 'DROP VIEW omnipath_utils.resolver_gene_protein_combined CASCADE';
  ELSIF k = 'm' THEN
    EXECUTE 'DROP MATERIALIZED VIEW omnipath_utils.resolver_gene_protein_combined CASCADE';
  END IF;
END $$;
CREATE MATERIALIZED VIEW omnipath_utils.resolver_gene_protein_combined AS
WITH sec_pri AS (
    SELECT m.source_id AS sec, m.target_id AS pri
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id AND st.name = 'uniprot-sec'
    JOIN omnipath_utils.id_type tt
      ON m.target_type_id = tt.id AND tt.name = 'uniprot-pri'
),
protein_source_raw AS (
    -- full-UniProt: native uniprot -> X, inverted to X -> uniprot
    SELECT
      m.ncbi_tax_id,
      tt.name AS source_type,
      m.target_id AS source_id,
      m.source_id AS ac
    FROM omnipath_utils.id_mapping_ftp m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id AND st.name = 'uniprot'
    JOIN omnipath_utils.id_type tt
      ON m.target_type_id = tt.id
     AND tt.name IN (
       'genesymbol', 'ensg', 'ensp', 'entrez', 'refseqp', 'uniprot_entry'
     )
    UNION ALL
    -- curated, native uniprot -> X (inverted)
    SELECT
      m.ncbi_tax_id,
      tt.name AS source_type,
      m.target_id AS source_id,
      m.source_id AS ac
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id AND st.name = 'uniprot'
    JOIN omnipath_utils.id_type tt
      ON m.target_type_id = tt.id
     AND tt.name IN ('genesymbol', 'ensg', 'ensp', 'entrez', 'refseqp')
    UNION ALL
    -- curated, forward X -> uniprot
    SELECT
      m.ncbi_tax_id,
      st.name AS source_type,
      m.source_id AS source_id,
      m.target_id AS ac
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id
     AND st.name IN ('genesymbol', 'ensg', 'ensp', 'entrez', 'refseqp')
    JOIN omnipath_utils.id_type tt
      ON m.target_type_id = tt.id AND tt.name = 'uniprot'
    UNION ALL
    -- primary UniProt accession identity
    SELECT DISTINCT
      m.ncbi_tax_id,
      'uniprot' AS source_type,
      m.source_id AS source_id,
      m.source_id AS ac
    FROM omnipath_utils.id_mapping_ftp m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id AND st.name = 'uniprot'
    WHERE m.source_id IS NOT NULL
    UNION ALL
    -- secondary UniProt accession -> primary UniProt accession
    SELECT
      NULLIF(m.ncbi_tax_id, 0) AS ncbi_tax_id,
      'uniprot' AS source_type,
      m.source_id AS source_id,
      m.target_id AS ac
    FROM omnipath_utils.id_mapping m
    JOIN omnipath_utils.id_type st
      ON m.source_type_id = st.id AND st.name = 'uniprot-sec'
    JOIN omnipath_utils.id_type tt
      ON m.target_type_id = tt.id AND tt.name = 'uniprot-pri'
    WHERE m.source_id IS NOT NULL
      AND m.target_id IS NOT NULL
),
protein_key_normalized AS (
    SELECT
      r.ncbi_tax_id,
      r.source_type,
      r.source_id,
      COALESCE(sp.pri, r.ac) AS primary_uniprot
    FROM protein_source_raw r
    LEFT JOIN sec_pri sp ON sp.sec = r.ac
    WHERE r.source_id IS NOT NULL
      AND r.ac IS NOT NULL
),
protein_key_flagged AS (
    SELECT
      pk.ncbi_tax_id,
      pk.source_type,
      pk.source_id,
      pk.primary_uniprot,
      (r.identifier IS NOT NULL) AS is_swissprot,
      bool_or(r.identifier IS NOT NULL) OVER (
        PARTITION BY pk.ncbi_tax_id, pk.source_type, pk.source_id
      ) AS key_has_swissprot
    FROM protein_key_normalized pk
    LEFT JOIN omnipath_utils.reflist r
      ON r.list_name = 'swissprot'
     AND r.identifier = pk.primary_uniprot
),
protein_key AS (
    SELECT DISTINCT
      ncbi_tax_id,
      source_type,
      source_id,
      primary_uniprot
    FROM protein_key_flagged
    WHERE ((key_has_swissprot AND is_swissprot) OR NOT key_has_swissprot)
      AND primary_uniprot ~ '^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
),
entrez_candidate AS (
    SELECT
      rg.ncbi_tax_id,
      rg.source_type,
      rg.source_id,
      rg.entrez
    FROM omnipath_utils.resolver_gene rg
    WHERE rg.source_id IS NOT NULL
      AND rg.entrez IS NOT NULL
    UNION
    -- If a source key has both a gene anchor and a primary UniProt, the primary
    -- UniProt should resolve to the same gene (e.g. Cngb1 -> A0A8I5ZN27 and
    -- Cngb1 -> Entrez).
    SELECT
      pk.ncbi_tax_id,
      'uniprot' AS source_type,
      pk.primary_uniprot AS source_id,
      rg.entrez
    FROM protein_key pk
    JOIN omnipath_utils.resolver_gene rg
      ON rg.ncbi_tax_id = pk.ncbi_tax_id
     AND rg.source_type = pk.source_type
     AND rg.source_id = pk.source_id
    WHERE pk.primary_uniprot IS NOT NULL
      AND rg.entrez IS NOT NULL
    UNION
    -- Secondary accessions inherit the primary accession's gene anchor.
    SELECT
      rg.ncbi_tax_id,
      'uniprot' AS source_type,
      sp.sec AS source_id,
      rg.entrez
    FROM sec_pri sp
    JOIN omnipath_utils.resolver_gene rg
      ON rg.source_type = 'uniprot'
     AND rg.source_id = sp.pri
    WHERE sp.sec IS NOT NULL
),
entrez_unique AS (
    SELECT
      ncbi_tax_id,
      source_type,
      source_id,
      min(entrez) AS canonical_id
    FROM entrez_candidate
    WHERE source_id IS NOT NULL
      AND entrez IS NOT NULL
    GROUP BY ncbi_tax_id, source_type, source_id
    HAVING count(DISTINCT entrez) = 1
),
uniprot_unique AS (
    SELECT
      pk.ncbi_tax_id,
      pk.source_type,
      pk.source_id,
      min(pk.primary_uniprot) AS canonical_id
    FROM protein_key pk
    LEFT JOIN entrez_unique eu
      ON eu.ncbi_tax_id IS NOT DISTINCT FROM pk.ncbi_tax_id
     AND eu.source_type = pk.source_type
     AND eu.source_id = pk.source_id
    WHERE eu.source_id IS NULL
      AND pk.source_id IS NOT NULL
      AND pk.primary_uniprot IS NOT NULL
    GROUP BY pk.ncbi_tax_id, pk.source_type, pk.source_id
    HAVING count(DISTINCT pk.primary_uniprot) = 1
)
SELECT
  ncbi_tax_id,
  source_type,
  source_id,
  'entrez' AS canonical_type,
  canonical_id
FROM entrez_unique
UNION ALL
SELECT
  ncbi_tax_id,
  source_type,
  source_id,
  'uniprot' AS canonical_type,
  canonical_id
FROM uniprot_unique;

CREATE INDEX IF NOT EXISTS resolver_gene_protein_combined_key_idx
    ON omnipath_utils.resolver_gene_protein_combined
    (ncbi_tax_id, source_type, source_id);
CREATE INDEX IF NOT EXISTS resolver_gene_protein_combined_target_idx
    ON omnipath_utils.resolver_gene_protein_combined
    (canonical_type, canonical_id, ncbi_tax_id);
