SET max_parallel_workers_per_gather = 0;
SET max_parallel_workers = 0;

-- ============================================================
-- metalinksdb_protein_annotations
-- One row per canonical human protein.
-- Sources: UniProt (10), GuideToPharma (21), TCDB (41)
-- 004-metalinksdb-view: added uniprot_pathway_participations (UniProt
-- 'Pathway Participation:OM:0607' -- closes the protein-pathway gap: the
-- category was named in the annotation contract but no live view covered it).
-- ============================================================
-- One row per canonical human protein with functional annotations aggregated from UniProt, GuideToPharma, and TCDB.
DROP MATERIALIZED VIEW IF EXISTS metalinksdb_protein_annotations;
CREATE MATERIALIZED VIEW metalinksdb_protein_annotations AS
WITH
-- Aggregate UniProt subcellular location, function, disease involvement, EC numbers, and protein family annotations per canonical protein entity.
uniprot_ann AS (
    SELECT eer.entity_id,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Subcellular Location:OM:0604') AS uniprot_subcellular_locations,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Function:OM:0603')             AS uniprot_functions,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Disease Involvement:OM:0606')  AS uniprot_disease_involvements,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Ec Number:OM:0611')            AS uniprot_ec_numbers,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Protein Family:OM:0610')       AS uniprot_protein_families,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Pathway Participation:OM:0607') AS uniprot_pathway_participations
    FROM entity_evidence_annotation eea
    JOIN annotation a ON a.annotation_key = eea.annotation_key
    JOIN entity_evidence_resolution eer ON eer.source_id = eea.source_id
        AND eer.entity_evidence_id = eea.entity_evidence_id AND eer.status_id = 1
    WHERE eea.source_id = get_source_id('uniprot')
      AND a.term IN (
          'Subcellular Location:OM:0604',
          'Function:OM:0603',
          'Disease Involvement:OM:0606',
          'Ec Number:OM:0611',
          'Protein Family:OM:0610',
          'Pathway Participation:OM:0607'
      )
    GROUP BY eer.entity_id
),
-- Aggregate GuideToPharma functional class (GPCR, Transporter, Enzyme, etc.) and protein subfamily annotations per canonical protein entity.
gtp_ann AS (
    SELECT eer.entity_id,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Protein Functional Class:OM:0637') AS gtp_functional_classes,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Protein Family:OM:0610')           AS gtp_families
    FROM entity_evidence_annotation eea
    JOIN annotation a ON a.annotation_key = eea.annotation_key
    JOIN entity_evidence_resolution eer ON eer.source_id = eea.source_id
        AND eer.entity_evidence_id = eea.entity_evidence_id AND eer.status_id = 1
    WHERE eea.source_id = get_source_id('guidetopharma')
      AND a.term IN ('Protein Functional Class:OM:0637', 'Protein Family:OM:0610')
    GROUP BY eer.entity_id
),
-- Aggregate TCDB transporter family classifications per canonical protein entity.
tcdb_ann AS (
    SELECT eer.entity_id,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Protein Family:OM:0610') AS tcdb_transporter_families
    FROM entity_evidence_annotation eea
    JOIN annotation a ON a.annotation_key = eea.annotation_key
    JOIN entity_evidence_resolution eer ON eer.source_id = eea.source_id
        AND eer.entity_evidence_id = eea.entity_evidence_id AND eer.status_id = 1
    WHERE eea.source_id = get_source_id('tcdb')
    GROUP BY eer.entity_id
)
SELECT
    e.entity_id                     AS protein_entity_id,
    e.canonical_identifier          AS protein_uniprot,
    vit.name                        AS protein_canonical_id_type,
    e.taxonomy_id,
    u.uniprot_subcellular_locations,
    u.uniprot_functions,
    u.uniprot_disease_involvements,
    u.uniprot_ec_numbers,
    u.uniprot_protein_families,
    u.uniprot_pathway_participations,
    g.gtp_functional_classes,
    g.gtp_families,
    t.tcdb_transporter_families
FROM entity e
JOIN vocab_identifier_type vit ON vit.identifier_type_id = e.canonical_identifier_type_id
LEFT JOIN uniprot_ann u ON u.entity_id = e.entity_id
LEFT JOIN gtp_ann     g ON g.entity_id = e.entity_id
LEFT JOIN tcdb_ann    t ON t.entity_id = e.entity_id
WHERE e.taxonomy_id = 9606
  AND (u.entity_id IS NOT NULL
    OR g.entity_id IS NOT NULL
    OR t.entity_id IS NOT NULL);

CREATE UNIQUE INDEX ON metalinksdb_protein_annotations (protein_entity_id);
CREATE INDEX ON metalinksdb_protein_annotations (protein_uniprot);
CREATE INDEX ON metalinksdb_protein_annotations USING gin (gtp_functional_classes);
CREATE INDEX ON metalinksdb_protein_annotations USING gin (uniprot_ec_numbers);

-- ============================================================
-- metalinksdb_compound_annotations
-- One row per canonical compound.
-- Sources: HMDB (22), LipidMaps (25)
-- 004-metalinksdb-view: metabolite PATHWAY, DISEASE, STITCH SCORE, and
-- ENDOGENOUS/EXOGENOUS were checked against the live build and found NOT
-- addable at this entity level:
--   - pathway: only source is KEGG's 'Pathway:OM:0874'/'Pathway Component:OM:0315',
--     both attached to Reaction entities (a multi-hop, non-MetaLinksDB-source
--     traversal) -- documented gap, not implemented.
--   - disease: 'Disease Involvement:OM:0606' only resolves to Gene/protein
--     entities (already surfaced on the protein side above) -- no compound-level
--     source exists -- documented gap.
--   - STITCH score and endogenous/exogenous are RELATION-level (per-interaction)
--     annotations, not per-compound ones -- they are surfaced as
--     `best_stitch_score` / `is_endogenous` columns on the combined
--     `metalinksdb_relations` view instead (see metalinksdb.sql), not here.
-- ============================================================
-- One row per canonical compound with HMDB location/biospecimen data, LipidMaps lipid classification, and an is_lipid flag.
DROP MATERIALIZED VIEW IF EXISTS metalinksdb_compound_annotations;
CREATE MATERIALIZED VIEW metalinksdb_compound_annotations AS
WITH
-- Aggregate HMDB subcellular location, biospecimen, and tissue annotations per canonical compound entity.
hmdb_ann AS (
    SELECT eer.entity_id,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Subcellular Location:OM:0604') AS hmdb_subcellular_locations,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Biospecimen:OM:0779')           AS hmdb_biospecimens,
        array_agg(DISTINCT a.value) FILTER (WHERE a.term = 'Tissue:OM:0764')                AS hmdb_tissues
    FROM entity_evidence_annotation eea
    JOIN annotation a ON a.annotation_key = eea.annotation_key
    JOIN entity_evidence_resolution eer ON eer.source_id = eea.source_id
        AND eer.entity_evidence_id = eea.entity_evidence_id AND eer.status_id = 1
    WHERE eea.source_id = get_source_id('hmdb')
      AND a.term IN (
          'Subcellular Location:OM:0604',
          'Biospecimen:OM:0779',
          'Tissue:OM:0764'
      )
    GROUP BY eer.entity_id
),
-- Collect LipidMaps lipid category, main class, and sub class per canonical compound entity.
lipidmaps_ann AS (
    SELECT eer.entity_id,
        MAX(a.value) FILTER (WHERE a.term = 'Lipid Category:OM:0614')   AS lipid_category,
        MAX(a.value) FILTER (WHERE a.term = 'Lipid Main Class:OM:0615') AS lipid_main_class,
        MAX(a.value) FILTER (WHERE a.term = 'Lipid Sub Class:OM:0616')  AS lipid_sub_class
    FROM entity_evidence_annotation eea
    JOIN annotation a ON a.annotation_key = eea.annotation_key
    JOIN entity_evidence_resolution eer ON eer.source_id = eea.source_id
        AND eer.entity_evidence_id = eea.entity_evidence_id AND eer.status_id = 1
    WHERE eea.source_id = get_source_id('lipidmaps')
      AND a.term IN (
          'Lipid Category:OM:0614',
          'Lipid Main Class:OM:0615',
          'Lipid Sub Class:OM:0616'
      )
    GROUP BY eer.entity_id
)
SELECT
    e.entity_id                     AS compound_entity_id,
    e.canonical_identifier          AS compound_canonical_id,
    vit.name                        AS compound_canonical_id_type,
    h.hmdb_subcellular_locations,
    h.hmdb_biospecimens,
    h.hmdb_tissues,
    l.lipid_category,
    l.lipid_main_class,
    l.lipid_sub_class,
    (l.entity_id IS NOT NULL)       AS is_lipid
FROM entity e
JOIN vocab_identifier_type vit ON vit.identifier_type_id = e.canonical_identifier_type_id
LEFT JOIN hmdb_ann      h ON h.entity_id = e.entity_id
LEFT JOIN lipidmaps_ann l ON l.entity_id = e.entity_id
WHERE (h.entity_id IS NOT NULL OR l.entity_id IS NOT NULL);

CREATE UNIQUE INDEX ON metalinksdb_compound_annotations (compound_entity_id);
CREATE INDEX ON metalinksdb_compound_annotations (compound_canonical_id);
CREATE INDEX ON metalinksdb_compound_annotations (lipid_category);
CREATE INDEX ON metalinksdb_compound_annotations (lipid_main_class);
CREATE INDEX ON metalinksdb_compound_annotations USING gin (hmdb_biospecimens);
CREATE INDEX ON metalinksdb_compound_annotations USING gin (hmdb_subcellular_locations);
