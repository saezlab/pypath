-- omnipath-build/custom_views/metalinksdb.sql
-- MetaLinksDB materialized views for each interaction source.
-- Apply with:  psql -U omnipath -d omnipath -f metalinksdb.sql
-- Refresh after each:  make load SOURCE=<source> && make derive
--   then: REFRESH MATERIALIZED VIEW CONCURRENTLY metalinksdb_<source>_relations;
--
-- Human-only filter: all views use a human_re CTE to pre-filter to human rows
-- before building identifier and annotation CTEs. This avoids processing the
-- full source dataset through expensive aggregations before discarding non-human rows.
-- For ChEMBL the human filter comes from annotation term 'Ncbi Tax Id:OM:0205'.
-- For all other sources it comes from entity_evidence_resolution -> entity.taxonomy_id = 9606.
--
-- Source IDs (as of omnipath-build c39dd5a rebuild 2026-05-28):
--   chembl=16  bindingdb=12  cellinker=14  guidetopharma=21
--   mrclinksdb=28  stitch=39  tcdb=41
--
-- Entity type IDs: SmallMolecule=2  Protein=3  Complex=4  Reaction=11  Transport=12
--
-- IMPORTANT (gene-anchored model, fixed 2026-06-15): the protein-side filter on
-- entity_evidence is `Protein:MI:0326` — entity_evidence carries the ORIGINAL
-- mention type; only the *resolved/canonical* entity becomes `Gene:MI:0250`. The
-- gene anchoring + human filter are enforced by the entity_evidence_resolution →
-- entity (taxonomy_id=9606) join, not by the evidence type. A previous version
-- filtered evidence on Gene:MI:0250, which matched 0 rows → every non-ChEMBL
-- per-source view was empty and the combined view was ChEMBL-only.

-- ────────────────────────────────────────────────────────────────────────────
-- ChEMBL  (source_id = 16)
-- Human-only via ncbi_tax_id annotation. human_re pre-filters ~3.4M rows to
-- ~300k human rows before all subsequent CTEs run.
-- ────────────────────────────────────────────────────────────────────────────

SET max_parallel_workers_per_gather = 4;
SET max_parallel_workers = 8;
SET work_mem = '256MB';

-- Drop combined view first so individual source views can be dropped without CASCADE

-- Helper function: resolves a source name to its source_id.
-- Declared STABLE so PostgreSQL evaluates it at plan time,
-- preserving partition pruning for all partitioned tables.
CREATE OR REPLACE FUNCTION get_source_id(source_name text)
RETURNS bigint LANGUAGE plpgsql STABLE AS $$
BEGIN
    RETURN (SELECT source_id FROM public.data_source WHERE name = source_name);
END;
$$;



CREATE OR REPLACE FUNCTION get_entity_type_id(type_name text)
RETURNS bigint LANGUAGE plpgsql STABLE AS $$
BEGIN
    RETURN (SELECT entity_type_id FROM public.vocab_entity_type WHERE name = type_name);
END;
$$;

CREATE OR REPLACE FUNCTION get_identifier_type_id(type_name text)
RETURNS bigint LANGUAGE plpgsql STABLE AS $$
BEGIN
    RETURN (SELECT identifier_type_id FROM public.vocab_identifier_type WHERE name = type_name);
END;
$$;

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_relations;

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_chembl_relations;
CREATE MATERIALIZED VIEW metalinksdb_chembl_relations AS
WITH

-- CURATION (decided 2026-06-15): keep only the CANONICAL drug-target set =
-- ChEMBL **mechanism-of-action** pairs (`Chembl Mechanism:OM:0227` annotation),
-- ~6,940 pairs. NON-HUMAN INCLUDED (no ncbi_tax_id gate).
-- Rationale: pChEMBL affinity thresholds flood the view (>6 = 1.6M, even >9 = 127k),
-- and MoA + pChEMBL sit in different ChEMBL tables (their AND = 0) — so the curated
-- contribution is the MoA mechanism set, not affinity-filtered assays. (CTE name
-- kept `human_re` for minimal-diff — it is now the curated MoA set, not human-only.)
human_re AS (
    SELECT DISTINCT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN relation_evidence_annotation rea_m
        ON  rea_m.relation_evidence_id = re.relation_evidence_id
        AND rea_m.source_id            = re.source_id
    JOIN annotation am ON am.annotation_key = rea_m.annotation_key
        AND am.term = 'Chembl Mechanism:OM:0227'
    WHERE re.source_id = get_source_id('chembl')
),

-- Collect raw ChEMBL compound identifiers (ChEMBL ID, InChIKey, SMILES, name) per entity_evidence row.
compound_ids AS (
    SELECT
        eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Chembl Compound:MI:0967') THEN ie.value END) AS chembl_compound_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Standard Inchi Key:MI:1101') THEN ie.value END) AS inchikey,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Smiles:MI:0239') THEN ie.value END) AS smiles,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Name:OM:0202') THEN ie.value END) AS compound_name
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('chembl')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw ChEMBL protein identifiers (UniProt, ChEMBL target ID) per entity_evidence row.
protein_ids AS (
    SELECT
        eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097')  THEN ie.value END) AS uniprot_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Chembl Target:MI:1348') THEN ie.value END) AS chembl_target_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('chembl')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Pivot ChEMBL activity annotations (pChEMBL, IC50, assay type, literature) into columns per relation_evidence row.
rel_annotations AS NOT MATERIALIZED (
    SELECT
        rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Chembl Activity:OM:0228'   THEN a.value END) AS chembl_activity_id,
        MAX(CASE WHEN a.term = 'Chembl Assay:OM:0225'      THEN a.value END) AS chembl_assay_id,
        MAX(CASE WHEN a.term = 'Chembl Document:OM:0226'   THEN a.value END) AS chembl_document_id,
        MAX(CASE WHEN a.term = 'Chembl Target:MI:1348'     THEN a.value END) AS chembl_target_ann,
        MAX(CASE WHEN a.term = 'Pchembl Value:OM:0708'     THEN a.value END) AS pchembl_value,
        MAX(CASE WHEN a.term = 'Confidence Score:OM:0762'  THEN a.value END) AS confidence_score,
        MAX(CASE WHEN a.term = 'Ic50:MI:0641'              THEN a.value END) AS ic50,
        MAX(CASE WHEN a.term = 'Ki:MI:0643'                THEN a.value END) AS ki,
        MAX(CASE WHEN a.term = 'Kd:MI:0646'                THEN a.value END) AS kd,
        MAX(CASE WHEN a.term = 'Ec50:MI:0642'              THEN a.value END) AS ec50,
        MAX(CASE WHEN a.term = 'Assay Category:OM:0761'    THEN a.value END) AS assay_category,
        MAX(CASE WHEN a.term = 'Cell Type:OM:0765'         THEN a.value END) AS cell_type,
        MAX(CASE WHEN a.term = 'Tissue:OM:0764'            THEN a.value END) AS tissue,
        MAX(CASE WHEN a.term = 'Ncbi Tax Id:OM:0205'       THEN a.value END) AS ncbi_tax_id,
        MAX(CASE WHEN a.term = 'Description:OM:0613'       THEN a.value END) AS description,
        MAX(CASE WHEN a.term = 'Chembl Mechanism:OM:0227'  THEN a.value END) AS chembl_mechanism,
        MAX(CASE WHEN a.term = 'Pubmed:MI:0446'            THEN a.value END) AS pubmed_id,
        MAX(CASE WHEN a.term = 'Doi:MI:0574'               THEN a.value END) AS doi,
        BOOL_OR(a.term = 'Binding:OM:0751')                                  AS is_binding_assay,
        BOOL_OR(a.term = 'Functional:OM:0752')                               AS is_functional_assay,
        BOOL_OR(a.term = 'Adme:OM:0753')                                     AS is_adme_assay,
        BOOL_OR(a.term = 'Toxicity:OM:0754')                                 AS is_toxicity_assay,
        BOOL_OR(a.term = 'Inhibition:OM:0931')                               AS is_inhibition,
        BOOL_OR(a.term = 'Activation:OM:0930')                               AS is_activation,
        BOOL_OR(a.term = 'Agonist:OM:0901')                                  AS is_agonist,
        BOOL_OR(a.term = 'Antagonist:OM:0920')                               AS is_antagonist
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('chembl')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'chembl'::text                              AS source,
    re.relation_evidence_id,

    e_compound.entity_id                        AS compound_entity_id,
    e_compound.canonical_identifier             AS compound_canonical_id,
    vit_c.name                                  AS compound_canonical_id_type,
    eer_compound.status_id                      AS compound_resolution_status,

    ci.chembl_compound_id,
    ci.inchikey                                 AS compound_inchikey,
    ci.smiles                                   AS compound_smiles,
    ci.compound_name,

    e_protein.entity_id                         AS protein_entity_id,
    e_protein.canonical_identifier              AS protein_canonical_id,
    vit_p.name                                  AS protein_canonical_id_type,
    eer_protein.status_id                       AS protein_resolution_status,

    pi.uniprot_id                               AS protein_uniprot,
    pi.chembl_target_id                         AS protein_chembl_target,

    vrp.name                                    AS predicate,
    vrc.name                                    AS relation_category,

    ann.chembl_activity_id,
    ann.chembl_assay_id,
    ann.chembl_document_id,
    ann.chembl_target_ann,

    ann.pchembl_value::numeric                  AS pchembl_value,
    ann.confidence_score::numeric               AS confidence_score,
    ann.ic50::numeric                           AS ic50,
    ann.ki::numeric                             AS ki,
    ann.kd::numeric                             AS kd,
    ann.ec50::numeric                           AS ec50,

    ann.assay_category,
    ann.is_binding_assay,
    ann.is_functional_assay,
    ann.is_adme_assay,
    ann.is_toxicity_assay,
    ann.is_inhibition,
    ann.is_activation,
    ann.is_agonist,
    ann.is_antagonist,

    ann.cell_type,
    ann.tissue,
    ann.ncbi_tax_id,
    ann.description,
    ann.chembl_mechanism,

    ann.pubmed_id,
    ann.doi

FROM human_re re

JOIN entity_evidence ee_compound
    ON  ee_compound.source_id = get_source_id('chembl')
    AND ee_compound.entity_evidence_id = re.subject_entity_evidence_id

JOIN entity_evidence ee_protein
    ON  ee_protein.source_id = get_source_id('chembl')
    AND ee_protein.entity_evidence_id = re.object_entity_evidence_id

LEFT JOIN compound_ids ci
    ON ci.entity_evidence_id = re.subject_entity_evidence_id

LEFT JOIN protein_ids pi
    ON pi.entity_evidence_id = re.object_entity_evidence_id

LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('chembl')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound
    ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id

LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('chembl')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein
    ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id

LEFT JOIN rel_annotations ann
    ON ann.relation_evidence_id = re.relation_evidence_id

JOIN vocab_relation_predicate vrp
    ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category vrc
    ON vrc.relation_category_id = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_chembl_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_chembl_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_chembl_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_chembl_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_chembl_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- BindingDB  (source_id = 12)
-- Human-only via entity resolution. human_re pre-filters ~2.4M rows
-- (BindingDB is mostly human, so the win here is smaller than ChEMBL).
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_bindingdb_relations;
CREATE MATERIALIZED VIEW metalinksdb_bindingdb_relations AS
WITH

-- CURATION: keep only high-affinity BindingDB pairs (pChEMBL >= 6) whose
-- protein resolves to a gene. NON-HUMAN INCLUDED (no taxonomy_id gate). Uses
-- the Pchembl Value annotation (present on ~707k rows);
-- TODO: also admit rows that lack pChEMBL but carry Ic50/Ki/Kd/Ec50 in
-- nM via the derived pChEMBL = 9 - log10(nM) (censored/non-positive excluded).
-- (CTE name kept `human_re` for minimal-diff — now the curated set.)
human_re AS (
    SELECT DISTINCT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence_resolution eer
        ON  eer.source_id          = re.source_id
        AND eer.entity_evidence_id = re.object_entity_evidence_id
        AND eer.status_id          = 1
    JOIN relation_evidence_annotation rea_p
        ON  rea_p.relation_evidence_id = re.relation_evidence_id
        AND rea_p.source_id            = re.source_id
    JOIN annotation ap ON ap.annotation_key = rea_p.annotation_key
        AND ap.term  = 'Pchembl Value:OM:0708'
        AND ap.value ~ '^[0-9]+(\.[0-9]+)?$'
        AND ap.value::numeric >= 6
    WHERE re.source_id = get_source_id('bindingdb')
),

-- Collect raw BindingDB compound identifiers (BindingDB ID, InChIKey, SMILES, PubChem CID, ChEMBL ID) per entity_evidence row.
compound_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Bindingdb:OM:0006') THEN ie.value END) AS bindingdb_compound_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Standard Inchi Key:MI:1101') THEN ie.value END) AS inchikey,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Smiles:MI:0239') THEN ie.value END) AS smiles,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Pubchem Compound:OM:0002') THEN ie.value END) AS pubchem_cid,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Chembl Compound:MI:0967') THEN ie.value END) AS chembl_compound_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('bindingdb')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw BindingDB protein identifier (UniProt) per entity_evidence row.
protein_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097') THEN ie.value END) AS uniprot_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('bindingdb')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Pivot BindingDB affinity measurements (pChEMBL, IC50, Ki, Kd, EC50) and experimental conditions into columns per relation_evidence row.
rel_annotations AS NOT MATERIALIZED (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Pchembl Value:OM:0708'        THEN a.value END) AS pchembl_value,
        MAX(CASE WHEN a.term = 'Ic50:MI:0641'                 THEN a.value END) AS ic50,
        MAX(CASE WHEN a.term = 'Ki:MI:0643'                   THEN a.value END) AS ki,
        MAX(CASE WHEN a.term = 'Kd:MI:0646'                   THEN a.value END) AS kd,
        MAX(CASE WHEN a.term = 'Ec50:MI:0642'                 THEN a.value END) AS ec50,
        MAX(CASE WHEN a.term = 'Pubmed:MI:0446'               THEN a.value END) AS pubmed_id,
        MAX(CASE WHEN a.term = 'Doi:MI:0574'                  THEN a.value END) AS doi,
        MAX(CASE WHEN a.term = 'Patent Number:OM:0206'        THEN a.value END) AS patent_number,
        MAX(CASE WHEN a.term = 'Comment:MI:0612'              THEN a.value END) AS comment,
        MAX(CASE WHEN a.term = 'Ph:MI:0837'                   THEN a.value END) AS ph,
        MAX(CASE WHEN a.term = 'Temperature Celsius:OM:0701'  THEN a.value END) AS temperature_celsius
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('bindingdb')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'bindingdb'::text               AS source,
    re.relation_evidence_id,
    e_compound.entity_id            AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name                      AS compound_canonical_id_type,
    eer_compound.status_id          AS compound_resolution_status,
    ci.bindingdb_compound_id,
    ci.inchikey                     AS compound_inchikey,
    ci.smiles                       AS compound_smiles,
    ci.pubchem_cid                  AS compound_pubchem_cid,
    ci.chembl_compound_id           AS compound_chembl_id,
    e_protein.entity_id             AS protein_entity_id,
    e_protein.canonical_identifier  AS protein_canonical_id,
    vit_p.name                      AS protein_canonical_id_type,
    eer_protein.status_id           AS protein_resolution_status,
    pi.uniprot_id                   AS protein_uniprot,
    vrp.name                        AS predicate,
    vrc.name                        AS relation_category,
    NULLIF(regexp_replace(ann.pchembl_value, '[^0-9.eE+-]', '', 'g'), '')::numeric AS pchembl_value,
    NULLIF(regexp_replace(ann.ic50,          '[^0-9.eE+-]', '', 'g'), '')::numeric AS ic50,
    NULLIF(regexp_replace(ann.ki,            '[^0-9.eE+-]', '', 'g'), '')::numeric AS ki,
    NULLIF(regexp_replace(ann.kd,            '[^0-9.eE+-]', '', 'g'), '')::numeric AS kd,
    NULLIF(regexp_replace(ann.ec50,          '[^0-9.eE+-]', '', 'g'), '')::numeric AS ec50,
    ann.pubmed_id,
    ann.doi,
    ann.patent_number,
    ann.comment,
    NULLIF(regexp_replace(ann.ph,            '[^0-9.eE+-]', '', 'g'), '')::numeric AS ph,
    NULLIF(split_part(ann.temperature_celsius, ' ', 1), '')::numeric AS temperature_celsius

FROM human_re re

LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('bindingdb')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('bindingdb')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = re.relation_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_bindingdb_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_bindingdb_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_bindingdb_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_bindingdb_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_bindingdb_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- CellLinker  (source_id = 14)  — metabolite-protein pairs only (SM=2→Prot=3)
-- Small source (~5884 human rows); human_re applied for consistency.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_cellinker_relations;
CREATE MATERIALIZED VIEW metalinksdb_cellinker_relations AS
WITH

-- Pre-filter to human SM→Protein rows: compound entity_type_id=2 checked to
-- exclude protein-protein rows; protein resolves to a human entity.
human_re AS (
    SELECT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence ee_compound
        ON  ee_compound.source_id          = re.source_id
        AND ee_compound.entity_evidence_id = re.subject_entity_evidence_id
        AND ee_compound.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_protein
        ON  ee_protein.source_id          = re.source_id
        AND ee_protein.entity_evidence_id = re.object_entity_evidence_id
        AND ee_protein.entity_type_id = get_entity_type_id('Protein:MI:0326')
    JOIN entity_evidence_resolution eer
        ON  eer.source_id          = re.source_id
        AND eer.entity_evidence_id = re.object_entity_evidence_id
        AND eer.status_id          = 1
    JOIN entity e ON e.entity_id = eer.entity_id AND e.taxonomy_id = 9606
    WHERE re.source_id = get_source_id('cellinker')
),

-- Collect raw CellLinker compound identifiers (HMDB ID, PubChem CID, SMILES) per entity_evidence row.
compound_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Hmdb:OM:0004')  THEN ie.value END) AS hmdb_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Pubchem Compound:OM:0002') THEN ie.value END) AS pubchem_cid,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Smiles:MI:0239') THEN ie.value END) AS smiles
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('cellinker')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw CellLinker protein identifiers (UniProt, Entrez, gene name) per entity_evidence row.
protein_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097') THEN ie.value END) AS uniprot_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Entrez:MI:0477') THEN ie.value END) AS entrez_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Gene Name Primary:OM:0200') THEN ie.value END) AS gene_name
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('cellinker')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect CellLinker literature references and interaction annotations per relation_evidence row.
rel_annotations AS (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Pubmed:MI:0446'                  THEN a.value END) AS pubmed_id,
        MAX(CASE WHEN a.term = 'Interaction Xref:OM:1206'        THEN a.value END) AS interaction_xref,
        MAX(CASE WHEN a.term = 'Interaction Annotation:OM:1207'  THEN a.value END) AS interaction_annotation
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('cellinker')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'cellinker'::text               AS source,
    re.relation_evidence_id,
    e_compound.entity_id            AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name                      AS compound_canonical_id_type,
    eer_compound.status_id          AS compound_resolution_status,
    ci.hmdb_id                      AS compound_hmdb_id,
    ci.pubchem_cid                  AS compound_pubchem_cid,
    ci.smiles                       AS compound_smiles,
    e_protein.entity_id             AS protein_entity_id,
    e_protein.canonical_identifier  AS protein_canonical_id,
    vit_p.name                      AS protein_canonical_id_type,
    eer_protein.status_id           AS protein_resolution_status,
    pi.uniprot_id                   AS protein_uniprot,
    pi.entrez_id                    AS protein_entrez,
    pi.gene_name                    AS protein_gene_name,
    vrp.name                        AS predicate,
    vrc.name                        AS relation_category,
    ann.pubmed_id,
    ann.interaction_xref,
    ann.interaction_annotation

FROM human_re re
LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('cellinker')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('cellinker')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = re.relation_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_cellinker_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_cellinker_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_cellinker_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_cellinker_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_cellinker_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- Guide to Pharmacology  (source_id = 21)
-- Small source (~23904 rows); human_re applied for consistency.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_guidetopharma_relations;
CREATE MATERIALIZED VIEW metalinksdb_guidetopharma_relations AS
WITH

-- Pre-filter to human-only relation_evidence rows via protein entity resolution.
human_re AS (
    SELECT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence ee_protein
        ON  ee_protein.source_id          = re.source_id
        AND ee_protein.entity_evidence_id = re.object_entity_evidence_id
        AND ee_protein.entity_type_id = get_entity_type_id('Protein:MI:0326')
    JOIN entity_evidence_resolution eer
        ON  eer.source_id          = re.source_id
        AND eer.entity_evidence_id = re.object_entity_evidence_id
        AND eer.status_id          = 1
    JOIN entity e ON e.entity_id = eer.entity_id AND e.taxonomy_id = 9606
    WHERE re.source_id = get_source_id('guidetopharma')
),

-- Collect raw GuideToPharma compound identifiers (GtP ligand ID, InChIKey, SMILES, PubChem CID, ChEMBL ID) per entity_evidence row.
compound_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Guidetopharma:OM:0008') THEN ie.value END) AS guidetopharma_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Standard Inchi Key:MI:1101') THEN ie.value END) AS inchikey,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Smiles:MI:0239') THEN ie.value END) AS smiles,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Pubchem Compound:OM:0002') THEN ie.value END) AS pubchem_cid,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Chembl Compound:MI:0967') THEN ie.value END) AS chembl_compound_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('guidetopharma')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw GuideToPharma protein identifiers (UniProt, HGNC, gene name) per entity_evidence row.
protein_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097') THEN ie.value END) AS uniprot_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Hgnc:MI:1095') THEN ie.value END) AS hgnc_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Gene Name Primary:OM:0200') THEN ie.value END) AS gene_name
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('guidetopharma')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Pivot GuideToPharma pharmacological annotations (affinity values, endogenous flag, agonist/antagonist/modulator flags) into columns per relation_evidence row.
rel_annotations AS (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Pubmed:MI:0446'                THEN a.value END) AS pubmed_id,
        MAX(CASE WHEN a.term = 'Affinity Median:OM:0628'       THEN a.value END) AS affinity_median,
        MAX(CASE WHEN a.term = 'Affinity Low:OM:0627'          THEN a.value END) AS affinity_low,
        MAX(CASE WHEN a.term = 'Affinity High:OM:0626'         THEN a.value END) AS affinity_high,
        MAX(CASE WHEN a.term = 'Ncbi Tax Id:OM:0205'           THEN a.value END) AS ncbi_tax_id,
        BOOL_OR(a.term = 'Endogenous:OM:0625')                                   AS is_endogenous,
        BOOL_OR(a.term = 'Inhibitor:OM:1004')                                    AS is_inhibitor,
        BOOL_OR(a.term = 'Inhibition:OM:0931')                                   AS is_inhibition,
        BOOL_OR(a.term = 'Agonist:OM:1001')                                      AS is_agonist,
        BOOL_OR(a.term = 'Antagonist:OM:1002')                                   AS is_antagonist,
        BOOL_OR(a.term = 'Full Agonist:OM:0902')                                 AS is_full_agonist,
        BOOL_OR(a.term = 'Partial Agonist:OM:0903')                              AS is_partial_agonist,
        BOOL_OR(a.term = 'Inverse Agonist:OM:0904')                              AS is_inverse_agonist,
        BOOL_OR(a.term = 'Activator:OM:1003')                                    AS is_activator,
        BOOL_OR(a.term = 'Binding:OM:0980')                                      AS is_binding,
        BOOL_OR(a.term = 'Allosteric Modulator:OM:1005')                         AS is_allosteric_modulator,
        BOOL_OR(a.term = 'Channel Blocker:OM:1020')                              AS is_channel_blocker
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('guidetopharma')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'guidetopharma'::text           AS source,
    re.relation_evidence_id,
    e_compound.entity_id            AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name                      AS compound_canonical_id_type,
    eer_compound.status_id          AS compound_resolution_status,
    ci.guidetopharma_id             AS compound_guidetopharma_id,
    ci.inchikey                     AS compound_inchikey,
    ci.smiles                       AS compound_smiles,
    ci.pubchem_cid                  AS compound_pubchem_cid,
    ci.chembl_compound_id           AS compound_chembl_id,
    e_protein.entity_id             AS protein_entity_id,
    e_protein.canonical_identifier  AS protein_canonical_id,
    vit_p.name                      AS protein_canonical_id_type,
    eer_protein.status_id           AS protein_resolution_status,
    pi.uniprot_id                   AS protein_uniprot,
    pi.hgnc_id                      AS protein_hgnc,
    pi.gene_name                    AS protein_gene_name,
    vrp.name                        AS predicate,
    vrc.name                        AS relation_category,
    ann.pubmed_id,
    ann.affinity_median::numeric    AS affinity_median,
    ann.affinity_low::numeric       AS affinity_low,
    ann.affinity_high::numeric      AS affinity_high,
    ann.ncbi_tax_id,
    ann.is_endogenous,
    ann.is_inhibitor,
    ann.is_inhibition,
    ann.is_agonist,
    ann.is_antagonist,
    ann.is_full_agonist,
    ann.is_partial_agonist,
    ann.is_inverse_agonist,
    ann.is_activator,
    ann.is_binding,
    ann.is_allosteric_modulator,
    ann.is_channel_blocker

FROM human_re re
LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('guidetopharma')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('guidetopharma')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = re.relation_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_guidetopharma_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_guidetopharma_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_guidetopharma_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_guidetopharma_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_guidetopharma_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- MRCLinksDB  (source_id = 28)  — inherently human, small source (~1468 rows).
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_mrclinksdb_relations;
CREATE MATERIALIZED VIEW metalinksdb_mrclinksdb_relations AS
WITH

-- Pre-filter to human SM→Protein rows: compound entity_type_id=2 checked to
-- exclude any non-metabolite rows.
human_re AS (
    SELECT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence ee_compound
        ON  ee_compound.source_id          = re.source_id
        AND ee_compound.entity_evidence_id = re.subject_entity_evidence_id
        AND ee_compound.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_protein
        ON  ee_protein.source_id          = re.source_id
        AND ee_protein.entity_evidence_id = re.object_entity_evidence_id
        AND ee_protein.entity_type_id = get_entity_type_id('Protein:MI:0326')
    JOIN entity_evidence_resolution eer
        ON  eer.source_id          = re.source_id
        AND eer.entity_evidence_id = re.object_entity_evidence_id
        AND eer.status_id          = 1
    JOIN entity e ON e.entity_id = eer.entity_id AND e.taxonomy_id = 9606
    WHERE re.source_id = get_source_id('mrclinksdb')
),

-- Collect raw MRCLinksDB compound identifiers (HMDB ID, PubChem CID, SMILES) per entity_evidence row.
compound_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Hmdb:OM:0004')  THEN ie.value END) AS hmdb_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Pubchem Compound:OM:0002') THEN ie.value END) AS pubchem_cid,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Smiles:MI:0239') THEN ie.value END) AS smiles
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('mrclinksdb')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw MRCLinksDB protein identifiers (UniProt, Entrez, gene name) per entity_evidence row.
protein_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097') THEN ie.value END) AS uniprot_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Entrez:MI:0477') THEN ie.value END) AS entrez_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Gene Name Primary:OM:0200') THEN ie.value END) AS gene_name
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('mrclinksdb')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect MRCLinksDB literature references and curator comments per relation_evidence row.
rel_annotations AS (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Pubmed:MI:0446'   THEN a.value END) AS pubmed_id,
        MAX(CASE WHEN a.term = 'Comment:MI:0612'  THEN a.value END) AS comment
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('mrclinksdb')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'mrclinksdb'::text              AS source,
    re.relation_evidence_id,
    e_compound.entity_id            AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name                      AS compound_canonical_id_type,
    eer_compound.status_id          AS compound_resolution_status,
    ci.hmdb_id                      AS compound_hmdb_id,
    ci.pubchem_cid                  AS compound_pubchem_cid,
    ci.smiles                       AS compound_smiles,
    e_protein.entity_id             AS protein_entity_id,
    e_protein.canonical_identifier  AS protein_canonical_id,
    vit_p.name                      AS protein_canonical_id_type,
    eer_protein.status_id           AS protein_resolution_status,
    pi.uniprot_id                   AS protein_uniprot,
    pi.entrez_id                    AS protein_entrez,
    pi.gene_name                    AS protein_gene_name,
    vrp.name                        AS predicate,
    vrc.name                        AS relation_category,
    ann.pubmed_id,
    ann.comment

FROM human_re re
LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('mrclinksdb')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('mrclinksdb')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = re.relation_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_mrclinksdb_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_mrclinksdb_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_mrclinksdb_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_mrclinksdb_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_mrclinksdb_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- STITCH  (source_id = 39)
-- Human filter uses entity_evidence.taxonomy_id = 9606 directly — simpler than
-- identifier-based approach and gives the planner correct cardinality estimates.
DROP MATERIALIZED VIEW IF EXISTS metalinksdb_stitch_relations;
CREATE MATERIALIZED VIEW metalinksdb_stitch_relations AS
WITH

-- Pre-filter to human SM→Protein rows using taxonomy_id directly on entity_evidence.
human_re AS (
    SELECT DISTINCT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence ee_compound
        ON  ee_compound.source_id          = re.source_id
        AND ee_compound.entity_evidence_id = re.subject_entity_evidence_id
        AND ee_compound.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_protein
        ON  ee_protein.source_id          = re.source_id
        AND ee_protein.entity_evidence_id = re.object_entity_evidence_id
        AND ee_protein.entity_type_id = get_entity_type_id('Protein:MI:0326')
        AND ee_protein.taxonomy_id        = 9606
    WHERE re.source_id = get_source_id('stitch')
),

-- Collect raw STITCH compound identifier (PubChem CID) per entity_evidence row.
compound_ids AS NOT MATERIALIZED (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Pubchem Compound:OM:0002') THEN ie.value END) AS pubchem_cid
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('stitch')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw STITCH protein identifier (Ensembl ID) per entity_evidence row; UniProt comes from canonical resolution.
protein_ids AS NOT MATERIALIZED (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Ensembl:MI:0476') THEN ie.value END) AS ensembl_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('stitch')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Pivot STITCH interaction scores and action type flags (binding, inhibition, activation, enzymatic) into columns per relation_evidence row.
rel_annotations AS NOT MATERIALIZED (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Stitch Action Score:OM:1213'    THEN a.value END) AS stitch_action_score,
        MAX(CASE WHEN a.term = 'Confidence Value:OM:1201'       THEN a.value END) AS confidence_value,
        MAX(CASE WHEN a.term = 'Control Type:OM:1212'           THEN a.value END) AS control_type,
        BOOL_OR(a.term = 'Stereospecific:OM:1214')                                AS is_stereospecific,
        BOOL_OR(a.term = 'Binding:OM:0980')                                       AS is_binding,
        BOOL_OR(a.term = 'Enzymatic Reaction:MI:0414')                            AS is_enzymatic,
        BOOL_OR(a.term = 'Inhibition:OM:0931')                                    AS is_inhibition,
        BOOL_OR(a.term = 'Activation:OM:0930')                                    AS is_activation
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('stitch')
      AND rea.relation_evidence_id IN (SELECT relation_evidence_id FROM human_re)
    GROUP BY rea.relation_evidence_id
)

SELECT
    'stitch'::text                  AS source,
    re.relation_evidence_id,
    e_compound.entity_id            AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name                      AS compound_canonical_id_type,
    eer_compound.status_id          AS compound_resolution_status,
    ci.pubchem_cid                  AS compound_pubchem_cid,
    e_protein.entity_id             AS protein_entity_id,
    e_protein.canonical_identifier  AS protein_canonical_id,
    vit_p.name                      AS protein_canonical_id_type,
    eer_protein.status_id           AS protein_resolution_status,
    pi.ensembl_id                   AS protein_ensembl,
    e_protein.canonical_identifier  AS protein_uniprot,
    vrp.name                        AS predicate,
    vrc.name                        AS relation_category,
    ann.stitch_action_score::numeric AS stitch_action_score,
    ann.confidence_value::numeric    AS confidence_value,
    ann.control_type,
    ann.is_stereospecific,
    ann.is_binding,
    ann.is_enzymatic,
    ann.is_inhibition,
    ann.is_activation

FROM human_re re

LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('stitch')
    AND eer_compound.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('stitch')
    AND eer_protein.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = re.relation_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_stitch_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_stitch_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_stitch_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_stitch_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_stitch_relations (protein_uniprot);


-- ────────────────────────────────────────────────────────────────────────────
-- TCDB  (source_id = 41)  — transporter-substrate pairs
-- Direction REVERSED: Protein(3, subject) → SmallMolecule(2, object)
-- Small source (20k rows); human_re applied for consistency.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_tcdb_relations;
CREATE MATERIALIZED VIEW metalinksdb_tcdb_relations AS
WITH

-- Pre-filter to human-only rows; TCDB direction is reversed — protein is the subject.
human_re AS (
    SELECT
        re.relation_evidence_id,
        re.subject_entity_evidence_id,
        re.object_entity_evidence_id,
        re.predicate_id,
        re.relation_category_id
    FROM relation_evidence re
    JOIN entity_evidence ee_compound
        ON  ee_compound.source_id          = re.source_id
        AND ee_compound.entity_evidence_id = re.object_entity_evidence_id   -- compound is OBJECT for TCDB
        AND ee_compound.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_protein
        ON  ee_protein.source_id          = re.source_id
        AND ee_protein.entity_evidence_id = re.subject_entity_evidence_id   -- protein is SUBJECT for TCDB
        AND ee_protein.entity_type_id = get_entity_type_id('Protein:MI:0326')
    JOIN entity_evidence_resolution eer
        ON  eer.source_id          = re.source_id
        AND eer.entity_evidence_id = re.subject_entity_evidence_id
        AND eer.status_id          = 1
    JOIN entity e ON e.entity_id = eer.entity_id AND e.taxonomy_id = 9606
    WHERE re.source_id = get_source_id('tcdb')
),

-- Collect raw TCDB compound identifiers (ChEBI ID, compound name) per entity_evidence row.
compound_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Chebi:MI:0474')  THEN ie.value END) AS chebi_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Name:OM:0202') THEN ie.value END) AS compound_name
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('tcdb')
      AND eei.entity_evidence_id IN (SELECT object_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
),

-- Collect raw TCDB protein identifiers (UniProt, TCDB family ID) per entity_evidence row.
protein_ids AS (
    SELECT eei.entity_evidence_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Uniprot:MI:1097')  THEN ie.value END) AS uniprot_id,
        MAX(CASE WHEN ie.identifier_type_id = get_identifier_type_id('Tcdb:OM:0238') THEN ie.value END) AS tcdb_id
    FROM entity_evidence_identifier eei
    JOIN identifier_evidence ie ON ie.identifier_id = eei.identifier_id
    WHERE eei.source_id = get_source_id('tcdb')
      AND eei.entity_evidence_id IN (SELECT subject_entity_evidence_id FROM human_re)
    GROUP BY eei.entity_evidence_id
)

SELECT
    'tcdb'::text                        AS source,
    re.relation_evidence_id,
    e_compound.entity_id                AS compound_entity_id,
    e_compound.canonical_identifier     AS compound_canonical_id,
    vit_c.name                          AS compound_canonical_id_type,
    eer_compound.status_id              AS compound_resolution_status,
    ci.chebi_id                         AS compound_chebi_id,
    ci.compound_name,
    e_protein.entity_id                 AS protein_entity_id,
    e_protein.canonical_identifier      AS protein_canonical_id,
    vit_p.name                          AS protein_canonical_id_type,
    eer_protein.status_id               AS protein_resolution_status,
    e_protein.canonical_identifier      AS protein_uniprot,
    pi.uniprot_id                       AS protein_raw_uniprot,
    pi.tcdb_id,
    vrp.name                            AS predicate,
    vrc.name                            AS relation_category

FROM human_re re
LEFT JOIN entity_evidence_resolution eer_protein
    ON  eer_protein.source_id = get_source_id('tcdb')
    AND eer_protein.entity_evidence_id = re.subject_entity_evidence_id
LEFT JOIN entity e_protein
    ON  e_protein.entity_id            = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p
    ON  vit_p.identifier_type_id       = e_protein.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON  eer_compound.source_id = get_source_id('tcdb')
    AND eer_compound.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN entity e_compound
    ON  e_compound.entity_id           = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c
    ON  vit_c.identifier_type_id       = e_compound.canonical_identifier_type_id
LEFT JOIN compound_ids ci ON ci.entity_evidence_id = re.object_entity_evidence_id
LEFT JOIN protein_ids  pi ON pi.entity_evidence_id = re.subject_entity_evidence_id
JOIN vocab_relation_predicate vrp ON vrp.relation_predicate_id = re.predicate_id
JOIN vocab_relation_category  vrc ON vrc.relation_category_id  = re.relation_category_id;

CREATE UNIQUE INDEX ON metalinksdb_tcdb_relations (relation_evidence_id);
CREATE INDEX ON metalinksdb_tcdb_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_tcdb_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_tcdb_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_tcdb_relations (protein_uniprot);

-- ────────────────────────────────────────────────────────────────────────────
-- New sources added by 004-metalinksdb-view: RECON3D, Rhea, Human-GEM (transport),
-- CellPhoneDB, NeuronChat (signaling). Their data shapes differ from the older
-- sources: transport edges arrive under the Transport:OM:0035 entity type, and
-- Human-GEM is loaded under the source name 'metatlas'.
-- ────────────────────────────────────────────────────────────────────────────

-- ────────────────────────────────────────────────────────────────────────────
-- RECON3D transport  (source_id varies; metabolite-transporter pairs)
-- Two-hop model: protein --controls--> Transport entity --has_participant--> Chemical.
-- compartment_from/compartment_to come from the participant's 'Subcellular Location'
-- annotation, split by its 'Reactant'/'Product' role within the same Transport event
-- (both are presence-only flags with NULL value -- detect via BOOL_OR, not MAX(value)).
-- No raw UniProt is carried by RECON3D itself (Entrez/gene-name only) -- protein_uniprot
-- is left NULL here and resolved at the combined-view level via gene_protein_representative.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_recon3d_relations;
CREATE MATERIALIZED VIEW metalinksdb_recon3d_relations AS
WITH
transport_control AS (
    SELECT re.relation_evidence_id AS control_re_id,
           re.subject_entity_evidence_id AS protein_entity_evidence_id,
           re.object_entity_evidence_id AS transport_entity_evidence_id
    FROM relation_evidence re
    WHERE re.source_id = get_source_id('recon3d')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'controls')
      AND re.object_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('recon3d') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
),
participant_raw AS (
    SELECT re.relation_evidence_id AS participant_re_id,
           re.subject_entity_evidence_id AS transport_entity_evidence_id,
           re.object_entity_evidence_id AS compound_entity_evidence_id,
           BOOL_OR(a.term = 'Reactant:OM:0310') AS is_reactant,
           BOOL_OR(a.term = 'Product:OM:0311') AS is_product,
           MAX(CASE WHEN a.term = 'Subcellular Location:OM:0604' THEN a.value END) AS compartment
    FROM relation_evidence re
    LEFT JOIN relation_evidence_annotation rea ON rea.source_id = re.source_id AND rea.relation_evidence_id = re.relation_evidence_id
    LEFT JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE re.source_id = get_source_id('recon3d')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'has_participant')
      AND re.subject_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('recon3d') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
    GROUP BY 1,2,3
),
participant_resolved AS (
    SELECT pr.transport_entity_evidence_id,
           eer.entity_id AS compound_entity_id,
           MAX(CASE WHEN pr.is_reactant THEN pr.compartment END) AS compartment_from,
           MAX(CASE WHEN pr.is_product  THEN pr.compartment END) AS compartment_to
    FROM participant_raw pr
    JOIN entity_evidence_resolution eer
        ON eer.source_id = get_source_id('recon3d')
       AND eer.entity_evidence_id = pr.compound_entity_evidence_id
       AND eer.status_id = 1
    GROUP BY 1,2
)
SELECT
    'recon3d'::text AS source,
    tc.control_re_id AS relation_evidence_id,
    e_compound.entity_id AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name AS compound_canonical_id_type,
    1 AS compound_resolution_status,
    e_protein.entity_id AS protein_entity_id,
    e_protein.canonical_identifier AS protein_canonical_id,
    vit_p.name AS protein_canonical_id_type,
    1 AS protein_resolution_status,
    NULL::text AS protein_uniprot,
    pr.compartment_from,
    pr.compartment_to,
    CASE WHEN pr.compartment_from IS NOT NULL AND pr.compartment_to IS NOT NULL THEN 'reactant_to_product'
         WHEN pr.compartment_from IS NOT NULL THEN 'reactant'
         WHEN pr.compartment_to IS NOT NULL THEN 'product' END AS reaction_direction,
    NULL::text AS interaction_type
FROM transport_control tc
JOIN participant_resolved pr ON pr.transport_entity_evidence_id = tc.transport_entity_evidence_id
JOIN entity_evidence_resolution eer_protein
    ON eer_protein.source_id = get_source_id('recon3d')
   AND eer_protein.entity_evidence_id = tc.protein_entity_evidence_id
   AND eer_protein.status_id = 1
JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
JOIN entity e_compound ON e_compound.entity_id = pr.compound_entity_id
LEFT JOIN vocab_identifier_type vit_c ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id;

CREATE INDEX ON metalinksdb_recon3d_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_recon3d_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_recon3d_relations (compound_canonical_id);

-- ────────────────────────────────────────────────────────────────────────────
-- Rhea transport  -- same two-hop model as RECON3D, but compartment comes from
-- 'Membrane Side' (in/out), which is sparser than RECON3D's Subcellular Location
-- (only present on a subset of participant rows) -- compartment_from/to will be
-- NULL more often here; that is expected, not a bug.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_rhea_relations;
CREATE MATERIALIZED VIEW metalinksdb_rhea_relations AS
WITH
transport_control AS (
    SELECT re.relation_evidence_id AS control_re_id,
           re.subject_entity_evidence_id AS protein_entity_evidence_id,
           re.object_entity_evidence_id AS transport_entity_evidence_id
    FROM relation_evidence re
    WHERE re.source_id = get_source_id('rhea')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'controls')
      AND re.object_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('rhea') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
),
participant_raw AS (
    SELECT re.relation_evidence_id AS participant_re_id,
           re.subject_entity_evidence_id AS transport_entity_evidence_id,
           re.object_entity_evidence_id AS compound_entity_evidence_id,
           BOOL_OR(a.term = 'Reactant:OM:0310') AS is_reactant,
           BOOL_OR(a.term = 'Product:OM:0311') AS is_product,
           MAX(CASE WHEN a.term = 'Membrane Side:OM:1231' THEN a.value END) AS compartment
    FROM relation_evidence re
    LEFT JOIN relation_evidence_annotation rea ON rea.source_id = re.source_id AND rea.relation_evidence_id = re.relation_evidence_id
    LEFT JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE re.source_id = get_source_id('rhea')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'has_participant')
      AND re.subject_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('rhea') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
    GROUP BY 1,2,3
),
participant_resolved AS (
    SELECT pr.transport_entity_evidence_id,
           eer.entity_id AS compound_entity_id,
           MAX(CASE WHEN pr.is_reactant THEN pr.compartment END) AS compartment_from,
           MAX(CASE WHEN pr.is_product  THEN pr.compartment END) AS compartment_to
    FROM participant_raw pr
    JOIN entity_evidence_resolution eer
        ON eer.source_id = get_source_id('rhea')
       AND eer.entity_evidence_id = pr.compound_entity_evidence_id
       AND eer.status_id = 1
    GROUP BY 1,2
)
SELECT
    'rhea'::text AS source,
    tc.control_re_id AS relation_evidence_id,
    e_compound.entity_id AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name AS compound_canonical_id_type,
    1 AS compound_resolution_status,
    e_protein.entity_id AS protein_entity_id,
    e_protein.canonical_identifier AS protein_canonical_id,
    vit_p.name AS protein_canonical_id_type,
    1 AS protein_resolution_status,
    NULL::text AS protein_uniprot,
    pr.compartment_from,
    pr.compartment_to,
    CASE WHEN pr.compartment_from IS NOT NULL AND pr.compartment_to IS NOT NULL THEN 'reactant_to_product'
         WHEN pr.compartment_from IS NOT NULL THEN 'reactant'
         WHEN pr.compartment_to IS NOT NULL THEN 'product' END AS reaction_direction,
    NULL::text AS interaction_type
FROM transport_control tc
JOIN participant_resolved pr ON pr.transport_entity_evidence_id = tc.transport_entity_evidence_id
JOIN entity_evidence_resolution eer_protein
    ON eer_protein.source_id = get_source_id('rhea')
   AND eer_protein.entity_evidence_id = tc.protein_entity_evidence_id
   AND eer_protein.status_id = 1
JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
JOIN entity e_compound ON e_compound.entity_id = pr.compound_entity_id
LEFT JOIN vocab_identifier_type vit_c ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id;

CREATE INDEX ON metalinksdb_rhea_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_rhea_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_rhea_relations (compound_canonical_id);

-- ────────────────────────────────────────────────────────────────────────────
-- Human-GEM transport. NOTE: Human-GEM is loaded in pypath/this build under the
-- data_source name 'metatlas' (Metabolic Atlas' unified MAR/MAM namespace; confirmed
-- Human-GEM-shaped -- same Transport/Reaction entity-type split as RECON3D,
-- taxonomy_id NULL throughout, single-organism model). All get_source_id() lookups
-- below read from 'metatlas'; the view itself is named/labelled 'humangem' to match
-- the onboarded source name and the 12-source list. Same two-hop model + Subcellular
-- Location compartment as RECON3D. Protein-side resolution is currently poor
-- (Ensembl-keyed identifiers mostly unresolved -- a resolver-layer gap, out of this
-- view-layer-only spec's scope) so this view is mechanically correct but near-empty
-- today; an empty or near-empty contribution from a source is logged, not a
-- build failure.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_humangem_relations;
CREATE MATERIALIZED VIEW metalinksdb_humangem_relations AS
WITH
transport_control AS (
    SELECT re.relation_evidence_id AS control_re_id,
           re.subject_entity_evidence_id AS protein_entity_evidence_id,
           re.object_entity_evidence_id AS transport_entity_evidence_id
    FROM relation_evidence re
    WHERE re.source_id = get_source_id('metatlas')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'controls')
      AND re.object_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('metatlas') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
),
participant_raw AS (
    SELECT re.relation_evidence_id AS participant_re_id,
           re.subject_entity_evidence_id AS transport_entity_evidence_id,
           re.object_entity_evidence_id AS compound_entity_evidence_id,
           BOOL_OR(a.term = 'Reactant:OM:0310') AS is_reactant,
           BOOL_OR(a.term = 'Product:OM:0311') AS is_product,
           MAX(CASE WHEN a.term = 'Subcellular Location:OM:0604' THEN a.value END) AS compartment
    FROM relation_evidence re
    LEFT JOIN relation_evidence_annotation rea ON rea.source_id = re.source_id AND rea.relation_evidence_id = re.relation_evidence_id
    LEFT JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE re.source_id = get_source_id('metatlas')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'has_participant')
      AND re.subject_entity_evidence_id IN (
          SELECT entity_evidence_id FROM entity_evidence
          WHERE source_id = get_source_id('metatlas') AND entity_type_id = get_entity_type_id('Transport:OM:0035')
      )
    GROUP BY 1,2,3
),
participant_resolved AS (
    SELECT pr.transport_entity_evidence_id,
           eer.entity_id AS compound_entity_id,
           MAX(CASE WHEN pr.is_reactant THEN pr.compartment END) AS compartment_from,
           MAX(CASE WHEN pr.is_product  THEN pr.compartment END) AS compartment_to
    FROM participant_raw pr
    JOIN entity_evidence_resolution eer
        ON eer.source_id = get_source_id('metatlas')
       AND eer.entity_evidence_id = pr.compound_entity_evidence_id
       AND eer.status_id = 1
    GROUP BY 1,2
)
SELECT
    'humangem'::text AS source,
    tc.control_re_id AS relation_evidence_id,
    e_compound.entity_id AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name AS compound_canonical_id_type,
    1 AS compound_resolution_status,
    e_protein.entity_id AS protein_entity_id,
    e_protein.canonical_identifier AS protein_canonical_id,
    vit_p.name AS protein_canonical_id_type,
    1 AS protein_resolution_status,
    NULL::text AS protein_uniprot,
    pr.compartment_from,
    pr.compartment_to,
    CASE WHEN pr.compartment_from IS NOT NULL AND pr.compartment_to IS NOT NULL THEN 'reactant_to_product'
         WHEN pr.compartment_from IS NOT NULL THEN 'reactant'
         WHEN pr.compartment_to IS NOT NULL THEN 'product' END AS reaction_direction,
    NULL::text AS interaction_type
FROM transport_control tc
JOIN participant_resolved pr ON pr.transport_entity_evidence_id = tc.transport_entity_evidence_id
JOIN entity_evidence_resolution eer_protein
    ON eer_protein.source_id = get_source_id('metatlas')
   AND eer_protein.entity_evidence_id = tc.protein_entity_evidence_id
   AND eer_protein.status_id = 1
JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
JOIN entity e_compound ON e_compound.entity_id = pr.compound_entity_id
LEFT JOIN vocab_identifier_type vit_c ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id;

CREATE INDEX ON metalinksdb_humangem_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_humangem_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_humangem_relations (compound_canonical_id);

-- ────────────────────────────────────────────────────────────────────────────
-- ────────────────────────────────────────────────────────────────────────────
-- CellPhoneDB signaling  -- direct Chemical --interacts_with--> Protein edges.
-- No distinct 'interaction type' term exists in the data (only Ligand/Receptor
-- presence flags); 'ligand_receptor' is used as a fixed label since that role
-- pairing is exactly what this predicate represents.
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_cellphonedb_relations;
CREATE MATERIALIZED VIEW metalinksdb_cellphonedb_relations AS
WITH
signaling_re AS (
    SELECT re.relation_evidence_id, re.subject_entity_evidence_id AS compound_entity_evidence_id,
           re.object_entity_evidence_id AS protein_entity_evidence_id
    FROM relation_evidence re
    JOIN entity_evidence ee_c ON ee_c.source_id = re.source_id AND ee_c.entity_evidence_id = re.subject_entity_evidence_id
        AND ee_c.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_p ON ee_p.source_id = re.source_id AND ee_p.entity_evidence_id = re.object_entity_evidence_id
        AND ee_p.entity_type_id = get_entity_type_id('Protein:MI:0326')
    WHERE re.source_id = get_source_id('cellphonedb')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'interacts_with')
)
SELECT
    'cellphonedb'::text AS source,
    sr.relation_evidence_id,
    e_compound.entity_id AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name AS compound_canonical_id_type,
    eer_compound.status_id AS compound_resolution_status,
    e_protein.entity_id AS protein_entity_id,
    e_protein.canonical_identifier AS protein_canonical_id,
    vit_p.name AS protein_canonical_id_type,
    eer_protein.status_id AS protein_resolution_status,
    NULL::text AS protein_uniprot,
    NULL::text AS compartment_from,
    NULL::text AS compartment_to,
    NULL::text AS reaction_direction,
    'ligand_receptor'::text AS interaction_type
FROM signaling_re sr
LEFT JOIN entity_evidence_resolution eer_compound
    ON eer_compound.source_id = get_source_id('cellphonedb') AND eer_compound.entity_evidence_id = sr.compound_entity_evidence_id AND eer_compound.status_id = 1
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON eer_protein.source_id = get_source_id('cellphonedb') AND eer_protein.entity_evidence_id = sr.protein_entity_evidence_id AND eer_protein.status_id = 1
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
WHERE eer_compound.entity_id IS NOT NULL AND eer_protein.entity_id IS NOT NULL;

CREATE INDEX ON metalinksdb_cellphonedb_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_cellphonedb_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_cellphonedb_relations (compound_canonical_id);

-- ────────────────────────────────────────────────────────────────────────────
-- NeuronChat signaling -- direct Chemical --interacts_with--> Protein edges,
-- with a real-valued 'Interaction Type' annotation per row (ligand-receptor,
-- gas-effector, ligand degradation, ligand uptake).
-- ────────────────────────────────────────────────────────────────────────────

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_neuronchat_relations;
CREATE MATERIALIZED VIEW metalinksdb_neuronchat_relations AS
WITH
signaling_re AS (
    SELECT re.relation_evidence_id, re.subject_entity_evidence_id AS compound_entity_evidence_id,
           re.object_entity_evidence_id AS protein_entity_evidence_id
    FROM relation_evidence re
    JOIN entity_evidence ee_c ON ee_c.source_id = re.source_id AND ee_c.entity_evidence_id = re.subject_entity_evidence_id
        AND ee_c.entity_type_id = get_entity_type_id('Chemical:OM:0037')
    JOIN entity_evidence ee_p ON ee_p.source_id = re.source_id AND ee_p.entity_evidence_id = re.object_entity_evidence_id
        AND ee_p.entity_type_id = get_entity_type_id('Protein:MI:0326')
    WHERE re.source_id = get_source_id('neuronchat')
      AND re.predicate_id = (SELECT relation_predicate_id FROM vocab_relation_predicate WHERE name = 'interacts_with')
),
rel_annotations AS (
    SELECT rea.relation_evidence_id,
        MAX(CASE WHEN a.term = 'Interaction Type:OM:1237' THEN a.value END) AS interaction_type
    FROM relation_evidence_annotation rea
    JOIN annotation a ON a.annotation_key = rea.annotation_key
    WHERE rea.source_id = get_source_id('neuronchat')
    GROUP BY 1
)
SELECT
    'neuronchat'::text AS source,
    sr.relation_evidence_id,
    e_compound.entity_id AS compound_entity_id,
    e_compound.canonical_identifier AS compound_canonical_id,
    vit_c.name AS compound_canonical_id_type,
    eer_compound.status_id AS compound_resolution_status,
    e_protein.entity_id AS protein_entity_id,
    e_protein.canonical_identifier AS protein_canonical_id,
    vit_p.name AS protein_canonical_id_type,
    eer_protein.status_id AS protein_resolution_status,
    NULL::text AS protein_uniprot,
    NULL::text AS compartment_from,
    NULL::text AS compartment_to,
    NULL::text AS reaction_direction,
    ann.interaction_type
FROM signaling_re sr
LEFT JOIN rel_annotations ann ON ann.relation_evidence_id = sr.relation_evidence_id
LEFT JOIN entity_evidence_resolution eer_compound
    ON eer_compound.source_id = get_source_id('neuronchat') AND eer_compound.entity_evidence_id = sr.compound_entity_evidence_id AND eer_compound.status_id = 1
LEFT JOIN entity e_compound ON e_compound.entity_id = eer_compound.entity_id
LEFT JOIN vocab_identifier_type vit_c ON vit_c.identifier_type_id = e_compound.canonical_identifier_type_id
LEFT JOIN entity_evidence_resolution eer_protein
    ON eer_protein.source_id = get_source_id('neuronchat') AND eer_protein.entity_evidence_id = sr.protein_entity_evidence_id AND eer_protein.status_id = 1
LEFT JOIN entity e_protein ON e_protein.entity_id = eer_protein.entity_id
LEFT JOIN vocab_identifier_type vit_p ON vit_p.identifier_type_id = e_protein.canonical_identifier_type_id
WHERE eer_compound.entity_id IS NOT NULL AND eer_protein.entity_id IS NOT NULL;

CREATE INDEX ON metalinksdb_neuronchat_relations (compound_entity_id);
CREATE INDEX ON metalinksdb_neuronchat_relations (protein_entity_id);
CREATE INDEX ON metalinksdb_neuronchat_relations (compound_canonical_id);

DROP MATERIALIZED VIEW IF EXISTS metalinksdb_relations;
CREATE MATERIALIZED VIEW metalinksdb_relations AS
SELECT
    base.compound_entity_id,
    base.compound_canonical_id,
    base.compound_canonical_id_type,
    base.protein_entity_id,
    base.protein_canonical_id,
    base.protein_canonical_id_type,
    base.protein_uniprot,
    base.sources,
    base.source_count,
    base.relation_types,
    base.best_pchembl_value,
    base.pubmed_ids,
    base.dois,
    base.compartment_from,
    base.compartment_to,
    base.reaction_direction,
    base.interaction_type,
    base.best_stitch_score,
    base.is_endogenous,
    mca.hmdb_subcellular_locations,
    mca.hmdb_biospecimens,
    mca.hmdb_tissues,
    mca.lipid_category,
    mca.lipid_main_class,
    mca.lipid_sub_class,
    mca.is_lipid,
    mpa.uniprot_subcellular_locations,
    mpa.uniprot_functions,
    mpa.uniprot_disease_involvements,
    mpa.uniprot_ec_numbers,
    mpa.uniprot_protein_families,
    mpa.uniprot_pathway_participations,
    mpa.gtp_functional_classes,
    mpa.gtp_families,
    mpa.tcdb_transporter_families
FROM (
    SELECT
        combined.compound_entity_id,
        MAX(combined.compound_canonical_id)                                 AS compound_canonical_id,
        MAX(combined.compound_canonical_id_type)                            AS compound_canonical_id_type,
        combined.protein_entity_id,
        MAX(combined.protein_canonical_id)                                  AS protein_canonical_id,
        MAX(combined.protein_canonical_id_type)                             AS protein_canonical_id_type,
        COALESCE(MAX(combined.protein_uniprot), MAX(gpr.representative_uniprot)) AS protein_uniprot,
        array_agg(DISTINCT combined.source ORDER BY combined.source)        AS sources,
        COUNT(DISTINCT combined.source)                                     AS source_count,
        array_agg(DISTINCT combined.relation_type ORDER BY combined.relation_type) AS relation_types,
        MAX(combined.pchembl_value)                                         AS best_pchembl_value,
        array_agg(DISTINCT combined.pubmed_id) FILTER (WHERE combined.pubmed_id IS NOT NULL) AS pubmed_ids,
        array_agg(DISTINCT combined.doi)       FILTER (WHERE combined.doi       IS NOT NULL) AS dois,
        MAX(combined.compartment_from)                                      AS compartment_from,
        MAX(combined.compartment_to)                                        AS compartment_to,
        MAX(combined.reaction_direction)                                    AS reaction_direction,
        MAX(combined.interaction_type)                                      AS interaction_type,
        MAX(combined.stitch_action_score)                                   AS best_stitch_score,
        BOOL_OR(combined.is_endogenous)                                     AS is_endogenous
    FROM (
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'interaction'::text AS relation_type, pchembl_value, pubmed_id, doi,
               NULL::text AS compartment_from, NULL::text AS compartment_to, NULL::text AS reaction_direction, NULL::text AS interaction_type,
               NULL::numeric AS stitch_action_score, NULL::boolean AS is_endogenous
        FROM metalinksdb_chembl_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        -- BindingDB DROPPED from the curated combined view (decided 2026-06-15)
        -- -- preserved as-is by 004-metalinksdb-view: per-source matview kept for
        -- provenance, never unioned into the combined contract.

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'interaction'::text, NULL::numeric, pubmed_id, NULL::text,
               NULL, NULL, NULL, NULL, NULL::numeric, NULL::boolean
        FROM metalinksdb_cellinker_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'interaction'::text, NULL::numeric, pubmed_id, NULL::text,
               NULL, NULL, NULL, NULL, NULL::numeric, is_endogenous
        FROM metalinksdb_guidetopharma_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'interaction'::text, NULL::numeric, pubmed_id, NULL::text,
               NULL, NULL, NULL, NULL, NULL::numeric, NULL::boolean
        FROM metalinksdb_mrclinksdb_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'interaction'::text, NULL::numeric, NULL::text, NULL::text,
               NULL, NULL, NULL, NULL, stitch_action_score, NULL::boolean
        FROM metalinksdb_stitch_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'transport'::text,
               NULL::numeric, NULL::text, NULL::text,
               NULL, NULL, NULL, NULL, NULL::numeric, NULL::boolean
        FROM metalinksdb_tcdb_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'transport'::text, NULL::numeric, NULL::text, NULL::text,
               compartment_from, compartment_to, reaction_direction, interaction_type, NULL::numeric, NULL::boolean
        FROM metalinksdb_recon3d_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'transport'::text, NULL::numeric, NULL::text, NULL::text,
               compartment_from, compartment_to, reaction_direction, interaction_type, NULL::numeric, NULL::boolean
        FROM metalinksdb_rhea_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'transport'::text, NULL::numeric, NULL::text, NULL::text,
               compartment_from, compartment_to, reaction_direction, interaction_type, NULL::numeric, NULL::boolean
        FROM metalinksdb_humangem_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'signaling'::text, NULL::numeric, NULL::text, NULL::text,
               compartment_from, compartment_to, reaction_direction, interaction_type, NULL::numeric, NULL::boolean
        FROM metalinksdb_cellphonedb_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1

        UNION ALL
        SELECT source, compound_entity_id, compound_canonical_id, compound_canonical_id_type,
               protein_entity_id, protein_canonical_id, protein_canonical_id_type,
               protein_uniprot, 'signaling'::text, NULL::numeric, NULL::text, NULL::text,
               compartment_from, compartment_to, reaction_direction, interaction_type, NULL::numeric, NULL::boolean
        FROM metalinksdb_neuronchat_relations
        WHERE compound_resolution_status = 1 AND protein_resolution_status = 1
    ) combined
    LEFT JOIN gene_protein_representative gpr ON gpr.entity_id = combined.protein_entity_id
    -- metabolite-class filter (net-new, applied to all 12 sources)
    WHERE combined.compound_entity_id IN (
        SELECT entity_id FROM entity
        WHERE chemical_class_id = (SELECT chemical_class_id FROM vocab_chemical_class WHERE name = 'metabolite')
    )
    GROUP BY combined.compound_entity_id, combined.protein_entity_id
) base
LEFT JOIN metalinksdb_compound_annotations mca ON mca.compound_entity_id = base.compound_entity_id
LEFT JOIN metalinksdb_protein_annotations mpa ON mpa.protein_entity_id = base.protein_entity_id;

CREATE UNIQUE INDEX ON metalinksdb_relations (compound_entity_id, protein_entity_id);
CREATE INDEX ON metalinksdb_relations (compound_canonical_id);
CREATE INDEX ON metalinksdb_relations (protein_uniprot);
CREATE INDEX ON metalinksdb_relations (source_count);
CREATE INDEX ON metalinksdb_relations USING gin (relation_types);

