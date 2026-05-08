"""
Parse MACdb data and emit Entity records.

This module converts annotations of metabolites-cancer associations into Entity
records using the declarative schema pattern.
"""

from functools import partial

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.macdb import iter_associations
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    OntologyDataset
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.internals.ontology_builder import (
    OntologyBuilder,
    RelationshipBuilder
)
from pypath.internals.ontology_schema import OntologyTypedef
from pypath.internals.cv_terms import (
    EntityTypeCv,
    InteractionTypeCv,
    OntologyCv,
    DiseaseAnnotationCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'https://ngdc.cncb.ac.cn/macdb/api/get_download_file?file_type=%s'
TABLES = ['metabolite', 'trait', 'study', 'publication']

config = ResourceConfig(
    id=ResourceCv.MACDB,
    name='MACdb',
    url='https://ngdc.cncb.ac.cn/macdb/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR, # Static?
    pubmed='37027007',
    primary_category='metabolism',
    description=(
        'MACdb is a comprehensive knowledgebase based totally on manual '
        'curation of cancer metabolic literature, providing high quality '
        'metabolite associations, online visualizing networks, and other '
        'effective tools for cancer or metabonomics researches.'
    ),
)

# ================================== DOWNLOAD ==================================

download = dict()

for t in TABLES:

    download[t] = Download(
            url=BASE_URL % t,
            filename=f'downloads.{t}.txt',
            subfolder='macdb',
            large=True,
            ext='.txt',
            default_mode='r',
        )

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'metacID': r'^(METAC_\d+)$',
        'year':  r'^(\d{4})',
    },
    map={},
    transform={},
)

trait_schema = OntologyBuilder(
    id='Trait_Ontology_ID',
    name='Trait_Ontology',
    relationships=[
        RelationshipBuilder(
            type='part_of',
            target=f('Trait_Type'),
        ),
        RelationshipBuilder(
            type='is_a',
            target=f('EFO_ID'),
            target_name=f('EFO_Ontology'),
        ),
    ]
)

association_schema = EntityBuilder(
    entity_type=EntityTypeCv.ASSOCIATION,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.METAC,
            value=f('Cohort_id', extract='metacID')
        ),
        CV(
            term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
            value=f('pubchem_CID')
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionTypeCv.PHENOTYPE_RESULT),
        CV(
            term=IdentifierNamespaceCv.METAC,
            value=f('Cohort_id', extract='metacID')
        ),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_id')),
        CV(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=f('pmc_id')),
        CV(
            term=MoleculeAnnotationsCv.CONCENTRATION_MEAN,
            value=f('case_concentration')
        ),
        CV(
            term=MoleculeAnnotationsCv.CONCENTRATION_MIN,
            value=f('case_concentration_low')
        ),
        CV(
            term=MoleculeAnnotationsCv.CONCENTRATION_MAX,
            value=f('case_concentration_high')
        ),
        CV(
            term=MoleculeAnnotationsCv.CONCENTRATION_SD,
            value=f('case_confidence_interval')
        ),
        CV(
            term=AssayAnnotationsCv.CONTRAST_P_VAL,
            value=f('case_control_p-value')
        ),
        CV(
            term=AssayAnnotationsCv.CONTRAST_LOGFC,
            value=f('log2FC')
        ),
        CV(term=DiseaseAnnotationCv.TYPE, value=f('study_Cancer_type')),
        CV(term=DiseaseAnnotationCv.SUBTYPE, value=f('study_Cancer_subtype')),
        CV(term=AssayAnnotationsCv.CASE_DESCRIPTION, value=f('study_Case_name')),
        CV(term=AssayAnnotationsCv.CASE_AGE, value=f('study_Case_age_group')),
        CV(term=AssayAnnotationsCv.CASE_SEX, value=f('study_Case_sex')),
        CV(term=AssayAnnotationsCv.CASE_SAMPLE_COUNT, value=f('study_Case_size')),
        CV(
            term=AssayAnnotationsCv.CONTROL_DESCRIPTION,
            value=f('study_Control_name')
        ),
        CV(term=AssayAnnotationsCv.CONTROL_AGE, value=f('study_Control_age_group')),
        CV(term=AssayAnnotationsCv.CONTROL_SEX, value=f('study_Control_sex')),
        CV(
            term=AssayAnnotationsCv.CONTROL_SAMPLE_COUNT,
            value=f('study_Control_size')
        ),
        CV(term=AssayAnnotationsCv.DESCRIPTION, value=f('study_Condition')),
        CV(term=AssayAnnotationsCv.CONCLUSION, value=f('study_Conclusion')),
        CV(term=MoleculeAnnotationsCv.EXPERIMENTAL_METHOD, value=f('study_Platform')),
        CV(term=AssayAnnotationsCv.TISSUE, value=f('study_Tissue')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.CV_TERM,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.CV_TERM_ACCESSION,
                        value=f('trait_id'),
                    ),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('trait_name')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=MoleculeSubtypeCv.METABOLITE,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                        value=f('pubchem_CID')
                    ),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('original_metabolite_name')),
                ),
            ),
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    trait=OntologyDataset(
        download=download['trait'],
        mapper=trait_schema,
        raw_parser=iter_tsv,
        ontology_id='macdb_traits',
        remark='MACdb trait ontology exported via pypath.',
        typedefs=[
            OntologyTypedef(id='part_of', name='part_of'),
            OntologyTypedef(id='is_a', name='is_a'),
        ],
        extension='obo',
        file_stem='macdb',
    ),
    associations=Dataset(
        download=download['metabolite'],
        mapper=association_schema,
        raw_parser=partial(
            iter_associations,
            study_download=download['study'],
            publication_download=download['publication'],
            trait_download=download['trait'],
        ),
    ),
)

# ================================= REFERENCE ==================================

# metabolite
#   original_metabolite_name
#   pubchem_CID
#   case_concentration
#   case_concentration_low
#   case_concentration_high
#   case_confidence_interval
#   control_concentration
#   control_concentration_low
#   control_concentration_high
#   control_confidence_interval
#   case_control_p-value
#   Delta_concentration
#   Delta_concentration
#   log2FC
#   Cohort_id
#   e.g. Glutamate	33032	0.802	0.46	1.302		1.175	0.814	1.223		0.403			-0.23	METAC_915

# trait
#   Trait_Ontology
#   Trait_Ontology_ID
#   Trait_Type
#   EFO_Ontology
#   EFO_ID
#   e.g. Acinar Cell Adenocarcinoma (ACA)	1	cancer	acinar cell carcinoma	EFO:0000216

# study
#   Cohort_id
#   Trait_onto_ID
#   Cancer_type
#   Cancer_subtype
#   Cancer_DOID
#   Case_name
#   Case_age_group
#   Case_sex
#   Case_size
#   Control_name
#   Control_age_group
#   Control_sex
#   Control_size
#   Condition
#   Conclusion
#   Platform
#   Tissue
#   Type
#   pubmed_id
#   e.g. METAC_11	128	breast cancer	early breast cancer patients (EBC)	1612	Pre-menopause	42.5 ± 5.5	F	173	Post-menopause	54.3 ± 6.8	F	28	Metabolic characteristics stratified by menopausal status	At diagnosis, compared with post?menopausal women, pre-menopausal patients were more likely to have lower glucose, HbA1c, and triglyceride levels and a lower HOMA score.	Laboratory assays	Blood	case vs.control group	26424165

# publication
#   PMID
#   PMC_ID
#   Title
#   Journal_Title
#   Date_of_Publication
#   Authors_Full_Name
#   Authors_Affiliation
#   Pages
#   Volume
#   Issue
#   e.g. 19845817	PMC3822498	Metabolic profiling reveals key metabolic features of renal cell carcinoma.	Journal of cellular and molecular medicine	2011 Jan	Catchpole, Gareth||Platzer, Alexander||Weikert, Cornelia||Kempkensteffen, Carsten||Johannsen, Manfred||Krause, Hans||Jung, Klaus||Miller, Kurt||Willmitzer, Lothar||Selbig, Joachim||Weikert, Steffen	Department of Central Metabolism, Max-Planck-Institute of Molecular Plant Physiology, Golm, Germany.	109-18	15	1
