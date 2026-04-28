"""
Parse MACdb data and emit Entity records.

This module converts annotations of metabolites-cancer associations into Entity
records using the declarative schema pattern.
"""

import os
import re
import requests
from functools import partial
from pathlib import Path

from bs4 import BeautifulSoup

from pypath.inputs_v2.parsers.base import iter_tsv
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
)
from pypath.internals.ontology_builder import (
    OntologyBuilder,
    RelationshipBuilder
)
from pypath.internals.ontology_schema import OntologyDocument, OntologyTypedef
from pypath.internals.cv_terms import (
    EntityTypeCv,
    CurationCv,
    InteractionMetadataCv,
    BiologicalRoleCv,
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
            filename='downloads.%s.txt',
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

metabolite_schema = EntityBuilder(
    entity_type=MoleculeSubtypeCv.METABOLITE,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
            value=f('pubchem_CID')
        ),
        CV(
            term=IdentifierNamespaceCv.METAC,
            value=f('Cohort_id', extract='metacID')
        ),
    ),
    annotations=AnnotationsBuilder(
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
        #CV
        # (term=MoleculeAnnotationsCv.CONCENTRATION_MEAN,
        # value=f('control_concentration')
        #),
        #CV
        # (term=MoleculeAnnotationsCv.CONCENTRATION_MIN,
        # value=f('control_concentration_low')
        #),
        #CV
        # (term=MoleculeAnnotationsCv.CONCENTRATION_MAX,
        # value=f('control_concentration_high')
        #),
        #CV
        # (term=MoleculeAnnotationsCv.CONCENTRATION_SD,
        # value=f('control_confidence_interval')
        #),
        CV(
            term=AssayAnnotationsCv.CONTRAST_P_VAL,
            value=f('case_control_p-value')
        ),
        CV(
            term=AssayAnnotationsCv.CONTRAST_LOGFC,
            value=f('log2FC')
        ),
    ),
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

trait_ontology_document = OntologyDocument(
    ontology='macdb_traits',
    default_namespace='macdb_traits',
    remark='MACdb trait ontology exported via pypath.',
    typedefs=[
        OntologyTypedef(id='part_of', name='part_of'),
        OntologyTypedef(id='is_a', name='is_a')
    ],
)

study_schema = EntityBuilder(
    entity_type=EntityTypeCv.ASSAY,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.METAC,
            value=f('Cohort_id', extract='metacID')
        ),
        CV(term=IdentifierNamespaceCv.DOID, value=f('Cancer_DOID')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            value=f('Trait_onto_ID', extract='metacID'),
        ),
        CV(term=DiseaseAnnotationCv.TYPE, value=f('Cancer_type')),
        CV(term=DiseaseAnnotationCv.SUBTYPE, value=f('Cancer_subtype')),
        CV(term=AssayAnnotationsCv.CASE_DESCRIPTION, value=f('Case_name')),
        CV(term=AssayAnnotationsCv.CASE_AGE, value=f('Case_age_group')),
        CV(term=AssayAnnotationsCv.CASE_SEX, value=f('Case_sex')),
        CV(term=AssayAnnotationsCv.CASE_SAMPLE_COUNT, value=f('Case_size')),
        CV(
            term=AssayAnnotationsCv.CONTROL_DESCRIPTION,
            value=f('Control_name')
        ),
        CV(term=AssayAnnotationsCv.CONTROL_AGE, value=f('Control_age_group')),
        CV(term=AssayAnnotationsCv.CONTROL_SEX, value=f('Control_sex')),
        CV(
            term=AssayAnnotationsCv.CONTROL_SAMPLE_COUNT,
            value=f('Control_size')
        ),
        CV(term=AssayAnnotationsCv.DESCRIPTION, value=f('Condition')),
        CV(term=AssayAnnotationsCv.CONCLUSION, value=f('Conclusion')),
        CV(term=MoleculeAnnotationsCv.EXPERIMENTAL_METHOD, value=f('Platform')),
        CV(term=AssayAnnotationsCv.TISSUE, value=f('Tissue')),
    ),
)

publication_schema = EntityBuilder(
    entity_type=EntityTypeCv.PUBLICATION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PMID')),
        CV(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=f('PMC_ID')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=CurationCv.TITLE, value=f('Title')),
        CV(term=CurationCv.JOURNAL, value=f('Journal_Title')),
        CV(
            term=CurationCv.YEAR,
            value=f('Date_of_Publication', extract='year')
        ),
        CV(
            term=CurationCv.AUTHORS,
            value=f('Authors_Full_Name', delimiter='||')
        ),
        CV(term=CurationCv.PAGES, value=f('Pages')),
        CV(term=CurationCv.VOLUME, value=f('Volume')),
        CV(term=CurationCv.ISSUE, value=f('Issue')),
    ),
)

# ================================= RESOURCE ===================================

kwargs = dict()

for t in TABLES:

    if t == 'trait':

        kwargs[t] = OntologyDataset(
            download=download[t],
            mapper=locals().get('%s_schema' % t),
            raw_parser=iter_tsv,
            document=trait_ontology_document,
            extension='obo',
            file_stem='macdb',
        )

    else:

        kwargs[t] = Dataset(
            download=download[t],
            mapper=locals().get('%s_schema' % t),
            raw_parser=iter_tsv,
        )

resource = Resource(config=config, **kwargs)

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
