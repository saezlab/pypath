"""MACdb chemical–trait associations with study evidence.

Traits retain source IDs; the EFO crosswalk is a topic, not asserted equivalence.
Contrast p-values are typed. Study conclusions and conditions are narrative;
scalar assay values retain source-field context as typed quantities; ambiguous
intervals and nonnumeric assay fields stay in the original evidence payload.
"""

from functools import partial
import math
import re

from omnipath_core.measurements import Measurement

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.macdb import iter_associations
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.internals.cv_terms import LicenseCV, UpdateCategoryCV, ResourceCv

# =================================== SET-UP ===================================

BASE_URL = 'https://ngdc.cncb.ac.cn/macdb/api/get_download_file?file_type=%s'
TABLES = ['metabolite', 'trait', 'study', 'publication']

config = ResourceConfig(
    id=ResourceCv.MACDB,
    name='MACdb',
    url='https://ngdc.cncb.ac.cn/macdb/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR,  # Static?
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
        'year': r'^(\d{4})',
    },
    map={},
    transform={},
)

# MACdb trait IDs are source concepts. The EFO crosswalk is preserved as a topic,
# not identity or subclass, because that stronger assertion is not given by the table.
trait_terms_schema = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MACDB_TRAIT, value=f('Trait_Ontology_ID')),
        CV(term=Namespace.NAME, value=f('Trait_Ontology')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.has_topic, value=f('EFO_ID'))),
)


_QUANTITATIVE_FIELDS = (
    'case_concentration',
    'case_concentration_low',
    'case_concentration_high',
    'case_confidence_interval',
    'control_concentration',
    'control_concentration_low',
    'control_concentration_high',
    'control_confidence_interval',
    'Delta_concentration',
    'log2FC',
    'study_Case_size',
    'study_Control_size',
)
_MEASUREMENT_RE = re.compile(
    r'^\s*(<=|>=|[<>=~≈])?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$'
)


def _measurement(row, field):
    raw = row.get(field)
    match = _MEASUREMENT_RE.fullmatch(str(raw)) if raw is not None else None
    if not match or not math.isfinite(float(match[2])):
        return None
    comparator = match[1]
    relation = {'<': 'less_than', '>': 'greater_than', '=': 'equal_to'}.get(
        comparator
    )
    return Measurement(
        model.QuantityValue(
            has_numeric_value=float(match[2]), has_binary_relation=relation
        ),
        source_field=field,
        comparator=comparator,
    )


association_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.PUBCHEM,
                value=f('pubchem_CID', extract='^(\\d+)$'),
            ),
            CV(term=Namespace.NAME, value=f('original_metabolite_name')),
        ),
        annotations=AnnotationsBuilder(),
    ),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=model.OntologyClass,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.MACDB_TRAIT, value=f('trait_id')),
            CV(term=Namespace.NAME, value=f('trait_name')),
        ),
        annotations=AnnotationsBuilder(),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('study_Conclusion')),
        CV(term=slots.description, value=f('study_Condition')),
        CV(
            term=slots.p_value,
            value=lambda row: _measurement(row, 'case_control_p-value'),
        ),
        CV(
            term=slots.has_quantitative_value,
            value=lambda row: [
                v
                for field in _QUANTITATIVE_FIELDS
                if (v := _measurement(row, field)) is not None
            ],
        ),
        CV(
            term=slots.publications,
            value=f(
                'pubmed_id',
                extract='^(\\d+)$',
                transform=lambda v: 'PMID:' + str(v),
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'pmc_id',
                extract='^(?:PMC)?(\\d+)$',
                transform=lambda v: 'PMC:' + str(v),
            ),
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    trait=Dataset(
        download=download['trait'],
        mapper=trait_terms_schema,
        raw_parser=iter_tsv,
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
