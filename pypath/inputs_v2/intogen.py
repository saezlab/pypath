"""
IntOGen compendium of mutational cancer driver genes.

IntOGen runs several driver discovery methods on the somatic mutations of
cohorts of tumors and combines them. This module reads the catalog of driver
genes of a release (``Compendium_Cancer_Genes.tsv``): a gene is a driver in a
cohort of a cancer type. The statement is the gene ``associated_with`` the
cancer type, with the cohort, the methods, the mode of action and the
significance on it. The cancer types of the cohorts (``cohorts.tsv``) are the
``terms`` dataset.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.intogen import iter_cancer_types, iter_drivers
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.INTOGEN,
    name='IntOGen',
    full='Integrative OncoGenomics',
    url='https://www.intogen.org/',
    # Each archive holds LICENSE.txt, the text of the Creative Commons CC0 1.0
    # Universal legal code (as does the registry of the legacy, resources.json).
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    # Martinez-Jimenez et al., "A compendium of mutational cancer driver genes",
    # Nature Reviews Cancer 2020.
    pubmed='32778778',
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.TRAIT_ASSOCIATIONS),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary='Cancer driver genes by cohort and cancer type, with the methods, the mode of action and the significance.',
    description=(
        'IntOGen is a pipeline and a catalog of the genes whose mutations drive '
        'tumorigenesis. It applies seven driver discovery methods to the somatic '
        'mutations of the tumors of 271 cohorts and combines their '
        'results. For each driver gene of a cohort the catalog gives the cancer '
        'type, the methods that find it, the number of mutations and samples, the '
        'significance of the combined output and the gene\'s mode of action '
        '(activating or loss of function).'
    ),
)

# The 2024.09.20 release, the newest on https://www.intogen.org/download (the
# page offers one archive per release and no alias for the latest). The legacy
# module read the release of 2014.12 (Drivers_type_role.tsv).
download_drivers = Download(
    url='https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip',
    filename='IntOGen-Drivers-20240920.zip',
    subfolder='intogen',
    ext='zip',
)
download_cohorts = Download(
    url='https://www.intogen.org/download?file=IntOGen-Cohorts-20240920.zip',
    filename='IntOGen-Cohorts-20240920.zip',
    subfolder='intogen',
    ext='zip',
)

# The cancer types are the nomenclature and acronyms of the MSKCC OncoTree,
# customized by IntOGen (https://intogen.readthedocs.io/en/latest/postprocessing.html:
# "we used MSKCC oncotree [...] We customized the oncotree"), so the codes are
# OncoTree codes. OncoTree is not imported as an ontology, and the
# long names are those of ``cohorts.tsv``.
CANCER_NAMESPACE = Namespace.ONCOTREE

# Every value of the ``ROLE`` column (release 2024.09.20), carried as written:
# the mode of action of the gene, "LoF (Loss of Function), Act (Activating) or
# amb (Ambiguous)" in the README, ``ambiguous`` in the data. A role that is not
# here raises.
ROLES = frozenset({'Act', 'LoF', 'ambiguous'})
# Every method of the ``METHODS`` column (the first six are the driver
# discovery methods; ``combination`` is the combination of their p-values).
METHODS = frozenset({
    'cbase', 'combination', 'dndscv', 'hotmaps', 'mutpanning', 'oncodriveclustl',
    'oncodrivefml', 'smregions',
})

_COHORT = 'intogen:cohort'
_METHOD = 'intogen:method'
_ROLE = 'intogen:role'
_TRANSCRIPT = 'intogen:transcript'
# Counts, and the q-value of the combined output of the methods (0 to 1).
_MUTATIONS = 'intogen:mutations'
_SAMPLES = 'intogen:samples'
_QVALUE = 'intogen:qvalue_combination'

f = FieldConfig()

# The cancer type: the same entity in the statements and in the ``terms`` dataset.
cancer_type_object = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(CV(term=CANCER_NAMESPACE, value=f('cancer_type'))),
)
cancer_type_terms = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=CANCER_NAMESPACE, value=f('cancer_type')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
)

# Gene symbols of IntOGen are HUGO symbols; the cohorts are of human tumors.
gene = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(CV(term=Namespace.GENESYMBOL, value=f('symbol'))),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value='NCBITaxon:9606')),
)


def _count(column: str):
    return lambda row: measurement(row[column], source_field=column)


driver_schema = RelationBuilder(
    subject=gene,
    predicate=slots.associated_with,
    object=cancer_type_object,
    annotations=AnnotationsBuilder(
        CV(term=_COHORT, value=f('cohort')),
        CV(term=_METHOD, value=f('methods')),
        CV(term=_ROLE, value=f('role')),
        CV(term=_TRANSCRIPT, value=f('transcript')),
        CV(term=_MUTATIONS, value=_count('mutations')),
        CV(term=_SAMPLES, value=_count('samples')),
        CV(term=_QVALUE, value=_count('qvalue')),
    ),
)


def map_driver(row):
    """The statement of a driver gene; an unknown role or method raises."""
    if row['role'] not in ROLES:
        raise ValueError(f'Unknown IntOGen role: {row["role"]!r}')
    unknown = [method for method in row['methods'] if method not in METHODS]
    if unknown:
        raise ValueError(f'Unknown IntOGen method: {unknown[0]!r}')
    if not row['symbol'] or not row['cancer_type']:
        raise ValueError(f'IntOGen driver without gene or cancer type: {row!r}')
    return driver_schema(row)


resource = Resource(
    config,
    terms=Dataset(
        download=download_cohorts,
        mapper=cancer_type_terms,
        raw_parser=iter_cancer_types,
    ),
    drivers=Dataset(
        download=download_drivers,
        mapper=map_driver,
        raw_parser=iter_drivers,
    ),
)
