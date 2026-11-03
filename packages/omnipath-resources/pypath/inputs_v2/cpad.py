"""
CPAD regulator-pathway-cancer associations.

CPAD collects, from the literature, entries that show how a regulator (a gene,
a microRNA, a lncRNA, a transcription factor, a drug) acts on a signaling
pathway in a cancer. Each entry has the regulator's expression pattern
(upregulated or downregulated), the pathway's state (activating or inhibiting),
the outcome for the cancer, the method, species, tissue, cell lines and the
publication. This module maps each entry to the pairs it relates (regulator
and pathway, regulator and cancer, pathway and cancer), each an association.

The file does not say what the regulator is, and the two columns of the
legacy module named ``effect_on_pathway`` and ``effect_on_cancer`` are, on
the website, the regulator's expression pattern and the pathway's state: the
regulator's effect on the pathway is not given, only the two observed states.
So no pair has a direction; the two values are attributes of the pairs they
describe.

The pathways and the cancers are the entities of the ``terms`` dataset, named as
the file names them; the pairs are the ``entries`` dataset.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cpad import iter_entries, iter_terms
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
    id=ResourceCv.CPAD,
    name='CPAD',
    full='Cancer Pathway Associated Database',
    url='http://bio-bigdata.hrbmu.edu.cn/CPAD/',
    # The site states no license ("freely available", news of 2018-06-29); the
    # legacy registry (resources.json) has "No license".
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
    # Zhang et al. 2020, "A comprehensive overview of oncogenic pathways in
    # human cancer", Briefings in Bioinformatics.
    pubmed='31155677',
    topics=(Topic.DISEASE_PHENOTYPE, Topic.SIGNALING),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.TRAIT_ASSOCIATIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-curated entries relating regulators, signaling pathways and cancers, with expression pattern, pathway state, outcome and publication.',
    description=(
        'CPAD (Cancer Pathway Associated Database) is a manually curated '
        'collection of experimentally supported relationships between cancers, '
        'biological pathways and their upstream regulators (microRNAs, lncRNAs, '
        'transcription factors, drugs and others). It documents 4709 entries of '
        '49 cancers and 1557 pathways with 2427 regulators in 7 species, from '
        'about 8500 papers.'
    ),
)

# The server answers 403 to the default user agents of download clients and
# serves the file to one that names itself.
download_entries = Download(
    url='https://bio-bigdata.hrbmu.edu.cn/CPAD/download/CPAD-data.txt',
    filename='CPAD-data.txt',
    subfolder='cpad',
    # The file is in the GBK encoding; the legacy module read it as ISO-8859-1,
    # which turns the Greek letters of the names into two Latin ones.
    encoding='gb18030',
    download_kwargs={'headers': ['User-Agent: pypath-omnipath (+https://omnipathdb.org/)']},
)

# Source-owned vocabularies of the pathways and the cancers, as the file names them.
PATHWAY_NAMESPACE = 'cpad_pathway'
CANCER_NAMESPACE = 'cpad_cancer'

# Every value of the ``Species`` column (4709 entries) with its taxon. The file
# names two genera only, for which no taxon is carried.
SPECIES = {
    'Homo sapiens': 'NCBITaxon:9606',
    'Mus musculus': 'NCBITaxon:10090',
    'Rattus norvegicus': 'NCBITaxon:10116',
    'Canis lupus familiaris': 'NCBITaxon:9615',
    'Danio rerio': 'NCBITaxon:7955',
    'Drosophila': None,
    'Cricetus': None,
}

_REGULATOR_EXPRESSION = 'cpad:regulator_expression'
_PATHWAY_STATE = 'cpad:pathway_state'
_OUTCOME = 'cpad:outcome'
_DISEASE_ID = 'cpad:disease_id'
_KEGG_PATHWAY = 'cpad:kegg_pathway'
_PATHWAY_CATEGORY = 'cpad:pathway_category'
_METHOD = 'cpad:method'
_CELL_LINE = 'cpad:cell_line'
_TISSUE = 'cpad:tissue'

f = FieldConfig()

# The regulator is known by its name only: the file does not say whether it
# is a gene, an RNA or a compound.
regulator = EntityBuilder(
    entity_type=model.NamedThing,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('regulator'))),
)
pathway = EntityBuilder(
    entity_type=model.Pathway,
    identifiers=IdentifiersBuilder(CV(term=PATHWAY_NAMESPACE, value=f('pathway'))),
)
cancer = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(CV(term=CANCER_NAMESPACE, value=f('cancer'))),
)

pathway_terms = EntityBuilder(
    entity_type=model.Pathway,
    identifiers=IdentifiersBuilder(
        CV(term=PATHWAY_NAMESPACE, value=f('name')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
)
cancer_terms = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=CANCER_NAMESPACE, value=f('name')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
)
TERM_SCHEMAS = {'pathway': pathway_terms, 'cancer': cancer_terms}


def map_term(row: dict[str, Any]):
    """The entity of a pathway or a cancer."""
    if row['kind'] not in TERM_SCHEMAS:
        raise ValueError(f'Unknown CPAD term kind: {row["kind"]!r}')
    return TERM_SCHEMAS[row['kind']](row)


def _species(row: dict[str, Any]) -> str | None:
    return SPECIES[row['species']] if row['species'] else None


# What every pair of an entry carries: the study.
_STUDY = (
    CV(term=slots.publications, value=lambda row: [f'PMID:{row["pmid"]}'] if row['pmid'] else []),
    CV(term=slots.in_taxon, value=_species),
    CV(term=_TISSUE, value=f('tissue')),
    CV(term=_CELL_LINE, value=f('cell_lines')),
    CV(term=_METHOD, value=f('methods')),
    CV(term=slots.description, value=f('description')),
)
_PATHWAY = (
    CV(term=_PATHWAY_STATE, value=f('pathway_state')),
    CV(term=_KEGG_PATHWAY, value=f('kegg_pathways')),
    CV(term=_PATHWAY_CATEGORY, value=f('pathway_categories')),
)
_CANCER = (
    CV(term=_OUTCOME, value=f('outcome')),
    CV(term=_DISEASE_ID, value=f('disease_id')),
)
_REGULATOR = (CV(term=_REGULATOR_EXPRESSION, value=f('regulator_expression')),)

SCHEMAS = {
    'regulator_pathway': RelationBuilder(
        subject=regulator,
        predicate=slots.associated_with,
        object=pathway,
        annotations=AnnotationsBuilder(*_REGULATOR, *_PATHWAY, *_STUDY),
    ),
    'regulator_cancer': RelationBuilder(
        subject=regulator,
        predicate=slots.associated_with,
        object=cancer,
        annotations=AnnotationsBuilder(*_REGULATOR, *_CANCER, *_STUDY),
    ),
    'pathway_cancer': RelationBuilder(
        subject=pathway,
        predicate=slots.associated_with,
        object=cancer,
        annotations=AnnotationsBuilder(*_PATHWAY, *_CANCER, *_STUDY),
    ),
}


def map_entry(row: dict[str, Any]):
    """The association of a pair of an entry; an unknown species or pair raises."""
    if row['species'] and row['species'] not in SPECIES:
        raise ValueError(f'Unknown CPAD species: {row["species"]!r}')
    if row['pair'] not in SCHEMAS:
        raise ValueError(f'Unknown CPAD pair: {row["pair"]!r}')
    return SCHEMAS[row['pair']](row)


resource = Resource(
    config,
    terms=Dataset(download=download_entries, mapper=map_term, raw_parser=iter_terms),
    entries=Dataset(download=download_entries, mapper=map_entry, raw_parser=iter_entries),
)
