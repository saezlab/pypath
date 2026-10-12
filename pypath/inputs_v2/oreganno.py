"""
Map the ORegAnno transcription factor binding records and emit Relation records.

ORegAnno (Open Regulatory Annotation) is a catalogue of regulatory regions,
transcription factor binding sites and regulatory variants of many species.
This module takes the records of the combined export (see
``parsers/oreganno.py`` for the file and for which lines are records) that
annotate a binding site of a named transcription factor, with a positive
outcome, to a named gene:

    OREG0000010  Homo sapiens  POSITIVE OUTCOME  TRANSCRIPTION FACTOR BINDING SITE
    MACROD1  ENSG00000133315  ENSEMBL  ESR1  ENSG00000091831  ENSEMBL  ...  15691879

What the record states is that the factor binds a site annotated to the gene,
with a positive result in the experiment of the article. The export has no
technique and no mode of regulation, and binding does not establish a change of
expression, so the relation is ``interacts_with`` without direction or aspect.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.oreganno import iter_binding_sites
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

# Species of the export -> NCBI Taxonomy ID. The file names only the species,
# so this is the exact crosswalk of its names; another name is an error.
SPECIES = {
    'Homo sapiens': '9606',
    'Mus musculus': '10090',
    'Rattus norvegicus': '10116',
    'Gallus gallus': '9031',
    'Danio rerio': '7955',
    'Xenopus tropicalis': '8364',
    'Takifugu rubripes': '31033',
    'Ciona intestinalis': '7719',
    'Ciona savignyi': '51511',
    'Drosophila melanogaster': '7227',
    'Caenorhabditis elegans': '6239',
    'Caenorhabditis briggsae': '6238',
    'Saccharomyces cerevisiae': '4932',
    'Escherichia coli K-12': '83333',
}
# Source of an identifier column -> identifier namespace. USER DEFINED
# identifiers are the curators' own, not accessions; they are not kept.
_ENSEMBL_HUMAN_GENE = re.compile(r'^ENSG\d{11}$')
_ID_SOURCES = frozenset({'NCBI', 'ENSEMBL', 'RegulonDB', 'USER DEFINED'})
# "Dataset" values of the records: the collection a record was imported from.
DATASETS = frozenset(
    {
        'FLYREG',
        'NFIRegulomeDB',
        'NRSF/REST ChIPSeq sites',
        'RegulonDB',
        'STAT1 literature-derived sites',
        'Yeast Regulatory Map',
    }
)


def _taxon(row: dict[str, Any]) -> str:
    species = row.get('species')
    if species not in SPECIES:
        raise ValueError(f'Unknown ORegAnno species: {species!r}')
    return f'NCBITaxon:{SPECIES[species]}'


def _identifier(row: dict[str, Any], role: str, source: str):
    """The identifier of ``role`` (tf or gene) if its source column is ``source``."""
    found = row.get(f'{role}_source')
    if found is not None and found not in _ID_SOURCES:
        raise ValueError(f'Unknown ORegAnno identifier source: {found!r}')
    identifier = row.get(f'{role}_id')
    if found != source or identifier is None:
        return None
    return identifier


def _ensembl_gene(row: dict[str, Any], role: str):
    identifier = _identifier(row, role, 'ENSEMBL')
    return identifier if identifier and _ENSEMBL_HUMAN_GENE.match(identifier) else None


def _ensembl_other(row: dict[str, Any], role: str):
    identifier = _identifier(row, role, 'ENSEMBL')
    return identifier if identifier and not _ENSEMBL_HUMAN_GENE.match(identifier) else None


def _data_source(row: dict[str, Any]):
    dataset = row.get('dataset')
    if dataset is not None and dataset not in DATASETS:
        raise ValueError(f'Unknown ORegAnno dataset: {dataset!r}')
    return dataset


config = ResourceConfig(
    id=ResourceCv.OREGANNO,
    name='ORegAnno',
    full='Open Regulatory Annotation',
    url='http://www.oreganno.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.LGPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='26578589',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Transcription factor binding sites annotated to a target gene, '
        'from the literature and curated data collections, for many species.'
    ),
    description=(
        'ORegAnno is a community-curated, open-access catalogue of regulatory '
        'annotation: regulatory regions, transcription factor binding sites, '
        'regulatory polymorphisms and haplotypes, with the experimental '
        'outcome, the species, the genomic coordinates in several genome '
        'builds and the PubMed ID of the article. Records come from curation '
        'of the literature and from imports of curated collections such as '
        'RegulonDB, FlyReg and the NFI Regulome database. This module '
        'imports the binding sites of a named transcription factor that have '
        'a positive outcome and a named target gene, from the 2016-01-19 '
        'combined export (the original site is offline).'
    ),
)

f = FieldConfig(
    extract={'pmid': PMID_EXTRACT},
    transform={'pmid_curie': lambda v: f'PMID:{v}'},
)


def _participant(entity_type, role: str) -> EntityBuilder:
    symbol = f'{role}_symbol'
    return EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(symbol)),
            CV(
                term=Namespace.ENTREZ,
                value=lambda row: _identifier(row, role, 'NCBI'),
            ),
            CV(
                term=Namespace.ENSG,
                value=lambda row: _ensembl_gene(row, role),
            ),
            CV(
                term=Namespace.ENSEMBL,
                value=lambda row: _ensembl_other(row, role),
            ),
            CV(
                term='regulondb',
                value=lambda row: _identifier(row, role, 'RegulonDB'),
            ),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
    )


# A row without a transcription factor or a target yields no relation: the
# participant has no identifier, so the builder returns None.
interactions_schema = RelationBuilder(
    subject=_participant(model.Protein, 'tf'),
    predicate=slots.interacts_with,
    object=_participant(model.Gene, 'gene'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f('pmid', extract='pmid', transform='pmid_curie'),
        ),
        CV(term=slots.supporting_data_source, value=_data_source),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/ORegAnno_Combined_2016.01.19.tsv',
            filename='ORegAnno_Combined_2016.01.19.tsv',
            subfolder='oreganno',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_binding_sites,
    ),
)
