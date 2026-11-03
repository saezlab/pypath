"""
DisGeNET curated gene-disease associations.

DisGeNET collects gene-disease associations from expert curated databases and
the literature. This module reads the file of the curated associations (CTD
human, UniProt, ClinGen, CGI, Orphanet, Genomics England and PsyGeNET). The
statement is the gene ``associated_with`` the disease; the disease is its own
entity, identified by its UMLS concept unique identifier and named by the
``terms`` dataset. The databases that report the association are its
supporting data sources.

Downloads from disgenet.org now need an account. The file is the copy that
OmniPath rescued of the last open release.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.disgenet import iter_associations, iter_diseases
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
    id=ResourceCv.DISGENET,
    name='DisGeNET',
    url='https://www.disgenet.com/',
    # The downloads page of the last open release (archived 2024-01-19,
    # https://www.disgenet.org/downloads): "The DisGeNET database is made
    # available under the Attribution-NonCommercial-ShareAlike 4.0 International
    # License". The site now sells plans and states no license for downloads;
    # the legacy registry (resources.json) records the same license.
    license=LicenseCV.CC_BY_NC_SA_4_0,
    update_category=UpdateCategoryCV.STATIC,
    # Pinero et al., "The DisGeNET knowledge platform for disease genomics: 2019
    # update", Nucleic Acids Research 2020.
    pubmed='31680165',
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.TRAIT_ASSOCIATIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Curated gene-disease associations with the DisGeNET score, the supporting databases and publication counts.',
    description=(
        'DisGeNET is a platform of gene-disease and variant-disease associations. '
        'The curated gene-disease associations come from CTD (human), UniProt, '
        'ClinGen, the Cancer Genome Interpreter, Orphanet, Genomics England and '
        'PsyGeNET. Each association has the DisGeNET score, the disease specificity '
        'and pleiotropy indices of the gene, the years of its first and last '
        'publication and the number of publications and variants that support it.'
    ),
)

# https://www.disgenet.org/static/disgenet_ap1/files/downloads/
# curated_gene_disease_associations.tsv.gz now redirects to the plans page.
# OmniPath's rescued copy of the file (dated 2025-08-12; the file names no
# release) is what the legacy module read (``annotations_rescued``); it holds
# the curated set only: the files of the literature, BeFree and all
# associations (the other values of the legacy ``dataset`` argument) are not
# in the rescue folder.
download_associations = Download(
    url='https://rescued.omnipathdb.org/disgenet_curated_gene_disease_associations.tsv.gz',
    filename='disgenet_curated_gene_disease_associations.tsv.gz',
    subfolder='disgenet',
    ext='gz',
)

# The diseases are UMLS concepts (``C0036341``). UMLS is not imported as an ontology.
DISEASE_NAMESPACE = Namespace.UMLS
# Every database of the ``source`` column, as the file writes it.
SOURCES = {
    name: name
    for name in ('CGI', 'CLINGEN', 'CTD_human', 'GENOMICS_ENGLAND', 'ORPHANET', 'PSYGENET', 'UNIPROT')
}
# Every value of the ``diseaseType`` column.
DISEASE_TYPES = frozenset({'disease', 'group', 'phenotype'})

_DISEASE_TYPE = 'disgenet:disease_type'
_SCORE = 'disgenet:score'
_EI = 'disgenet:evidence_index'
_DSI = 'disgenet:dsi'
_DPI = 'disgenet:dpi'
_YEAR_INITIAL = 'disgenet:year_initial'
_YEAR_FINAL = 'disgenet:year_final'
_N_PMIDS = 'disgenet:publication_count'
_N_SNPS = 'disgenet:variant_count'

f = FieldConfig(map={'source': SOURCES})

disease_terms = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=DISEASE_NAMESPACE, value=f('disease_id')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(CV(term=_DISEASE_TYPE, value=f('type'))),
)


def _number(column: str):
    return lambda row: measurement(row[column], source_field=column)


# The gene: the Entrez Gene identifier and the symbol, human.
association_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.ENTREZ, value=f('entrez_id')),
            CV(term=Namespace.GENESYMBOL, value=f('symbol')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value='NCBITaxon:9606')),
    ),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=model.OntologyClass,
        identifiers=IdentifiersBuilder(CV(term=DISEASE_NAMESPACE, value=f('disease_id'))),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.supporting_data_source, value=f('sources', map='source')),
        # The score ranks the associations by the number and the type of sources and of
        # publications (0 to 1); it is not a confidence in the association.
        CV(term=_SCORE, value=_number('score')),
        CV(term=_EI, value=_number('ei')),
        CV(term=_DSI, value=_number('dsi')),
        CV(term=_DPI, value=_number('dpi')),
        CV(term=_YEAR_INITIAL, value=f('year_initial')),
        CV(term=_YEAR_FINAL, value=f('year_final')),
        CV(term=_N_PMIDS, value=_number('n_pmids')),
        CV(term=_N_SNPS, value=_number('n_snps')),
    ),
)


def map_association(row: dict[str, Any]):
    """The statement of an association; an unknown database or disease type raises."""
    unknown = [source for source in row['sources'] if source not in SOURCES]
    if unknown:
        raise ValueError(f'Unknown DisGeNET source database: {unknown[0]!r}')
    if row['disease_type'] not in DISEASE_TYPES:
        raise ValueError(f'Unknown DisGeNET disease type: {row["disease_type"]!r}')
    if not (row['entrez_id'] and row['symbol'] and row['disease_id']):
        raise ValueError(f'DisGeNET association without gene or disease: {row!r}')
    return association_schema(row)


resource = Resource(
    config,
    terms=Dataset(
        download=download_associations,
        mapper=disease_terms,
        raw_parser=iter_diseases,
    ),
    associations=Dataset(
        download=download_associations,
        mapper=map_association,
        raw_parser=iter_associations,
    ),
)
